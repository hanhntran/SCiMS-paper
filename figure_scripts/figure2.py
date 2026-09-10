"""
CAMI simulation benchmark figure — Figure 2.

Layout (three rows, each independent):
  Row a (4 panels):  Stacked bars per tool (correct/incorrect/uncertain) by depth
  Row b (3 panels):  Accuracy curves — pooled | male | female
  Row c (3 panels):  Precision | Recall | F1 (per class, grouped by method)

Recall denominator uses the standard convention:
    recall = TP / (TP + FN + uncertain)
i.e. uncertain calls on the target class count against recall.

Rx and Ry calls are computed from the 95% CI bounds (not the point
estimate), matching the convention of the original Rx/Ry tools.
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec


# ---------------------------------------------------------------------------
# Typography  (single source of truth for all text sizing)
# ---------------------------------------------------------------------------
FONT_BODY    = 7     # ticks, default text
FONT_LABEL   = 9     # axis labels
FONT_TITLE   = 10    # panel titles
FONT_LEGEND  = 8
FONT_CAPTION = 11
FONT_SMALL   = 5     # in-bar fraction labels (must stay small; many per panel)

plt.rcParams.update({
    "font.family":           "sans-serif",
    "font.sans-serif":       ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size":             FONT_BODY,
    "axes.titlesize":        FONT_TITLE,
    "axes.labelsize":        FONT_LABEL,
    "xtick.labelsize":       FONT_BODY,
    "ytick.labelsize":       FONT_BODY,
    "legend.fontsize":       FONT_LEGEND,
    "legend.title_fontsize": FONT_LABEL,
    "figure.titlesize":      FONT_CAPTION,
    "svg.fonttype":          "none",
})


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

PATHS = {
    "scims": "hg38_cami_scims_out/hg38_cami_simulation_metadata_scims_updated.txt",
    "bexy":  "bexy/hg38_cami_simulation_bexy_out.txt",
    "rxry":  "rxry/results/hg38_cami_simulation_rxry_output.txt",
}

OUTPUT_DIR = "./figures"

METHODS = ["SCiMS", "BeXY", "Rx", "Ry"]
CLASSES = ["male", "female"]

DEPTHS = [150, 250, 350, 450, 1000, 10000]
DEPTH_LABELS = [f"{d:,}" for d in DEPTHS]

METHOD_COLORS = {
    "SCiMS": "#58508f",
    "BeXY":  "#bd5090",
    "Rx":    "#ac9546",
    "Ry":    "#eddca5",
}
METHOD_MARKERS = {
    "SCiMS": "D", "BeXY": "h", "Rx": "o", "Ry": "s",
}
OUTCOME_COLORS = {
    "correct":   "#0099CC",
    "incorrect": "#FF6633",
    "uncertain": "#c3c3c3",
}

# label fraction segments at least this tall (avoids unreadable slivers)
LABEL_MIN_FRAC = 0.05


# ---------------------------------------------------------------------------
# Prediction parsing (uses 95% CI columns)
# ---------------------------------------------------------------------------

def _parse_ci(ci_str):
    if not isinstance(ci_str, str):
        return None
    s = ci_str.strip("()")
    if "nan" in s.lower() or "inf" in s.lower():
        return None
    try:
        parts = [p.strip() for p in s.split(",")]
        if len(parts) < 2:
            return None
        return float(parts[0]), float(parts[1])
    except (TypeError, ValueError):
        return None


def get_rx_prediction(ci_str):
    parsed = _parse_ci(ci_str)
    if parsed is None: return "uncertain"
    low, high = parsed
    if low  > 0.8: return "female"
    if high < 0.6: return "male"
    return "uncertain"


def get_ry_prediction(ci_str):
    parsed = _parse_ci(ci_str)
    if parsed is None: return "uncertain"
    low, high = parsed
    if low  > 0.077: return "male"
    if high < 0.016: return "female"
    return "uncertain"


def wilson_sem(p, n):
    if n == 0: return 0.0
    return np.sqrt(p * (1 - p) / n)


# ---------------------------------------------------------------------------
# Load + merge
# ---------------------------------------------------------------------------

def load_and_merge():
    scims = pd.read_csv(PATHS["scims"], sep="\t")
    base = scims[["Sample", "sex", "host_depth", "host_fraction", "SCiMS_sex"]].copy()
    base.columns = ["Sample", "actual_sex", "host_depth", "host_fraction", "SCiMS"]
    base["SCiMS"]      = base["SCiMS"].astype(str).str.strip().str.lower()
    base["actual_sex"] = base["actual_sex"].astype(str).str.strip().str.lower()

    rxry = pd.read_csv(PATHS["rxry"], sep="\t")
    rxry.rename(columns={"Sample": "Sample"}, inplace=True)
    rxry["Rx_call"] = rxry["Rx 95% CI"].apply(get_rx_prediction)
    rxry["Ry_call"] = rxry["Ry 95% CI"].apply(get_ry_prediction)
    rxry = rxry[["Sample", "Rx_call", "Ry_call"]]

    bexy = pd.read_csv(PATHS["bexy"], sep="\t")
    bexy["sample"] = bexy["sample"].astype(str).str.replace(".sorted", "", regex=False)
    bexy.rename(columns={"sample": "Sample"}, inplace=True)
    bexy_map = {"XX": "female", "XY": "male"}
    bexy["BeXY"] = bexy["sex_karyotype"].map(bexy_map).fillna("uncertain")
    bexy = bexy[["Sample", "BeXY"]]

    df = base.merge(rxry, on="Sample", how="left")
    df = df.merge(bexy, on="Sample", how="left")
    df.rename(columns={"Rx_call": "Rx", "Ry_call": "Ry"}, inplace=True)
    df[["BeXY", "Rx", "Ry"]] = df[["BeXY", "Rx", "Ry"]].fillna("uncertain")
    return df


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------

def compute_accuracy_by_depth(df, sex_filter=None):
    sub_df = df if sex_filter is None else df[df["actual_sex"] == sex_filter]
    rows = []
    for d in DEPTHS:
        sub = sub_df[sub_df["host_depth"] == d]
        for method in METHODS:
            n = len(sub)
            if n == 0:
                acc, sem = np.nan, np.nan
            else:
                correct = (sub[method] == sub["actual_sex"]).sum()
                acc = correct / n
                sem = wilson_sem(acc, n)
            rows.append({"depth": d, "method": method,
                         "accuracy": acc, "sem": sem, "n": n})
    return pd.DataFrame(rows)


def compute_per_class_metrics(df):
    rows = []
    actual = df["actual_sex"]
    for method in METHODS:
        preds = df[method]
        for cls in CLASSES:
            is_cls = actual == cls
            n_correct   = (is_cls & (preds == cls)).sum()
            n_incorrect = (is_cls & preds.isin(CLASSES) & (preds != cls)).sum()
            n_uncertain = (is_cls & ~preds.isin(CLASSES)).sum()
            n_fp        = (~is_cls & actual.isin(CLASSES) & (preds == cls)).sum()

            prec_den = n_correct + n_fp
            rec_den  = n_correct + n_incorrect + n_uncertain

            precision = n_correct / prec_den if prec_den > 0 else np.nan
            recall    = n_correct / rec_den  if rec_den  > 0 else np.nan
            if (np.isnan(precision) or np.isnan(recall) or
                (precision + recall) == 0):
                f1 = np.nan
            else:
                f1 = 2 * precision * recall / (precision + recall)

            rows.append({
                "method": method, "class": cls,
                "precision": precision, "recall": recall, "f1-score": f1,
                "n_correct":   int(n_correct),
                "n_incorrect": int(n_incorrect),
                "n_uncertain": int(n_uncertain),
                "n_fp":        int(n_fp),
            })
    return pd.DataFrame(rows)


def compute_stacked_breakdown(df):
    rows = []
    for d in DEPTHS:
        sub = df[df["host_depth"] == d]
        n = len(sub)
        for method in METHODS:
            if n == 0:
                rows.append({"depth": d, "method": method,
                             "correct_frac": np.nan, "incorrect_frac": np.nan,
                             "uncertain_frac": np.nan, "n": 0})
                continue
            preds = sub[method]
            actual = sub["actual_sex"]
            n_uncertain = (~preds.isin(CLASSES)).sum()
            n_correct   = (preds.isin(CLASSES) & (preds == actual)).sum()
            n_incorrect = n - n_uncertain - n_correct
            rows.append({
                "depth": d, "method": method,
                "correct_frac":   n_correct   / n,
                "incorrect_frac": n_incorrect / n,
                "uncertain_frac": n_uncertain / n,
                "n": n,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_accuracy_curve(ax, accuracy_df, title, show_ylabel=True,
                         show_yticklabels=True):
    x_pos = np.arange(len(DEPTHS))
    for method in METHODS:
        sub = accuracy_df[accuracy_df["method"] == method].set_index("depth")
        sub = sub.reindex(DEPTHS)
        y    = sub["accuracy"].values
        yerr = sub["sem"].values
        valid = ~np.isnan(y)

        ax.plot(x_pos[valid], y[valid],
                color=METHOD_COLORS[method], marker=METHOD_MARKERS[method],
                linewidth=1.5, markersize=5, label=method)
        ax.fill_between(x_pos[valid],
                        (y - yerr)[valid], (y + yerr)[valid],
                        alpha=0.25, color=METHOD_COLORS[method])

    ax.set_xticks(x_pos)
    ax.set_xticklabels(DEPTH_LABELS, rotation=45, ha="right")
    ax.set_xlabel("Host reads", weight="bold")
    ax.set_ylim(0, 1.03)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
    ax.set_title(title, weight="bold")
    if show_ylabel:
        ax.set_ylabel("Accuracy", weight="bold")
    if not show_yticklabels:
        ax.set_yticklabels([])
    if ax.get_legend() is not None:
        ax.get_legend().remove()


def _plot_metric_bars(ax, long_df, metric_name, show_ylabel=True,
                      show_yticklabels=True):
    sns.barplot(
        data=long_df[long_df["metric"] == metric_name],
        x="class", y="score", hue="method",
        palette=METHOD_COLORS, edgecolor="black", linewidth=0.6, ax=ax,
    )
    ax.set_title(metric_name.capitalize(), weight="bold")
    ax.set_xlabel("")
    ax.set_ylim(0, 1.03)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5)
    if show_ylabel:
        ax.set_ylabel("Score", weight="bold")
    else:
        ax.set_ylabel("")
    if not show_yticklabels:
        ax.set_yticklabels([])
    if ax.get_legend() is not None:
        ax.get_legend().remove()


def plot_figure(accuracy_pooled, accuracy_male, accuracy_female,
                metrics_df, stacked_df, out_path):
    # Three-row layout. Each row is independent — different column counts:
    #   Row a: 4 panels (one per tool)
    #   Row b: 3 panels (pooled / male / female accuracy)
    #   Row c: 3 panels (precision / recall / F1)
    #
    # We use a 12-column underlying grid:
    #   Row a panels span 3 columns each (4 panels x 3 = 12)
    #   Row b panels span 4 columns each (3 panels x 4 = 12)
    #   Row c panels span 4 columns each (3 panels x 4 = 12)
    s     = 6.69 / 4
    fig_w_in = 180 / 25.4    # 180 mm = ~7.087 inches
    fig_h_in = fig_w_in * 0.8  # adjust this ratio to taste
    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=300)
    gs = GridSpec(
        3, 12, figure=fig,
        hspace=0.55, wspace=0.40,
        left=0.07, right=0.985, top=0.93, bottom=0.08,
    )

    # ------------------------------------------------------------------
    # ROW A: stacked bars per tool
    # ------------------------------------------------------------------
    x_pos = np.arange(len(DEPTHS))
    bar_w = 0.8
    ax_a_list = []

    for i, method in enumerate(METHODS):
        ax = fig.add_subplot(gs[0, i*3:(i+1)*3])
        ax_a_list.append(ax)
        sub = stacked_df[stacked_df["method"] == method].set_index("depth")
        sub = sub.reindex(DEPTHS)

        correct   = np.nan_to_num(sub["correct_frac"].values,   nan=0.0)
        incorrect = np.nan_to_num(sub["incorrect_frac"].values, nan=0.0)
        uncertain = np.nan_to_num(sub["uncertain_frac"].values, nan=0.0)

        ax.bar(x_pos, correct, width=bar_w,
               color=OUTCOME_COLORS["correct"],
               edgecolor="black", linewidth=0.4, label="Correct")
        ax.bar(x_pos, incorrect, width=bar_w, bottom=correct,
               color=OUTCOME_COLORS["incorrect"],
               edgecolor="black", linewidth=0.4, label="Incorrect")
        ax.bar(x_pos, uncertain, width=bar_w, bottom=correct + incorrect,
               color=OUTCOME_COLORS["uncertain"],
               edgecolor="black", linewidth=0.4, label="Uncertain")

        # --- fraction labels centered in each segment ---
        for xi in range(len(x_pos)):
            segments = [
                (correct[xi],   0.0,                         OUTCOME_COLORS["correct"]),
                (incorrect[xi], correct[xi],                 OUTCOME_COLORS["incorrect"]),
                (uncertain[xi], correct[xi] + incorrect[xi], OUTCOME_COLORS["uncertain"]),
            ]
            for frac, base, color in segments:
                if frac > LABEL_MIN_FRAC:
                    ax.text(
                        x_pos[xi], base + frac / 2,
                        f"{frac:.2f}",
                        ha="center", va="center", fontsize=FONT_SMALL,
                        color="black" if color == OUTCOME_COLORS["uncertain"] else "white",
                    )

        ax.set_xticks(x_pos)
        ax.set_xticklabels(DEPTH_LABELS, rotation=45, ha="right")
        ax.set_ylim(0, 1.03)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_title(method, weight="bold")
        ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5)
        if i == 0:
            ax.set_ylabel("Fraction of samples", weight="bold")
            ax.legend(loc="lower right", frameon=False)
        else:
            ax.set_yticklabels([])

    # ------------------------------------------------------------------
    # ROW B: accuracy curves — pooled | male | female
    # ------------------------------------------------------------------
    ax_pooled = fig.add_subplot(gs[1, 0:4])
    ax_male   = fig.add_subplot(gs[1, 4:8])
    ax_female = fig.add_subplot(gs[1, 8:12])

    _plot_accuracy_curve(ax_pooled, accuracy_pooled, "Pooled accuracy",
                         show_ylabel=True,  show_yticklabels=True)
    _plot_accuracy_curve(ax_male,   accuracy_male,   "Male accuracy",
                         show_ylabel=False, show_yticklabels=False)
    _plot_accuracy_curve(ax_female, accuracy_female, "Female accuracy",
                         show_ylabel=False, show_yticklabels=False)

    # ------------------------------------------------------------------
    # ROW C: precision | recall | F1
    # ------------------------------------------------------------------
    ax_prec = fig.add_subplot(gs[2, 0:4])
    ax_rec  = fig.add_subplot(gs[2, 4:8])
    ax_f1   = fig.add_subplot(gs[2, 8:12])

    long = metrics_df.melt(
        id_vars=["method", "class"],
        value_vars=["precision", "recall", "f1-score"],
        var_name="metric", value_name="score",
    )
    _plot_metric_bars(ax_prec, long, "precision",
                      show_ylabel=True,  show_yticklabels=True)
    _plot_metric_bars(ax_rec,  long, "recall",
                      show_ylabel=False, show_yticklabels=False)
    _plot_metric_bars(ax_f1,   long, "f1-score",
                      show_ylabel=False, show_yticklabels=False)

    # ------------------------------------------------------------------
    # Cosmetic + global legend
    # ------------------------------------------------------------------
    all_axes = [*ax_a_list, ax_pooled, ax_male, ax_female,
                ax_prec, ax_rec, ax_f1]
    for ax in all_axes:
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)

    handles, labels = ax_pooled.get_legend_handles_labels()
    fig.legend(handles, labels, title="Method",
               ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.00),
               frameon=False)

    fig.savefig(out_path, bbox_inches="tight",
                transparent=True, format="eps", dpi=300)
    png_path = out_path.replace(".eps", ".png")
    fig.savefig(png_path, bbox_inches="tight", dpi=200)
    print(f"  wrote {out_path} (+ {png_path})")
    plt.show()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading + merging…")
    df = load_and_merge()
    print(f"  {len(df)} samples merged across tools")

    print("\nComputing pooled accuracy by depth…")
    accuracy_pooled = compute_accuracy_by_depth(df)

    print("\nComputing male-only accuracy by depth…")
    accuracy_male = compute_accuracy_by_depth(df, sex_filter="male")

    print("\nComputing female-only accuracy by depth…")
    accuracy_female = compute_accuracy_by_depth(df, sex_filter="female")

    print("\nComputing per-class metrics (uncertain hurts recall)…")
    metrics_df = compute_per_class_metrics(df)
    print(metrics_df.round(3).to_string(index=False))

    print("\nComputing stacked outcome breakdown…")
    stacked_df = compute_stacked_breakdown(df)

    print("\nPlotting…")
    out_path = f"{OUTPUT_DIR}/fig2.svg"
    plot_figure(accuracy_pooled, accuracy_male, accuracy_female,
                metrics_df, stacked_df, out_path)


if __name__ == "__main__":
    main()