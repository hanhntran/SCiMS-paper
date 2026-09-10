"""
Animal benchmark figure: one row per species, four panels per row.

  Col 0: classification outcome per tool (stacked: correct/incorrect/uncertain)
  Col 1: Precision (per sex, grouped by tool)
  Col 2: Recall    (per sex, grouped by tool)
  Col 3: F1-score  (per sex, grouped by tool)

Rows: one per species (Mouse, Cow, Baboon, Black Rhino, Pig,
      Mesquite Lizard, Chicken).

Recall convention (full denominator, matching all other scripts):
    recall = n_correct / (n_correct + n_incorrect + n_uncertain)
    precision = n_correct / (n_correct + n_incorrect)

NOTE: ground-truth sex is taken from host_sex (metadata). For ZW species
(chicken) and the lizard, ensure the tool-call columns in the input files
already use the correct (inverted where appropriate) Male/Female labels,
or accuracy for those species will be wrong.
"""

import re
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


# ---------------------------------------------------------------------------
# Typography
# ---------------------------------------------------------------------------
FONT_BODY    = 7
FONT_CAPTION = 11
mpl.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size":        FONT_BODY,
    "axes.titlesize":   FONT_CAPTION,
    "axes.labelsize":   FONT_BODY,
    "xtick.labelsize":  FONT_BODY,
    "ytick.labelsize":  FONT_BODY,
    "figure.titlesize": FONT_CAPTION,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "svg.fonttype":     "none",
})


# ---------------------------------------------------------------------------
# Config: one entry per species
# ---------------------------------------------------------------------------
ORGANISMS = {
    "Mouse":           {"scims": "mouse/mouse_scims_out/mouse_metadata_scims_updated.txt",
                        "bexy":  "mouse/bexy/mouse_bexy_output.txt",
                        "rxry":  "mouse/rxry/results/mouse_rxry_output.txt"},
    "Cow":             {"scims": "cow/cow_scims_out/cow_metadata_scims_updated.txt",
                        "bexy":  "cow/bexy/cow_bexy_out.txt",
                        "rxry":  "cow/rxry/results/cow_rxry_output.txt"},
    "Baboon":          {"scims": "baboon/baboon_scims_out/baboon_metagenomics_scims_updated.txt",
                        "bexy":  "baboon/bexy/baboon_bexy_out.txt",
                        "rxry":  "baboon/rxry/results/baboon_rxry_output.txt"},
    "Black Rhino":     {"scims": "black_rhino/black_rhino_scims_out/black_rhino_metadata_scims_updated.txt",
                        "bexy":  "black_rhino/bexy/black_rhino_bexy_out.txt",
                        "rxry":  "black_rhino/rxry/results/black_rhino_rxry_output.txt"},
    "Pig":             {"scims": "pig/pig_scims_out/pig_metadata_scims_updated.txt",
                        "bexy":  "pig/bexy/pig_bexy_out.txt",
                        "rxry":  "pig/rxry/results/pig_rxry_output.txt"},
    "Mesquite Lizard": {"scims": "mesquite_lizard/mesquite_lizard_scims_out/mesquite_lizard_metadata_scims_updated.txt",
                        "bexy":  "mesquite_lizard/bexy/mesquite_lizard_bexy_out.txt",
                        "rxry":  "mesquite_lizard/rxry/results/mesquite_lizard_rxry_output.txt"},
    "Chicken":         {"scims": "chicken/chicken_scims_out/chicken_metadata_scims_updated.txt",
                        "bexy":  "chicken/bexy/all_chicken_bexy_output.txt",
                        "rxry":  "chicken/rxry/results/all_chicken_rxry_output.txt"},
}

CLASSES = ["Male", "Female"]
TOOLS   = ["SCiMS", "BeXY", "Rx", "Ry"]
CATEGORY_ORDER = ["correct", "incorrect", "uncertain"]

OUTCOME_COLORS = {"correct": "#0099cc", "incorrect": "#f08c2f", "uncertain": "#c3c3c3"}
TOOL_COLORS    = {"SCiMS": "#58508f", "BeXY": "#bd5090", "Rx": "#ac9546", "Ry": "#eddca5"}

OUTPUT_DIR = Path("./figures")


# ---------------------------------------------------------------------------
# Call extraction
# ---------------------------------------------------------------------------
def _norm_sex(x):
    if not isinstance(x, str): return np.nan
    x = x.strip().lower()
    if x in {"male", "m"}:   return "Male"
    if x in {"female", "f"}: return "Female"
    return np.nan

def _parse_ci(ci_str):
    if not isinstance(ci_str, str): return None
    m = re.findall(r"-?\d+\.?\d*", ci_str)
    if len(m) < 2: return None
    return float(m[0]), float(m[1])

def call_bexy(k):
    if k == "XX": return "Female"
    if k == "XY": return "Male"
    return "uncertain"

def call_rx(ci):
    p = _parse_ci(ci)
    if p is None: return "uncertain"
    low, high = p
    if low  > 0.8: return "Female"
    if high < 0.6: return "Male"
    return "uncertain"

def call_ry(ci):
    p = _parse_ci(ci)
    if p is None: return "uncertain"
    low, high = p
    if low  > 0.077: return "Male"
    if high < 0.016: return "Female"
    return "uncertain"


def load_organism(paths):
    scims = pd.read_csv(paths["scims"], sep="\t").rename(columns={"Run": "Sample"})
    scims["host_sex"] = scims["host_sex"].map(_norm_sex)
    scims["SCiMS"]    = scims["SCiMS_sex"].map(_norm_sex).fillna("uncertain")
    base = scims[["Sample", "host_sex", "SCiMS"]].copy()

    bexy = pd.read_csv(paths["bexy"], sep="\t").rename(columns={"sample": "Sample"})
    bexy["Sample"] = bexy["Sample"].astype(str).str.replace(".sorted", "", regex=False)
    bexy["BeXY"]   = bexy["sex_karyotype"].map(call_bexy)

    rxry = pd.read_csv(paths["rxry"], sep="\t").rename(columns={"Run": "Sample"})
    rxry["Rx"] = rxry["Rx 95% CI"].map(call_rx)
    rxry["Ry"] = rxry["Ry 95% CI"].map(call_ry)

    n_scims = len(base)
    # INNER join: only keep samples that ALL tools were actually run on.
    # (A left join would turn samples a tool never processed into spurious
    #  'uncertain' calls — e.g. cow/pig, where BeXY/Rx-Ry cover a subset.)
    df = base.merge(bexy[["Sample", "BeXY"]], on="Sample", how="inner")
    df = df.merge(rxry[["Sample", "Rx", "Ry"]], on="Sample", how="inner")
    # After an inner join there should be no missing tool calls. Any remaining
    # NaN would be a genuine parsing gap; map it to 'uncertain' explicitly.
    for t in ["BeXY", "Rx", "Ry"]:
        df[t] = df[t].fillna("uncertain")
    df = df.dropna(subset=["host_sex"])
    n_common = len(df)
    if n_common < n_scims:
        print(f"    [note] {n_scims} SCiMS samples -> {n_common} common to all "
              f"tools ({n_scims - n_common} dropped: not run by BeXY/Rx/Ry)")
    return df


# ---------------------------------------------------------------------------
# Per-species computations
# ---------------------------------------------------------------------------
def outcome_fractions(df):
    """Fraction correct/incorrect/uncertain per tool for one species."""
    out = {}
    for tool in TOOLS:
        preds = df[tool]
        total = len(df)
        unc = (preds == "uncertain").sum()
        cor = ((preds == df["host_sex"]) & (preds != "uncertain")).sum()
        inc = total - cor - unc
        out[tool] = {"correct": cor / total if total else 0,
                     "incorrect": inc / total if total else 0,
                     "uncertain": unc / total if total else 0}
    return out


def per_class_metrics(df):
    rows = []
    actual = df["host_sex"]
    for method in TOOLS:
        preds = df[method]
        for cls in CLASSES:
            is_cls = actual == cls
            n_correct   = (is_cls & (preds == cls)).sum()
            n_incorrect = (is_cls & preds.isin(CLASSES) & (preds != cls)).sum()
            n_uncertain = (is_cls & ~preds.isin(CLASSES)).sum()
            n_fp        = (~is_cls & actual.isin(CLASSES) & (preds == cls)).sum()

            prec_den = n_correct + n_fp
            rec_den  = n_correct + n_incorrect + n_uncertain      # FULL denominator
            precision = n_correct / prec_den if prec_den > 0 else np.nan
            recall    = n_correct / rec_den  if rec_den  > 0 else np.nan
            if np.isnan(precision) or np.isnan(recall) or (precision + recall) == 0:
                f1 = np.nan
            else:
                f1 = 2 * precision * recall / (precision + recall)
            rows.append({"method": method, "class": cls, "precision": precision,
                         "recall": recall, "f1-score": f1})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def make_figure(species_data, outpath):
    species = list(species_data.keys())
    n_sp = len(species)

    # Fit one journal page: 180 mm wide, height capped so 7 rows fit on a page.
    fig_w = 180 / 25.4                      # 180 mm  ≈ 7.09 in
    row_h = 1.18                            # in per species row
    fig_h = min(row_h * n_sp + 0.6, 9.6)    # cap ~9.6 in to stay on one page
    fig, axes = plt.subplots(n_sp, 4, figsize=(fig_w, fig_h), dpi=300)
    if n_sp == 1:
        axes = axes[np.newaxis, :]

    x = np.arange(len(TOOLS))

    for r, sp in enumerate(species):
        is_bottom = (r == n_sp - 1)
        df = species_data[sp]
        outcomes = outcome_fractions(df)
        metrics  = per_class_metrics(df)
        long = metrics.melt(id_vars=["method", "class"],
                            value_vars=["precision", "recall", "f1-score"],
                            var_name="metric", value_name="score")

        # --- Col 0: stacked outcome bars per tool ---
        ax0 = axes[r, 0]
        bottoms = np.zeros(len(TOOLS))
        for cat in CATEGORY_ORDER:
            heights = [outcomes[t][cat] for t in TOOLS]
            ax0.bar(x, heights, bottom=bottoms, color=OUTCOME_COLORS[cat],
                    width=0.6, edgecolor="black", linewidth=0.4,
                    label=cat.capitalize() if r == 0 else "")
            bottoms += np.array(heights)
        ax0.set_xticks(x)
        if is_bottom:
            ax0.set_xticklabels(TOOLS, fontsize=6, rotation=45, ha="right")
        else:
            ax0.set_xticklabels([])
        ax0.set_ylim(0, 1.05)
        ax0.set_ylabel(sp, fontsize=8, fontweight="bold")
        ax0.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5)
        ax0.tick_params(axis="y", labelsize=6)
        if r == 0:
            ax0.set_title("Outcome", fontsize=8, fontweight="bold")

        # --- Cols 1-3: precision / recall / F1 ---
        for j, met in enumerate(["precision", "recall", "f1-score"], start=1):
            ax = axes[r, j]
            sns.barplot(data=long[long["metric"] == met], x="class", y="score",
                        hue="method", palette=TOOL_COLORS, edgecolor="black",
                        linewidth=0.4, ax=ax)
            ax.set_ylim(0, 1.05)
            ax.set_xlabel(""); ax.set_ylabel("")
            ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5)
            if is_bottom:
                ax.tick_params(axis="x", labelsize=6)
            else:
                ax.set_xticklabels([])
            ax.tick_params(axis="y", labelsize=6)
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            if r == 0:
                ax.set_title(met.capitalize().replace("-score", " score"),
                             fontsize=8, fontweight="bold")

    # global legends: outcome (from col0 row0) + tool (build manually)
    out_handles = [plt.Rectangle((0, 0), 1, 1, facecolor=OUTCOME_COLORS[c],
                                 edgecolor="black", linewidth=0.4,
                                 label=c.capitalize()) for c in CATEGORY_ORDER]
    tool_handles = [plt.Rectangle((0, 0), 1, 1, facecolor=TOOL_COLORS[t],
                                  edgecolor="black", linewidth=0.4, label=t)
                    for t in TOOLS]
    fig.legend(handles=out_handles, title="Outcome", ncol=3,
               loc="upper left", bbox_to_anchor=(0.02, 1.005),
               frameon=False, fontsize=7, title_fontsize=8)
    fig.legend(handles=tool_handles, title="Tool", ncol=4,
               loc="upper right", bbox_to_anchor=(0.98, 1.005),
               frameon=False, fontsize=7, title_fontsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.97], h_pad=0.6, w_pad=0.5)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".svg"), bbox_inches="tight")
    #fig.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
    #fig.savefig(outpath.with_suffix(".eps"), bbox_inches="tight", format="eps")
    plt.close(fig)
    print(f"  wrote {outpath} (+ .svg)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    print("Loading species…")
    data = {}
    for sp, paths in ORGANISMS.items():
        df = load_organism(paths)
        data[sp] = df
        print(f"  {sp}: {len(df)} samples "
              f"({int((df.host_sex=='Male').sum())}M / "
              f"{int((df.host_sex=='Female').sum())}F)")
              # --- save metric tables ---
    metric_frames  = []
    outcome_frames = []
    for sp, df in data.items():
        m = per_class_metrics(df)
        m.insert(0, "species", sp)
        metric_frames.append(m)

        oc = outcome_fractions(df)
        oc_df = (pd.DataFrame(oc).T
                   .rename_axis("method")
                   .reset_index())
        oc_df.insert(0, "species", sp)
        outcome_frames.append(oc_df)

    metrics_all  = pd.concat(metric_frames,  ignore_index=True)
    outcomes_all = pd.concat(outcome_frames, ignore_index=True)

    metrics_all.to_csv(OUTPUT_DIR / "fig4_per_class_metrics.csv", index=False)
    outcomes_all.to_csv(OUTPUT_DIR / "fig4_outcome_fractions.csv", index=False)
    print(f"  wrote {OUTPUT_DIR/'fig4_per_class_metrics.csv'} "
          f"and {OUTPUT_DIR/'fig4_outcome_fractions.csv'}")

    print("\nRendering multi-species figure…")
    make_figure(data, OUTPUT_DIR / "fig4_multispecies_animal.svg")
    print("\nDone.")


if __name__ == "__main__":
    main()
