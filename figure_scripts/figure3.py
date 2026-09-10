"""
SCiMS benchmarking: three-section human figure.

  Panel A: SCiMS classification outcome per cohort
           (stacked horizontal bars: correct / incorrect / uncertain)
  Panel B: Accuracy  heatmap = correct / confident calls
           (cohort x depth bin, one column block per tool; orange-white-blue)
  Panel C: Call rate heatmap = confident calls / total samples
           (cohort x depth bin, one column block per tool; viridis_r)

Cells with no samples are gray.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch


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
# Config
# ---------------------------------------------------------------------------
ORGANISMS = {
    "HMP_oral":           {"scims": "../../human_data/hmp_oral/dbGap_metadata_scims_updated_oral.txt",
                           "bexy":  "../../human_data/hmp_oral/hmp_bexy_output_0.95.txt",
                           "rxry":  "../../human_data/hmp_oral/hmp_rxry_output.txt"},
    "HMP_anterior_nares": {"scims": "../../human_data/hmp_anterior_nares/dbGap_metadata_scims_updated_anterior_nares.txt",
                           "bexy":  "../../human_data/hmp_anterior_nares/hmp_bexy_output_0.95.txt",
                           "rxry":  "../../human_data/hmp_anterior_nares/hmp_rxry_output.txt"},
    "HMP_vaginal":        {"scims": "../../human_data/hmp_vaginal/dbGap_metadata_scims_updated_vaginal.txt",
                           "bexy":  "../../human_data/hmp_vaginal/hmp_bexy_output_0.95.txt",
                           "rxry":  "../../human_data/hmp_vaginal/hmp_rxry_output.txt"},
    "HMP_fecal":          {"scims": "../../human_data/hmp_fecal/dbGap_metadata_scims_updated_fecal.txt",
                           "bexy":  "../../human_data/hmp_fecal/hmp_bexy_output_0.95.txt",
                           "rxry":  "../../human_data/hmp_fecal/hmp_rxry_output.txt"},
    #"Hadza_fecal":        {"scims": "human_data/hadza_fecal/hadza_PRJEB49206_metadata_scims_updated_99.txt",
    #                       "bexy":  "human_data/hadza_fecal/hadza_PRJEB49206_bexy_out.txt",
   #                        "rxry":  "human_data/hadza_fecal/hadza_PRJEB49206_rxry_output.txt"},
    "India_fecal":        {"scims": "../../human_data/india_fecal/indidan_metadata_scims_updated_filt.txt",
                           "bexy":  "../../human_data/india_fecal/indian_metagenomic_PRJNA397112_bexy_out.txt",
                           "rxry":  "../../human_data/india_fecal/indian_metagenomic_PRJNA397112_rxry_output.txt"},
}

BIN_EDGES  = [0, 500, 1000, 10000, np.inf]
BIN_LABELS = ["<500", "500–1k", "1k–10k", ">10k"]

TOOLS   = ["SCiMS", "BeXY", "RX", "RY"]
CLASSES = ["Male", "Female"]

OUTPUT_DIR = Path("./figures_human")

BLUE   = "#0099CC"
ORANGE = "#FF6633"
CMAP_ACC  = LinearSegmentedColormap.from_list("acc_div", [ORANGE, "#ffffff", BLUE])
CMAP_RATE = plt.get_cmap("viridis_r")

# Panel A outcome colors
C_CORRECT   = "#0099cc"
C_INCORRECT = "#f08c2f"
C_UNCERTAIN = "#c3c3c3"


# ---------------------------------------------------------------------------
# Per-tool call extraction
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

def call_bexy(karyotype):
    if karyotype == "XX": return "Female"
    if karyotype == "XY": return "Male"
    return "uncertain"

def call_rx(ci_str):
    parsed = _parse_ci(ci_str)
    if parsed is None: return "uncertain"
    low, high = parsed
    if low  > 0.8: return "Female"
    if high < 0.6: return "Male"
    return "uncertain"

def call_ry(ci_str):
    parsed = _parse_ci(ci_str)
    if parsed is None: return "uncertain"
    low, high = parsed
    if low  > 0.077: return "Male"
    if high < 0.016: return "Female"
    return "uncertain"


# ---------------------------------------------------------------------------
# Load + merge per organism
# ---------------------------------------------------------------------------
def load_organism(paths):
    scims = pd.read_csv(paths["scims"], sep="\t").rename(columns={"Run": "Sample"})
    scims["host_sex"] = scims["host_sex"].map(_norm_sex)
    scims["SCiMS"]    = scims["SCiMS_sex"].map(_norm_sex).fillna("uncertain")
    scims["depth"]    = pd.to_numeric(scims["SCiMS_reads_mapped"], errors="coerce")
    base = scims[["Sample", "host_sex", "depth", "SCiMS"]].copy()

    bexy = pd.read_csv(paths["bexy"], sep="\t").rename(columns={"sample": "Sample"})
    bexy["Sample"] = bexy["Sample"].astype(str).str.replace(".sorted", "", regex=False)
    bexy["BeXY"]   = bexy["sex_karyotype"].map(call_bexy)

    rxry = pd.read_csv(paths["rxry"], sep="\t").rename(columns={"Run": "Sample"})
    rxry["RX"] = rxry["Rx 95% CI"].map(call_rx)
    rxry["RY"] = rxry["Ry 95% CI"].map(call_ry)

    n_scims = len(base)
    # INNER join: keep only samples ALL tools were run on, so no tool is
    # penalized as 'uncertain' for samples it never processed.
    df = base.merge(bexy[["Sample", "BeXY"]], on="Sample", how="inner")
    df = df.merge(rxry[["Sample", "RX", "RY"]], on="Sample", how="inner")
    for t in ["BeXY", "RX", "RY"]:
        df[t] = df[t].fillna("uncertain")   # genuine parsing gaps only
    df = df.dropna(subset=["host_sex", "depth"])
    if len(df) < n_scims:
        print(f"    [note] {n_scims} SCiMS samples -> {len(df)} common to all "
              f"tools ({n_scims - len(df)} dropped: not run by all tools)")
    return df


# ---------------------------------------------------------------------------
# Panel A: outcome fractions per cohort, for ANY tool
# ---------------------------------------------------------------------------
def tool_outcomes(organisms_data, orgs, tool):
    """Correct / incorrect / uncertain counts+fractions per cohort for `tool`."""
    rows = []
    for org in orgs:
        df = organisms_data[org]
        total = len(df)
        unc = int((df[tool] == "uncertain").sum())
        cor = int(((df[tool] == df["host_sex"]) & (df[tool] != "uncertain")).sum())
        inc = total - cor - unc
        rows.append({"cohort": org, "correct": cor, "incorrect": inc,
                     "uncertain": unc, "total": total})
    return pd.DataFrame(rows).set_index("cohort").reindex(orgs)


# ---------------------------------------------------------------------------
# Panels B/C: accuracy + call rate per (organism, depth_bin, tool)
# ---------------------------------------------------------------------------
def build_matrices(organisms_data):
    orgs   = list(organisms_data.keys())
    n_org  = len(orgs)
    n_bin  = len(BIN_LABELS)

    acc      = {t: np.full((n_org, n_bin), np.nan) for t in TOOLS}
    rate     = {t: np.full((n_org, n_bin), np.nan) for t in TOOLS}
    n_called = {t: np.zeros((n_org, n_bin), dtype=int) for t in TOOLS}
    n_total  = {t: np.zeros((n_org, n_bin), dtype=int) for t in TOOLS}

    for i, org in enumerate(orgs):
        df = organisms_data[org].copy()
        df["bin"] = pd.cut(df["depth"], bins=BIN_EDGES, labels=BIN_LABELS,
                           right=False, include_lowest=True)
        for j, b in enumerate(BIN_LABELS):
            sub = df[df["bin"] == b]
            if len(sub) == 0:
                continue
            for t in TOOLS:
                total   = len(sub)
                called  = (sub[t] != "uncertain").sum()
                correct = ((sub[t] == sub["host_sex"]) & (sub[t] != "uncertain")).sum()
                n_total[t][i, j]  = total
                n_called[t][i, j] = called
                if called > 0:
                    acc[t][i, j] = correct / called
                rate[t][i, j] = called / total
    return orgs, acc, rate, n_called, n_total


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
def _text_color_for(val, norm, cmap):
    if np.isnan(val): return "#555"
    rgba = cmap(norm(val))
    r, g, b = rgba[:3]
    lum = 0.299 * r + 0.587 * g + 0.114 * b
    return "white" if lum < 0.55 else "#222"


def _plot_facet(ax, matrix, n_anno, organisms, title,
                cmap, norm, show_y, show_x, show_xlabel):
    cmap_with_bad = cmap.copy()
    cmap_with_bad.set_bad(color="#d9d9d9")
    masked = np.ma.masked_invalid(matrix)
    ax.imshow(masked, aspect="auto", cmap=cmap_with_bad, norm=norm,
              interpolation="nearest")

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]; n = n_anno[i, j]
            if np.isnan(val):
                txt = f"n={n}" if n > 0 else ""; color = "#555"
            else:
                txt = f"{val:.2f}\n(n={n})"; color = _text_color_for(val, norm, cmap)
            ax.text(j, i, txt, ha="center", va="center",
                    color=color, linespacing=1.0, fontsize=6)

    ax.set_xticks(range(len(BIN_LABELS)))
    ax.set_xticklabels(BIN_LABELS if show_x else [], rotation=0)
    if show_xlabel:
        ax.set_xlabel("Host reads")
    ax.set_yticks(range(len(organisms)))
    ax.set_yticklabels(organisms if show_y else [])
    if title is not None:
        ax.set_title(title, pad=6, fontweight="medium")
    ax.set_xticks(np.arange(-.5, len(BIN_LABELS), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(organisms), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)


def _plot_panel_a(ax, outcomes, orgs, title, show_y):
    """Stacked horizontal bars of one tool's outcome fractions per cohort."""
    bar_h = 0.62
    for i, org in enumerate(orgs):
        row = outcomes.loc[org]
        total = row["total"]
        cor = row["correct"] / total if total else 0
        inc = row["incorrect"] / total if total else 0
        unc = row["uncertain"] / total if total else 0
        ax.barh(i, cor, color=C_CORRECT, edgecolor="black", height=bar_h, linewidth=0.4)
        ax.barh(i, inc, left=cor, color=C_INCORRECT, edgecolor="black", height=bar_h, linewidth=0.4)
        ax.barh(i, unc, left=cor + inc, color=C_UNCERTAIN, edgecolor="black", height=bar_h, linewidth=0.4)
        # only annotate the larger segments to avoid clutter across 4 panels
        if cor > 0.10:
            ax.text(cor / 2, i, f"{cor*100:.0f}", ha="center", va="center", fontsize=5)
        if inc > 0.10:
            ax.text(cor + inc / 2, i, f"{inc*100:.0f}", ha="center", va="center", fontsize=5)
        if unc > 0.10:
            ax.text(cor + inc + unc / 2, i, f"{unc*100:.0f}", ha="center", va="center", fontsize=5)

    ax.set_yticks(range(len(orgs)))
    ax.set_yticklabels(orgs if show_y else [])
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["0", "1"])
    ax.set_xlabel("Fraction")
    ax.set_title(title, pad=6, fontweight="medium")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linestyle="--", alpha=0.4)


# ---------------------------------------------------------------------------
# Combined 3-section figure
# ---------------------------------------------------------------------------
def make_combined_figure(orgs, outcomes_by_tool, acc, rate, n_called, n_total, outpath):
    n_org  = len(orgs)
    n_tool = len(TOOLS)

    fig_w_in = 180 / 25.4
    row_h_in = 0.40 * n_org + 0.50
    # height: Panel A (slightly shorter) + two heatmap row blocks; cap for one page
    fig_h_in = min(row_h_in * 3 + 0.9, 9.4)

    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=300)
    # 3 row-blocks: A is a per-tool outcome grid; B and C are heatmap grids
    gs = GridSpec(
        3, n_tool, figure=fig,
        height_ratios=[0.85, 1.0, 1.0],
        hspace=0.38, wspace=0.08,
        left=0.13, right=0.85, top=0.94, bottom=0.07,
    )

    norm_acc  = TwoSlopeNorm(vmin=0.0, vcenter=0.5, vmax=1.0)
    norm_rate = Normalize(vmin=0.0, vmax=1.0)

    # ---- Panel A: one outcome sub-panel per tool ----
    panel_a_axes = []
    for k, tool in enumerate(TOOLS):
        ax = fig.add_subplot(gs[0, k])
        _plot_panel_a(ax, outcomes_by_tool[tool], orgs,
                      title=tool, show_y=(k == 0))
        panel_a_axes.append(ax)

    # one shared outcome legend, centered just below the Panel A row
    pa_left  = panel_a_axes[0].get_position(fig)
    pa_right = panel_a_axes[-1].get_position(fig)
    legend_handles = [Patch(facecolor=C_CORRECT, edgecolor="black", label="Correct"),
                      Patch(facecolor=C_INCORRECT, edgecolor="black", label="Incorrect"),
                      Patch(facecolor=C_UNCERTAIN, edgecolor="black", label="Uncertain")]
    fig.legend(handles=legend_handles, ncol=3, frameon=False, fontsize=7,
               loc="center",
               bbox_to_anchor=(((pa_left.x0 + pa_right.x1) / 2), pa_left.y0 - 0.025))

    # ---- Panel B: accuracy heatmaps ----
    for k, tool in enumerate(TOOLS):
        ax = fig.add_subplot(gs[1, k])
        _plot_facet(ax, acc[tool], n_called[tool], orgs,
                    title=tool, cmap=CMAP_ACC, norm=norm_acc,
                    show_y=(k == 0), show_x=False, show_xlabel=False)

    # ---- Panel C: call-rate heatmaps ----
    for k, tool in enumerate(TOOLS):
        ax = fig.add_subplot(gs[2, k])
        _plot_facet(ax, rate[tool], n_total[tool], orgs,
                    title=None, cmap=CMAP_RATE, norm=norm_rate,
                    show_y=(k == 0), show_x=True, show_xlabel=True)

    # ---- Colorbars (for the two heatmap blocks) ----
    cbar_w = 0.022
    p_b = gs[1, 0].get_position(fig)
    cax_acc = fig.add_axes([0.87, p_b.y0, cbar_w, p_b.y1 - p_b.y0])
    sm_acc = ScalarMappable(norm=norm_acc, cmap=CMAP_ACC); sm_acc.set_array([])
    cb_acc = fig.colorbar(sm_acc, cax=cax_acc, orientation="vertical")
    cb_acc.set_label("Accuracy\n(correct / confident calls)")
    cb_acc.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

    p_c = gs[2, 0].get_position(fig)
    cax_rate = fig.add_axes([0.87, p_c.y0, cbar_w, p_c.y1 - p_c.y0])
    sm_rate = ScalarMappable(norm=norm_rate, cmap=CMAP_RATE); sm_rate.set_array([])
    cb_rate = fig.colorbar(sm_rate, cax=cax_rate, orientation="vertical")
    cb_rate.set_label("Call rate\n(confident calls / total)")
    cb_rate.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

    # Panel letters
    for letter, gs_cell in zip("abc", [gs[0, 0], gs[1, 0], gs[2, 0]]):
        pos = gs_cell.get_position(fig)
        fig.text(0.02, pos.y1 - 0.01, letter, fontsize=13, fontweight="bold",
                 va="top", ha="left")

    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".svg"), bbox_inches="tight")
    #fig.savefig(outpath.with_suffix(".eps"), bbox_inches="tight")
    #fig.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {outpath} (+ .svg)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUTPUT_DIR.mkdir(exist_ok=True)

    print("Loading organisms…")
    data = {}
    for org, paths in ORGANISMS.items():
        df = load_organism(paths)
        data[org] = df
        n_m = int((df["host_sex"] == "Male").sum())
        n_f = int((df["host_sex"] == "Female").sum())
        print(f"  {org}: {len(df)} samples ({n_m} male, {n_f} female)")

    orgs = list(data.keys())
    outcomes_by_tool = {t: tool_outcomes(data, orgs, t) for t in TOOLS}
    print("\nPer-tool outcome counts per cohort:")
    for t in TOOLS:
        print(f"\n[{t}]")
        print(outcomes_by_tool[t].to_string())

    orgs, acc, rate, n_called, n_total = build_matrices(data)

    print("\nRendering 3-section figure…")
    make_combined_figure(
        orgs, outcomes_by_tool, acc, rate, n_called, n_total,
        outpath=OUTPUT_DIR / "outcome_accuracy_callrate_human_all_tools.svg",
    )
    print("\nDone.")


if __name__ == "__main__":
    main()
