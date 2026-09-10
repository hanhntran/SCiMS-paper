"""
Human metagenomic split-panel figure, extended to six cohorts.

Panel A (left):  SCiMS classification outcome (correct / incorrect / uncertain)
                 per cohort  (stacked horizontal bars)
Panel A (right): host read-depth distribution per cohort (stacked bars)
Panel B (bottom):per-sex precision / recall / F1 for SCiMS, BeXY, Rx, Ry,
                 pooled across all six cohorts

Cohorts: HMP_oral, HMP_anterior_nares, HMP_fecal, HMP_vaginal,
         Hadza_fecal, India_fecal.

Recall convention (matches all other scripts in the project):
    recall = TP / (TP + FN + uncertain)        # uncertain hurts recall
    precision = TP / (TP + FP)
"""

import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns


# ---------------------------------------------------------------------------
# Config: one entry per cohort, each with its own three files
# ---------------------------------------------------------------------------
ORGANISMS = {
    "HMP_oral":           {"scims": "human_data/hmp_oral/dbGap_metadata_scims_updated_oral.txt",
                           "bexy":  "human_data/hmp_oral/hmp_bexy_output_0.95.txt",
                           "rxry":  "human_data/hmp_oral/hmp_rxry_output.txt"},
    "HMP_anterior_nares": {"scims": "human_data/hmp_anterior_nares/dbGap_metadata_scims_updated_anterior_nares.txt",
                           "bexy":  "human_data/hmp_anterior_nares/hmp_bexy_output_0.95.txt",
                           "rxry":  "human_data/hmp_anterior_nares/hmp_rxry_output.txt"},
    "HMP_fecal":          {"scims": "human_data/hmp_fecal/dbGap_metadata_scims_updated_fecal.txt",
                           "bexy":  "human_data/hmp_fecal/hmp_bexy_output_0.95.txt",
                           "rxry":  "human_data/hmp_fecal/hmp_rxry_output.txt"},
    "HMP_vaginal":        {"scims": "human_data/hmp_vaginal/dbGap_metadata_scims_updated_vaginal.txt",
                           "bexy":  "human_data/hmp_vaginal/hmp_bexy_output_0.95.txt",
                           "rxry":  "human_data/hmp_vaginal/hmp_rxry_output.txt"},
    "Hadza_fecal":        {"scims": "human_data/hadza_fecal/hadza_PRJEB49206_metadata_scims_updated_99.txt",
                           "bexy":  "human_data/hadza_fecal/hadza_PRJEB49206_bexy_out.txt",
                           "rxry":  "human_data/hadza_fecal/hadza_PRJEB49206_rxry_output.txt"},
    "India_fecal":        {"scims": "human_data/india_fecal/indidan_metadata_scims_updated_filt.txt",
                           "bexy":  "human_data/india_fecal/indian_metagenomic_PRJNA397112_bexy_out.txt",
                           "rxry":  "human_data/india_fecal/indian_metagenomic_PRJNA397112_rxry_output.txt"},
}

# Display order (top -> bottom), roughly high host DNA -> low host DNA
COHORT_ORDER = ["HMP_oral", "HMP_anterior_nares", "HMP_vaginal",
                "HMP_fecal", "Hadza_fecal", "India_fecal"]

READ_DEPTH_COL = "SCiMS_reads_mapped"
DEPTH_BINS   = [0, 1000, 10000, np.inf]
DEPTH_LABELS = ["Low (<1k)", "Med (1-10k)", "High (>10k)"]

CLASSES = ["Male", "Female"]
METHODS = ["SCiMS", "BeXY", "Rx", "Ry"]

if not os.path.isdir("./figures"):
    os.mkdir("./figures")


# ---------------------------------------------------------------------------
# Call extraction (consistent with the other scripts)
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

def call_bexy(k):
    if k == "XX": return "Female"
    if k == "XY": return "Male"
    return "uncertain"


def load_organism(paths, cohort):
    scims = pd.read_csv(paths["scims"], sep="\t").rename(columns={"Run": "Sample"})
    scims["host_sex"]  = scims["host_sex"].map(_norm_sex)
    scims["SCiMS"]     = scims["SCiMS_sex"].map(_norm_sex).fillna("uncertain")
    scims["depth"]     = pd.to_numeric(scims[READ_DEPTH_COL], errors="coerce")
    base = scims[["Sample", "host_sex", "depth", "SCiMS"]].copy()

    bexy = pd.read_csv(paths["bexy"], sep="\t").rename(columns={"sample": "Sample"})
    bexy["Sample"] = bexy["Sample"].astype(str).str.replace(".sorted", "", regex=False)
    bexy["BeXY"]   = bexy["sex_karyotype"].map(call_bexy)

    rxry = pd.read_csv(paths["rxry"], sep="\t").rename(columns={"Run": "Sample"})
    rxry["Rx"] = rxry["Rx 95% CI"].map(call_rx)
    rxry["Ry"] = rxry["Ry 95% CI"].map(call_ry)

    df = base.merge(bexy[["Sample", "BeXY"]], on="Sample", how="left")
    df = df.merge(rxry[["Sample", "Rx", "Ry"]], on="Sample", how="left")
    for t in ["BeXY", "Rx", "Ry"]:
        df[t] = df[t].fillna("uncertain")
    df["cohort"] = cohort
    return df.dropna(subset=["host_sex", "depth"])


# ---------------------------------------------------------------------------
# Load all cohorts
# ---------------------------------------------------------------------------
frames = []
for cohort in COHORT_ORDER:
    d = load_organism(ORGANISMS[cohort], cohort)
    frames.append(d)
    print(f"  {cohort}: {len(d)} samples")
merged_df = pd.concat(frames, ignore_index=True)

merged_df["depth_group"] = pd.cut(merged_df["depth"], bins=DEPTH_BINS,
                                  labels=DEPTH_LABELS)


# ---------------------------------------------------------------------------
# Panel A data: SCiMS outcome counts per cohort + depth distribution
# ---------------------------------------------------------------------------
def scims_outcomes(df):
    rows = []
    for cohort in COHORT_ORDER:
        grp = df[df["cohort"] == cohort]
        total = len(grp)
        unc = (grp["SCiMS"] == "uncertain").sum()
        cor = ((grp["SCiMS"] == grp["host_sex"]) & (grp["SCiMS"] != "uncertain")).sum()
        inc = total - cor - unc
        rows.append({"cohort": cohort, "correct": cor,
                     "incorrect": inc, "uncertain": unc, "total": total})
    return pd.DataFrame(rows).set_index("cohort").reindex(COHORT_ORDER)

acc_scims = scims_outcomes(merged_df)

depth_stats = (merged_df.groupby(["cohort", "depth_group"], observed=False)
               .size().unstack(fill_value=0))
depth_pct = depth_stats.div(depth_stats.sum(axis=1), axis=0).reindex(COHORT_ORDER)


# ---------------------------------------------------------------------------
# Panel B data: per-sex precision/recall/f1 pooled across cohorts
# ---------------------------------------------------------------------------
def per_class_metrics(df):
    rows = []
    actual = df["host_sex"]
    for method in METHODS:
        preds = df[method]
        for cls in CLASSES:
            is_cls = actual == cls
            n_correct   = (is_cls & (preds == cls)).sum()
            n_incorrect = (is_cls & preds.isin(CLASSES) & (preds != cls)).sum()
            n_uncertain = (is_cls & ~preds.isin(CLASSES)).sum()
            n_fp        = (~is_cls & actual.isin(CLASSES) & (preds == cls)).sum()

            prec_den = n_correct + n_fp
            rec_den  = n_correct + n_incorrect + n_uncertain     # FULL denominator
            precision = n_correct / prec_den if prec_den > 0 else np.nan
            recall    = n_correct / rec_den  if rec_den  > 0 else np.nan
            if np.isnan(precision) or np.isnan(recall) or (precision + recall) == 0:
                f1 = np.nan
            else:
                f1 = 2 * precision * recall / (precision + recall)
            rows.append({"method": method, "class": cls,
                         "precision": precision, "recall": recall,
                         "f1-score": f1})
    return pd.DataFrame(rows)

metrics = per_class_metrics(merged_df)

# Save tables
merged_df.to_csv("./figures/merged_df_human_allcohorts.csv", index=False)
depth_pct.to_csv("./figures/depth_pct_human_allcohorts.csv")
metrics.to_csv("./figures/metrics_human_allcohorts.csv", index=False)
acc_scims.to_csv("./figures/acc_scims_human_allcohorts.csv")


# ---------------------------------------------------------------------------
# PLOTTING
# ---------------------------------------------------------------------------
c_perf = {"Accuracy": "#39c0c8", "Misclassification": "#f08c2f",
          "Uncertainty": "#c3c3c3"}
c_read = {"High (>10k)": "#245668", "Med (1-10k)": "#0d8f81",
          "Low (<1k)": "#6ec574"}
line_colors = {"SCiMS": "#58508f", "BeXY": "#bd5090",
               "Rx": "#ac9546", "Ry": "#eddca5"}

labels = COHORT_ORDER
y_pos  = range(len(labels))

fig = plt.figure(figsize=(12, 6.5))
gs = gridspec.GridSpec(2, 6, height_ratios=[1.1, 1.2], wspace=1.0, hspace=0.5)

bar_height = 0.6

# ---- Panel A left: classification outcome ----
ax_perf = plt.subplot(gs[0, :3])
for i, cohort in enumerate(labels):
    row = acc_scims.loc[cohort]
    total = row["total"]
    acc_rate = row["correct"]   / total if total else 0
    inc_rate = row["incorrect"] / total if total else 0
    unc_rate = row["uncertain"] / total if total else 0

    ax_perf.barh(i, acc_rate, color=c_perf["Accuracy"], edgecolor="black", height=bar_height)
    ax_perf.barh(i, inc_rate, left=acc_rate, color=c_perf["Misclassification"], edgecolor="black", height=bar_height)
    ax_perf.barh(i, unc_rate, left=acc_rate + inc_rate, color=c_perf["Uncertainty"], edgecolor="black", height=bar_height)

    ax_perf.text(acc_rate / 2, i, f"{acc_rate*100:.1f}%", ha="center", va="center", fontsize=8)
    if inc_rate > 0.05:
        ax_perf.text(acc_rate + inc_rate / 2, i, f"{inc_rate*100:.1f}%", ha="center", va="center", fontsize=8)
    if unc_rate > 0.05:
        ax_perf.text(acc_rate + inc_rate + unc_rate / 2, i, f"{unc_rate*100:.1f}%", ha="center", va="center", fontsize=8)

ax_perf.set_yticks(y_pos)
ax_perf.set_yticklabels(labels, fontsize=10, fontweight="bold")
ax_perf.set_xlim(0, 1.05)
ax_perf.set_title("SCiMS classification outcome", fontsize=12, pad=10)
ax_perf.invert_yaxis()
ax_perf.spines[["top", "right", "bottom"]].set_visible(False)
ax_perf.grid(axis="x", linestyle="--", alpha=0.5)
# legend for outcome colors
from matplotlib.patches import Patch
outcome_handles = [Patch(facecolor=c_perf[k], edgecolor="black", label=lbl)
                   for k, lbl in [("Accuracy", "Correct"),
                                  ("Misclassification", "Incorrect"),
                                  ("Uncertainty", "Uncertain")]]
ax_perf.legend(handles=outcome_handles, loc="upper center",
               bbox_to_anchor=(0.5, -0.08), ncol=3, frameon=False, fontsize=9)

# ---- Panel A right: read-depth distribution ----
ax_read = plt.subplot(gs[0, 3:], sharey=ax_perf)
for i, cohort in enumerate(labels):
    vals = depth_pct.loc[cohort]
    v_high = vals.get("High (>10k)", 0)
    v_med  = vals.get("Med (1-10k)", 0)
    v_low  = vals.get("Low (<1k)", 0)
    ax_read.barh(i, v_high, color=c_read["High (>10k)"], edgecolor="black",
                 height=bar_height, label="High (>10k)" if i == 0 else "")
    ax_read.barh(i, v_med, left=v_high, color=c_read["Med (1-10k)"], edgecolor="black",
                 height=bar_height, label="Med (1-10k)" if i == 0 else "")
    ax_read.barh(i, v_low, left=v_high + v_med, color=c_read["Low (<1k)"], edgecolor="black",
                 height=bar_height, label="Low (<1k)" if i == 0 else "")

ax_read.set_xlim(0, 1.05)
ax_read.set_title("Host read depth distribution", fontsize=12, pad=10)
ax_read.spines[["top", "right", "left", "bottom"]].set_visible(False)
ax_read.tick_params(left=False, labelleft=False, bottom=True)
ax_read.grid(axis="x", linestyle="--", alpha=0.5)
ax_read.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=3,
               frameon=False, fontsize=9)

# ---- Panel B: per-sex metrics ----
long = metrics.melt(id_vars=["method", "class"],
                    value_vars=["precision", "recall", "f1-score"],
                    var_name="metric", value_name="score")
plot_metrics = ["precision", "recall", "f1-score"]
for i, met in enumerate(plot_metrics):
    ax = plt.subplot(gs[1, (i*2):(i*2)+2])
    sns.barplot(data=long[long["metric"] == met], x="class", y="score",
                hue="method", palette=line_colors, edgecolor="black",
                linewidth=0.6, ax=ax)
    ax.set_title(met.capitalize().replace("-score", " score"), fontsize=11, weight="bold")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.5)
    ax.set_xlabel(""); ax.set_ylabel("")
    if i != 1:
        if ax.get_legend(): ax.get_legend().remove()
    else:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=4,
                  frameon=False, title="")

plt.tight_layout()
plt.savefig("./figures/fig3_split_panel_allcohorts.eps", dpi=300,
            bbox_inches="tight", transparent=True, format="eps")
plt.savefig("./figures/fig3_split_panel_allcohorts.png", dpi=300,
            bbox_inches="tight")
print("wrote ./figures/fig3_split_panel_allcohorts.{eps,png}")
plt.show()
