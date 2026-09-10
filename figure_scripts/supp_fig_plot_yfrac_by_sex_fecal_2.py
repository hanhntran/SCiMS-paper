#!/usr/bin/env python3
"""
Supplementary figure: effect of host depletion on SCiMS sex inference
across three fecal cohorts — HMP (raw, dbGaP), India/PRJNA397112 (raw
deposit), and Hadza/PRJEB49206 (host-depleted deposit).

  Panel (a)  Y/(X+Y) read fraction by labeled host sex, one subpanel per
             cohort, with reference bands from the pooled raw cohorts and
             points colored by SCiMS outcome.
  Panel (b)  Classification accuracy vs. host read depth, cohorts overlaid
             as lines with Jeffreys binomial confidence intervals.

Caption statistics for both panels are printed to stdout.

Usage:
    python supp_fig_plot_yfrac_by_sex_fecal.py \
        --hmp dbGap_metadata_scims_updated_fecal.txt \
        --india indidan_metadata_scims_updated.txt \
        --hadza hadza_PRJEB49206_metadata_scims_updated_filt.txt \
        --out supp_fig_host_depletion

Options:
    --panel-b-cohorts india hadza     # restrict panel (b); default all three
    --no-panel-b                      # panel (a) only
"""

import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------- config
MIN_SEXCHROM_READS = 20     # panel (a): exclude samples with <= this many X+Y reads
JITTER_SEED = 7             # reproducible point jitter
JITTER_WIDTH = 0.16
BAND_Q = (0.05, 0.95)       # reference band quantiles, same for both sexes
BAND_CONCORDANT_ONLY = True  # build bands from samples whose SCiMS call
                             # matched the recorded label, so the reference
                             # reflects verified karyotype, not metadata error
INFLATION_THRESH = {"male": 0.20, "female": 0.10}  # caption statistics

# panel (b)
BINS = [0, 500, 1000, 10000, np.inf]
BIN_LABELS = ["<500", "500–1k", "1k–10k", ">10k"]
CI_LEVEL = 0.95
MIN_N = 3                   # don't draw a point estimated from fewer than this
DODGE = 0.055               # horizontal offset so error bars don't overlap

COLORS = {
    "correct": "#2b8cbe",
    "wrong": "#e6550d",
    "uncertain": "#999999",
    "male_band": "#2b8cbe",
    "female_band": "#de77ae",
}

# keys must match the cohort keys built in main()
COHORT_STYLE = {
    "hmp":   dict(color="#1b7837", marker="o", ls="-",  label="HMP fecal (raw)"),
    "india": dict(color="#2b8cbe", marker="s", ls="-",  label="India fecal (raw)"),
    "hadza": dict(color="#e6550d", marker="^", ls="--", label="Hadza fecal (depleted)"),
}

PANEL_A_TITLES = {
    "hmp":   "HMP fecal\n(raw, dbGaP)",
    "india": "India fecal\n(raw deposit)",
    "hadza": "Hadza fecal\n(host-depleted deposit)",
}

PLOT_STYLE = {
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",   # editable text in SVG
    "pdf.fonttype": 42,       # editable text in PDF (TrueType)
}


# ---------------------------------------------------------------- loading
def load_cohort(path, label_col):
    """Read a SCiMS results table; return per-sample dataframe.

    Keeps all confidently-labeled samples (used by panel b) and adds the
    Y-fraction plus a sex-chromosome-read filter flag (used by panel a).
    """
    df = pd.read_csv(path, sep="\t")
    d = df.rename(
        columns={
            label_col: "host_sex",
            "SCiMS_sex": "call",
            "SCiMS_reads_mapped": "host",
            "SCiMS_reads_mapped_to_X": "X",
            "SCiMS_reads_mapped_to_Y": "Y",
        }
    )[["Run", "host_sex", "call", "host", "X", "Y"]]
    d["host_sex"] = d["host_sex"].astype(str).str.lower().str.strip()
    d["call"] = d["call"].astype(str).str.lower().str.strip()
    d = d[d["host_sex"].isin(["male", "female"])]
    d = d[d["call"].isin(["male", "female", "uncertain"])].copy()
    d["sexchrom_ok"] = (d["X"] + d["Y"]) > MIN_SEXCHROM_READS
    with np.errstate(invalid="ignore", divide="ignore"):
        d["yfrac"] = d["Y"] / (d["X"] + d["Y"])
    d["status"] = np.where(
        d["call"] == "uncertain",
        "uncertain",
        np.where(d["call"] == d["host_sex"], "correct", "wrong"),
    )
    d["bin"] = pd.cut(d["host"], BINS, labels=BIN_LABELS, right=False)
    return d


def yfrac_subset(d):
    """Panel (a) subset: samples with enough sex-chromosome reads."""
    return d[d["sexchrom_ok"]]


def accuracy_table(d):
    """Panel (b): per-bin accuracy among confident calls, with Jeffreys CIs."""
    rows = []
    a = (1 - CI_LEVEL) / 2
    for b in BIN_LABELS:
        sub = d[(d["bin"] == b) & d["call"].isin(["male", "female"])]
        n = len(sub)
        k = int((sub["call"] == sub["host_sex"]).sum())
        if n == 0:
            rows.append((b, 0, 0, np.nan, np.nan, np.nan))
            continue
        lo, hi = stats.beta.ppf([a, 1 - a], k + 0.5, n - k + 0.5)
        rows.append((b, n, k, k / n, lo, hi))
    return pd.DataFrame(rows, columns=["bin", "n", "k", "acc", "lo", "hi"])


# ---------------------------------------------------------------- stats
def caption_stats(cohorts, order, m_band, f_band, panel_b_keys):
    """Print the numbers the figure caption / response letter needs."""
    print("=" * 68)
    print("PANEL (a): sex-chromosome read composition")
    print("=" * 68)
    print(f"Filter: X+Y > {MIN_SEXCHROM_READS} sex-chromosome reads")
    src = "concordant raw samples" if BAND_CONCORDANT_ONLY else "all labeled raw samples"
    print(f"Reference bands ({BAND_Q[0]:.0%}-{BAND_Q[1]:.0%} of {src}): "
          f"male ({m_band[0]:.3f}, {m_band[1]:.3f}), "
          f"female ({f_band[0]:.3f}, {f_band[1]:.3f}); "
          f"gap = {m_band[0] - f_band[1]:+.3f}\n")
    for key in order:
        d = yfrac_subset(cohorts[key])
        name = COHORT_STYLE[key]["label"]
        for sex in ["male", "female"]:
            sub = d[d["host_sex"] == sex]
            yf = sub["yfrac"]
            thr = INFLATION_THRESH[sex]
            print(
                f"{name:<24} {sex:<6} n={len(sub):>3}  "
                f"median={yf.median():.3f}  "
                f"IQR=({yf.quantile(.25):.3f},{yf.quantile(.75):.3f})  "
                f">{thr:.2f}: {(yf > thr).sum()}/{len(sub)} "
                f"({(yf > thr).mean():.0%})  "
                f"misclassified: {(sub['status'] == 'wrong').sum()}"
            )
        print()

    raw_keys = [k for k in order if k != "hadza"]
    raw_m = pd.concat(
        [yfrac_subset(cohorts[k]).query("host_sex == 'male'")["yfrac"]
         for k in raw_keys])
    had_m = yfrac_subset(cohorts["hadza"]).query("host_sex == 'male'")["yfrac"]
    print(f"Hadza males vs pooled raw males (n={len(had_m)} vs {len(raw_m)}):")
    print(f"  Mann-Whitney P = {stats.mannwhitneyu(raw_m, had_m).pvalue:.2e}")
    print(f"  Kolmogorov-Smirnov P = {stats.ks_2samp(raw_m, had_m).pvalue:.2e}")
    inband = int(((had_m >= m_band[0]) & (had_m <= m_band[1])).sum())
    print(f"  Hadza males inside pooled raw male band: {inband}/{len(had_m)}")

    if not panel_b_keys:
        return
    print("\n" + "=" * 68)
    print("PANEL (b): accuracy vs. host read depth")
    print("=" * 68)
    tabs = {}
    for key in panel_b_keys:
        tab = accuracy_table(cohorts[key])
        tabs[key] = tab
        print(f"\n{COHORT_STYLE[key]['label']}")
        print(tab.to_string(index=False,
                            float_format=lambda v: f"{v:.3f}"))
    # Fisher test in the best-powered shared bin
    if {"india", "hadza"} <= set(panel_b_keys):
        b = "1k–10k"
        def kn(key):
            s = cohorts[key]
            s = s[(s["bin"] == b) & s["call"].isin(["male", "female"])]
            return int((s["call"] == s["host_sex"]).sum()), len(s)
        hk, hn = kn("hadza")
        ik, inn = kn("india")
        if hn and inn:
            p = stats.fisher_exact([[hk, hn - hk], [ik, inn - ik]]).pvalue
            print(f"\nFisher ({b}): Hadza {hk}/{hn} vs India {ik}/{inn}, "
                  f"P = {p:.4f}")


# ---------------------------------------------------------------- figure
def make_figure(cohorts, order, m_band, f_band, panel_b_keys, outstem):
    rng = np.random.default_rng(JITTER_SEED)
    n_a = len(order)
    two_panel = bool(panel_b_keys)

    if two_panel:
        fig = plt.figure(figsize=(4.2 * n_a, 8.4))
        gs = fig.add_gridspec(2, n_a, height_ratios=[1, 0.95], hspace=0.42)
        axes_a = [fig.add_subplot(gs[0, i]) for i in range(n_a)]
        ax_b = fig.add_subplot(gs[1, :])
    else:
        fig, axes_a = plt.subplots(1, n_a, figsize=(4.2 * n_a, 4.4),
                                   sharey=True)
        axes_a = list(np.atleast_1d(axes_a))
        ax_b = None

    # ---- panel (a) -------------------------------------------------
    for i, (ax, key) in enumerate(zip(axes_a, order)):
        d = yfrac_subset(cohorts[key])
        ax.axhspan(*m_band, color=COLORS["male_band"], alpha=0.10, lw=0)
        ax.axhspan(*f_band, color=COLORS["female_band"], alpha=0.12, lw=0)
        for j, sex in enumerate(["female", "male"]):
            sub = d[d["host_sex"] == sex]
            x = j + rng.uniform(-JITTER_WIDTH, JITTER_WIDTH, len(sub))
            for st in ["correct", "uncertain", "wrong"]:  # wrong drawn on top
                m = (sub["status"] == st).values
                ax.scatter(
                    x[m], sub["yfrac"].values[m],
                    s=20, c=COLORS[st],
                    alpha=0.95 if st == "wrong" else 0.75,
                    edgecolors="white", linewidths=0.3,
                    zorder=3 if st == "wrong" else 2,
                )
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Labeled\nfemale", "Labeled\nmale"])
        ax.set_title(PANEL_A_TITLES[key], fontsize=11)
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylim(-0.02, 0.75)
        if i:
            ax.tick_params(labelleft=False)   # safe with sharey
    axes_a[0].set_ylabel("Y / (X + Y) read fraction")

    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=COLORS[k], markersize=7, label=lbl)
        for k, lbl in [("correct", "Correct"),
                       ("wrong", "Misclassified"),
                       ("uncertain", "Uncertain")]
    ]
    handles += [
        plt.Rectangle((0, 0), 1, 1, fc=COLORS["male_band"], alpha=0.10,
                      label=f"Male range (raw, {BAND_Q[0]:.0%}-{BAND_Q[1]:.0%})"),
        plt.Rectangle((0, 0), 1, 1, fc=COLORS["female_band"], alpha=0.12,
                      label=f"Female range (raw, {BAND_Q[0]:.0%}-{BAND_Q[1]:.0%})"),
    ]
    axes_a[-1].legend(handles=handles, frameon=False, fontsize=8.5,
                      loc="upper left")

    # ---- panel (b) -------------------------------------------------
    if two_panel:
        xpos = np.arange(len(BIN_LABELS))
        offsets = np.linspace(-DODGE, DODGE, len(panel_b_keys))
        for key, off in zip(panel_b_keys, offsets):
            tab = accuracy_table(cohorts[key])
            style = COHORT_STYLE[key]
            ok = (tab["n"] >= MIN_N).values
            x = xpos[ok] + off
            y = tab.loc[ok, "acc"].values
            lo = np.clip(y - tab.loc[ok, "lo"].values, 0, None)
            hi = np.clip(tab.loc[ok, "hi"].values - y, 0, None)
            ax_b.errorbar(
                x, y, yerr=[lo, hi],
                color=style["color"], marker=style["marker"], ls=style["ls"],
                lw=1.8, ms=6, capsize=3, elinewidth=1.2,
                label=style["label"], zorder=3,
            )
            for xi, yi, ni in zip(x, y, tab.loc[ok, "n"].values):
                ax_b.annotate(f"n={ni}", (xi, yi), textcoords="offset points",
                              xytext=(0, -14), ha="center", fontsize=7,
                              color=style["color"])
        ax_b.axhline(0.5, color="#999999", lw=0.8, ls=":", zorder=1)
        ax_b.annotate("chance", (len(BIN_LABELS) - 0.55, 0.505), fontsize=8,
                      color="#777777", va="bottom", ha="right")
        ax_b.set_xticks(xpos)
        ax_b.set_xticklabels(BIN_LABELS)
        ax_b.set_xlabel("Host read count")
        ax_b.set_ylabel("Accuracy (correct / confident calls)")
        ax_b.set_ylim(0.3, 1.05)
        ax_b.set_xlim(-0.45, len(BIN_LABELS) - 0.55)
        ax_b.legend(frameon=False, fontsize=9, loc="lower right")

    # ---- panel letters --------------------------------------------
    axes_a[0].text(-0.18, 1.10, "a", transform=axes_a[0].transAxes,
                   fontsize=15, fontweight="bold", va="top", ha="right")
    if two_panel:
        ax_b.text(-0.06, 1.08, "b", transform=ax_b.transAxes,
                  fontsize=15, fontweight="bold", va="top", ha="right")

    if not two_panel:
        fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{outstem}.{ext}", dpi=300, bbox_inches="tight")
        print(f"\nwrote {outstem}.{ext}")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hmp", required=True, help="HMP fecal SCiMS table (tsv)")
    ap.add_argument("--india", required=True, help="India cohort SCiMS table (tsv)")
    ap.add_argument("--hadza", required=True, help="Hadza cohort SCiMS table (tsv)")
    ap.add_argument("--out", default="supp_fig_host_depletion", help="output stem")
    ap.add_argument("--panel-b-cohorts", nargs="+", default=["hmp", "india", "hadza"],
                    choices=["hmp", "india", "hadza"],
                    help="which cohorts to draw in panel (b)")
    ap.add_argument("--no-panel-b", action="store_true",
                    help="draw panel (a) only")
    args = ap.parse_args()

    mpl.rcParams.update(PLOT_STYLE)

    cohorts = {
        "hmp": load_cohort(args.hmp, "host_sex"),
        "india": load_cohort(args.india, "Sex"),
        "hadza": load_cohort(args.hadza, "sex"),
    }
    order = ["hmp", "india", "hadza"]          # raw cohorts first, depleted last
    panel_b_keys = [] if args.no_panel_b else [
        k for k in order if k in args.panel_b_cohorts]

    # reference bands from the pooled raw cohorts
    raw = pd.concat([yfrac_subset(cohorts[k]) for k in order if k != "hadza"])
    if BAND_CONCORDANT_ONLY:
        raw = raw[raw["status"] == "correct"]
    m_band = tuple(raw.loc[raw["host_sex"] == "male", "yfrac"].quantile(BAND_Q))
    f_band = tuple(raw.loc[raw["host_sex"] == "female", "yfrac"].quantile(BAND_Q))
    print(f"Reference n: {(raw['host_sex'] == 'male').sum()} males, "
          f"{(raw['host_sex'] == 'female').sum()} females"
          f"{' (concordant only)' if BAND_CONCORDANT_ONLY else ''}\n")

    caption_stats(cohorts, order, m_band, f_band, panel_b_keys)
    make_figure(cohorts, order, m_band, f_band, panel_b_keys, args.out)


if __name__ == "__main__":
    main()
