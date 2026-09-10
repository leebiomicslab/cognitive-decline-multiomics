"""
Purpose:          CISS group boxplots (C1/C2/C3) and Spearman correlation with MoCA score.
                  Companion visualization script for feature_ratio_modeling_cont_clin.py.
Manuscript:       Results — CISS group comparison and MoCA association
Figure/Table:     Figure 5d (CISS boxplots), Figure 5e (CISS vs MoCA correlation)
Input:            sample_ratio_features.csv (from feature_ratio_modeling_cont_clin.py),
                  vascular.csv (contains moca_1 column)
Output:           ratio_boxplot_groups.png, ratio_moca_correlation.png, ratio_moca_spearman.csv
Main dependencies: pandas, numpy, matplotlib, seaborn, scipy

Status: COMPANION VISUALIZATION — Figure 5d/5e
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from pathlib import Path

# ── Paths (from config/paths.py) ────────────────────────────────────────────────────────
RATIO_CSV = CISS_RATIO_CSV
OUT_DIR   = MODELING_DIR / "classification" / "feature" / "logistic" / "ratio_modeling_cont_clin"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

GROUP_ORDER  = ["C1", "C2", "C3"]
GROUP_COLORS = {"C1": "#ee7a5b", "C2": "#eca362", "C3": "#6fcbbc"}
RATIO_COLS   = ["Micro_Ratio", "Clin_Ratio", "Prot_Ratio"]
RATIO_LABELS = {"Micro_Ratio": "Microbiome CISS", "Clin_Ratio": "Clinical CISS", "Prot_Ratio": "Proteomics CISS"}

# ─── Load ─────────────────────────────────────────────────────────────────────
df_ratio = pd.read_csv(RATIO_CSV).rename(columns={"Unnamed: 0": "SampleID"})
df_ratio["SampleID"] = df_ratio["SampleID"].astype(str)

df_vasc = pd.read_csv(VASCULAR_CSV)
df_vasc = df_vasc.rename(columns={"id": "SampleID"})
df_vasc["SampleID"] = df_vasc["SampleID"].astype(str)
df_vasc = df_vasc[["SampleID", "moca_1"]].dropna()

df = pd.merge(df_ratio, df_vasc, on="SampleID", how="left")
print(f"Samples with MoCA data: {df['moca_1'].notna().sum()} / {len(df)}")
print(f"Group distribution:\n{df['Group'].value_counts()}")

# ─── 1. Boxplot C1 / C2 / C3 ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 6), sharey=False)

for ax, col in zip(axes, RATIO_COLS):
    data_by_grp = [df.loc[df["Group"] == g, col].dropna().values for g in GROUP_ORDER]
    bp = ax.boxplot(data_by_grp, patch_artist=True, notch=False,
                    medianprops=dict(color="black", linewidth=2),
                    whiskerprops=dict(linewidth=1.5),
                    capprops=dict(linewidth=1.5),
                    flierprops=dict(marker="o", markersize=3, alpha=0.4, linestyle="none"))

    for patch, grp in zip(bp["boxes"], GROUP_ORDER):
        patch.set_facecolor(GROUP_COLORS[grp])
        patch.set_alpha(0.8)

    # Overlay raw dots with jitter
    for i, (grp, vals) in enumerate(zip(GROUP_ORDER, data_by_grp)):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                   color=GROUP_COLORS[grp], edgecolor="k", linewidth=0.4,
                   s=18, alpha=0.5, zorder=3)

    # Pairwise Wilcoxon
    pairs = [("C1", "C2", 1, 2), ("C1", "C3", 1, 3), ("C2", "C3", 2, 3)]
    y_top = max(max(v) for v in data_by_grp if len(v) > 0)
    tick_up = (y_top - min(min(v) for v in data_by_grp if len(v) > 0)) * 0.07
    y_sig = y_top + tick_up

    for g1, g2, x1, x2 in pairs:
        v1 = df.loc[df["Group"] == g1, col].dropna().values
        v2 = df.loc[df["Group"] == g2, col].dropna().values
        _, p = stats.mannwhitneyu(v1, v2, alternative="two-sided")
        if p < 0.001:   sig = "***"
        elif p < 0.01:  sig = "**"
        elif p < 0.05:  sig = "*"
        else:           sig = "ns"
        ax.plot([x1, x1, x2, x2], [y_sig, y_sig + tick_up * 0.3, y_sig + tick_up * 0.3, y_sig], color="black", lw=1.2)
        ax.text((x1 + x2) / 2, y_sig + tick_up * 0.35, sig, ha="center", va="bottom", fontsize=10)
        y_sig += tick_up * 1.5

    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(GROUP_ORDER, fontsize=11)
    ax.set_title(RATIO_LABELS[col], fontsize=12, pad=10)
    ax.set_ylabel("CISS Score", fontsize=10)
    ax.axhline(0, color="grey", lw=0.8, ls="--", alpha=0.6)

legend_patches = [mpatches.Patch(facecolor=GROUP_COLORS[g], edgecolor="k", label=g) for g in GROUP_ORDER]
fig.legend(handles=legend_patches, title="Group", loc="upper right", frameon=True, fontsize=10)
fig.suptitle("CISS Feature Distribution by Cognitive Group", fontsize=15, y=1.02, fontweight="bold")
plt.tight_layout()
plt.savefig(OUT_DIR / "ratio_boxplot_groups.png", dpi=600, bbox_inches="tight")
plt.close()
print("Saved: ratio_boxplot_groups.png")

# ─── 2. MoCA Correlation Scatter (3 panels) ───────────────────────────────────
df_corr = df.dropna(subset=["moca_1"])

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

corr_rows = []
for ax, col in zip(axes, RATIO_COLS):
    x = df_corr[col].values
    y = df_corr["moca_1"].values

    # Spearman correlation
    rho, pval = stats.spearmanr(x, y)
    if pval < 0.001:   pstr = "p < 0.001"
    elif pval < 0.01:  pstr = f"p = {pval:.3f}"
    else:              pstr = f"p = {pval:.3f}"
    corr_rows.append({"Ratio": col, "Spearman_rho": rho, "p_value": pval})

    # Color by group
    for grp in GROUP_ORDER:
        mask = df_corr["Group"] == grp
        ax.scatter(x[mask], y[mask], c=GROUP_COLORS[grp], edgecolor="k",
                   linewidth=0.3, s=30, alpha=0.7, label=grp, zorder=3)

    # Trend line (all samples)
    m, b = np.polyfit(x, y, 1)
    x_line = np.linspace(x.min(), x.max(), 200)
    ax.plot(x_line, m * x_line + b, color="black", lw=1.8, ls="--", alpha=0.8)

    ax.set_xlabel(RATIO_LABELS[col], fontsize=10)
    ax.set_ylabel("MoCA Score", fontsize=10)
    ax.set_title(f"{RATIO_LABELS[col]}\nρ = {rho:.3f}, {pstr}", fontsize=11, pad=8)
    ax.grid(alpha=0.25)

# Shared legend
handles = [mpatches.Patch(facecolor=GROUP_COLORS[g], edgecolor="k", label=g) for g in GROUP_ORDER]
fig.legend(handles=handles, title="Group", loc="upper right", fontsize=10)
fig.suptitle("CISS Features vs. MoCA Score (Spearman Correlation)", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(OUT_DIR / "ratio_moca_correlation.png", dpi=600, bbox_inches="tight")
plt.close()
print("Saved: ratio_moca_correlation.png")

# ─── 3. Save correlation table ────────────────────────────────────────────────
corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(OUT_DIR / "ratio_moca_spearman.csv", index=False)
print("\nSpearman Correlation Table:")
print(corr_df.to_string(index=False))
