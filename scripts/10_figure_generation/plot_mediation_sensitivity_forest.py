import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Forest plot of mediation sensitivity across metabolites
Manuscript:       Results - Mediation sensitivity
Figure/Table:     Supplementary Fig. 9
"""


"""
Forest plot for serial mediation sensitivity analysis (LOG2 strategy).
Shows TMAO_pre_2 and TMAO_pre_avg across 4 covariate models ? 2 directions.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# ?? Paths ?????????????????????????????????????????????????????????????????????
LOG2_CSV = str(MODELING_DIR / "sequential_mediation_all_metabolites_log2.csv")
OUT_DIR  = str(FIGURES_DIR / "mediation")
os.makedirs(OUT_DIR, exist_ok=True)

# ?? Data ??????????????????????????????????????????????????????????????????????
df = pd.read_csv(LOG2_CSV)
targets = ["TMAO_pre_avg"]
df = df[df["Metabolite"].isin(targets)].copy()

# Labels
MODEL_LABELS = {
    "Model_1_Core":              "Model 1\nCore",
    "Model_2_Sens_Inflam_Renal": "Model 2\n+CRP, eGFR",
    "Model_3_Sens_Micro_Drug":   "Model 3\n+Statin, PPI, Abx",
    "Model_4_Sens_Chronic_Med":  "Model 4\n+dm_med, htn_med,\nchol_med",
}
MODEL_ORDER = list(MODEL_LABELS.keys())
N_MODELS = len(MODEL_ORDER)

DIR_LABELS = {
    "Forward (Micro->Met->Prot->MoCA)": "Forward\nMicro?MAO?rot?oCA",
    "Reverse (Prot->Met->Micro->MoCA)": "Reverse\nProt?MAO?icro?oCA",
}
DIR_SHORT = {
    "Forward (Micro->Met->Prot->MoCA)": "Forward",
    "Reverse (Prot->Met->Micro->MoCA)": "Reverse",
}

# Colors
MET_COLORS = {
    "TMAO_pre_avg": {"Forward": "#c0392b", "Reverse": "#2b7bb9"},
}

# ?? Figure layout ?????????????????????????????????????????????????????????????
fig, ax_single = plt.subplots(1, 1, figsize=(8, 6))
axes = [ax_single]
fig.patch.set_facecolor("#f8f9fa")

PANEL_TITLES = {
    "TMAO_pre_avg": "帠-Butyrobetaine",
}

for ax, met in zip(axes, targets):
    ax.set_facecolor("#f8f9fa")

    sub = df[df["Metabolite"] == met].copy()

    # y-axis: model ? direction  (interleaved, Forward on top)
    # We place 2 dots per model, separated by small gap
    y_ticks = []
    y_labels = []

    directions = ["Forward (Micro->Met->Prot->MoCA)", "Reverse (Prot->Met->Micro->MoCA)"]

    for i, mod in enumerate(reversed(MODEL_ORDER)):
        for j, direction in enumerate(directions):
            row = sub[(sub["Covariate_Model"] == mod) & (sub["Direction"] == direction)]
            if row.empty:
                continue

            row = row.iloc[0]
            y = i * 2.8 + j * 1.0          # spacing: 2.8 between models, 1.0 between dirs
            ind = row["indirect_serial"]
            lo  = row["CI95_lower"]
            hi  = row["CI95_upper"]
            p   = row["p_value"]

            col = MET_COLORS[met][DIR_SHORT[direction]]

            # CI bar
            ax.hlines(y, lo, hi, colors=col, linewidth=2.5, alpha=0.85)
            # Point
            ax.scatter(ind, y, color=col, s=80, zorder=5, edgecolors="white", linewidth=0.8)
            # p-value annotation
            sig_str = f"p={p:.3f}" if p >= 0.001 else "p<0.001"
            sig_star = ""
            if p < 0.001: sig_star = "***"
            elif p < 0.01: sig_star = "**"
            elif p < 0.05: sig_star = "*"

            ax.text(hi + 0.00005, y, f" {sig_star}  {sig_str}",
                    va="center", ha="left", fontsize=7.5, color="#333333",
                    fontstyle="italic" if sig_star else "normal")

            y_ticks.append(y)
            y_labels.append(f"  {DIR_SHORT[direction]}")

        # Model label at midpoint
        mid_y = i * 2.8 + 0.5
        ax.text(-0.0098, mid_y,
                MODEL_LABELS[mod],
                va="center", ha="left",
                fontsize=8.5, fontweight="bold", color="#222222",
                transform=ax.get_yaxis_transform() if False else ax.transData)

    # Zero line
    ax.axvline(0, color="#555555", linewidth=1.0, linestyle="--", alpha=0.7)

    # Shaded null region
    ax.axvspan(-0.0005, 0.0005, color="#888888", alpha=0.07)

    # Axis formatting
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=8)
    # No inversion needed ??reversed MODEL_ORDER puts Model 1 at top naturally

    ax.set_xlabel("Standardised Indirect Effect (95% Bootstrap CI)",
                  fontsize=10, labelpad=8)
    ax.set_title(PANEL_TITLES[met], fontsize=11, fontweight="bold", pad=10, color="#1a1a2e")
    ax.spines[["top","right","left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=9)
    ax.set_xlim(-0.013, 0.006)

    # Horizontal separator lines between models
    for i in range(1, N_MODELS):
        ax.axhline(i * 2.8 - 0.65, color="#cccccc", linewidth=0.7, linestyle=":")

# Legend
fwd_patch = mpatches.Patch(color="#c0392b", label="Forward  (Micro ??帠-BB ??Prot ??MoCA)")
rev_patch = mpatches.Patch(color="#2b7bb9", label="Reverse  (Prot ??帠-BB ??Micro ??MoCA)")

fig.legend(
    handles=[fwd_patch, rev_patch],
    loc="lower center", bbox_to_anchor=(0.5, 0.0),
    fontsize=9, title_fontsize=9,
    framealpha=0.9, edgecolor="#cccccc",
    ncol=2
)

fig.suptitle(
    "Serial Mediation: 帠-Butyrobetaine as Mediator\n"
    "Microbiome Ratio ??Proteomics Ratio ??MoCA\n"
    r"(log$_2$-transformed; Bootstrap 95% CI, $n$=369; * $p$<0.05)",
    fontsize=12, fontweight="bold", y=1.02, color="#1a1a2e"
)

plt.tight_layout(rect=[0, 0.1, 1, 1])

out_path = os.path.join(OUT_DIR, "mediation_TMAO_sensitivity_forest.png")
fig.savefig(out_path, dpi=600, bbox_inches="tight", facecolor=fig.get_facecolor())
print(f"Saved: {out_path}")
plt.close()

