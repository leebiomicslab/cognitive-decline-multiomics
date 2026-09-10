import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Module eigengene PCA loadings for Leiden microbial clusters
Manuscript:       Methods - Module eigengenes
Figure/Table:     Supplementary Fig. 5
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from pathlib import Path
import os

ASV_FILE = MICRO_GENUS_COUNTS
CLINICAL_FILE = CLINICAL_CSV
CLUSTER_FILE = FASTSPAR_MICRO_DIR / "network_topology" / "leiden_clusters_24_taxa_G_all.csv"
OUTPUT_DIR = FASTSPAR_MICRO_DIR / "network_topology"

def clr_transformation(df, epsilon=None):
    if epsilon is None:
        values = df.values.flatten()
        min_pos = np.min(values[values > 0])
        epsilon = min_pos * 0.65
    df_imp = df.replace(0, epsilon)
    log_df = np.log(df_imp)
    return log_df.sub(log_df.mean(axis=1), axis=0)

# 1. Load clinical & CLR data
clinical = pd.read_csv(CLINICAL_FILE)
clinical.rename(columns={'id': 'SampleID'}, inplace=True)
group_map = dict(zip(clinical['SampleID'].astype(str), clinical['group']))

expr_raw = pd.read_csv(ASV_FILE)
expr_raw.set_index(expr_raw.columns[0], inplace=True)
expr_raw = expr_raw.loc[:, ~expr_raw.columns.duplicated()]
clr_df = clr_transformation(expr_raw.T.astype(float))

# 2. Load Clusters
clusters = pd.read_csv(CLUSTER_FILE)

# 3. Compute MEs and Loadings per Cluster
me_dict = {}
loadings_res = []

for cid in sorted(clusters['Leiden_Cluster'].unique()):
    taxa = clusters[clusters['Leiden_Cluster'] == cid]['Taxon'].tolist()
    # Intersect with available CLR data
    valid_taxa = [t for t in taxa if t in clr_df.columns]
    
    sub_df = clr_df[valid_taxa]
    pca = PCA(n_components=1)
    # PCA score (Eigengene)
    me_scores = pca.fit_transform(sub_df).flatten()
    loadings = pca.components_[0]
    
    # Standardize sign: Make the ME positively correlated with the majority of its taxa
    if np.mean(loadings) < 0:
        loadings = -loadings
        me_scores = -me_scores
        
    me_dict[f"Cluster {cid}"] = me_scores
    
    for t, l in zip(valid_taxa, loadings):
        loadings_res.append({"Leiden_Cluster": cid, "Taxon": t, "PC1_Loading": l})

df_loadings = pd.DataFrame(loadings_res)
df_loadings.to_csv(OUTPUT_DIR / "leiden_clusters_loadings.csv", index=False)

# Format loadings to friendly markdown using print
print("\n=== PC1 Loadings by Cluster ===")
for cid in sorted(df_loadings['Leiden_Cluster'].unique()):
    print(f"\n>> Cluster {cid} Loadings:")
    sub = df_loadings[df_loadings['Leiden_Cluster'] == cid].sort_values(by="PC1_Loading", ascending=False)
    for _, row in sub.iterrows():
        print(f"  - {row['Taxon']:<35} : {row['PC1_Loading']:.4f}")

# 4. ME Heatmap
me_df = pd.DataFrame(me_dict, index=clr_df.index)
me_df['Group'] = me_df.index.astype(str).map(group_map)

# Drop NA groups if any, group by cognitive group C1, C2, C3
me_df.dropna(subset=['Group'], inplace=True)
me_mean = me_df.groupby('Group').mean().T
# Ensure order
me_mean = me_mean[["C1", "C2", "C3"]]

# Plot Heatmap
sns.set_theme(style="white")
# Scale robustly for diverging colormap
vmax = np.abs(me_mean.values).max() * 1.1

plt.figure(figsize=(8, 5))
ax = sns.heatmap(me_mean, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                 annot=True, fmt=".2f", annot_kws={"size": 14, "weight": "bold"},
                 linewidths=2.0, linecolor='white')

ax.set_title("Mean Module Eigengene (ME) Score\nof Leiden Clusters Across Cognitive Groups", fontsize=16, fontweight='bold', pad=20)
ax.set_ylabel("")
ax.set_xlabel("")
plt.setp(ax.get_xticklabels(), fontsize=14, fontweight='bold')
plt.setp(ax.get_yticklabels(), fontsize=13, fontweight='bold', rotation=0)

out_png = OUTPUT_DIR / "leiden_clusters_ME_heatmap.png"
plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"\n[DONE] Heatmap saved to {out_png}")

# ?????????????????????????????????????????????????????????????????????????????
# 5. Boxplot: each module ? three groups  (C1 / C2 / C3)
# ?????????????????????????????????????????????????????????????????????????????
from scipy import stats
from itertools import combinations

COLORS = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'}
MODULE_LABEL = {cid: f"M{i+1}" for i, cid in
                enumerate(sorted(clusters['Leiden_Cluster'].unique()))}
GROUP_ORDER = ['C1', 'C2', 'C3']
PAIRS = list(combinations(GROUP_ORDER, 2))  # (C1,C2),(C1,C3),(C2,C3)

def pairwise_bh(group_vals, group_order, pairs):
    """
    Pairwise Welch's t-tests with Benjamini-Hochberg FDR correction.
    """
    from statsmodels.stats.multitest import multipletests
    raw_ps = []
    for g1, g2 in pairs:
        a = group_vals[group_order.index(g1)]
        b = group_vals[group_order.index(g2)]
        _, p = stats.ttest_ind(a, b, equal_var=False)
        raw_ps.append(p)
    reject, adj_ps, _, _ = multipletests(raw_ps, method='fdr_bh')
    return {pair: adj_p for pair, adj_p in zip(pairs, adj_ps)}

def sig_label(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

def draw_bracket(ax, x1, x2, y, h, label, color='black'):
    """Draw a significance bracket between x1 and x2 at height y."""
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.2, color=color)
    ax.text((x1+x2)/2, y+h, label, ha='center', va='bottom',
            fontsize=9, fontweight='bold', color=color)

n_modules = len(MODULE_LABEL)
fig, axes = plt.subplots(1, n_modules, figsize=(4.5 * n_modules, 5.5),
                         sharey=False)
if n_modules == 1:
    axes = [axes]

for ax, (cid, mlabel) in zip(axes, MODULE_LABEL.items()):
    col_name = f"Cluster {cid}"
    plot_data = me_df[[col_name, 'Group']].copy()
    plot_data.rename(columns={col_name: 'ME'}, inplace=True)
    plot_data = plot_data[plot_data['Group'].isin(GROUP_ORDER)]

    # Main ANOVA test
    group_vals = [plot_data[plot_data['Group'] == g]['ME'].values for g in GROUP_ORDER]
    from statsmodels.stats.oneway import anova_oneway
    res = anova_oneway(group_vals, use_var='unequal')
    kw_p = res.pvalue
    anova_type = "Welch ANOVA"
    pairwise_type = "Welch t-test"

    # Boxplot
    bp = ax.boxplot(
        [group_vals[i] for i in range(len(GROUP_ORDER))],
        positions=range(len(GROUP_ORDER)),
        widths=0.5,
        patch_artist=True,
        medianprops=dict(color='black', linewidth=2),
        whiskerprops=dict(linewidth=1.2),
        capprops=dict(linewidth=1.2),
        flierprops=dict(marker='o', markersize=3, alpha=0.5),
        boxprops=dict(linewidth=1.2),
    )
    for patch, grp in zip(bp['boxes'], GROUP_ORDER):
        patch.set_facecolor(COLORS[grp])
        patch.set_alpha(0.75)

    # Strip (jitter)
    for xi, grp in enumerate(GROUP_ORDER):
        vals = group_vals[xi]
        jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(vals))
        ax.scatter(xi + jitter, vals, color=COLORS[grp],
                   s=18, alpha=0.6, zorder=3, edgecolors='white', linewidths=0.4)

    # Pairwise brackets
    y_max = max(v.max() for v in group_vals if len(v) > 0)
    y_min = min(v.min() for v in group_vals if len(v) > 0)
    y_range = y_max - y_min
    h_unit  = y_range * 0.06   # bracket arm height
    bracket_starts = [y_max + y_range * 0.12,
                      y_max + y_range * 0.27,
                      y_max + y_range * 0.42]

    adj_p_dict = pairwise_bh(group_vals, GROUP_ORDER, PAIRS)
    x_pos = {g: i for i, g in enumerate(GROUP_ORDER)}
    for bi, (g1, g2) in enumerate(PAIRS):
        p_val = adj_p_dict[(g1, g2)]
        label = sig_label(p_val)
        draw_bracket(ax, x_pos[g1], x_pos[g2],
                     bracket_starts[bi], h_unit, label)

    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels(GROUP_ORDER, fontsize=12, fontweight='bold')
    ax.set_ylabel("ME Score", fontsize=11, fontweight='bold')
    ax.set_title(
        f"{mlabel}\n{anova_type} p={kw_p:.3f}{'*' if kw_p < 0.05 else ''}\npairwise: {pairwise_type} BH",
        fontsize=10, fontweight='bold')
    ax.axhline(0, color='grey', lw=0.8, ls='--', alpha=0.5)
    ax.spines[['top', 'right']].set_visible(False)
    # expand y-axis to fit brackets
    ax.set_ylim(y_min - y_range * 0.05,
                bracket_starts[-1] + h_unit * 2.5)

fig.suptitle("Module Eigengene Score by Cognitive Group",
             fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
boxplot_png = OUTPUT_DIR / "leiden_clusters_ME_boxplot.png"
plt.savefig(boxplot_png, dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print(f"[DONE] Boxplot saved to {boxplot_png}")



