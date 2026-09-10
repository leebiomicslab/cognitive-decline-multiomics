import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Protein module eigengene boxplots per group
Manuscript:       Results - Protein modules
Figure/Table:     Fig. 4a
"""

from pathlib import Path
"""
Proteomics Module Visualization
1. ME Heatmap (C1/C2/C3)
2. PC1 Loadings Lollipop per module
3. Microbiome ? Proteomics ME Correlation Heatmap
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.stats import pearsonr
from pathlib import Path
import os

# ??? Paths ???????????????????????????????????????????????????????????????????

# Proteomics
PROT_DA      = PROT_ASSOC_DIR / "Main_Model" / "proteomics_ancova_detailed_results.csv"
PROT_EXPR    = PROTEOMICS_LOG2
PROT_MODULE  = FASTSPAR_PROT_DIR / "pearson" / "module_protein_assignment.csv"
PROT_OUT     = FASTSPAR_PROT_DIR / "pearson"

# Microbiota
MICRO_ASV    = MICRO_GENUS_COUNTS
MICRO_MODULE = FASTSPAR_MICRO_DIR / "network_topology" / "leiden_clusters_24_taxa_G_all.csv"

# Shared
CLINICAL     = CLINICAL_CSV

# ??? Helpers ?????????????????????????????????????????????????????????????????
def group_map_from_clinical():
    clin = pd.read_csv(CLINICAL)
    clin.rename(columns={'id': 'SampleID'}, inplace=True, errors='ignore')
    if 'SampleID' not in clin.columns:
        clin.rename(columns={clin.columns[0]: 'SampleID'}, inplace=True)
    if 'group' in clin.columns:
        clin.rename(columns={'group': 'Group'}, inplace=True)
    return dict(zip(clin['SampleID'].astype(str), clin['Group']))

def clr_transformation(expr_T: pd.DataFrame) -> pd.DataFrame:
    """Input: samples ? taxa (raw counts). Returns CLR."""
    values = expr_T.values.flatten()
    eps = np.min(values[values > 0]) * 0.65
    df_imp = expr_T.replace(0, eps).astype(float)
    log_df = np.log(df_imp)
    return log_df.sub(log_df.mean(axis=1), axis=0)

# ????????????????????????????????????????????????????????????????????????????????# 1. PROTEOMICS ME HEATMAP + LOADINGS
# ????????????????????????????????????????????????????????????????????????????????
print("=== Loading proteomics data ===")
gmap = group_map_from_clinical()

# Load significant proteins
da = pd.read_csv(PROT_DA)
sig_proteins = set(da[da['P_Value'] < 0.05]['Protein'].tolist())

# Load expression
prot_raw = pd.read_csv(PROT_EXPR)
prot_raw.rename(columns={prot_raw.columns[0]: 'SampleID'}, inplace=True)
valid_cols = [c for c in prot_raw.columns if c in sig_proteins]
prot_df = prot_raw[['SampleID'] + valid_cols].copy()
prot_df['Group'] = prot_df['SampleID'].astype(str).map(gmap)
prot_df.dropna(subset=['Group'], inplace=True)
prot_df.set_index('SampleID', inplace=True)

# Load module assignments
mod_assign = pd.read_csv(PROT_MODULE)
clusters = sorted(mod_assign['Module'].unique())

# ?? Compute ME & Loadings ??
me_dict_prot   = {}
loadings_prot  = []

for cid in clusters:
    proteins = mod_assign[mod_assign['Module'] == cid]['Protein'].tolist()
    valid_p  = [p for p in proteins if p in prot_df.columns]
    if not valid_p:
        continue
    sub = prot_df[valid_p].astype(float)
    pca = PCA(n_components=1)
    me_scores = pca.fit_transform(sub).flatten()
    loads     = pca.components_[0]
    if np.mean(loads) < 0:
        loads     = -loads
        me_scores = -me_scores
    me_dict_prot[f"P{cid+1}"] = pd.Series(me_scores, index=prot_df.index)
    for p, l in zip(valid_p, loads):
        loadings_prot.append({"Module": cid, "Protein": p, "PC1_Loading": l})

df_me_prot = pd.DataFrame(me_dict_prot)
df_me_prot['Group'] = df_me_prot.index.astype(str).map(gmap)
df_me_prot.dropna(subset=['Group'], inplace=True)
me_mean_prot = df_me_prot.groupby('Group').mean().T[["C1", "C2", "C3"]]

df_loads_prot = pd.DataFrame(loadings_prot)
df_loads_prot.to_csv(PROT_OUT / "module_me_loadings.csv", index=False)

print(f"  Computed ME for {len(me_dict_prot)} proteomics modules")

# ?? Plot 1a: ME Heatmap ??
vmax = np.abs(me_mean_prot.values).max() * 1.1
plt.figure(figsize=(8, max(4, len(clusters) * 1.2)))
ax = sns.heatmap(me_mean_prot, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                 annot=True, fmt=".2f", annot_kws={"size": 14, "weight": "bold"},
                 linewidths=2.0, linecolor='white')
ax.set_title("Mean Module Eigengene (ME) Score\nProteomics Leiden Clusters ? Cognitive Groups",
             fontsize=15, fontweight='bold', pad=15)
ax.set_ylabel("")
ax.set_xlabel("")
plt.setp(ax.get_xticklabels(), fontsize=13, fontweight='bold')
plt.setp(ax.get_yticklabels(), fontsize=12, fontweight='bold', rotation=0)
plt.tight_layout()
plt.savefig(PROT_OUT / "proteomics_ME_heatmap.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: proteomics_ME_heatmap.png")

# ?? Plot 1c: Boxplot per proteomics module (C1/C2/C3) ??
from scipy import stats
from itertools import combinations

COLORS = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'}
MODULE_LABEL_PROT = {f"P{cid+1}": f"P{cid+1}" for cid in clusters if f"P{cid+1}" in df_me_prot.columns}
GROUP_ORDER = ['C1', 'C2', 'C3']
PAIRS = list(combinations(GROUP_ORDER, 2))

def pairwise_bh(group_vals, group_order, pairs):
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
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.2, color=color)
    ax.text((x1+x2)/2, y+h, label, ha='center', va='bottom',
            fontsize=9, fontweight='bold', color=color)

n_modules_prot = len(MODULE_LABEL_PROT)
fig_bp, axes_bp = plt.subplots(1, n_modules_prot, figsize=(4.5 * n_modules_prot, 5.5), sharey=False)
if n_modules_prot == 1:
    axes_bp = [axes_bp]

for ax, (col_name, mlabel) in zip(axes_bp, MODULE_LABEL_PROT.items()):
    plot_data = df_me_prot[[col_name, 'Group']].copy()
    plot_data.rename(columns={col_name: 'ME'}, inplace=True)
    plot_data = plot_data[plot_data['Group'].isin(GROUP_ORDER)]

    group_vals = [plot_data[plot_data['Group'] == g]['ME'].values for g in GROUP_ORDER]
    from statsmodels.stats.oneway import anova_oneway
    res = anova_oneway(group_vals, use_var='unequal')
    kw_p = res.pvalue
    anova_type = "Welch ANOVA"
    pairwise_type = "Welch t-test"

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

    for xi, grp in enumerate(GROUP_ORDER):
        vals = group_vals[xi]
        jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(vals))
        ax.scatter(xi + jitter, vals, color=COLORS[grp],
                   s=18, alpha=0.6, zorder=3, edgecolors='white', linewidths=0.4)

    y_max = max(v.max() for v in group_vals if len(v) > 0)
    y_min = min(v.min() for v in group_vals if len(v) > 0)
    y_range = y_max - y_min
    h_unit  = y_range * 0.06
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
    ax.set_ylim(y_min - y_range * 0.05, bracket_starts[-1] + h_unit * 2.5)

fig_bp.suptitle("Proteomics Module Eigengene Score by Cognitive Group", fontsize=15, fontweight='bold', y=1.01)
plt.tight_layout()
boxplot_prot_png = PROT_OUT / "proteomics_ME_boxplot.png"
plt.savefig(boxplot_prot_png, dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print("  [DONE] Saved: proteomics_ME_boxplot.png\n")
# ?? Plot 1b: Loadings Lollipop (one panel per module) ??
n_mods = len(clusters)
fig, axes = plt.subplots(n_mods, 1, figsize=(11, max(12, n_mods * 4.5)),
                          gridspec_kw={'hspace': 0.45})
if n_mods == 1:
    axes = [axes]

for i, cid in enumerate(clusters):
    ax  = axes[i]
    sub = df_loads_prot[df_loads_prot['Module'] == cid].sort_values('PC1_Loading', ascending=True)
    colors = ['#E41A1C' if v >= 0 else '#377EB8' for v in sub['PC1_Loading']]
    ax.hlines(sub['Protein'], 0, sub['PC1_Loading'], color=colors, linewidth=3.0, alpha=0.8)
    ax.scatter(sub['PC1_Loading'], sub['Protein'], color=colors, s=150,
               edgecolors='black', linewidth=1.5, zorder=3)
    ax.axvline(0, color='gray', linestyle='--', linewidth=1.5)
    ax.set_title(f"P{cid+1} ??PC1 Loadings", fontsize=14, fontweight='bold', pad=12)
    ax.set_xlabel("PC1 Loading", fontsize=11, fontweight='bold')
    ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='x', labelsize=10)
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    for spine in ['right','top','left']: ax.spines[spine].set_visible(False)
    xmin, xmax = ax.get_xlim()
    xr = xmax - xmin
    for _, row in sub.iterrows():
        off = xr * 0.02
        ha  = 'left' if row['PC1_Loading'] >= 0 else 'right'
        ax.text(row['PC1_Loading'] + (off if ha=='left' else -off),
                row['Protein'], f"{row['PC1_Loading']:.2f}",
                va='center', ha=ha, fontsize=9, fontweight='bold')

plt.savefig(PROT_OUT / "proteomics_module_loadings_lollipop.png", dpi=300,
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: proteomics_module_loadings_lollipop.png")


# ????????????????????????????????????????????????????????????????????????????????# 2. MICROBIOME ME (recompute from CLR)
# ????????????????????????????????????????????????????????????????????????????????print("\n=== Computing microbiome MEs ===")
asv_raw = pd.read_csv(MICRO_ASV)
asv_raw.set_index(asv_raw.columns[0], inplace=True)
asv_raw = asv_raw.loc[:, ~asv_raw.columns.duplicated()]
clr_df  = clr_transformation(asv_raw.T.astype(float))

micro_assign = pd.read_csv(MICRO_MODULE)
micro_clusters = sorted(micro_assign['Leiden_Cluster'].unique())

me_dict_micro = {}
for cid in micro_clusters:
    taxa   = micro_assign[micro_assign['Leiden_Cluster'] == cid]['Taxon'].tolist()
    v_taxa = [t for t in taxa if t in clr_df.columns]
    if not v_taxa:
        continue
    sub = clr_df[v_taxa]
    pca = PCA(n_components=1)
    me_scores = pca.fit_transform(sub).flatten()
    loads     = pca.components_[0]
    if np.mean(loads) < 0:
        me_scores = -me_scores
    me_dict_micro[f"M{cid+1}"] = pd.Series(me_scores, index=clr_df.index)

df_me_micro = pd.DataFrame(me_dict_micro)
print(f"  Computed ME for {len(me_dict_micro)} microbiome modules")


# ????????????????????????????????????????????????????????????????????????????????# 3. MICROBIOME ? PROTEOMICS ME CORRELATION HEATMAP (per group: C1, C2, C3)
# ????????????????????????????????????????????????????????????????????????????????print("\n=== Computing Microbiome ? Proteomics ME correlation (per group) ===")

# Merge group info into microbiome ME
df_me_micro_grp = df_me_micro.copy()
df_me_micro_grp['Group'] = df_me_micro_grp.index.astype(str).map(gmap)
df_me_micro_grp.dropna(subset=['Group'], inplace=True)

# Align samples
common_idx = df_me_prot.index.intersection(df_me_micro_grp.index)
prot_cols  = [c for c in df_me_prot.columns if c != 'Group']
micro_cols = [c for c in df_me_micro_grp.columns if c != 'Group']

prot_me_sub  = df_me_prot.loc[common_idx, prot_cols]
micro_me_sub = df_me_micro_grp.loc[common_idx, micro_cols]
group_labels  = df_me_prot.loc[common_idx, 'Group']

print(f"  Common samples: {len(common_idx)}")

GROUPS = ["C1", "C2", "C3"]

# Annotation helper
def make_annot(r_df, p_df):
    annot = r_df.copy().astype(str)
    for row in r_df.index:
        for col in r_df.columns:
            r = r_df.loc[row, col]
            p = p_df.loc[row, col]
            if pd.isna(r):
                annot.loc[row, col] = "NA"
            elif p < 0.01:
                annot.loc[row, col] = f"{r:.2f}**"
            elif p < 0.05:
                annot.loc[row, col] = f"{r:.2f}*"
            else:
                annot.loc[row, col] = f"{r:.2f}"
    return annot

# Compute per-group r / p matrices
group_r = {}
group_p = {}
for grp in GROUPS:
    idx_g = group_labels[group_labels == grp].index
    pm  = prot_me_sub.loc[idx_g]
    mm  = micro_me_sub.loc[idx_g]
    r_m = pd.DataFrame(index=pm.columns, columns=mm.columns, dtype=float)
    p_m = pd.DataFrame(index=pm.columns, columns=mm.columns, dtype=float)
    for pmod in pm.columns:
        for mmod in mm.columns:
            x = pm[pmod].values.astype(float)
            y = mm[mmod].values.astype(float)
            mask = np.isfinite(x) & np.isfinite(y)
            if mask.sum() < 5:
                r_m.loc[pmod, mmod] = np.nan
                p_m.loc[pmod, mmod] = np.nan
            else:
                rv, pv = pearsonr(x[mask], y[mask])
                r_m.loc[pmod, mmod] = rv
                p_m.loc[pmod, mmod] = pv
    group_r[grp] = r_m.astype(float)
    group_p[grp] = p_m.astype(float)
    print(f"  {grp}: n={len(idx_g)} samples")

# Global colour scale across all three groups
all_r_vals = np.concatenate([group_r[g].values.flatten() for g in GROUPS])
vmax_r = max(0.5, np.nanmax(np.abs(all_r_vals)))

n_micro = len(micro_cols)
n_prot  = len(prot_cols)

cell_w  = 1.8   # width per microbiome module column
cell_h  = 1.4   # height per proteomics module row
cbar_w  = 0.6   # colourbar
gap     = 0.5   # gap between panels

total_w = 3 * (n_micro * cell_w) + 2 * gap + cbar_w + 1.0
total_h = max(4, n_prot * cell_h + 2.5)

fig, axes = plt.subplots(1, 3,
                          figsize=(total_w, total_h),
                          gridspec_kw={'wspace': gap / (n_micro * cell_w)})

GROUP_COLORS = {"C1": "#4393C3", "C2": "#66C2A5", "C3": "#FC8D62"}

for i, grp in enumerate(GROUPS):
    ax  = axes[i]
    r_m = group_r[grp]
    p_m = group_p[grp]
    ann = make_annot(r_m, p_m)
    n_g = (group_labels == grp).sum()

    # Only show colourbar on last panel
    cbar = (i == 2)
    sns.heatmap(r_m, annot=ann, fmt='', cmap='RdBu_r',
                vmin=-vmax_r, vmax=vmax_r, center=0,
                linewidths=1.5, linecolor='white', ax=ax,
                annot_kws={"size": 12, "weight": "bold"},
                cbar=cbar,
                cbar_kws={'shrink': 0.8, 'label': 'Pearson r'} if cbar else {})

    ax.set_title(f"{grp}  (n={n_g})",
                 fontsize=15, fontweight='bold', pad=10,
                 color=GROUP_COLORS[grp])
    ax.set_ylabel("Proteomics Module" if i == 0 else "",
                  fontsize=12, fontweight='bold')
    ax.set_xlabel("Microbiome Module", fontsize=12, fontweight='bold')
    plt.setp(ax.get_xticklabels(), fontsize=12, fontweight='bold',
             rotation=0, ha='center')
    plt.setp(ax.get_yticklabels(), fontsize=12, fontweight='bold', rotation=0)
    if i > 0:
        ax.set_yticklabels([])

fig.suptitle("Microbiome ? Proteomics ME Correlation (Pearson r)\n* p<0.05  ** p<0.01",
             fontsize=16, fontweight='bold', y=1.02)

CROSS_OUT = RESULTS_DIR / "fastspar"
plt.tight_layout()
plt.savefig(CROSS_OUT / "microbiome_x_proteomics_ME_corr.png", dpi=300,
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: microbiome_x_proteomics_ME_corr.png  (C1 / C2 / C3 panels)")

# Save per-group correlation tables
for grp in GROUPS:
    group_r[grp].to_csv(CROSS_OUT / f"microbiome_x_proteomics_ME_r_{grp}.csv")
    group_p[grp].to_csv(CROSS_OUT / f"microbiome_x_proteomics_ME_pval_{grp}.csv")

print("\n[DONE] All visualizations complete!")


