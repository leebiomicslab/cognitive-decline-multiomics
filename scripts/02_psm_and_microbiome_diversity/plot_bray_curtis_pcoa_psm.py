import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Bray-Curtis PCoA after PSM
Manuscript:       Results - Microbiome diversity (PSM)
Figure/Table:     Supplementary Fig. 7
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
from matplotlib.patches import Ellipse

# ?????????????????????????????????????????????
# Define paths
# ?????????????????????????????????????????????

# Inputs
group_file  = str(PSM_MATCHED_CSV)
asv_file    = os.path.join(project_root, "data", "processed", "MICROBIOTA", "asv_genus_table_reads.csv")
union_file  = str(PROTEOMICS_LOG2)

# Output
output_dir  = str(RESULTS_DIR)
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, "BrayCurtis_PCoA_PSM.png")

# ?????????????????????????????????????????????
# Load data
# ?????????????????????????????????????????????
print("Loading data...")
try:
    # Metadata / group labels
    meta = pd.read_csv(group_file)
    if 'id' in meta.columns:
        meta['id'] = meta['id'].astype(str)
        meta.set_index('id', inplace=True)
    elif 'SampleID' in meta.columns:
        meta['SampleID'] = meta['SampleID'].astype(str)
        meta.set_index('SampleID', inplace=True)
        
    # ASV table: rows = taxa, columns = samples  ?? transpose to samples ? taxa
    asv_df = pd.read_csv(asv_file, index_col=0)
    asv_df = asv_df.T                          # now samples ? taxa
    asv_df.index = asv_df.index.astype(str)

except Exception as e:
    print(f"Error loading files: {e}")
    exit(1)

# ?????????????????????????????????????????????
# Load union genera and filter ASV table
# ?????????????????????????????????????????????
print("Loading union genus set...")
union_df = pd.read_csv(union_file)
target_genera = union_df['taxon'].str.strip().tolist()
print(f"  Target genera in union set: {len(target_genera)}")

asv_cols_lower  = {col: col.lower() for col in asv_df.columns}
target_lower    = {g.lower(): g for g in target_genera}

keep_cols = [col for col, col_l in asv_cols_lower.items() if col_l in target_lower]
print(f"  ASV columns matching target genera: {len(keep_cols)}")

if len(keep_cols) == 0:
    print("No ASV columns matched the target genera. Exiting.")
    exit(1)

asv_filtered = asv_df[keep_cols].copy()

# ?????????????????????????????????????????????
# Restrict to common samples
# ?????????????????????????????????????????????
common_samples = asv_filtered.index.intersection(meta.index)

print(f"Samples in ASV (filtered): {len(asv_filtered)}")
print(f"Samples in Meta (PSM):     {len(meta)}")
print(f"Final common samples:      {len(common_samples)}")

if len(common_samples) == 0:
    print("No common samples found. Exiting.")
    exit(1)

asv_filtered = asv_filtered.loc[common_samples]
meta         = meta.loc[common_samples]
groups       = meta['group']

# Remove zero-sum samples
sample_sums = asv_filtered.sum(axis=1)
valid_mask  = sample_sums > 0
asv_filtered = asv_filtered[valid_mask]
groups       = groups[valid_mask]
print(f"Samples after removing zero-sum: {len(asv_filtered)}")

# ?????????????????????????????????????????????
# Bray-Curtis distance + PCoA
# ?????????????????????????????????????????????
print("Calculating Bray-Curtis distance (union genera only)...")
bc_dm   = beta_diversity('braycurtis', asv_filtered.values, asv_filtered.index)

print("Calculating PCoA...")
bc_pcoa = pcoa(bc_dm)
pc_df   = bc_pcoa.samples[['PC1', 'PC2']].copy()
pc_df['group'] = groups.values

pc1_var = bc_pcoa.proportion_explained['PC1'] * 100
pc2_var = bc_pcoa.proportion_explained['PC2'] * 100

# ?????????????????????????????????????????????
# PERMANOVA
# ?????????????????????????????????????????????
print("Running PERMANOVA...")
perm_res = permanova(bc_dm, groups, permutations=999)
pval     = perm_res['p-value']
fval     = perm_res['test statistic']
n_samples = perm_res['sample size']
n_groups = len(groups.unique())
df_a = n_groups - 1
df_w = n_samples - n_groups
r2 = (fval * df_a) / (fval * df_a + df_w)
print(f"PERMANOVA: F={fval:.2f}, R2={r2:.3f}, p={pval:.3f}")

# ?????????????????????????????????????????????
# Plotting settings
# ?????????????????????????????????????????????
colors    = {"C1": "#ee7a5b", "C2": "#eca362", "C3": "#6fcbbc"}
hue_order = ["C1", "C2", "C3"]

plt.rcParams['font.weight']       = 'bold'
plt.rcParams['axes.labelweight']  = 'bold'
plt.rcParams['axes.titleweight']  = 'bold'
plt.rcParams['font.sans-serif']   = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams.update({
    'font.size':        14,
    'axes.titlesize':   18,
    'axes.labelsize':   16,
    'xtick.labelsize':  12,
    'ytick.labelsize':  12,
    'legend.fontsize':  12,
    'figure.dpi':       600,
})

def plot_confidence_ellipse(ax, x, y, n_std=2, facecolor='none', **kwargs):
    if len(x) < 3:
        return
    cov  = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta  = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)
    ellipse = Ellipse(
        xy=(np.mean(x), np.mean(y)),
        width=width, height=height,
        angle=theta,
        facecolor=facecolor, **kwargs
    )
    ax.add_patch(ellipse)

# ?????????????????????????????????????????????
# Draw plot
# ?????????????????????????????????????????????
print("Plotting...")
plt.figure(figsize=(8, 6))
ax = plt.gca()

sns.scatterplot(
    data=pc_df,
    x='PC1', y='PC2',
    hue='group',
    hue_order=hue_order,
    palette=colors,
    s=70, edgecolor='k', lw=0.5,
    alpha=0.9,
    ax=ax
)

for grp in hue_order:
    if grp in pc_df['group'].unique():
        sub   = pc_df[pc_df['group'] == grp]
        color = colors[grp]
        plot_confidence_ellipse(
            ax,
            sub['PC1'].values, sub['PC2'].values,
            n_std=2,
            edgecolor=color,
            facecolor=color,
            alpha=0.15,
            lw=2
        )

ax.set_xlabel(f"PC1 ({pc1_var:.1f}% variance)")
ax.set_ylabel(f"PC2 ({pc2_var:.1f}% variance)")
ax.set_title(
    f"Bray?urtis PCoA (PSM Matched)\n(PERMANOVA: F={fval:.2f}, R2={r2:.3f}, p={pval:.3g})",
    fontweight='bold'
)

plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.tight_layout()

print(f"Saving plot to {output_file}...")
plt.savefig(output_file, dpi=600, bbox_inches='tight')
print("Done.")



