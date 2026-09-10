import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Bray-Curtis PCoA of microbiome beta diversity
Manuscript:       Results - Microbiome diversity
Figure/Table:     Fig. 2a
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

# Define paths

# Inputs
group_file = str(GROUP_CSV)
# Using Species level table as it's the most granular available (CSV)
asv_file = os.path.join(project_root, "data", "processed", "MICROBIOTA", "asv_species_table_reads.csv")

# Output
output_dir = str(RESULTS_DIR)
output_file = os.path.join(output_dir, "BrayCurtis_PCoA_withCI_stylish.png")

os.makedirs(output_dir, exist_ok=True)

# Load Data
print("Loading data...")
try:
    # Load Groups
    meta = pd.read_csv(group_file)
    # Ensure columns are correct: assuming metadata file has 'id' and 'group'
    # Check header from previous steps: ['id', 'group']
    if 'group' not in meta.columns:
        # Fallback if header is different or missing
        # Based on previous `head` output: it has headers 'id', 'group'
        pass
    meta['id'] = meta['id'].astype(str)
    meta.set_index('id', inplace=True)

    # Load ASV Table
    # Rows = Features, Cols = Samples (based on previous check)
    asv_df = pd.read_csv(asv_file, index_col=0)
    # Transpose to Samples x Features for skbio
    asv_df = asv_df.T
    asv_df.index = asv_df.index.astype(str)

except Exception as e:
    print(f"Error loading files: {e}")
    exit(1)

# Load common samples
common_samples_file = str(COMMON_SAMPLES_CSV)
if os.path.exists(common_samples_file):
    print(f"Loading common samples from {common_samples_file}...")
    common_samples_list = pd.read_csv(common_samples_file)['id'].astype(str).tolist()
else:
    print("Common samples file not found. Using all samples.")
    common_samples_list = []

# Intersect samples
common_samples = asv_df.index.intersection(meta.index)
if common_samples_list:
    common_samples = common_samples.intersection(common_samples_list)

print(f"Samples in ASV: {len(asv_df)}")
print(f"Samples in Meta: {len(meta)}")
print(f"Final Common samples (intersected with all omics): {len(common_samples)}")

if len(common_samples) == 0:
    print("No common samples found!")
    exit(1)

figures_dir = os.path.join(output_dir, "figures")
os.makedirs(figures_dir, exist_ok=True)
output_file = os.path.join(figures_dir, "BrayCurtis_PCoA_withCI_stylish.png")

asv_df = asv_df.loc[common_samples]
meta = meta.loc[common_samples]
groups = meta['group']

# Filter zero sum samples if any
# (Although intersection usually handles valid IDs, skbio needs non-zero vectors)
sample_sums = asv_df.sum(axis=1)
valid_mask = sample_sums > 0
asv_df = asv_df[valid_mask]
groups = groups[valid_mask]
print(f"Samples after removing zero-sum: {len(asv_df)}")

# Beta Diversity (Bray-Curtis)
print("Calculating Bray-Curtis distance...")
bc_dm = beta_diversity('braycurtis', asv_df.values, asv_df.index)

# PCoA
print("Calculating PCoA...")
bc_pcoa = pcoa(bc_dm)
pc_df = bc_pcoa.samples[['PC1', 'PC2']]
pc_df['group'] = groups

pc1_var = bc_pcoa.proportion_explained['PC1'] * 100
pc2_var = bc_pcoa.proportion_explained['PC2'] * 100

# PERMANOVA
print("Running PERMANOVA...")
perm_res = permanova(bc_dm, groups, permutations=999)
pval = perm_res['p-value']
fval = perm_res['test statistic']
n_samples = perm_res['sample size']
n_groups = len(groups.unique())
df_a = n_groups - 1
df_w = n_samples - n_groups
r2 = (fval * df_a) / (fval * df_a + df_w)
print(f"PERMANOVA: F={fval:.2f}, R2={r2:.3f}, p={pval:.3f}")

# Plotting Settings
colors = {
    "C1": "#ee7a5b",
    "C2": "#eca362",
    "C3": "#6fcbbc"
}
hue_order = ["C1", "C2", "C3"]

plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.dpi': 300
})

def plot_confidence_ellipse(ax, x, y, n_std=2, facecolor='none', **kwargs):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)
    ellipse = Ellipse(
        xy=(np.mean(x), np.mean(y)),
        width=width, height=height,
        angle=theta,
        facecolor=facecolor, **kwargs
    )
    ax.add_patch(ellipse)

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

# Ellipses
for grp in hue_order:
    if grp in pc_df['group'].unique():
        sub = pc_df[pc_df['group'] == grp]
        color = colors[grp]
        plot_confidence_ellipse(
            ax,
            sub['PC1'], sub['PC2'],
            n_std=2,
            edgecolor=color,
            facecolor=color,
            alpha=0.15,
            lw=2
        )

ax.set_xlabel(f"PC1 ({pc1_var:.1f}% variance)")
ax.set_ylabel(f"PC2 ({pc2_var:.1f}% variance)")
ax.set_title(f"Bray?urtis PCoA (PERMANOVA: F={fval:.2f}, R2={r2:.3f}, p={pval:.3g})", fontweight='bold')

plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.tight_layout()

print(f"Saving plot to {output_file}...")
plt.savefig(output_file, dpi=600, bbox_inches='tight')
print("Done.")



