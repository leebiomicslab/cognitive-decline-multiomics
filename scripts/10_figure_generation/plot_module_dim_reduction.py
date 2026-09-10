import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          PCA/UMAP of module eigengene space
Manuscript:       Methods - Module structure
Figure/Table:     Supplementary Fig. 6
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import umap
from pathlib import Path

# ??? Config ??????????????????????????????????????????????????????????????????
OUT_DIR = RESULTS_DIR / "fastspar" / "module_visualizations"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CLINICAL     = CLINICAL_CSV
MICRO_ASV    = MICRO_GENUS_COUNTS
MICRO_MODULE = FASTSPAR_MICRO_DIR / "network_topology" / "leiden_clusters_24_taxa_G_all.csv"
PROT_EXPR    = PROTEOMICS_LOG2
PROT_MODULE  = FASTSPAR_PROT_DIR / "pearson" / "module_protein_assignment.csv"

COLORS = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'}
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']

# ??? Helpers ?????????????????????????????????????????????????????????????????
def get_clr(df_raw, eps):
    df_imp = df_raw.replace(0, eps).astype(float)
    log_df = np.log(df_imp)
    return log_df.sub(log_df.mean(axis=1), axis=0)

def plot_module_viz(data, labels, module_name, output_path):
    """Generates a 1x3 grid containing PCA, t-SNE, and UMAP."""
    print(f"  - Processing {module_name}...")
    
    # Standardize for DR
    X_scaled = StandardScaler().fit_transform(data)
    
    # 1. PCA
    pca = PCA(n_components=2, random_state=42).fit_transform(X_scaled)
    # 2. t-SNE
    tsne = TSNE(n_components=2, perplexity=min(30, len(data)-1), random_state=42).fit_transform(X_scaled)
    # 3. UMAP
    reducer = umap.UMAP(n_neighbors=min(15, len(data)-1), random_state=42)
    u_map = reducer.fit_transform(X_scaled)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    methods = [
        ('PCA', pca, 'PC1', 'PC2'),
        ('t-SNE', tsne, 't-SNE 1', 't-SNE 2'),
        ('UMAP', u_map, 'UMAP 1', 'UMAP 2')
    ]
    
    for ax, (title, coords, xl, yl) in zip(axes, methods):
        sns.scatterplot(
            x=coords[:, 0], y=coords[:, 1], hue=labels, 
            palette=COLORS, ax=ax, s=60, alpha=0.8, edgecolor='w'
        )
        ax.set_title(f"{title}", fontsize=14, fontweight='bold')
        ax.set_xlabel(xl, fontweight='bold')
        ax.set_ylabel(yl, fontweight='bold')
        ax.grid(alpha=0.3)
        ax.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.suptitle(f"Dimensionality Reduction for Module: {module_name}", fontsize=18, fontweight='bold', y=1.05)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

# ??? Data Loading ????????????????????????????????????????????????????????????
print("=== Loading Metadata and Expression Data ===")
clin = pd.read_csv(CLINICAL).rename(columns={'id': 'SampleID', 'group': 'Group'})
clin['SampleID'] = clin['SampleID'].astype(str)
y_all = clin[['SampleID', 'Group']].dropna()
y_all = y_all[y_all['Group'].isin(['C1', 'C2', 'C3'])].set_index('SampleID')

# Microbiota
asv_csv = pd.read_csv(MICRO_ASV)
asv_raw = asv_csv.set_index(asv_csv.columns[0]).loc[:, ~asv_csv.columns[1:].duplicated()].T
asv_raw.index = asv_raw.index.astype(str)
m_micro = pd.read_csv(MICRO_MODULE)

# Proteomics
prot_raw = pd.read_csv(PROT_EXPR).rename(columns={pd.read_csv(PROT_EXPR).columns[0]: 'SampleID'})
prot_raw['SampleID'] = prot_raw['SampleID'].astype(str)
prot_raw.set_index('SampleID', inplace=True)
m_prot  = pd.read_csv(PROT_MODULE)

# Find common samples
common = y_all.index.intersection(asv_raw.index).intersection(prot_raw.index)
y_sub = y_all.loc[common, 'Group']

# ??? Execution ???????????????????????????????????????????????????????????????
print("\n>>> Generating Microbiota Module Plots (M1-M3)...")
eps = np.min(asv_raw.values[asv_raw.values > 0]) * 0.65
asv_clr = get_clr(asv_raw.loc[common], eps)

for cid in sorted(m_micro['Leiden_Cluster'].unique()):
    m_name = f"M{cid+1}"
    taxa = m_micro[m_micro['Leiden_Cluster']==cid]['Taxon'].tolist()
    data_sub = asv_clr[[t for t in taxa if t in asv_clr.columns]]
    plot_module_viz(data_sub, y_sub, m_name, OUT_DIR / f"{m_name}_dim_reduction.png")

print("\n>>> Generating Proteomics Module Plots (P1-P5)...")
prot_sub = prot_raw.loc[common]

for cid in sorted(m_prot['Module'].unique()):
    m_name = f"P{cid+1}"
    proteins = m_prot[m_prot['Module']==cid]['Protein'].tolist()
    data_sub = prot_sub[[p for p in proteins if p in prot_sub.columns]]
    plot_module_viz(data_sub, y_sub, m_name, OUT_DIR / f"{m_name}_dim_reduction.png")

print(f"\n??All 8 module visualizations saved to: {OUT_DIR}")


