import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Heatmap of differentially abundant microbiome genera
Manuscript:       Results - Microbiome DA
Figure/Table:     Fig. 5b
"""

from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Paths
DA_RESULTS = MICRO_ASSOC_DIR / "Model_A_Main" / "microbiota_da_results_final.csv"
CLINICAL_FILE = CLINICAL_CSV
ASV_FILE = MICRO_GENUS_COUNTS
OUTPUT_DIR = FASTSPAR_MICRO_DIR / "network_topology"

def clr_transformation(df, epsilon=None):
    if epsilon is None:
        values = df.values.flatten()
        min_pos = np.min(values[values > 0])
        epsilon = min_pos * 0.65
    df_imp = df.replace(0, epsilon)
    log_df = np.log(df_imp)
    clr_df = log_df.sub(log_df.mean(axis=1), axis=0)
    return clr_df

# 1. Compute Mean CLR for each Group
print("Loading data for CLR calculation...")
clinical = pd.read_csv(CLINICAL_FILE)
clinical.rename(columns={'id': 'SampleID'}, inplace=True)

expr_raw = pd.read_csv(ASV_FILE)
genus_col = expr_raw.columns[0]
expr_raw.set_index(genus_col, inplace=True)
expr_raw = expr_raw.loc[:, ~expr_raw.columns.duplicated()]

expr = expr_raw.T
print("Applying CLR transformation...")
clr_df_full = clr_transformation(expr.astype(float))
clr_df_full['SampleID'] = clr_df_full.index.astype(str)

# 2. Get Significant DA Taxa & Modules
da = pd.read_csv(DA_RESULTS)
sig_taxa = da[da["Union_Sig"] == True]["taxon"].str.strip().tolist()

# Load mapping for module sidebar
LEIDEN_LIST = FASTSPAR_MICRO_DIR / "network_topology" / "leiden_clusters_24_taxa_G_all.csv"
taxa_mod = pd.read_csv(LEIDEN_LIST)
taxa_mod['Taxon'] = taxa_mod['Taxon'].str.strip()
taxa_mod['ModuleName'] = taxa_mod['Leiden_Cluster'].map({0: 'M1', 1: 'M2', 2: 'M3'})
mod_dict = dict(zip(taxa_mod['Taxon'], taxa_mod['ModuleName']))

# Prepare High-Contrast Module Sidebar Palette
mod_colors_map = {
    'M1': '#d62728', # Red
    'M2': '#1f77b4', # Blue
    'M3': '#2ca02c'  # Green
}

def generate_heatmap(clinical_subset, suffix, title_suffix):
    group_map = dict(zip(clinical_subset['SampleID'].astype(str), clinical_subset['group']))
    
    clr_df = clr_df_full.copy()
    clr_df['Group'] = clr_df['SampleID'].map(group_map)
    # Filter out samples not in the clinical_subset
    clr_df = clr_df.dropna(subset=['Group'])
    
    # 3. Create Heatmap Matrix
    actual_cols = [tax for tax in sig_taxa if tax in clr_df.columns]
    heatmap_df = clr_df.groupby("Group")[actual_cols].mean().T
    heatmap_df = heatmap_df[["C1", "C2", "C3"]]
    
    # Calculate group counts and rename columns
    group_counts = clr_df['Group'].value_counts()
    heatmap_df.columns = [
        f"C1\n(n={group_counts.get('C1', 0)})", 
        f"C2\n(n={group_counts.get('C2', 0)})", 
        f"C3\n(n={group_counts.get('C3', 0)})"
    ]
    
    row_colors = pd.Series(heatmap_df.index).map(mod_dict).map(mod_colors_map)
    row_colors.index = heatmap_df.index
    row_colors.name = "Module"
    
    # 5. Visualization
    sns.set_theme(style="white")
    vmax = np.abs(heatmap_df.values).max() * 0.9
    
    g = sns.clustermap(
        heatmap_df,
        cmap="RdBu_r",
        col_cluster=False, 
        row_cluster=True,  
        figsize=(8, 12),
        vmin=-vmax, vmax=vmax,
        row_colors=row_colors,
        cbar_pos=(0.02, 0.82, 0.04, 0.15),
        annot=True, fmt=".2f", 
        annot_kws={"size": 10, "weight": "bold"},
        linewidths=1.0,
        linecolor='white'
    )
    
    title = f"Mean CLR Abundance (Taxa)\n{title_suffix}" if title_suffix else "Mean CLR Abundance (Taxa)"
    g.ax_heatmap.set_title(title, fontsize=16, fontweight='bold', pad=15)
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.set_xlabel("")
    
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=0, fontsize=14, fontweight='bold')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=11, fontstyle='italic')
    
    # Adjust Colorbar label
    g.ax_cbar.set_title("CLR", fontsize=10, fontweight='bold', pad=10)
    
    # Add Legend for Modules
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=l) for l, c in mod_colors_map.items()]
    plt.legend(handles=legend_elements, title="Modules", bbox_to_anchor=(1.5, 1), loc='upper left', fontsize=12)
    
    out_path = OUTPUT_DIR / f"da_taxa_mean_clr_heatmap{suffix}.png"
    g.savefig(out_path, dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(g.fig)
    print(f"Heatmap saved to {out_path}")

print("\nGenerating Heatmap for All Samples...")
generate_heatmap(clinical, "", "(All Samples)")

print("\nGenerating Heatmap excluding Statin users...")
generate_heatmap(clinical[clinical['Statin'] == 0], "_no_statin", "(excluding Statin users)")

print("\nGenerating Heatmap excluding PPI users...")
generate_heatmap(clinical[clinical['PPI'] == 0], "_no_ppi", "(excluding PPI users)")

print("\nGenerating Heatmap excluding Antibiotics users...")
generate_heatmap(clinical[clinical['Antibiotics'] == 0], "_no_antibiotics", "(excluding Antibiotics users)")


