import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Correlation of module eigengenes with MoCA scores
Manuscript:       Results - Module-cognition link
Figure/Table:     Fig. 5c
"""

from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def geometric_mean(series, epsilon=1e-6):
    """Calculates the geometric mean of a series with pseudo-count."""
    return np.exp(np.log(series + epsilon).mean())

def plot_me_moca_correlation():
    # 1. Setup paths
    vascular_path = base_dir / "data" / "vascular.csv"
    micro_abundance_path = MICRO_GENUS_COUNTS
    leiden_list_path = FASTSPAR_MICRO_DIR / "network_topology" / "leiden_clusters_24_taxa_G_all.csv"
    prot_me_path = FASTSPAR_PROT_DIR / "pearson" / "ME_proteomics_ALL_pearson.csv"
    output_dir = RESULTS_DIR / "fastspar" / "module_visualizations"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data...")
    # 2. Load MoCA scores
    df_vascular = pd.read_csv(vascular_path)
    df_moca = df_vascular[['id', 'moca_1']].dropna()
    df_moca.columns = ['id', 'MoCA']

    # 3. Calculate Microbiota Geometric Means (M1, M2, M3)
    taxa_info = pd.read_csv(leiden_list_path)
    abund = pd.read_csv(micro_abundance_path).set_index("Genus").T
    abund.index = abund.index.astype(str)

    module_map = {0: 'M1', 1: 'M2', 2: 'M3'}
    taxa_info['Module'] = taxa_info['Leiden_Cluster'].map(module_map)

    gm_results = pd.DataFrame(index=abund.index)
    for mod_name in ['M1', 'M2', 'M3']:
        mod_taxa = taxa_info[taxa_info['Module'] == mod_name]['Taxon'].tolist()
        valid_taxa = [t for t in mod_taxa if t in abund.columns]
        if valid_taxa:
            gm_results[mod_name] = abund[valid_taxa].apply(lambda x: geometric_mean(x, 1e-6), axis=1)
    
    gm_results = gm_results.reset_index().rename(columns={'index': 'id'})
    micro_cols = ['M1', 'M2', 'M3']

    # 4. Load Proteomics MEs
    df_prot = pd.read_csv(prot_me_path)
    df_prot = df_prot.rename(columns={
        'SampleID': 'id',
        'M0_prot': 'P1',
        'M1_prot': 'P2',
        'M2_prot': 'P3',
        'M3_prot': 'P4',
        'M4_prot': 'P5'
    })
    prot_cols = ['P1', 'P2', 'P3', 'P4', 'P5']
    df_prot = df_prot[['id'] + prot_cols]

    # 5. Merge
    df_merged = df_moca.merge(gm_results, on='id').merge(df_prot, on='id')
    print(f"Number of common samples: {len(df_merged)}")

    # 6. Calculate Correlations
    me_cols = micro_cols + prot_cols
    results = []
    
    for col in me_cols:
        # Spearman
        r_s, p_s = spearmanr(df_merged[col], df_merged['MoCA'])
        # Pearson
        r_p, p_p = pearsonr(df_merged[col], df_merged['MoCA'])
        
        results.append({
            'Module': col,
            'Spearman_r': r_s,
            'Spearman_p': p_s,
            'Pearson_r': r_p,
            'Pearson_p': p_p,
            'Type': 'Microbiota (GM)' if col.startswith('M') else 'Proteomics (ME)'
        })
    
    df_results = pd.DataFrame(results)
    
    # 7. Visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=600, sharey=True)
    sns.set_style("whitegrid")
    
    colors = {'Microbiota (GM)': '#6fcbbc', 'Proteomics (ME)': '#ee7a5b'}
    
    metrics = [('Pearson_r', 'Pearson_p', 'Pearson Correlation'), 
               ('Spearman_r', 'Spearman_p', 'Spearman Correlation')]
    
    # Consistent sorting by Spearman
    plot_order = df_results.sort_values('Spearman_r', ascending=False)['Module'].tolist()

    for i, (r_col, p_col, title) in enumerate(metrics):
        ax = axes[i]
        sns.barplot(data=df_results, x=r_col, y='Module', hue='Type', 
                    palette=colors, dodge=False, ax=ax, order=plot_order)
        
        # Add p-value annotations
        for j, mod in enumerate(plot_order):
            row = df_results[df_results['Module'] == mod].iloc[0]
            p = row[p_col]
            r_val = row[r_col]
            
            if p < 0.001: annot = "***"
            elif p < 0.01: annot = "**"
            elif p < 0.05: annot = "*"
            else: annot = ""
            
            if annot:
                x_pos = r_val + (0.01 if r_val > 0 else -0.04)
                ax.text(x_pos, j, annot, va='center', fontweight='bold', fontsize=12)

        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel("Correlation Coefficient", fontsize=12)
        ax.set_ylabel("Module" if i==0 else "", fontsize=12)
        ax.axvline(0, color='black', linewidth=1)
        ax.set_xlim(-0.25, 0.25)
        
        if i > 0: ax.get_legend().remove()

    plt.suptitle("Microbiota (Geometric Mean) & Proteomics (ME) vs. MoCA", fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "me_gm_moca_correlation.png", bbox_inches='tight')
    df_results.to_csv(output_dir / "me_gm_moca_correlation_stats.csv", index=False)
    
    print(f"Results saved to {output_dir}")
    print(df_results)

if __name__ == "__main__":
    plot_me_moca_correlation()



if __name__ == "__main__":
    plot_me_moca_correlation()


