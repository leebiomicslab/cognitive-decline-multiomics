import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Network topology comparison across C1/C2/C3
Manuscript:       Results - Network topology
Figure/Table:     Fig. 4d
"""


import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

colors = {"C1": "#ee7a5b", "C2": "#eca362", "C3": "#6fcbbc"}

def plot_network_stats(files, out_dir, prefix):
    edge_weights = []
    edge_counts = []
    groups = []
    
    for group in ["C1", "C2", "C3"]:
        path = files.get(group)
        if not path or not os.path.exists(path):
            print(f"File not found: {path}")
            continue
            
        G = nx.read_graphml(path)
        weights = [float(data.get('weight', 0)) for u, v, data in G.edges(data=True)]
        
        edge_weights.extend(weights)
        groups.extend([group] * len(weights))
        edge_counts.append({'Group': group, 'Count': len(weights)})
        
    if not edge_weights:
        print(f"No valid data found for {prefix}")
        return
        
    df_weights = pd.DataFrame({'Weight': edge_weights, 'Group': groups})
    df_counts = pd.DataFrame(edge_counts)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Bar plot for Edge Counts
    sns.barplot(data=df_counts, x='Group', y='Count', palette=colors, ax=axes[0])
    axes[0].set_title(f"{prefix} Edge Counts", fontweight='bold', fontsize=14)
    axes[0].set_ylabel("Number of Edges", fontsize=12)
    axes[0].set_xlabel("")
    for i, p in enumerate(axes[0].patches):
        axes[0].annotate(f"{int(p.get_height())}", (p.get_x() + p.get_width() / 2., p.get_height()), 
                         ha='center', va='bottom', fontsize=11, fontweight='bold')
                         
    # Box plot for Edge Weights
    sns.boxplot(data=df_weights, x='Group', y='Weight', palette=colors, ax=axes[1])
    axes[1].set_title(f"{prefix} Edge Weights", fontweight='bold', fontsize=14)
    axes[1].set_ylabel("Edge Weight", fontsize=12)
    axes[1].set_xlabel("")
    
    plt.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{prefix}_network_comparison.png")
    fig.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_path}")

# Microbiota
microbiota_files = {
    "C1": str(FASTSPAR_MICRO_DIR / "microbiome_G_C1.graphml"),
    "C2": str(FASTSPAR_MICRO_DIR / "microbiome_G_C2.graphml"),
    "C3": str(FASTSPAR_MICRO_DIR / "microbiome_G_C3.graphml")
}
plot_network_stats(microbiota_files, str(FIGURES_DIR / "microbiota"), "Microbiota")

# Proteomics
proteomics_files = {
    "C1": str(FASTSPAR_PROT_DIR / "pearson" / "GraphML_export" / "Proteomics_C1_pearson.graphml"),
    "C2": str(FASTSPAR_PROT_DIR / "pearson" / "GraphML_export" / "Proteomics_C2_pearson.graphml"),
    "C3": str(FASTSPAR_PROT_DIR / "pearson" / "GraphML_export" / "Proteomics_C3_pearson.graphml")
}
plot_network_stats(proteomics_files, str(FIGURES_DIR / "proteomics"), "Proteomics")

