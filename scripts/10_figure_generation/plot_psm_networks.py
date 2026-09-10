import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Microbial network topology after PSM
Manuscript:       Methods - PSM network sensitivity
Figure/Table:     Supplementary Fig. 8a
"""

from pathlib import Path
"""
Plot FastSpar Microbiome PSM Networks (C1, C2, C3)
Generates both individual high-res network plots and a combined 1x3 panel with shared legend.
"""
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
from pathlib import Path
import os

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42

# ?? Paths ??????????????????????????????????????????????????????????????????????
DA_RESULTS   = MICRO_ASSOC_DIR / "Model_A_Main" / "microbiota_da_results_final.csv"
GRAPH_DIR    = FASTSPAR_MICRO_DIR
CLUSTER_FILE = FASTSPAR_MICRO_DIR / "network_topology" / "leiden_clusters_24_taxa_G_all.csv"
CLINICAL_FILE= CLINICAL_CSV
ASV_FILE     = MICRO_GENUS_COUNTS

OUT_DIR      = GRAPH_DIR / "network_psm_plots"
os.makedirs(OUT_DIR, exist_ok=True)

GROUPS       = ("C1", "C2", "C3")
THRESH       = 0.1
P_VALUE_THRESH = 0.05
FONT_SIZE    = 24
GROUP_COLORS = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'}

# ?? Helper for CLR ?????????????????????????????????????????????????????????????
def clr_transformation(df, epsilon=None):
    if epsilon is None:
        values = df.values.flatten()
        min_pos = np.min(values[values > 0])
        epsilon = min_pos * 0.65
    df_imp = df.replace(0, epsilon)
    log_df = np.log(df_imp)
    clr_df = log_df.sub(log_df.mean(axis=1), axis=0)
    return clr_df

def clean_taxon(name):
    return name.replace('"', '').strip()

# ?? Load CLR data ??????????????????????????????????????????????????????????????
print("Loading clinical & expression data for CLR calculation...")
clinical = pd.read_csv(CLINICAL_FILE)
clinical.rename(columns={'id': 'SampleID'}, inplace=True)
group_map = dict(zip(clinical['SampleID'].astype(str), clinical['group']))

expr_raw = pd.read_csv(ASV_FILE)
genus_col = expr_raw.columns[0]
expr_raw.set_index(genus_col, inplace=True)
expr_raw = expr_raw.loc[:, ~expr_raw.columns.duplicated()]

expr = expr_raw.T
clr_df = clr_transformation(expr.astype(float))
clr_df['SampleID'] = clr_df.index.astype(str)
clr_df['Group'] = clr_df['SampleID'].map(group_map)

mean_clr_grp = {}
for g in GROUPS:
    mean_clr_grp[g] = clr_df[clr_df['Group'] == g].drop(columns=['SampleID', 'Group']).mean().to_dict()

# ?? Load Significant Taxa and Clusters ?????????????????????????????????????????
da = pd.read_csv(DA_RESULTS)
sig_taxa = set([clean_taxon(t) for t in da[da["Union_Sig"] == True]["taxon"].tolist()])

leiden_df = pd.read_csv(CLUSTER_FILE)
leiden_map = dict(zip(leiden_df['Taxon'].str.replace('"', '').str.strip(), leiden_df['Leiden_Cluster']))

# ?? Load PSM Graphs ????????????????????????????????????????????????????????????
raw_graphs = {}
union_nodes = set()
for g in GROUPS:
    gml_path = GRAPH_DIR / f"microbiome_G_{g}_psm.graphml"
    if gml_path.exists():
        G = nx.read_graphml(gml_path)
        # Filter nodes to significant taxa
        keep = [n for n in G.nodes() if clean_taxon(n) in sig_taxa]
        subG = G.subgraph(keep).copy()
        
        edges_to_keep = []
        for u, v, d in subG.edges(data=True):
            w = float(d.get("weight", 0))
            p = float(d.get("p_value", 1.0))
            if abs(w) >= THRESH and p < P_VALUE_THRESH:
                edges_to_keep.append((u, v, {'weight': w}))
        
        fG = nx.Graph()
        fG.add_nodes_from(keep)
        fG.add_edges_from(edges_to_keep)
        raw_graphs[g] = fG
        union_nodes.update([clean_taxon(n) for n in keep])
    else:
        print(f"Warning: {gml_path.name} not found.")

print(f"Total valid PSM nodes across groups: {len(union_nodes)}")

# Compute CLR range for color mapping
all_means_sig = []
for g in GROUPS:
    for tax in union_nodes:
        if tax in mean_clr_grp[g]:
            all_means_sig.append(mean_clr_grp[g][tax])

if all_means_sig:
    vmin_clr, vmax_clr = min(all_means_sig), max(all_means_sig)
else:
    vmin_clr, vmax_clr = -2.0, 2.0

# ?? Sorted Circle Layout ???????????????????????????????????????????????????????
def make_sorted_circle_layout(taxa, leiden_map):
    from collections import defaultdict
    cluster_members = defaultdict(list)
    for t in sorted(taxa):
        cluster_members[leiden_map.get(t, 0)].append(t)
    ordered = []
    for cid in sorted(cluster_members):
        ordered.extend(cluster_members[cid])
    n = len(ordered)
    pos = {}
    for i, t in enumerate(ordered):
        angle = 2 * np.pi * i / n
        pos[t] = np.array([np.cos(angle), np.sin(angle)])
    return pos

pos_all = make_sorted_circle_layout(list(union_nodes), leiden_map)

# ?? Color Palettes ?????????????????????????????????????????????????????????????
norm = mcolors.Normalize(vmin=vmin_clr, vmax=vmax_clr)
cmap_fill = plt.cm.RdBu_r
unique_clusters = sorted(set(leiden_map.values()))
cmap_cluster = plt.get_cmap('Set2' if len(unique_clusters) <= 8 else 'tab20')

# ?? Drawing Function for Panel ?????????????????????????????????????????????????
def draw_panel(ax, G, group, pos):
    edges = G.edges(data=True)
    deg_abs = {n: 0.0 for n in G.nodes()}
    for u, v, d in edges:
        w = abs(float(d.get('weight', 0)))
        deg_abs[u] += w
        deg_abs[v] += w
        
    max_deg = max(deg_abs.values()) if deg_abs and max(deg_abs.values()) > 0 else 1.0
    node_sizes = [500 + 1400 * (deg_abs.get(n, 0) / max_deg) for n in G.nodes()]
    
    node_colors = [cmap_fill(norm(mean_clr_grp[group].get(clean_taxon(n), 0))) for n in G.nodes()]
    node_border_colors = [cmap_cluster(leiden_map.get(clean_taxon(n), 0) % (8 if len(unique_clusters)<=8 else 20)) for n in G.nodes()]
    
    edge_colors = ['#D55E00' if float(d.get('weight', 0)) >= 0 else '#0072B2' for _, _, d in edges]
    edge_widths = [0.5 + 4.0 * (abs(float(d.get('weight', 0))) - THRESH) / (1.0 - THRESH)
                   if (1.0 - THRESH) > 0 else 1.0 for _, _, d in edges]
                   
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(u, v) for u, v, _ in edges],
                           edge_color=edge_colors, width=edge_widths, alpha=0.55,
                           connectionstyle="arc3,rad=0.15", arrows=True, arrowstyle='-')
                           
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color='white', linewidths=0)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=[s * 0.75 for s in node_sizes],
                           node_color=node_colors, linewidths=0)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color='none',
                           edgecolors=node_border_colors, linewidths=3.5)
                           
    nodes_list = list(G.nodes())
    for ni, node in enumerate(nodes_list):
        if clean_taxon(node) not in pos: continue
        x, y = pos[clean_taxon(node)]
        radius_offset = 1.05 + (node_sizes[ni] / 20000)
        text_x = x * radius_offset
        text_y = y * radius_offset
        ha = 'left' if x > 0.01 else ('right' if x < -0.01 else 'center')
        va = 'bottom' if y > 0.01 else ('top' if y < -0.01 else 'center')
        ax.text(text_x, text_y, clean_taxon(node), fontsize=FONT_SIZE, fontweight='bold',
                fontfamily='Arial', ha=ha, va=va, color='#1a1a1a', zorder=10)
                
    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.65, 1.65)
    ax.axis('off')
    ax.set_title(f"PSM - {group}  (|r| ??{THRESH}, p < {P_VALUE_THRESH})",
                 fontsize=32, fontweight='bold', pad=14, color=GROUP_COLORS[group])

# ?? 1. Create Combined 1x3 Figure ??????????????????????????????????????????????
print("\nPlotting combined 1x3 PSM network figure...")
fig, axes = plt.subplots(1, 3, figsize=(42, 15))
fig.patch.set_facecolor('white')

for i, g in enumerate(GROUPS):
    if g in raw_graphs:
        draw_panel(axes[i], raw_graphs[g], g, pos_all)

# Shared Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='#D55E00', lw=5, label='Positive correlation'),
    Line2D([0], [0], color='#0072B2', lw=5, label='Negative correlation'),
]
for i in unique_clusters:
    col = cmap_cluster(i % (8 if len(unique_clusters)<=8 else 20))
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor='w',
               markeredgecolor=col, markersize=22, markeredgewidth=4,
               label=f'Module M{i + 1}')
    )

fig.legend(handles=legend_elements, loc='upper right', ncol=1, fontsize=24,
           frameon=True, framealpha=0.9, edgecolor='#cccccc', bbox_to_anchor=(0.995, 0.97),
           title='Network Features', title_fontproperties={'weight': 'bold', 'size': 26})

# Shared Colorbar
sm = plt.cm.ScalarMappable(cmap=cmap_fill, norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.65])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.outline.set_visible(False)
cbar.ax.tick_params(labelsize=22, width=0, length=0)
cbar.set_label("Mean CLR Abundance", size=24, weight='bold', labelpad=18)

fig.suptitle(f"Microbiome PSM Co-abundance Networks (|r| ??{THRESH}, p < {P_VALUE_THRESH})",
             fontsize=36, fontweight='bold', y=1.01)

plt.tight_layout(rect=[0, 0.04, 0.91, 1.0])
combo_path = OUT_DIR / f"Network_C1C2C3_psm_thresh{THRESH}_combined.png"
plt.savefig(combo_path, dpi=600, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f"[OK] Combined figure saved -> {combo_path}")

# ?? 2. Create Individual Figures (just in case) ????????????????????????????????
print("\nPlotting individual PSM network figures...")
for g in GROUPS:
    if g not in raw_graphs: continue
    fig, ax = plt.subplots(figsize=(16, 14))
    fig.patch.set_facecolor('white')
    draw_panel(ax, raw_graphs[g], g, pos_all)
    
    # Legend
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.25, 1.0),
              fontsize=16, frameon=False, title="Network Features", title_fontproperties={'weight':'bold', 'size':18})
              
    # Colorbar
    cbar_ax_ind = fig.add_axes([0.95, 0.25, 0.02, 0.5])
    cbar_ind = fig.colorbar(sm, cax=cbar_ax_ind)
    cbar_ind.outline.set_visible(False)
    cbar_ind.ax.tick_params(labelsize=14, width=0, length=0)
    cbar_ind.set_label("Mean CLR Abundance", size=18, weight='bold', labelpad=15)
    
    ind_path = OUT_DIR / f"Network_{g}_psm_thresh{THRESH}.png"
    plt.savefig(ind_path, dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  [OK] Saved -> {ind_path.name}")

print("\n??All PSM network visualizations completed!")


