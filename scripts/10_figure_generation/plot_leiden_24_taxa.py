import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Leiden microbiome modules visualisation (24 DA taxa)
Manuscript:       Methods - Microbial network modules
Figure/Table:     Fig. 4b
"""

from pathlib import Path
import pandas as pd
import numpy as np
import networkx as nx
import igraph as ig
import leidenalg as la
import matplotlib.pyplot as plt
from pathlib import Path
import os
import matplotlib.colors as mcolors

# Paths
DA_RESULTS = MICRO_ASSOC_DIR / "Model_A_Main" / "microbiota_da_results_final.csv"
GRAPH_ALL = FASTSPAR_MICRO_DIR / "microbiome_G_all.graphml"
OUTPUT_DIR = FASTSPAR_MICRO_DIR / "network_topology"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# 1. Identify 24 DA Taxa
da = pd.read_csv(DA_RESULTS)
sig_taxa = da[da["Union_Sig"] == True]["taxon"].str.strip().tolist()
valid_taxa = set(sig_taxa)

# 2. Load ALL Network
G_all = nx.read_graphml(GRAPH_ALL)
nodes_to_keep = [n for n in G_all.nodes() if n.strip() in valid_taxa]
subG = G_all.subgraph(nodes_to_keep).copy()

print(f"Extracted subgraph with {subG.number_of_nodes()} nodes and {subG.number_of_edges()} edges.")

# 3. Convert to igraph for Leiden
nodes = list(subG.nodes())
node_idx = {n: i for i, n in enumerate(nodes)}

edges = []
weights = []
for u, v, d in subG.edges(data=True):
    w = float(d.get("weight", 0))
    # Leiden operates effectively on absolute weights for finding dense correlated communities
    abs_w = abs(w)
    if abs_w > 0: 
        edges.append((node_idx[u], node_idx[v]))
        weights.append(abs_w)

g_ig = ig.Graph()
g_ig.add_vertices(len(nodes))
g_ig.vs["name"] = nodes
g_ig.add_edges(edges)
g_ig.es["weight"] = weights

# 4. Run Leiden Alg (resolution=1.0 is standard)
RESOLUTION = 1.0
SEED = 42

part = la.find_partition(
    g_ig, la.RBConfigurationVertexPartition,
    weights=g_ig.es["weight"], resolution_parameter=RESOLUTION, seed=SEED
)

leiden_clusters = {}
cluster_sizes = {}
for cid, comm in enumerate(part):
    taxa_in_comm = [g_ig.vs[vid]['name'] for vid in comm]
    cluster_sizes[cid] = len(taxa_in_comm)
    for tax in taxa_in_comm:
        leiden_clusters[tax] = cid

print(f"\nLeiden algorithm detected {len(cluster_sizes)} distinctive clusters within the 24 DA taxa.")

df_cluster = pd.DataFrame([{"Taxon": k, "Leiden_Cluster": v} for k, v in leiden_clusters.items()])
df_cluster.to_csv(OUTPUT_DIR / "leiden_clusters_24_taxa_G_all.csv", index=False)

# 5. Plot the Network
# Edges for visualization threshold
DRAW_THRESH = 0.2
edges_to_draw = []
for u, v, d in subG.edges(data=True):
    if abs(float(d.get("weight", 0))) >= DRAW_THRESH:
        edges_to_draw.append((u, v))

# Create drawing graph
G_draw = nx.Graph()
G_draw.add_nodes_from(nodes)
G_draw.add_edges_from(edges_to_draw)

# Calculate sizes
deg_abs = {n: 0.0 for n in G_draw.nodes()}
for u, v in G_draw.edges():
    w = abs(float(subG[u][v].get("weight", 0)))
    deg_abs[u] += w
    deg_abs[v] += w
max_deg = max(deg_abs.values()) if deg_abs and max(deg_abs.values()) > 0 else 1.0
node_sizes = [300 + 1500 * (deg_abs[n] / max_deg) for n in nodes]

num_clusters = len(cluster_sizes)
cmap = plt.get_cmap('Set2') if num_clusters <= 8 else plt.get_cmap('tab20')
node_colors = [cmap(leiden_clusters[n]) for n in nodes]

edge_colors = []
edge_widths = []
for u, v in G_draw.edges():
    w = float(subG[u][v]["weight"])
    edge_colors.append('#E41A1C' if w >= 0 else '#377EB8')
    w_abs = abs(w)
    norm_w = 1.0 + 3.0 * (w_abs - DRAW_THRESH) / (1.0 - DRAW_THRESH) if (1.0 - DRAW_THRESH) > 0 else 1.0
    edge_widths.append(norm_w)

# Position calculation
# Using spring layout but highly clustering nodes naturally via weights
pos = nx.spring_layout(G_draw, k=0.9, seed=SEED)

fig, ax = plt.subplots(figsize=(10, 10))

nx.draw_networkx_nodes(G_draw, pos, ax=ax, node_size=node_sizes, node_color=node_colors, edgecolors='black', linewidths=1.5)
nx.draw_networkx_edges(G_draw, pos, ax=ax, edge_color=edge_colors, width=edge_widths, alpha=0.6)
nx.draw_networkx_labels(G_draw, pos, ax=ax, font_size=10, font_weight="bold", font_family="sans-serif")

# Legend
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(cid), markeredgecolor='black', 
                          markersize=12, label=f'Cluster {cid} (n={cluster_sizes[cid]})') 
                   for cid in sorted(cluster_sizes.keys())]
# Add edge color legend
legend_elements.append(Line2D([0], [0], color='#E41A1C', lw=3, label='Positive Correl'))
legend_elements.append(Line2D([0], [0], color='#377EB8', lw=3, label='Negative Correl'))

ax.legend(handles=legend_elements, loc='upper left', title="Leiden Modules & Edges", fontsize=11, title_fontsize=13)
ax.set_title(f"Leiden Clustering of 24 DA Taxa\n(Trained on ALL network, plotted edges |r| >= {DRAW_THRESH})", fontsize=16, fontweight='bold', pad=15)
ax.axis("off")

out_png = OUTPUT_DIR / "leiden_network_24_taxa_G_all.png"
plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"\nNetwork plot successfully saved to {out_png}")

with open(OUTPUT_DIR / "leiden_clusters_summary.txt", "w", encoding="utf-8") as f:
    summary_hdr = "Leiden Clustering Results for the 24 DA Taxa (ALL Network):\n"
    print(summary_hdr)
    f.write(summary_hdr)
    for cid in sorted(cluster_sizes.keys()):
        taxa = [n for n in nodes if leiden_clusters[n] == cid]
        summary = f"- Cluster {cid} (n={len(taxa)}):\n    " + ", ".join(taxa) + "\n"
        print(summary)
        f.write(summary + "\n")


