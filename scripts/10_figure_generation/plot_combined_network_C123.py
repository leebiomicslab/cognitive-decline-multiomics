import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Combined network visualisation for C1/C2/C3
Manuscript:       Results - Network comparison
Figure/Table:     Fig. 4d
"""

from pathlib import Path
"""
Combined Network Plot: C1, C2, C3 (Pearson r >= 0.3) ??one figure, shared legend.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import igraph as ig
import leidenalg as la
from pathlib import Path
import os

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42

# ?? Paths ??????????????????????????????????????????????????????????????????????
DA_RESULTS    = PROT_ASSOC_DIR / "Main_Model" / "proteomics_ancova_detailed_results.csv"
EXPRESSION_FILE = PROTEOMICS_LOG2
GROUP_FILE    = CLINICAL_CSV
OUT_DIR       = FASTSPAR_PROT_DIR / "pearson"
SAVE_PATH     = OUT_DIR / "Network_C1C2C3_pearson_r03_combined.png"

GROUPS        = ("C1", "C2", "C3")
R_THRESH      = 0.3
RESOLUTION    = 1.0
SEED          = 42
P_THRESH      = 0.05
FONT_SIZE     = 26   # node label font size

# ?? Protein name map ????????????????????????????????????????????????????????????
PROTEIN_RENAME = {
    'KIF2A&KIF2B': 'KIF2A',
    'DLG5&SLC49A3': 'DLG5',
    'TRAF3&Testis-specific serine/threonine-protein kinase 6': 'TRAF3',
    'Mannan-binding lectin serine protease 1': 'MASP1',
    'Mitogen-activated protein kinase kinase kinase 21': 'MAP3K21',
    'CFAP410&CPXM1&PHB2': 'PHB2',
    'Mannose-binding protein C': 'MBP-C',
    'FOXO1&GPC1&NEMF&PEX5L&SINHCAF&SPPL3': 'FOXO1',
    'DSG1&DSG4': 'DSG1',
    'Probable non-functional immunoglobulin lambda variable 1-50': 'IGLV1-50',
    'Probable non-functional immunoglobulin lamda variable 1-50': 'IGLV1-50',
    'NFX1-type zinc finger-containing protein 1': 'ZNFX1',
    'Putative uncharacterized protein encoded by LINC00469': 'LINC00469',
}

# ?? Load data ??????????????????????????????????????????????????????????????????
print("Loading data...")
da = pd.read_csv(DA_RESULTS)
da['Protein'] = da['Protein'].replace(PROTEIN_RENAME)
sig_proteins = set(da[da['P_Value'] < P_THRESH]['Protein'].tolist())
print(f"  Significant proteins: {len(sig_proteins)}")

expr_raw = pd.read_csv(EXPRESSION_FILE)
expr_raw.rename(columns={expr_raw.columns[0]: 'SampleID'}, inplace=True)
expr_raw.rename(columns=PROTEIN_RENAME, inplace=True)

clin = pd.read_csv(GROUP_FILE)
clin.rename(columns={'id': 'SampleID'}, inplace=True, errors='ignore')
if 'group' in clin.columns:
    clin.rename(columns={'group': 'Group'}, inplace=True)

valid_cols = [col for col in expr_raw.columns if col in sig_proteins]
print(f"  Matched proteins: {len(valid_cols)}")
df = expr_raw[['SampleID'] + valid_cols].copy()
df = df.merge(clin[['SampleID', 'Group']], on='SampleID', how='inner')
print(f"  Samples: {len(df)}")

# ?? Build networks ?????????????????????????????????????????????????????????????
def build_network(df, group):
    gdf = df[df['Group'] == group].drop(columns=['Group', 'SampleID'])
    gdf = gdf.loc[:, gdf.nunique() > 1]
    corr = gdf.corr(method='pearson', min_periods=3)
    G = nx.Graph()
    G.add_nodes_from(corr.columns)
    cols = corr.columns
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            r = corr.iloc[i, j]
            if np.isfinite(r):
                G.add_edge(cols[i], cols[j], weight=float(r))
    return G

print("Building networks...")
graphs = {g: build_network(df, g) for g in GROUPS}
G_ALL = nx.compose_all(graphs.values())

# ?? Leiden clustering on G_ALL ?????????????????????????????????????????????????
nodes = list(G_ALL.nodes())
node_idx = {n: i for i, n in enumerate(nodes)}
edges_ig = [(node_idx[u], node_idx[v]) for u, v in G_ALL.edges()]
weights_ig = [abs(d.get("weight", 0)) for _, _, d in G_ALL.edges(data=True)]

g_ig = ig.Graph()
g_ig.add_vertices(len(nodes))
g_ig.vs["name"] = nodes
g_ig.add_edges(edges_ig)
g_ig.es["weight"] = weights_ig

part = la.find_partition(g_ig, la.RBConfigurationVertexPartition,
                         weights=g_ig.es["weight"],
                         resolution_parameter=RESOLUTION, seed=SEED)
module_map = {}
for cid, comm in enumerate(part):
    for vid in comm:
        module_map[g_ig.vs[vid]['name']] = cid

n_mods = len(set(module_map.values()))
print(f"  Leiden modules: {n_mods}")

# ?? Layout (shared position from G_ALL) ????????????????????????????????????????
def abs_weighted_degree(G):
    deg = {n: 0.0 for n in G.nodes()}
    for u, v, d in G.edges(data=True):
        w = abs(d.get("weight", 0))
        deg[u] += w
        deg[v] += w
    return deg

def compute_positions(G, gap=0.02, radius=1.0):
    nodes = list(G.nodes())
    if not nodes:
        return {}
    idx = {n: i for i, n in enumerate(nodes)}
    edges, weights = [], []
    for u, v, d in G.edges(data=True):
        edges.append((idx[u], idx[v]))
        weights.append(abs(d.get("weight", 0)))
    gi = ig.Graph()
    gi.add_vertices(len(nodes))
    gi.vs["name"] = nodes
    gi.add_edges(edges)
    gi.es["weight"] = weights
    part = la.find_partition(gi, la.RBConfigurationVertexPartition,
                             weights=gi.es["weight"],
                             resolution_parameter=RESOLUTION, seed=SEED)
    communities = [[gi.vs[i]["name"] for i in comm] for comm in part]
    communities = sorted(communities, key=len, reverse=True)
    pos = {}
    total = sum(len(c) for c in communities)
    start = 0.0
    deg_abs = abs_weighted_degree(G)
    for comm in communities:
        length = len(comm) / total
        if length == 0:
            continue
        seg_start = start + gap / 2
        seg_end   = start + length - gap / 2
        order = sorted(comm, key=lambda n: deg_abs.get(n, 0), reverse=True)
        if len(order) == 1:
            thetas = np.array([(seg_start + seg_end) / 2]) * 2 * np.pi
        else:
            thetas = np.linspace(seg_start, seg_end, len(order), endpoint=False) * 2 * np.pi
        for n, t in zip(order, thetas):
            pos[n] = (radius * np.cos(t), radius * np.sin(t))
        start += length
    return pos

pos_all = compute_positions(G_ALL)

# ?? Global expression range ????????????????????????????????????????????????????
all_means = []
for g in GROUPS:
    gdf = df[df['Group'] == g].drop(columns=['Group', 'SampleID'])
    all_means.append(gdf.mean().values)
expr_vmin = float(np.nanmin(np.concatenate(all_means)))
expr_vmax = float(np.nanmax(np.concatenate(all_means)))

# ?? Module color setup ????????????????????????????????????????????????????????
unique_modules = sorted(set(module_map.values()))
cmap_border = plt.get_cmap('Set2' if len(unique_modules) <= 8 else 'tab20')
cmap_fill   = plt.cm.Reds
norm        = plt.Normalize(vmin=expr_vmin, vmax=expr_vmax)

# ?? Combined Figure ????????????????????????????????????????????????????????????
print("Plotting combined figure...")
fig, axes = plt.subplots(1, 3, figsize=(42, 15))
fig.patch.set_facecolor('white')

GROUP_COLORS = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'}
GROUP_LABELS = {'C1': 'C1 (MCI-free)', 'C2': 'C2 (Mild Cognitive Impairment)', 'C3': 'C3 (Moderate-Severe)'}

def draw_panel(ax, G, group, pos):
    gdf = df[df['Group'] == group].drop(columns=['Group', 'SampleID'])
    mean_expr = gdf.mean().to_dict()

    # Edges
    edges = [(u, v, d) for u, v, d in G.edges(data=True) if abs(d.get("weight", 0)) >= R_THRESH]
    edge_colors = ['#D55E00' if d['weight'] >= 0 else '#0072B2' for _, _, d in edges]
    edge_widths = [0.5 + 4 * (abs(d['weight']) - R_THRESH) / (1 - R_THRESH)
                   if (1 - R_THRESH) > 0 else 1.0 for _, _, d in edges]

    # Node sizes
    deg_abs = {n: 0 for n in G.nodes()}
    for u, v, d in G.edges(data=True):
        if abs(d.get("weight", 0)) >= R_THRESH:
            w = abs(d.get("weight", 0))
            deg_abs[u] += w
            deg_abs[v] += w
    max_deg = max(deg_abs.values()) if deg_abs and max(deg_abs.values()) > 0 else 1.0
    node_sizes = [500 + 1400 * (deg_abs.get(n, 0) / max_deg) for n in G.nodes()]

    # Node fill (expression)
    node_fill_colors = [cmap_fill(norm(mean_expr.get(n, 0))) for n in G.nodes()]

    # Node border (module)
    node_border_colors = [cmap_border(module_map.get(n, 0) % (8 if len(unique_modules) <= 8 else 20))
                          for n in G.nodes()]

    # Draw edges
    nx.draw_networkx_edges(G, pos, ax=ax,
                           edgelist=[(u, v) for u, v, _ in edges],
                           edge_color=edge_colors, width=edge_widths,
                           alpha=0.55, connectionstyle="arc3,rad=0.15",
                           arrows=True, arrowstyle='-')

    # Draw nodes (layered)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                           node_color='white', linewidths=0)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=[s * 0.75 for s in node_sizes],
                           node_color=node_fill_colors, linewidths=0)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                           node_color='none', edgecolors=node_border_colors, linewidths=3.5)

    # Labels
    nodes_list = list(G.nodes())
    for ni, node in enumerate(nodes_list):
        if node not in pos:
            continue
        x, y = pos[node]
        radius_offset = 1.05 + (node_sizes[ni] / 20000)
        text_x = x * radius_offset
        text_y = y * radius_offset
        ha = 'left' if x > 0.01 else ('right' if x < -0.01 else 'center')
        va = 'bottom' if y > 0.01 else ('top' if y < -0.01 else 'center')
        ax.text(text_x, text_y, node, fontsize=FONT_SIZE, fontweight='bold',
                fontfamily='Arial', ha=ha, va=va, color='#1a1a1a', zorder=10)

    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.65, 1.65)
    ax.axis('off')
    ax.set_title(f"{group}  (|r| ??{R_THRESH})",
                 fontsize=32, fontweight='bold', pad=14, color=GROUP_COLORS[group])

for i, g in enumerate(GROUPS):
    draw_panel(axes[i], graphs[g], g, pos_all)

# ?? Shared Legend ??????????????????????????????????????????????????????????????
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

legend_elements = [
    Line2D([0], [0], color='#D55E00', lw=5, label='Positive correlation'),
    Line2D([0], [0], color='#0072B2', lw=5, label='Negative correlation'),
]
for i in unique_modules:
    col = cmap_border(i % (8 if len(unique_modules) <= 8 else 20))
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor='w',
               markeredgecolor=col, markersize=22, markeredgewidth=4,
               label=f'Module P{i + 1}')
    )

fig.legend(handles=legend_elements, loc='upper right',
           ncol=1, fontsize=24,
           frameon=True, framealpha=0.9, edgecolor='#cccccc',
           bbox_to_anchor=(0.995, 0.97),
           title='Network Features', title_fontproperties={'weight': 'bold', 'size': 26})

# Shared colorbar
sm = plt.cm.ScalarMappable(cmap=cmap_fill, norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.65])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.outline.set_visible(False)
cbar.ax.tick_params(labelsize=22, width=0, length=0)
cbar.set_label("Average Expression (Log2)", size=24, weight='bold', labelpad=18)

fig.suptitle("Proteomics Co-expression Networks ??Pearson Correlation (|r| ??0.3)",
             fontsize=36, fontweight='bold', y=1.01)

plt.tight_layout(rect=[0, 0.04, 0.91, 1.0])
plt.savefig(SAVE_PATH, dpi=600, bbox_inches='tight', facecolor='white')
print(f"\n[OK] Saved -> {SAVE_PATH}")
plt.close(fig)


