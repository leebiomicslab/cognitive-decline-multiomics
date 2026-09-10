import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Protein network topology after PSM
Manuscript:       Methods - PSM network sensitivity
Figure/Table:     Supplementary Fig. 8b
"""

from pathlib import Path
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from pathlib import Path

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42

# ?? Paths ??????????????????????????????????????????????????????????????????????
DA_RESULTS = PROT_ASSOC_DIR / "Main_Model" / "proteomics_ancova_detailed_results.csv"
EXPRESSION_FILE = PROTEOMICS_LOG2
PSM_GROUP_FILE = PSM_MATCHED_CSV
CLUSTER_FILE = FASTSPAR_PROT_DIR / "pearson" / "module_protein_assignment.csv"

OUT_DIR = FASTSPAR_PROT_DIR / "network_psm_plots"
os.makedirs(OUT_DIR, exist_ok=True)

GROUPS = ("C1", "C2", "C3")
R_THRESH = 0.3
P_VALUE_THRESH = 0.05
FONT_SIZE = 24
GROUP_COLORS = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'}

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

# ?? 1. Load Data ??????????????????????????????????????????????????????????????
print("Loading clinical PSM & expression data...")
clin_psm = pd.read_csv(PSM_GROUP_FILE)
clin_psm.rename(columns={'id': 'SampleID'}, inplace=True, errors='ignore')
if 'SampleID' not in clin_psm.columns:
    clin_psm.rename(columns={clin_psm.columns[0]: 'SampleID'}, inplace=True)
if 'group' in clin_psm.columns:
    clin_psm.rename(columns={'group': 'Group'}, inplace=True)

group_map = dict(zip(clin_psm['SampleID'].astype(str), clin_psm['Group']))

da = pd.read_csv(DA_RESULTS)
da['Protein'] = da['Protein'].replace(PROTEIN_RENAME)
sig_proteins = set(da[da['P_Value'] < P_VALUE_THRESH]['Protein'].tolist())

expr_raw = pd.read_csv(EXPRESSION_FILE)
expr_raw.rename(columns={expr_raw.columns[0]: 'SampleID'}, inplace=True)
expr_raw.rename(columns=PROTEIN_RENAME, inplace=True)
expr_raw['SampleID'] = expr_raw['SampleID'].astype(str)

# Filter expression to only PSM samples
expr_psm = expr_raw[expr_raw['SampleID'].isin(group_map.keys())].copy()
expr_psm['Group'] = expr_psm['SampleID'].map(group_map)

# Keep only significant proteins
valid_cols = [col for col in expr_psm.columns if col in sig_proteins]
df = expr_psm[['SampleID', 'Group'] + valid_cols].copy()

# Load Cluster Map (from Main Model)
mod_assign_df = pd.read_csv(CLUSTER_FILE)
module_map = dict(zip(mod_assign_df['Protein'], mod_assign_df['Module']))
unique_modules = sorted(set(module_map.values()))

# ?? 2. Build PSM Pearson Networks ??????????????????????????????????????????????
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

print("Building PSM networks...")
raw_graphs = {}
union_nodes = set()
for g in GROUPS:
    G = build_network(df, g)
    
    # Filter by threshold
    edges_to_keep = []
    for u, v, d in G.edges(data=True):
        w = float(d.get("weight", 0))
        if abs(w) >= R_THRESH:
            edges_to_keep.append((u, v, {'weight': w}))
            
    fG = nx.Graph()
    fG.add_nodes_from(G.nodes())
    fG.add_edges_from(edges_to_keep)
    
    raw_graphs[g] = fG
    union_nodes.update(fG.nodes())

# Calculate Expression ranges
all_means = []
mean_expr_grp = {}
for g in GROUPS:
    gdf = df[df['Group'] == g].drop(columns=['Group', 'SampleID'])
    mean_expr = gdf.mean()
    mean_expr_grp[g] = mean_expr.to_dict()
    all_means.append(mean_expr.values)

expr_vmin = float(np.nanmin(np.concatenate(all_means)))
expr_vmax = float(np.nanmax(np.concatenate(all_means)))

# ?? 3. Sorted Circular Layout (Grouped by Module) ?????????????????????????????
def make_sorted_circle_layout(taxa, mod_map):
    from collections import defaultdict
    cluster_members = defaultdict(list)
    for t in sorted(taxa):
        cluster_members[mod_map.get(t, 0)].append(t)
    ordered = []
    for cid in sorted(cluster_members):
        ordered.extend(cluster_members[cid])
    n = len(ordered)
    pos = {}
    for i, t in enumerate(ordered):
        angle = 2 * np.pi * i / n
        pos[t] = np.array([np.cos(angle), np.sin(angle)])
    return pos

pos_all = make_sorted_circle_layout(list(union_nodes), module_map)

# ?? 4. Drawing Functions ???????????????????????????????????????????????????????
cmap_border = plt.get_cmap('Set2' if len(unique_modules) <= 8 else 'tab20')
cmap_fill = plt.cm.Reds
norm = plt.Normalize(vmin=expr_vmin, vmax=expr_vmax)

def draw_panel(ax, G, group, pos):
    edges = [(u, v, d) for u, v, d in G.edges(data=True)]
    
    deg_abs = {n: 0.0 for n in G.nodes()}
    for u, v, d in edges:
        w = abs(float(d.get('weight', 0)))
        deg_abs[u] += w
        deg_abs[v] += w
        
    max_deg = max(deg_abs.values()) if deg_abs and max(deg_abs.values()) > 0 else 1.0
    node_sizes = [500 + 1400 * (deg_abs.get(n, 0) / max_deg) for n in G.nodes()]
    
    node_colors = [cmap_fill(norm(mean_expr_grp[group].get(n, 0))) for n in G.nodes()]
    node_border_colors = [cmap_border(module_map.get(n, 0) % (8 if len(unique_modules)<=8 else 20)) for n in G.nodes()]
    
    edge_colors = ['#D55E00' if float(d.get('weight', 0)) >= 0 else '#0072B2' for _, _, d in edges]
    edge_widths = [0.5 + 4.0 * (abs(float(d.get('weight', 0))) - R_THRESH) / (1.0 - R_THRESH)
                   if (1.0 - R_THRESH) > 0 else 1.0 for _, _, d in edges]
                   
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(u, v) for u, v, _ in edges],
                           edge_color=edge_colors, width=edge_widths, alpha=0.55,
                           connectionstyle="arc3,rad=0.15", arrows=True, arrowstyle='-')
                           
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color='white', linewidths=0)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=[s * 0.75 for s in node_sizes],
                           node_color=node_colors, linewidths=0)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color='none',
                           edgecolors=node_border_colors, linewidths=3.5)
                           
    def format_label(name):
        if '&' in name and len(name) > 20:
            parts = name.split('&')
            mid = len(parts) // 2 + (len(parts) % 2)
            return '&'.join(parts[:mid]) + '&\n' + '&'.join(parts[mid:])
        return name

    nodes_list = list(G.nodes())
    for ni, node in enumerate(nodes_list):
        if node not in pos: continue
        x, y = pos[node]
        radius_offset = 1.05 + (node_sizes[ni] / 20000)
        text_x = x * radius_offset
        text_y = y * radius_offset
        ha = 'left' if x > 0.01 else ('right' if x < -0.01 else 'center')
        va = 'bottom' if y > 0.01 else ('top' if y < -0.01 else 'center')
        ax.text(text_x, text_y, format_label(node), fontsize=FONT_SIZE, fontweight='bold',
                fontfamily='Arial', ha=ha, va=va, color='#1a1a1a', zorder=10)
                
    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.65, 1.65)
    ax.axis('off')
    ax.set_title(f"PSM - {group}  (|r| ??{R_THRESH})",
                 fontsize=32, fontweight='bold', pad=14, color=GROUP_COLORS[group])

# ?? 5. Create Combined 1x3 Figure ??????????????????????????????????????????????
print("\nPlotting combined 1x3 PSM network figure...")
fig, axes = plt.subplots(1, 3, figsize=(42, 15))
fig.patch.set_facecolor('white')

for i, g in enumerate(GROUPS):
    draw_panel(axes[i], raw_graphs[g], g, pos_all)

# Shared Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='#D55E00', lw=5, label='Positive correlation'),
    Line2D([0], [0], color='#0072B2', lw=5, label='Negative correlation'),
]
for i in unique_modules:
    col = cmap_border(i % (8 if len(unique_modules)<=8 else 20))
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor='w',
               markeredgecolor=col, markersize=22, markeredgewidth=4,
               label=f'Module P{i + 1}')
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
cbar.set_label("Average Expression (Log2)", size=24, weight='bold', labelpad=18)

fig.suptitle(f"Proteomics PSM Co-expression Networks (|r| ??{R_THRESH})",
             fontsize=36, fontweight='bold', y=1.01)

plt.tight_layout(rect=[0, 0.04, 0.91, 1.0])
combo_path = OUT_DIR / f"Network_C1C2C3_psm_pearson_r03_combined.png"
plt.savefig(combo_path, dpi=600, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f"[OK] Combined figure saved -> {combo_path}")

# ?? 6. Create Individual Figures ???????????????????????????????????????????????
print("\nPlotting individual PSM network figures...")
for g in GROUPS:
    fig, ax = plt.subplots(figsize=(16, 14))
    fig.patch.set_facecolor('white')
    draw_panel(ax, raw_graphs[g], g, pos_all)
    
    # Colorbar ONLY (No legend)
    cbar_ax_ind = fig.add_axes([0.95, 0.25, 0.02, 0.5])
    cbar_ind = fig.colorbar(sm, cax=cbar_ax_ind)
    cbar_ind.outline.set_visible(False)
    cbar_ind.ax.tick_params(labelsize=14, width=0, length=0)
    cbar_ind.set_label("Average Expression (Log2)", size=18, weight='bold', labelpad=15)
    
    ind_path = OUT_DIR / f"Network_{g}_psm_pearson_r03.png"
    plt.savefig(ind_path, dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  [OK] Saved -> {ind_path.name}")

# Generate standalone legend
fig_leg = plt.figure(figsize=(6, 5))
ax_leg = fig_leg.add_subplot(111)
ax_leg.axis('off')
ax_leg.legend(handles=legend_elements, loc='center', fontsize=20,
              frameon=True, framealpha=0.9, edgecolor='#cccccc',
              title='Network Features', title_fontproperties={'weight': 'bold', 'size': 24})
leg_path = OUT_DIR / "Network_psm_pearson_legend.png"
fig_leg.savefig(leg_path, dpi=600, bbox_inches='tight', facecolor='white')
plt.close(fig_leg)
print(f"  [OK] Standalone legend saved -> {leg_path.name}")

print("\nFinished generating all PSM network variants!")


