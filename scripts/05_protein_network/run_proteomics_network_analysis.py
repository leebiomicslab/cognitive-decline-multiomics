"""
Purpose:          Build Pearson and Spearman protein co-expression networks per
                  cognitive phenotype group. Apply Leiden community detection
                  (resolution=1.0, seed=42) to identify protein modules.
                  Compute module eigengenes via PCA. Export GraphML files.
Manuscript:       Methods — Protein co-expression network analysis
Figure/Table:     Fig. 4a-c; Supplementary Fig. 6
Input:            proteomics_ancova_detailed_results.csv (from 03_proteomics_association),
                  PROTEOMICS_LOG2, CLINICAL_CSV
Output:           results/fastspar/proteomics/{pearson,spearman}/
Main dependencies: pandas, numpy, networkx, igraph, leidenalg, scikit-learn, seaborn
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import igraph as ig
import leidenalg as la
from sklearn.decomposition import PCA
from pathlib import Path
import seaborn as sns

# -------------------------
# Configuration
# -------------------------
BASE_DIR = FASTSPAR_PROT_DIR.parent.parent  # resolves to results/
OUTPUT_DIR = FASTSPAR_PROT_DIR
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Input files
DA_RESULTS = PROT_ASSOC_DIR / "Main_Model" / "proteomics_ancova_detailed_results.csv"
EXPRESSION_FILE = PROTEOMICS_LOG2
GROUP_FILE = CLINICAL_CSV


# Parameters
GROUPS = ("C1", "C2", "C3")
RESOLUTION = 1.0
SEED = 42
P_THRESH = 0.05

# -------------------------
# Helper Functions (identical logic to reference script)
# -------------------------

def abs_weighted_degree(G):
    degree = {n: 0.0 for n in G.nodes()}
    for u, v, d in G.edges(data=True):
        w = abs(d.get("weight", 0))
        degree[u] += w
        degree[v] += w
    return degree

def louvain_positions(G, gap=0.02, radius=1.0, resolution=RESOLUTION, seed=SEED):
    """Leiden clustering + sorted circular layout (same as reference)."""
    nodes = list(G.nodes())
    if not nodes:
        return {}, {}

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

    part = la.find_partition(
        gi, la.RBConfigurationVertexPartition,
        weights=gi.es["weight"], resolution_parameter=resolution, seed=seed
    )

    module_map = {}
    for cid, comm in enumerate(part):
        for vid in comm:
            module_map[gi.vs[vid]['name']] = cid

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
        seg_end = start + length - gap / 2

        order = sorted(comm, key=lambda n: deg_abs.get(n, 0), reverse=True)
        if len(order) == 1:
            thetas = np.array([(seg_start + seg_end) / 2]) * 2 * np.pi
        else:
            thetas = np.linspace(seg_start, seg_end, len(order), endpoint=False) * 2 * np.pi

        for n, t in zip(order, thetas):
            pos[n] = (radius * np.cos(t), radius * np.sin(t))
        start += length

    return pos, module_map


def build_network(df, group, method='pearson', label_col='Group', id_col='SampleID'):
    print(f"  Building {group} network ({method})...")
    if group == "ALL":
        gdf = df.drop(columns=[label_col, id_col])
    else:
        gdf = df[df[label_col] == group].drop(columns=[label_col, id_col])
    gdf = gdf.loc[:, gdf.nunique() > 1]
    corr = gdf.corr(method=method, min_periods=3)
    G = nx.Graph()
    G.add_nodes_from(corr.columns)
    cols = corr.columns
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            r = corr.iloc[i, j]
            if np.isfinite(r):
                G.add_edge(cols[i], cols[j], weight=float(r))
    print(f"  [OK] {group}: nodes={G.number_of_nodes()}, edges={G.number_of_edges()}")
    return G


def global_expr_minmax(df, groups, label_col='Group', id_col='SampleID'):
    vals = []
    for g in groups:
        gdf = df[df[label_col] == g].drop(columns=[label_col, id_col])
        if not gdf.empty:
            vals.append(gdf.mean().values)
    if vals:
        all_vals = np.concatenate(vals)
        return float(np.nanmin(all_vals)), float(np.nanmax(all_vals))
    return 0.0, 1.0


def plot_network(G, df, group, pos, title, module_map=None,
                 label_col='Group', id_col='SampleID',
                 r_thresh=0.3, expr_vmin=None, expr_vmax=None,
                 save_path=None, font_size=11):
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    if group != "ALL":
        gdf = df[df[label_col] == group].drop(columns=[label_col, id_col])
    else:
        gdf = df.drop(columns=[label_col, id_col])

    mean_expr = gdf.mean().to_dict()

    # Edges
    edges = [(u, v, d) for u, v, d in G.edges(data=True) if abs(d.get("weight", 0)) >= r_thresh]
    # Premium edge colors (Colorblind safe: #D55E00 vermillion, #0072B2 blue)
    edge_colors = ['#D55E00' if d['weight'] >= 0 else '#0072B2' for _, _, d in edges]
    edge_widths = [0.5 + 4 * (abs(d['weight']) - r_thresh) / (1 - r_thresh)
                   if (1 - r_thresh) > 0 else 1.0 for _, _, d in edges]

    # Node fill color: Reds (expression)
    expr_vals = np.array([mean_expr.get(n, 0) for n in G.nodes()])
    if expr_vmin is None or expr_vmax is None:
        expr_vmin, expr_vmax = np.min(expr_vals), np.max(expr_vals)

    cmap_fill = plt.cm.Reds
    norm = plt.Normalize(vmin=expr_vmin, vmax=expr_vmax)
    node_fill_colors = [cmap_fill(norm(mean_expr.get(n, 0))) for n in G.nodes()]

    # Node border: module color
    if module_map:
        unique_modules = sorted(list(set(module_map.values())))
        # Use Set2 or tab20 for categorical
        cmap_border = plt.get_cmap('Set2' if len(unique_modules) <= 8 else 'tab20')
        node_border_colors = [cmap_border(module_map.get(n, 0) % (8 if len(unique_modules)<=8 else 20)) for n in G.nodes()]
    else:
        node_border_colors = ['#333333'] * G.number_of_nodes()
        unique_modules = []
        cmap_border = None

    # Node size: weighted degree
    deg_abs = {n: 0 for n in G.nodes()}
    for u, v, d in G.edges(data=True):
        if abs(d.get("weight", 0)) >= r_thresh:
            w = abs(d.get("weight", 0))
            deg_abs[u] += w
            deg_abs[v] += w
    max_deg = max(deg_abs.values()) if deg_abs and max(deg_abs.values()) > 0 else 1.0
    node_sizes = [400 + 1200 * (deg_abs.get(n, 0) / max_deg) for n in G.nodes()]

    fig, ax = plt.subplots(figsize=(16, 14))
    
    # Curved edges
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(u, v) for u, v, _ in edges],
                           edge_color=edge_colors, width=edge_widths, alpha=0.6,
                           connectionstyle="arc3,rad=0.15", arrows=True, arrowstyle='-')

    # Premium Nodes
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color='white', linewidths=0)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=[s * 0.75 for s in node_sizes], 
                           node_color=node_fill_colors, linewidths=0)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color='none',
                           edgecolors=node_border_colors, linewidths=3.5)

    # Radial Labels
    def format_label(name):
        if '&' in name and len(name) > 20:
            parts = name.split('&')
            mid = len(parts) // 2 + (len(parts) % 2) # Slightly favor top line
            return '&'.join(parts[:mid]) + '&\n' + '&'.join(parts[mid:])
        return name

    for node, (x, y) in pos.items():
        angle = np.arctan2(y, x)
        idx = list(G.nodes()).index(node)
        radius_offset = 1.05 + (node_sizes[idx] / 20000)
        
        text_x = x * radius_offset
        text_y = y * radius_offset
        
        ha = 'left' if x > 0.01 else ('right' if x < -0.01 else 'center')
        va = 'bottom' if y > 0.01 else ('top' if y < -0.01 else 'center')
            
        ax.text(text_x, text_y, format_label(node), fontsize=font_size, fontweight='bold', 
                fontfamily='Arial', ha=ha, va=va, color='#1a1a1a', zorder=10)

    # Minimalist Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_fill, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.95, 0.25, 0.02, 0.5]) 
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=14, width=0, length=0)
    cbar.set_label("Average Expression (Log2)", size=18, weight='bold', labelpad=15)

    fig.suptitle(f"{title} | Edges |r| >= {r_thresh}", fontsize=18, fontweight="bold", fontfamily='Arial', y=0.98)
    
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.axis("off")
    
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight', facecolor='white')
        print(f"  Saved -> {save_path.name}")
    plt.close(fig)


def export_graphml(G, df, group, module_map, filename, label_col='Group', id_col='SampleID'):
    G_copy = G.copy()
    if group != "ALL":
        gdf = df[df[label_col] == group].drop(columns=[label_col, id_col])
    else:
        gdf = df.drop(columns=[label_col, id_col])
    mean_expr = gdf.mean().to_dict()
    for u, v, d in G_copy.edges(data=True):
        w = float(d.get("weight", 0))
        d["weight"] = w
        d["sign"] = "positive" if w >= 0 else "negative"
    for n, d in G_copy.nodes(data=True):
        d["module"] = int(module_map.get(n, -1))
        d["mean_abundance"] = float(mean_expr.get(n, 0))
    nx.write_graphml(G_copy, str(filename))
    print(f"  [OK] GraphML exported: {filename.name}")


def compute_module_eigengenes(df, module_map, output_csv, label_col='Group', id_col='SampleID'):
    gdf = df.drop(columns=[label_col, id_col])
    gdf.index = df[id_col]
    modules = {}
    for feat, mod in module_map.items():
        modules.setdefault(mod, []).append(feat)

    eig_df = pd.DataFrame(index=gdf.index)
    loadings_list = []
    for m, features in modules.items():
        feats = [f for f in features if f in gdf.columns]
        if len(feats) < 2:
            continue
        pca = PCA(n_components=1)
        eig = pca.fit_transform(gdf[feats])
        eig_df[f"M{m}_prot"] = eig.ravel()
        loadings = pca.components_[0]
        explained_var = pca.explained_variance_ratio_[0]
        for feature, loading in zip(feats, loadings):
            loadings_list.append({
                "Module": m, "Feature": feature,
                "Loading": loading, "ExplainedVariance_PC1": explained_var
            })

    eig_df.to_csv(output_csv)
    print(f"  [OK] Module eigengenes: {output_csv.name}")
    if loadings_list:
        loadings_df = pd.DataFrame(loadings_list)
        loadings_out = Path(str(output_csv).replace(".csv", "_loadings.csv"))
        loadings_df.to_csv(loadings_out, index=False)
        print(f"  [OK] PCA loadings: {loadings_out.name}")
    return eig_df


# -------------------------
# Main Pipeline
# -------------------------
# Protein name standardization map
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

def run_analysis(method):
    print(f"\n{'='*50}\nMethod: {method.upper()}\n{'='*50}")
    method_dir = OUTPUT_DIR / method
    os.makedirs(method_dir, exist_ok=True)

    # 1. Significant proteins
    da = pd.read_csv(DA_RESULTS)
    da['Protein'] = da['Protein'].replace(PROTEIN_RENAME)
    sig_proteins = set(da[da['P_Value'] < P_THRESH]['Protein'].tolist())
    print(f"Significant proteins (p<{P_THRESH}): {len(sig_proteins)}")

    # 2. Load expression data
    expr_raw = pd.read_csv(EXPRESSION_FILE)
    expr_raw.rename(columns={expr_raw.columns[0]: 'SampleID'}, inplace=True)
    # Apply protein renaming to expression columns
    expr_raw.rename(columns=PROTEIN_RENAME, inplace=True)

    # 3. Group map from clinical file
    clin = pd.read_csv(GROUP_FILE)
    # Handle column naming variations
    clin.rename(columns={'id': 'SampleID'}, inplace=True, errors='ignore')
    if 'SampleID' not in clin.columns:
        clin.rename(columns={clin.columns[0]: 'SampleID'}, inplace=True)
    if 'group' in clin.columns:
        clin.rename(columns={'group': 'Group'}, inplace=True)

    # 4. Filter to valid sig proteins that exist in expression file
    valid_cols = [col for col in expr_raw.columns if col in sig_proteins]
    print(f"Matched in expression file: {len(valid_cols)}")

    df = expr_raw[['SampleID'] + valid_cols].copy()
    df = df.merge(clin[['SampleID', 'Group']], on='SampleID', how='inner')
    print(f"Samples after merge: {len(df)}")

    # 5. Build per-group networks
    graphs = [build_network(df, g, method=method) for g in GROUPS]
    G_C1, G_C2, G_C3 = graphs

    # 6. ALL network (union)
    G_ALL = nx.compose_all(graphs)
    print(f"ALL network: nodes={G_ALL.number_of_nodes()}, edges={G_ALL.number_of_edges()}")

    # 7. Leiden clustering on G_ALL
    nodes = list(G_ALL.nodes())
    node_idx = {n: i for i, n in enumerate(nodes)}
    edges_ig, weights_ig = [], []
    for u, v, d in G_ALL.edges(data=True):
        edges_ig.append((node_idx[u], node_idx[v]))
        weights_ig.append(abs(d.get("weight", 0)))

    g_ig = ig.Graph()
    g_ig.add_vertices(len(nodes))
    g_ig.vs["name"] = nodes
    g_ig.add_edges(edges_ig)
    g_ig.es["weight"] = weights_ig

    print("Running Leiden clustering...")
    part = la.find_partition(
        g_ig, la.RBConfigurationVertexPartition,
        weights=g_ig.es["weight"], resolution_parameter=RESOLUTION, seed=SEED
    )

    module_map = {}
    for cid, comm in enumerate(part):
        for vid in comm:
            module_map[g_ig.vs[vid]['name']] = cid

    n_mods = len(set(module_map.values()))
    print(f"Leiden detected {n_mods} modules.")

    # Save module assignments
    pd.DataFrame([(f, m) for f, m in module_map.items()],
                 columns=["Protein", "Module"]).to_csv(
        method_dir / "module_protein_assignment.csv", index=False
    )

    # 8. Layout
    pos_all, _ = louvain_positions(G_ALL, gap=0.02, resolution=RESOLUTION, seed=SEED)

    # 9. Module eigengenes
    eig_file = method_dir / f"ME_proteomics_ALL_{method}.csv"
    compute_module_eigengenes(df, module_map, eig_file)

    # 10. Plot for multiple thresholds
    thresholds = [0.3, 0.5]
    expr_vmin, expr_vmax = global_expr_minmax(df, GROUPS)
    graphml_dir = method_dir / "GraphML_export"
    os.makedirs(graphml_dir, exist_ok=True)

    print("Plotting networks...")
    for thresh in thresholds:
        thresh_str = str(thresh).replace('.', '')
        plot_network(G_ALL, df, "ALL", pos_all,
                     f"Proteomics - ALL ({method}) |r|>={thresh}",
                     module_map=module_map, r_thresh=thresh,
                     expr_vmin=expr_vmin, expr_vmax=expr_vmax,
                     save_path=method_dir / f"Network_ALL_{method}_r{thresh_str}.png",
                     font_size=24)
        for i, g in enumerate(GROUPS):
            plot_network(graphs[i], df, g, pos_all,
                         f"Proteomics - {g} ({method}) |r|>={thresh}",
                         module_map=module_map, r_thresh=thresh,
                         expr_vmin=expr_vmin, expr_vmax=expr_vmax,
                         save_path=method_dir / f"Network_{g}_{method}_r{thresh_str}.png",
                         font_size=24)

    # 11. Export GraphML
    print("Exporting GraphML...")
    export_graphml(G_ALL, df, "ALL", module_map,
                   graphml_dir / f"Proteomics_ALL_{method}.graphml")
    for i, g in enumerate(GROUPS):
        export_graphml(graphs[i], df, g, module_map,
                       graphml_dir / f"Proteomics_{g}_{method}.graphml")

    print(f"\n[OK] {method} analysis complete. Output: {method_dir}")

    # Generate standalone legend
    fig_leg = plt.figure(figsize=(6, 5))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis('off')

    from matplotlib.lines import Line2D
    unique_modules = sorted(set(module_map.values()))
    cmap_border = plt.get_cmap('Set2' if len(unique_modules) <= 8 else 'tab20')
    legend_elements = [
        Line2D([0], [0], color='#D55E00', lw=5, label='Positive correlation'),
        Line2D([0], [0], color='#0072B2', lw=5, label='Negative correlation')
    ]
    for cid in unique_modules:
        legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='w',
                                      markeredgecolor=cmap_border(cid % (8 if len(unique_modules)<=8 else 20)),
                                      markersize=22, markeredgewidth=4,
                                      label=f'Module P{cid+1}'))

    ax_leg.legend(handles=legend_elements, loc='center', fontsize=20,
                  frameon=True, framealpha=0.9, edgecolor='#cccccc',
                  title='Network Features', title_fontproperties={'weight': 'bold', 'size': 24})

    leg_path = method_dir / f"Network_{method}_legend.png"
    fig_leg.savefig(leg_path, dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig_leg)
    print(f"Standalone legend generated at: {leg_path.name}")


if __name__ == "__main__":
    run_analysis("pearson")
    run_analysis("spearman")
    print("\n[OK] All done!")
