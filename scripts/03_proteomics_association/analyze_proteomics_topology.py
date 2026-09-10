import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Protein network topology metrics
Manuscript:       Methods - Protein network
Figure/Table:     Supplementary Table 3b
"""

from pathlib import Path
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
import networkx.algorithms.community as nx_comm

def analyze_proteomics_topology():
    # ??? 1. Setup Paths ????????????????????????????????????????????????????????
    NOMINAL_PROTEINS = PROT_ASSOC_DIR / "Main_Model" / "nominal_significant_proteins.csv"
    GRAPH_DIR = FASTSPAR_PROT_DIR / "pearson" / "GraphML_export"
    OUTPUT_DIR = FASTSPAR_PROT_DIR / "network_topology"
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    THRESH = 0.3 # As requested, matching the PNGs r03
    
    print("Loading significant protein list...")
    # These are the entries (25 rows) that are officially significant
    nom_df = pd.read_csv(NOMINAL_PROTEINS)
    valid_features = set(nom_df["Protein"].str.strip().tolist())

    groups = ["C1", "C2", "C3"]
    graphs = {}

    # ??? 2. Build and Filter Graphs ??????????????????????????????????????????
    for g in groups:
        gml_path = GRAPH_DIR / f"Proteomics_{g}_pearson.graphml"
        if not gml_path.exists():
            print(f"Warning: {gml_path} not found.")
            continue
            
        G = nx.read_graphml(gml_path)
        # Keep nodes that are in our nominal significant list
        nodes_to_keep = [n for n in G.nodes() if n.strip() in valid_features]
        subG = G.subgraph(nodes_to_keep).copy()
        
        # Filter edges by threshold r >= 0.3
        # GraphML usually stores weight as 'weight' or 'r'
        edges_to_keep = []
        for u, v, d in subG.edges(data=True):
            # Try getting weight from 'weight' first, then 'correlation'
            w = float(d.get("weight", d.get("correlation", 0)))
            if abs(w) >= THRESH:
                edges_to_keep.append((u, v, {'weight': w, 'abs_weight': abs(w)}))
                
        f_subG = nx.Graph()
        f_subG.add_nodes_from(nodes_to_keep)
        f_subG.add_edges_from(edges_to_keep)
        graphs[g] = f_subG
        print(f"Group {g}: Built filtered graph with {f_subG.number_of_nodes()} nodes and {f_subG.number_of_edges()} edges.")

    # ??? 3. Topological Analysis ??????????????????????????????????????????????
    results = {}

    for g, G in graphs.items():
        print(f"Analyzing Group {g} metrics...")
        res = {}
        
        isolated = list(nx.isolates(G))
        G_active = G.copy()
        G_active.remove_nodes_from(isolated)
        
        # Edges
        pos_e = sum(1 for u, v, d in G.edges(data=True) if d.get('weight', 0) > 0)
        neg_e = sum(1 for u, v, d in G.edges(data=True) if d.get('weight', 0) < 0)
        
        res['Nodes_Total'] = G.number_of_nodes()
        res['Nodes_Isolated'] = len(isolated)
        res['Edges_Total'] = G.number_of_edges()
        res['Edges_Positive'] = pos_e
        res['Edges_Negative'] = neg_e
        res['Pos_to_Neg_Ratio'] = round(pos_e / neg_e, 2) if neg_e > 0 else float('inf') if pos_e > 0 else 0
        
        # Density
        res['Density'] = round(nx.density(G), 4)
        
        # Modularity
        if G_active.number_of_edges() > 0:
            comm = nx_comm.greedy_modularity_communities(G_active, weight='abs_weight')
            mod_score = nx_comm.modularity(G_active, comm, weight='abs_weight')
            res['Modularity_Q'] = round(mod_score, 4)
            res['Num_Communities'] = len(comm)
        else:
            res['Modularity_Q'] = 0.0
            res['Num_Communities'] = 0

        # Hubs
        if G.number_of_edges() > 0:
            deg_cent = nx.degree_centrality(G)
            bet_cent = nx.betweenness_centrality(G, weight='abs_weight')
            
            hubs = sorted(deg_cent.items(), key=lambda x: x[1], reverse=True)[:3]
            res['Top3_Hubs_Degree'] = ", ".join([f"{k} ({v:.2f})" for k, v in hubs if v > 0])
            
            bet_hubs = sorted(bet_cent.items(), key=lambda x: x[1], reverse=True)[:3]
            res['Top3_Hubs_Betweenness'] = ", ".join([f"{k} ({v:.2f})" for k, v in bet_hubs if v > 0])
        else:
            res['Top3_Hubs_Degree'] = "None"
            res['Top3_Hubs_Betweenness'] = "None"
            
        results[g] = res

    # ??? 4. Robustness Simulation ?????????????????????????????????????????????
    def simulate_attack(G):
        nodes_by_deg = sorted(dict(G.degree()).items(), key=lambda x: x[1], reverse=True)
        nodes_to_remove = [n[0] for n in nodes_by_deg]
        
        lcc_sizes = []
        G_temp = G.copy()
        initial_nodes = G.number_of_nodes()
        
        if initial_nodes == 0: return [0]
             
        def get_lcc_size(gr):
             if gr.number_of_nodes() == 0: return 0
             comps = list(nx.connected_components(gr))
             if not comps: return 0
             return len(max(comps, key=len)) / initial_nodes

        lcc_sizes.append(get_lcc_size(G_temp))
        for n in nodes_to_remove:
            G_temp.remove_node(n)
            lcc_sizes.append(get_lcc_size(G_temp))
        return lcc_sizes

    plt.figure(figsize=(10, 7))
    colors = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'} # Consistent colors
    for g in groups:
        if g not in graphs: continue
        G_eval = graphs[g]
        if G_eval.number_of_edges() == 0: continue
        
        lcc_trace = simulate_attack(G_eval)
        x_val = np.linspace(0, 1, len(lcc_trace))
        plt.plot(x_val, lcc_trace, label=f"Group {g}", color=colors[g], marker='o', markersize=4, linewidth=2.5)

    plt.title(f"Proteomics Network Robustness (Targeted Attack) | |r| >= {THRESH}", fontsize=15, weight='bold')
    plt.xlabel("Fraction of Hubs Removed", fontsize=13)
    plt.ylabel("Relative Size of LCC", fontsize=13)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.ylim(-0.05, 1.05)
    
    robust_path = OUTPUT_DIR / f"robustness_attack_{THRESH}.png"
    plt.savefig(robust_path, dpi=600, bbox_inches='tight')
    plt.close()

    # ??? 5. Save Results ??????????????????????????????????????????????????????
    res_df = pd.DataFrame(results).T
    metric_path = OUTPUT_DIR / f"network_metrics_{THRESH}.csv"
    res_df.to_csv(metric_path)
    
    print("\n=== Proteomics Network Topology Summary ===")
    print(res_df.to_string())
    print(f"\n??Analysis complete. Metrics saved to {metric_path}")
    print(f"??Robustness plot saved to {robust_path}")

if __name__ == "__main__":
    analyze_proteomics_topology()


