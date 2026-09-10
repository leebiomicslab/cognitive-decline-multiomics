"""
Purpose:          Compute microbial co-abundance network topology metrics
                  (degree, betweenness, modularity, robustness) for the
                  DA-taxa subgraph per cognitive phenotype group.
Manuscript:       Methods — Microbial network topology
Figure/Table:     Supplementary Table 3
Input:            microbiome_G_{C1,C2,C3}.graphml (from FastSpar),
                  microbiota_da_results_final.csv (from 04_microbiome_differential_abundance)
Output:           results/fastspar/microbiota/network_topology/network_metrics_*.csv,
                  robustness_attack_*.png
Main dependencies: networkx, pandas, numpy, matplotlib

# NOTE: correlation threshold THRESH=0.1 applied here
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import networkx.algorithms.community as nx_comm

DA_RESULTS  = MICRO_ASSOC_DIR / "Model_A_Main" / "microbiota_da_results_final.csv"
GRAPH_DIR   = FASTSPAR_MICRO_DIR
OUTPUT_DIR  = FASTSPAR_MICRO_DIR / "network_topology"
os.makedirs(OUTPUT_DIR, exist_ok=True)


THRESH = 0.1
P_VALUE_THRESH = 0.05

da = pd.read_csv(DA_RESULTS)
sig_taxa = da[da["Union_Sig"] == True]["taxon"].str.strip().tolist()
valid_taxa = set(sig_taxa)

groups = ["C1", "C2", "C3"]
graphs = {}

# 1. Load graphs and build subgraphs
for g in groups:
    gml_path = GRAPH_DIR / f"microbiome_G_{g}.graphml"
    G = nx.read_graphml(gml_path)
    nodes_to_keep = [n for n in G.nodes() if n.strip() in valid_taxa]
    subG = G.subgraph(nodes_to_keep).copy()
    
    edges_to_keep = []
    for u, v, d in subG.edges(data=True):
        w = float(d.get("weight", 0))
        p = float(d.get("p_value", 1.0))
        if abs(w) >= THRESH and p < P_VALUE_THRESH:
            # Need strict float since GraphML parses as string sometimes
            edges_to_keep.append((u, v, {'weight': w, 'abs_weight': abs(w)}))
            
    f_subG = nx.Graph()
    f_subG.add_nodes_from(nodes_to_keep)
    f_subG.add_edges_from(edges_to_keep)
    graphs[g] = f_subG

# Analysis Results Dictionary
results = {}

for g, G in graphs.items():
    print(f"\nAnalyzing Group {g}...")
    res = {}
    
    isolated = list(nx.isolates(G))
    G_no_iso = G.copy()
    G_no_iso.remove_nodes_from(isolated)
    
    # 1. Pos/Neg Ratio
    pos_e = sum(1 for u, v, d in G.edges(data=True) if float(d.get('weight', 0)) > 0)
    neg_e = sum(1 for u, v, d in G.edges(data=True) if float(d.get('weight', 0)) < 0)
    res['Nodes_Total'] = G.number_of_nodes()
    res['Nodes_Isolated'] = len(isolated)
    res['Edges_Total'] = pos_e + neg_e
    res['Edges_Positive'] = pos_e
    res['Edges_Negative'] = neg_e
    res['Pos_to_Neg_Ratio'] = round(pos_e / neg_e, 2) if neg_e > 0 else float('inf') if pos_e > 0 else 0
    
    # 2. Modularity
    if G_no_iso.number_of_edges() > 0:
        comm = nx_comm.greedy_modularity_communities(G_no_iso, weight='abs_weight')
        mod_score = nx_comm.modularity(G_no_iso, comm, weight='abs_weight')
        res['Modularity_Q'] = round(mod_score, 4)
        res['Num_Communities'] = len(comm)
    else:
        res['Modularity_Q'] = 0.0
        res['Num_Communities'] = 0

    # 3. Key Nodes (Centrality)
    if G.number_of_edges() > 0:
        deg_cent = nx.degree_centrality(G)
        bet_cent = nx.betweenness_centrality(G, weight='abs_weight')
        
        # Sort by degree centrality
        hubs = sorted(deg_cent.items(), key=lambda x: x[1], reverse=True)[:3]
        res['Top3_Hubs_Degree'] = ", ".join([f"{k.replace('\"','')} ({v:.2f})" for k, v in hubs if v > 0])
        
        bet_hubs = sorted(bet_cent.items(), key=lambda x: x[1], reverse=True)[:3]
        res['Top3_Hubs_Betweenness'] = ", ".join([f"{k.replace('\"','')} ({v:.2f})" for k, v in bet_hubs if v > 0])
    else:
        res['Top3_Hubs_Degree'] = "None"
        res['Top3_Hubs_Betweenness'] = "None"
        
    results[g] = res

# 4. Robustness Simulation (Targeted Attack by Degree)
def simulate_attack(G):
    # Sort nodes by strict degree count explicitly
    nodes_by_deg = sorted(dict(G.degree()).items(), key=lambda x: x[1], reverse=True)
    nodes_to_remove = [n[0] for n in nodes_by_deg]
    
    lcc_sizes = []
    G_temp = G.copy()
    initial_nodes = G.number_of_nodes()
    
    if initial_nodes == 0:
         return [0], [0]
         
    def get_lcc_size(gr):
         if gr.number_of_nodes() == 0: return 0
         if gr.number_of_edges() == 0: return 1 / initial_nodes  # 1 isolated node is technically a connected component of size 1
         comps = list(nx.connected_components(gr))
         if not comps: return 0
         return len(max(comps, key=len)) / initial_nodes

    lcc_sizes.append(get_lcc_size(G_temp))
    
    for n in nodes_to_remove:
        G_temp.remove_node(n)
        lcc_sizes.append(get_lcc_size(G_temp))
            
    return lcc_sizes

plt.figure(figsize=(8, 6))
colors = {'C1': '#377EB8', 'C2': '#4DAF4A', 'C3': '#E41A1C'}
for g in groups:
    G_eval = graphs[g]
    if G_eval.number_of_edges() == 0:
         continue
    lcc_trace = simulate_attack(G_eval)
    x_val = np.linspace(0, 1, len(lcc_trace))
    plt.plot(x_val, lcc_trace, label=f"Group {g}", color=colors[g], marker='o', markersize=5, linewidth=2)

plt.title(f"Network Robustness (Targeted Hub Removal) | r >= {THRESH}", fontsize=14, weight='bold')
plt.xlabel("Fraction of nodes removed", fontsize=12)
plt.ylabel("Relative size of Largest Connected Component (LCC)", fontsize=12)
plt.legend(fontsize=10)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
robust_path = OUTPUT_DIR / f"robustness_attack_{THRESH}.png"
plt.savefig(robust_path, dpi=300)
plt.close()

# Print / Write results
res_df = pd.DataFrame(results).T

# Output as markdown text directly print
print("\n=== Network Topology Metrics ===")
print("Metrics generated successfully.")
res_df.to_csv(OUTPUT_DIR / f"network_metrics_{THRESH}.csv")
print(res_df.to_string())
print(f"\n✅ Network analysis complete. Robustness plot saved to: {robust_path}")
