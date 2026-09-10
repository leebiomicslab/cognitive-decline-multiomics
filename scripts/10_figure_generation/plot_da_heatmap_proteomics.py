import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Heatmap of differentially expressed proteins
Manuscript:       Results - Proteomics DA
Figure/Table:     Fig. 3b
"""

from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def plot_da_heatmap_proteomics():
    # ??? 1. Setup Paths ????????????????????????????????????????????????????????
    CLINICAL_FILE = CLINICAL_CSV
    PROTEOMICS_DATA = PROTEOMICS_LOG2

    MODULE_MAP_FILE = FASTSPAR_PROT_DIR / "pearson" / "module_protein_assignment.csv"
    OUTPUT_DIR = FASTSPAR_PROT_DIR
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    # ??? 2. Load Clinical & Group Mapping ?????????????????????????????????????
    clinical = pd.read_csv(CLINICAL_FILE)
    clinical.rename(columns={'id': 'SampleID'}, inplace=True)
    group_map = dict(zip(clinical['SampleID'].astype(str), clinical['group']))

    # ??? 3. Load Nominal Significant Proteins & Module Mapping ???????????????
    NOMINAL_FILE = PROT_ASSOC_DIR / "Main_Model" / "nominal_significant_proteins.csv"
    nom_df = pd.read_csv(NOMINAL_FILE)
    
    # Load Module Assignments for Sidebar mapping
    MODULE_MAP_FILE = FASTSPAR_PROT_DIR / "pearson" / "module_protein_assignment.csv"
    mod_assign_df = pd.read_csv(MODULE_MAP_FILE)
    mod_assign_df['Module_Name'] = mod_assign_df['Module'].apply(lambda x: f"P{x+1}")
    mod_lookup = dict(zip(mod_assign_df['Protein'], mod_assign_df['Module_Name']))

    # ??? 4. Load Proteomics Expression ????????????????????????????????????????
    expr_raw = pd.read_csv(PROTEOMICS_DATA)
    if 'Protein Name' in expr_raw.columns:
        expr_raw = expr_raw.rename(columns={'Protein Name': 'SampleID'})
        
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
    expr_raw.rename(columns=PROTEIN_RENAME, inplace=True)
    
    expr_raw['SampleID'] = expr_raw['SampleID'].astype(str)
    
    # Filter for common samples
    expr_raw = expr_raw[expr_raw['SampleID'].isin(group_map.keys())].copy()
    expr_raw['Group'] = expr_raw['SampleID'].map(group_map)

    # ??? 5. Process Entries (Aggregate '&' groups) ??????????????????????????
    # We want exactly one row in heatmap for each entry in nominal_significant_proteins.csv
    # If the entry name is 'A&B', we take the mean(expr[A], expr[B])
    heatmap_rows = []
    row_module_map = {}
    
    print(f"Processing {len(nom_df)} distinct entries...")
    for label in nom_df['Protein']:
        label_str = str(label).strip()
        
        # Priority 1: Exact match with the label (e.g. "DSG1&DSG4")
        if label_str in expr_raw.columns:
            exist_parts = [label_str]
        else:
            # Priority 2: Split and find individual parts
            parts = [p.strip() for p in label_str.split('&')]
            exist_parts = [p for p in parts if p in expr_raw.columns]
        
        if len(exist_parts) > 0:
            # Mean of intensity across these features (if multiple matched)
            entry_expr = expr_raw[exist_parts].mean(axis=1)
            heatmap_rows.append(pd.Series(entry_expr, name=label_str))
            # Map to module
            row_module_map[label_str] = mod_lookup.get(label_str, "Other")
        else:
            print(f"Warning: No match found for '{label_str}' in data.")

    if not heatmap_rows:
        print("Error: No data rows found for heatmap.")
        return

    heatmap_data_all = pd.concat(heatmap_rows, axis=1)

    # ??? Z-score first, then Group Mean ??????????????????????????????????????
    # Standardize (Z-score) each protein across all samples
    z_scored_expr = heatmap_data_all.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
    z_scored_expr['Group'] = expr_raw['Group'].values
    
    # Aggregate by group (calculate mean of Z-scores)
    heatmap_df = z_scored_expr.groupby("Group").mean().T
    heatmap_df = heatmap_df[["C1", "C2", "C3"]]

    # ??? 6. Prepare Sidebar (Module Assignments) ????????????????????????????
    # row_module_map already contains the label -> module mapping
    # Palette: High contrast
    module_colors_map = {
        'P1': '#e377c2', 'P2': '#17becf', 'P3': '#bcbd22', 'P4': '#ff7f0e', 'P5': '#9467bd', 'Other': '#cccccc'
    }

    row_colors = pd.Series(heatmap_df.index).map(row_module_map).map(module_colors_map)
    row_colors.index = heatmap_df.index
    row_colors.name = "Module"


    # ??? 7. Visualization ????????????????????????????????????????????????????
    sns.set_theme(style="white")
    
    df_plot = heatmap_df
    
    # Clustering will group proteins by abundance pattern
    g = sns.clustermap(
        df_plot,
        cmap="RdBu_r",
        col_cluster=False,
        row_cluster=True,
        figsize=(10, 14),
        vmin=-1.5, vmax=1.5,
        row_colors=row_colors,
        cbar_pos=(0.02, 0.82, 0.04, 0.15),
        annot=True, fmt=".2f",
        annot_kws={"size": 9, "weight": "bold"},
        linewidths=0.5,
        linecolor='white'
    )

    g.ax_heatmap.set_title("Mean Protein Intensity (Z-score)", fontsize=16, fontweight='bold', pad=25)
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.set_xlabel("")

    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=0, fontsize=14, fontweight='bold')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)

    # Legend for row colors
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=l) for l, c in module_colors_map.items()]
    g.ax_heatmap.legend(handles=legend_elements, title="Modules", bbox_to_anchor=(1.25, 1), loc='upper left', fontsize=12)

    out_path = OUTPUT_DIR / "da_proteins_mean_intensity_heatmap.png"
    g.savefig(out_path, dpi=600, bbox_inches='tight', facecolor='white')
    print(f"Heatmap saved to {out_path}")



if __name__ == "__main__":
    plot_da_heatmap_proteomics()


