import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Venn diagram of microbiome/proteomics DA overlap
Manuscript:       Results - Cross-omics overlap
Figure/Table:     Supplementary Fig. 4
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os

# Define paths

# Inputs
microbiota_file = str(PROTEOMICS_LOG2)
proteomics_file = str(PROTEOMICS_LOG2)
group_file = str(GROUP_CSV)

# Outputs
output_dir = str(RESULTS_DIR)
figures_dir = os.path.join(output_dir, "figures")
common_samples_file = str(COMMON_SAMPLES_CSV)
venn_file = os.path.join(figures_dir, "venn_intersection.png")

os.makedirs(figures_dir, exist_ok=True)

print("Loading data for intersection...")

# 1. Microbiota (Rows=Features, Cols=Samples)
try:
    df_mic = pd.read_csv(microbiota_file, index_col=0)
    mic_samples = set(df_mic.columns.astype(str))
    print(f"Microbiota samples: {len(mic_samples)}")
except Exception as e:
    print(f"Error loading Microbiota: {e}")
    mic_samples = set()

# 2. Proteomics (Assuming 'id' column exists or index is ID)
try:
    df_prot = pd.read_csv(proteomics_file)
    # Check for ID column
    if 'id' in df_prot.columns:
        prot_samples = set(df_prot['id'].astype(str))
    elif 'ID' in df_prot.columns:
        prot_samples = set(df_prot['ID'].astype(str))
    else:
        # Assume index is ID? Or first column?
        # Let's try first column if 'id' not found
        prot_samples = set(df_prot.iloc[:, 0].astype(str))
    print(f"Proteomics samples: {len(prot_samples)}")
except Exception as e:
    print(f"Error loading Proteomics: {e}")
    prot_samples = set()

# 3. Groups (needed for classification)
try:
    df_group = pd.read_csv(group_file)
    if 'id' in df_group.columns:
        group_samples = set(df_group['id'].astype(str))
    else:
        group_samples = set(df_group.iloc[:, 0].astype(str))
    print(f"Grouped samples: {len(group_samples)}")
except Exception as e:
    print(f"Error loading Groups: {e}")
    group_samples = set()

# Intersection
common_samples = mic_samples & prot_samples & group_samples
print(f"Common samples (Intersection): {len(common_samples)}")

# Save common samples
pd.DataFrame(list(common_samples), columns=['id']).to_csv(common_samples_file, index=False)
print(f"Saved common samples to {common_samples_file}")

# Plot Venn Diagram (Publication Quality)
plt.figure(figsize=(10, 8))
v = venn3([mic_samples, prot_samples, group_samples], ('Microbiota', 'Proteomics', 'Groups'), 
          set_colors=('#dc856d', '#7bc0b4', '#eca668'), alpha=0.8)

# Enhance font sizes for readability in journals
for text in v.set_labels:
    if text:
        text.set_fontsize(22)
        text.set_fontweight('bold')
        
for text in v.subset_labels:
    if text:
        text.set_fontsize(20)
        text.set_fontweight('bold')

plt.title("Sample Intersection Across Datasets", fontsize=26, fontweight='bold', pad=25)

plt.savefig(venn_file, dpi=600, bbox_inches='tight', transparent=True)
plt.close()
print(f"Saved Venn diagram to {venn_file}")



