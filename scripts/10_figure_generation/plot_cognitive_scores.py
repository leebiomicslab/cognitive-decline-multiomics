import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Cognitive subtest score distributions across C1/C2/C3
Manuscript:       Results - Cognitive phenotyping
Figure/Table:     Fig. 1b
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Define paths (Relative to where the script is run: TVGH_congnitive_3group/script)
# We assume the script is run from inside 'TVGH_congnitive_3group/script' or we adjust relative paths accordingly.
# Ideally, we use absolute paths or careful relative paths.
# Let's use os.path.dirname to be safe.
# base_dir is up 3 levels from script (script -> TVGH_cognitive_3group -> TVGH -> Desktop)
# actually data is in TVGH/data
# structure:
# TVGH/
#   data/
#   TVGH_congnitive_3group/
#     script/
#       plot_cognitive_scores.py
#     results/

group_file = str(GROUP_CSV)
cognitive_file = str(COGNITIVE_XLSX)
output_dir = str(RESULTS_DIR)
output_file = os.path.join(output_dir, "cognitive_scores_3group_comparison.png")

print(f"Project root: {project_root}")
print(f"Output file: {output_file}")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load data
print("Loading data...")
try:
    df_group = pd.read_csv(group_file)
    df_cog = pd.read_excel(cognitive_file)
except Exception as e:
    print(f"Error loading files: {e}")
    # Fallback for debugging if paths are weird
    print(f"Tried loading: {group_file} and {cognitive_file}")
    exit(1)

# Inspect columns
print("Group columns:", df_group.columns.tolist())
print("Cognitive columns:", df_cog.columns.tolist())

# Load common samples
common_samples_file = str(COMMON_SAMPLES_CSV)
if os.path.exists(common_samples_file):
    print(f"Loading common samples from {common_samples_file}...")
    common_samples = pd.read_csv(common_samples_file)['id'].astype(str).tolist()
else:
    print("Common samples file not found. Using all samples.")
    common_samples = []

# Merge data
if 'id' not in df_cog.columns and 'ID' in df_cog.columns:
    df_cog.rename(columns={'ID': 'id'}, inplace=True)

print("Merging data...")
df_merged = pd.merge(df_group[['id', 'group']], df_cog, on='id', how='inner')

# Filter for common samples
if common_samples:
    original_len = len(df_merged)
    df_merged = df_merged[df_merged['id'].astype(str).isin(common_samples)]
    print(f"Filtered to common samples: {len(df_merged)} (from {original_len})")

print(f"Merged entries: {len(df_merged)}")

figures_dir = os.path.join(output_dir, "figures")
os.makedirs(figures_dir, exist_ok=True)
output_file = os.path.join(figures_dir, "cognitive_scores_3group_comparison.png")

# Define items to plot
items_to_plot = [
    "trail", "cube", "clock_shape", "clock_number", "clock_time",
    "naming", "con1", "con2", "con3", "lang1", "lang2",
    "abstract", "delay", "orientation"
]

valid_items = [col for col in items_to_plot if col in df_merged.columns]
if not valid_items:
    print("No valid items to plot.")
    exit(1)

# Melt to long format
df_long = df_merged.melt(
    id_vars=["group"],
    value_vars=valid_items,
    var_name="Cognitive_Item",
    value_name="Score"
)

colors = {
    "C1": "#ee7a5b",
    "C2": "#eca362",
    "C3": "#6fcbbc"
}
hue_order = ["C1", "C2", "C3"]

# Plotting
print("Generating plot...")

# Set publication style preferences
# Bold fonts, larger size
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
# plt.rcParams['figure.titleweight'] = 'bold' # not always available in all mpl versions

sns.set_style("whitegrid")
# Update rcParams for font sizes - EVEN BIGGER
plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 24,
    'axes.labelsize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'legend.title_fontsize': 18
})

plt.figure(figsize=(18, 8))

ax = sns.barplot(
    data=df_long,
    x="Cognitive_Item",
    y="Score",
    hue="group",
    hue_order=hue_order,
    palette=colors,
    errorbar="sd",
    capsize=0.15,
    # edgecolor="black", # REMOVED edge color as requested
    # linewidth=1.5
)

# Keep axis spines thicker but maybe less intrusive if needed
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)

plt.xticks(rotation=45, ha="right")
plt.xlabel("Cognitive Subtest")
plt.ylabel("Mean Score")
plt.title("Cognitive Subtest Scores by Group")
plt.legend(title="Group", loc='upper left') # Moved to UPPER LEFT

# Make tick labels bold (explicitly if rcParams doesn't catch them all)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

plt.tight_layout()

# Save plot
print(f"Saving plot to {output_file}...")
plt.savefig(output_file, dpi=600) # Higher DPI for publication
print("Done.")



