import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          MoCA score boxplot across C1/C2/C3 with Mann-Whitney U
Manuscript:       Results - Cognitive phenotype
Figure/Table:     Fig. 1c
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import mannwhitneyu

# Define paths
group_file = str(GROUP_CSV)
cognitive_file = str(COGNITIVE_XLSX)
output_dir = str(RESULTS_DIR)
output_file = os.path.join(output_dir, "moca_1_boxplot.png")

os.makedirs(output_dir, exist_ok=True)

# Load data
print("Loading data...")
try:
    df_group = pd.read_csv(group_file)
    df_cog = pd.read_excel(cognitive_file)
except Exception as e:
    print(f"Error loading files: {e}")
    exit(1)

# Load common samples
common_samples_file = str(COMMON_SAMPLES_CSV)
if os.path.exists(common_samples_file):
    print(f"Loading common samples from {common_samples_file}...")
    common_samples = pd.read_csv(common_samples_file)['id'].astype(str).tolist()
else:
    print("Common samples file not found. Using all samples.")
    common_samples = []

# Merge
if 'id' not in df_cog.columns and 'ID' in df_cog.columns:
    df_cog.rename(columns={'ID': 'id'}, inplace=True)

df = pd.merge(df_group[['id', 'group']], df_cog, on='id', how='inner')

# Filter for common samples
if common_samples:
    original_len = len(df)
    df = df[df['id'].astype(str).isin(common_samples)]
    print(f"Filtered to common samples: {len(df)} (from {original_len})")

print(f"Merged entries: {len(df)}")

figures_dir = os.path.join(output_dir, "figures")
os.makedirs(figures_dir, exist_ok=True)
output_file = os.path.join(figures_dir, "moca_1_boxplot.png")

target_col = 'moca_1'
if target_col not in df.columns:
    print(f"Column '{target_col}' not found. Available columns: {df.columns.tolist()}")
    exit(1)

# Style
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 24,
    'axes.labelsize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16
})

colors = {
    "C1": "#ee7a5b",
    "C2": "#eca362",
    "C3": "#6fcbbc"
}
order = ["C1", "C2", "C3"]

# Plot
plt.figure(figsize=(10, 8))
ax = sns.boxplot(
    data=df,
    x="group",
    y=target_col,
    order=order,
    palette=colors,
    width=0.6,
    linewidth=2
)
# Add points
sns.stripplot(
    data=df,
    x="group",
    y=target_col,
    order=order,
    color='black',
    alpha=0.3,
    jitter=True,
    size=5,
    ax=ax
)

# Statistics
print("Calculating statistics (Mann-Whitney U)...")
pairs = [("C1", "C2"), ("C1", "C3"), ("C2", "C3")]
p_values = {}

for g1, g2 in pairs:
    d1 = df[df["group"] == g1][target_col].dropna()
    d2 = df[df["group"] == g2][target_col].dropna()
    stat, p = mannwhitneyu(d1, d2, alternative='two-sided')
    p_values[(g1, g2)] = p
    print(f"{g1} vs {g2}: p={p:.4f}")

# Annotate significance
def label_diff(ax, data, column, group_col, order, p_values):
    # Find max value to determine y-position
    y_max = data[column].max()
    y_range = data[column].max() - data[column].min()

    # Define vertical offsets
    h = y_range * 0.05
    y_offset = y_max + h
    step = y_range * 0.1

    # Sort pairs by distance or just fixed order to stack them
    # Simple stacking
    current_y = y_offset

    for i, (g1, g2) in enumerate(pairs):
        p = p_values[(g1, g2)]

        # Determine significance marker
        if p >= 0.05:
            sig = "ns"
        elif p < 0.001:
            sig = "***"
        elif p < 0.01:
            sig = "**"
        elif p < 0.05:
            sig = "*"

        # Only plot significant or requested? User asked to "mark significance".
        # Usually checking all is good.

        # x coordinates
        x1 = order.index(g1)
        x2 = order.index(g2)

        # Draw line
        ax.plot([x1, x1, x2, x2], [current_y, current_y+h/2, current_y+h/2, current_y], lw=1.5, c='k')

        # Add text
        ax.text((x1+x2)*.5, current_y+h/2, sig, ha='center', va='bottom', color='k', fontsize=14, fontweight='bold')

        current_y += step

label_diff(ax, df, target_col, "group", order, p_values)

# Final adjustments
ax.set_xlabel("Group")
ax.set_ylabel("MoCA Score")
ax.set_title("MoCA Score Comparison")
# Make spines thick
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
print(f"Saving to {output_file}...")
plt.savefig(output_file, dpi=600)
print("Done.")



