import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          PCA of cognitive indicator scores across groups
Manuscript:       Methods - Cognitive phenotyping
Figure/Table:     Fig. 1a
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Define paths
group_file = str(GROUP_CSV)
cognitive_file = str(COGNITIVE_XLSX)
output_dir = str(RESULTS_DIR)
tables_dir = os.path.join(output_dir, "tables")
plot_file = os.path.join(output_dir, "pca_cognitive_scores.png")
loadings_file = os.path.join(tables_dir, "pca_loadings.csv")

os.makedirs(output_dir, exist_ok=True)
os.makedirs(tables_dir, exist_ok=True)

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
plot_file = os.path.join(figures_dir, "pca_cognitive_scores.png")

# Select columns
# "trail", "cube", "clock_shape", "clock_number", "clock_time",
# "naming", "con1", "con2", "con3", "lang1", "lang2",
# "abstract","delay", "moca_1"
features = [
    "trail", "cube", "clock_shape", "clock_number", "clock_time",
    "naming", "con1", "con2", "con3", "lang1", "lang2",
    "abstract", "delay", "moca_1"
]

# Check existing columns
valid_features = [f for f in features if f in df.columns]
missing = [f for f in features if f not in df.columns]
if missing:
    print(f"Warning: Missing features: {missing}")

if not valid_features:
    print("No valid features for PCA.")
    exit(1)

# Prepare data for PCA
df_pca = df[valid_features].dropna()
print(f"Entries after dropping NaNs: {len(df_pca)}")

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df_pca)

# PCA
pca = PCA(n_components=2)
principal_components = pca.fit_transform(X_scaled)
pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

# Add group info back
# We need to map back to original indices to get correct groups because we dropped NaNs
pca_df['group'] = df.loc[df_pca.index, 'group'].values

# Calculate loadings
loadings = pd.DataFrame(
    pca.components_.T,
    columns=['PC1', 'PC2'],
    index=valid_features
)
print("Loadings:")
print(loadings)

# Save loadings
print(f"Saving loadings to {loadings_file}...")
loadings.to_csv(loadings_file)

# Plotting
# Colors: C1 #ee7a5b, C2 #eca362, C3 #6fcbbc
colors = {
    "C1": "#ee7a5b",
    "C2": "#eca362",
    "C3": "#6fcbbc"
}
hue_order = ["C1", "C2", "C3"]

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
    'legend.fontsize': 16,
    'legend.title_fontsize': 18
})

plt.figure(figsize=(10, 8))
ax = sns.scatterplot(
    data=pca_df,
    x='PC1',
    y='PC2',
    hue='group',
    hue_order=hue_order,
    palette=colors,
    s=100, # size
    alpha=0.8,
    edgecolor='w',
    linewidth=0.5
)

# Explained variance
var_explained = pca.explained_variance_ratio_
plt.xlabel(f"PC1 ({var_explained[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({var_explained[1]*100:.1f}%)")
plt.title("PCA of Cognitive Indicators")
plt.legend(title="Group", loc='best')

# Spines
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
print(f"Saving plot to {plot_file}...")
plt.savefig(plot_file, dpi=600)
print("Done.")



