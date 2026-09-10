import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
from scipy import stats

# Set style
# Set style
sns.set(style="whitegrid")
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 14
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

def rarefy_data(df, depth, seed=42):
    """
    Rarefy data to a specified depth.
    df: samples x features
    """
    np.random.seed(seed)
    rarefied_data = []

    for idx, row in df.iterrows():
        total_counts = row.sum()
        if total_counts < depth:
            # Drop samples with less reads than depth
            continue

        # Create a prob distribution
        prob = row / total_counts
        # Sample
        sampled = np.random.multinomial(depth, prob)
        rarefied_data.append(pd.Series(sampled, index=df.columns, name=idx))

    return pd.DataFrame(rarefied_data)

def calculate_alpha_diversity(df):
    """
    df: samples x features (counts)
    Returns DataFrame with alpha metrics
    """
    metrics = pd.DataFrame(index=df.index)

    # Observed OTUs (Richness)
    metrics['Observed'] = (df > 0).sum(axis=1)

    # Shannon Index (H)
    # H = -sum(p * ln(p))
    total_counts = df.sum(axis=1)
    props = df.div(total_counts, axis=0)
    # Avoid log(0)
    props_masked = props.replace(0, np.nan)
    metrics['Shannon'] = - (props_masked * np.log(props_masked)).sum(axis=1)

    # Chao1
    # S_chao1 = S_obs + n1^2 / (2*n2)
    # n1: number of singletons (count == 1)
    # n2: number of doubletons (count == 2)

    chao1_list = []
    for idx, row in df.iterrows():
        n1 = (row == 1).sum()
        n2 = (row == 2).sum()
        s_obs = (row > 0).sum()

        if n2 == 0:
            # Bias-corrected form or just n1*(n1-1)/(2*(n2+1))?
            # If n2=0, standard formula divides by zero.
            # Use bias-corrected: S_chao1 = S_obs + n1*(n1-1) / (2*(n2+1))
            chao1 = s_obs + (n1 * (n1 - 1)) / (2 * (n2 + 1))
        else:
            chao1 = s_obs + (n1 ** 2) / (2 * n2)
        chao1_list.append(chao1)
    metrics['Chao1'] = chao1_list

    # Pielou's Evenness (J)
    # J = H / ln(S)
    metrics['Pielou'] = metrics['Shannon'] / np.log(metrics['Observed'])
    # Handle cases where S=1 => ln(S)=0 => J=nan (should be 0 or 1? Usually undefined, but fit to 0 if H is 0)
    metrics.loc[metrics['Observed'] <= 1, 'Pielou'] = 0

    return metrics

def main():
    parser = argparse.ArgumentParser(description="Plot Alpha Diversity")
    parser.add_argument("--level", required=True, help="Taxonomic level (Class, Family, Genus, Phylum)")
    parser.add_argument("--input", required=True, help="Path to reads CSV file")
    parser.add_argument("--group", required=True, help="Path to grouping CSV")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    args = parser.parse_args()

    # Load Data
    print(f"Loading data from {args.input}")
    try:
        df = pd.read_csv(args.input, index_col=0)
    except Exception as e:
        print(f"Error loading {args.input}: {e}")
        return

    # Check if data needs transposition
    # Usually samples are columns in 'asv_X_table.csv'. transform to samples as rows.
    # Check if index contains taxonomic names (strings) and columns are Sample IDs (E...)
    if df.shape[1] > 0 and df.columns[0].startswith('E'):
        df = df.T

    # Load Grouping
    groups = pd.read_csv(args.group)
    # Expect columns: SampleID, Group
    # Map SampleID to index if needed
    if 'SampleID' not in groups.columns:
        # Assuming first column is SampleID
        groups.columns = ['SampleID', 'Group'] + list(groups.columns[2:])

    # Filter for samples in grouping
    common_samples =  list(set(df.index) & set(groups['SampleID']))
    print(f"Found {len(common_samples)} common samples.")

    df_filtered = df.loc[common_samples]
    groups_filtered = groups[groups['SampleID'].isin(common_samples)].set_index('SampleID')
    groups_filtered = groups_filtered.loc[common_samples] # Align order

    # Ensure data is numeric
    df_filtered = df_filtered.apply(pd.to_numeric, errors='coerce').fillna(0)

    # Filter out samples with 0 reads (total count)
    sample_counts = df_filtered.sum(axis=1)
    valid_samples = sample_counts[sample_counts > 0].index

    n_removed = len(common_samples) - len(valid_samples)
    if n_removed > 0:
        print(f"Removed {n_removed} samples with 0 reads.")

    df_filtered = df_filtered.loc[valid_samples]
    groups_filtered = groups_filtered.loc[valid_samples]

    if df_filtered.empty:
        print("Error: No valid samples remaining after filtering 0-read samples.")
        return

    # Rarefaction
    # Find min depth (at least some threshold, e.g. 1000, or min of samples)
    min_observed_depth = df_filtered.sum(axis=1).min()
    print(f"Minimum rarefaction depth: {min_observed_depth}")

    # If min depth is very low (e.g. < 500), warn and maybe drop those samples?
    # For now, just rarefy to min_depth.
    df_rarefied = rarefy_data(df_filtered, int(min_observed_depth))

    # Calculate Metrics
    metrics = calculate_alpha_diversity(df_rarefied)

    # Merge with groups
    plot_data = metrics.join(groups_filtered)

    # Plotting
    metrics_list = ['Observed', 'Shannon', 'Chao1', 'Pielou']
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    axes = axes.flatten()

    ㄏ
    group_order = ['C1', 'C2', 'C3']

    # Kruskal-Wallis Test
    stats_results = []

    for i, metric in enumerate(metrics_list):
        ax = axes[i]

        # Violin Plot with Box
        sns.violinplot(x='Group', y=metric, data=plot_data, order=group_order,
                       palette=colors, ax=ax, inner=None, cut=0)
        sns.boxplot(x='Group', y=metric, data=plot_data, order=group_order,
                    width=0.2, color='white', showfliers=False, ax=ax)
        # Strip plot for points
        sns.stripplot(x='Group', y=metric, data=plot_data, order=group_order,
                      color='black', alpha=0.3, size=2, ax=ax, jitter=True)

        ax.set_title(metric)
        ax.set_xlabel('')

        # Stats
        groups_data = [plot_data[plot_data['Group'] == g][metric].values for g in group_order]
        # Check if samples sufficient
        if all(len(g) > 0 for g in groups_data):
            h_val, p_val = stats.kruskal(*groups_data)
            stats_results.append(f"{metric}: p={p_val:.4e}")

            # Annotate p-value
            ax.text(0.5, 1.05, f"Kruskal-Wallis p={p_val:.4g}",
                    transform=ax.transAxes, ha='center', va='bottom')
        else:
            stats_results.append(f"{metric}: Insufficient data")

    plt.tight_layout()

    # Save
    out_file = os.path.join(args.output_dir, f"{args.level}_alpha_violins.png")
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {out_file}")

    # Save CSV
    # Re-route from results/figures back to results/tables
    tables_dir = os.path.abspath(os.path.join(args.output_dir, "..", "tables"))
    os.makedirs(tables_dir, exist_ok=True)

    csv_file = os.path.join(tables_dir, f"{args.level}_alpha_metrics.csv")
    plot_data.to_csv(csv_file)
    print(f"Saved metrics to {csv_file}")

if __name__ == "__main__":
    main()
