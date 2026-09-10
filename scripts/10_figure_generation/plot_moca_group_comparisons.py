import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          MoCA subtest pairwise group comparisons
Manuscript:       Results - Cognitive phenotype
Figure/Table:     Supplementary Fig. 2
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from scipy import stats

cognitive_file = str(COGNITIVE_XLSX)
clinical_file = str(CLINICAL_CSV)
out_2g_dir = str(RESULTS_DIR)
out_3g_dir = str(RESULTS_DIR)

os.makedirs(out_2g_dir, exist_ok=True)
os.makedirs(out_3g_dir, exist_ok=True)

print("Loading data...")
df_cog = pd.read_excel(cognitive_file)
df_clin = pd.read_csv(clinical_file)

if 'ID' in df_cog.columns and 'id' not in df_cog.columns:
    df_cog = df_cog.rename(columns={'ID': 'id'})

print("Merging data...")
df_merged = pd.merge(df_cog, df_clin[['id', 'group']], on='id', how='inner')

items_to_plot = [
    "trail", "cube", "clock_shape", "clock_number", "clock_time",
    "naming", "con1", "con2", "con3", "lang1", "lang2",
    "abstract", "delay", "orientation"
]

moca_col = "moca_1"

# Prepare 2 groups (< 26 vs >= 26)
df_merged[moca_col] = pd.to_numeric(df_merged[moca_col], errors='coerce')
df_2g = df_merged.dropna(subset=[moca_col]).copy()
df_2g['MoCA_Group'] = np.where(df_2g[moca_col] >= 26, '>= 26', '< 26')

# Function to draw p-value bracket
def add_stat_annotation(ax, p_val, x1, x2, y, h):
    text = f'p = {p_val:.3f}' if p_val >= 0.001 else 'p < 0.001'
    if p_val < 0.05:
        fontweight = 'bold'
        color = '#D32F2F' # Red for significant
        if p_val < 0.001:
            text += ' ***'
        elif p_val < 0.01:
            text += ' **'
        else:
            text += ' *'
    else:
        fontweight = 'normal'
        color = 'black'
        text += ' (ns)'
        
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='black')
    ax.text((x1+x2)*.5, y+h + (h*0.3), text, ha='center', va='bottom', color=color, fontweight=fontweight, fontsize=14)

# Publication style
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 2.0
sns.set_style("ticks")
plt.rcParams.update({
    'font.size': 14, 'axes.titlesize': 18, 'axes.labelsize': 16,
    'xtick.labelsize': 14, 'ytick.labelsize': 14,
    'lines.linewidth': 2.5, 'lines.markersize': 8
})

# To handle sns warnings gracefully
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

print("Generating 2-group Bar plots (< 26 vs >= 26)...")
for item in items_to_plot:
    if item not in df_2g.columns: continue
    df_plot = df_2g[['MoCA_Group', item]].dropna()
    df_plot[item] = pd.to_numeric(df_plot[item], errors='coerce')
    df_plot = df_plot.dropna()
    
    g1 = df_plot[df_plot['MoCA_Group'] == '< 26'][item]
    g2 = df_plot[df_plot['MoCA_Group'] == '>= 26'][item]
    
    # Mann-Whitney U test
    if len(g1) > 0 and len(g2) > 0:
        stat, p_val = stats.mannwhitneyu(g1, g2, alternative='two-sided')
    else:
        p_val = 1.0
        
    plt.figure(figsize=(8, 6))
    
    palette = {'< 26': '#F44336', '>= 26': '#4CAF50'}
    order = ['< 26', '>= 26']
    
    # Check what kwargs are supported by sns.barplot based on version
    # Using errorbar='sd' instead of ci='sd' which is modern seaborn, capsize for standard deviation whiskers
    ax = sns.barplot(
        data=df_plot, x='MoCA_Group', y=item, order=order,
        palette=palette, width=0.5, errorbar="sd", capsize=0.1,
        edgecolor="black", linewidth=2.0, err_kws={'linewidth': 2.0, 'color': 'black'}
    )
    
    # Set xticklabels to include N
    counts = df_plot['MoCA_Group'].value_counts()
    labels = [f"{grp}\n(N={counts.get(grp, 0)})" for grp in order]
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, fontweight='bold')
    
    # Calculate max height including error bars (mean + std)
    group_stats = df_plot.groupby('MoCA_Group')[item].agg(['mean', 'std']).fillna(0)
    max_err_bar = (group_stats['mean'] + group_stats['std']).max()
    raw_max = df_plot[item].max()
    y_max = max(max_err_bar, raw_max)
    
    y_min = 0 # Bar plots should typically start at 0
    y_range = y_max if y_max > 0 else 1
    
    y_bracket = y_max + y_range * 0.08
    h = y_range * 0.05
    
    add_stat_annotation(ax, p_val, 0, 1, y_bracket, h)
    
    sns.despine(top=True, right=True)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.tick_params(width=2, length=6)
    
    ax.set_ylim(0, y_bracket + h + y_range * 0.2)
    
    plt.xlabel("Total MoCA Score Category")
    plt.ylabel(f"{item.capitalize()} Subtest Score")
    plt.title(f"{item.capitalize()} Score: MoCA < 26 vs >= 26", pad=20)
    
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
        
    plt.tight_layout()
    plt.savefig(os.path.join(out_2g_dir, f"{item}_2groups.png"), dpi=600, bbox_inches='tight')
    plt.close()

print("Generating 3-group Bar plots (C1, C2, C3)...")
for item in items_to_plot:
    if item not in df_merged.columns: continue
    df_plot = df_merged[['group', item]].dropna()
    df_plot = df_plot[df_plot['group'].isin(['C1', 'C2', 'C3'])]
    df_plot[item] = pd.to_numeric(df_plot[item], errors='coerce')
    df_plot = df_plot.dropna()
    
    g1 = df_plot[df_plot['group'] == 'C1'][item]
    g2 = df_plot[df_plot['group'] == 'C2'][item]
    g3 = df_plot[df_plot['group'] == 'C3'][item]
    
    # Kruskal-Wallis H test
    if len(g1) > 0 and len(g2) > 0 and len(g3) > 0:
        stat, p_val = stats.kruskal(g1, g2, g3)
    else:
        p_val = 1.0
        
    plt.figure(figsize=(9, 6))
    
    palette = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'}
    order = ['C1', 'C2', 'C3']
    
    ax = sns.barplot(
        data=df_plot, x='group', y=item, order=order,
        palette=palette, width=0.5, errorbar="sd", capsize=0.1,
        edgecolor="black", linewidth=2.0, err_kws={'linewidth': 2.0, 'color': 'black'}
    )
    
    counts = df_plot['group'].value_counts()
    labels = [f"{grp}\n(N={counts.get(grp, 0)})" for grp in order]
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, fontweight='bold')
    
    group_stats = df_plot.groupby('group')[item].agg(['mean', 'std']).fillna(0)
    max_err_bar = (group_stats['mean'] + group_stats['std']).max()
    raw_max = df_plot[item].max()
    y_max = max(max_err_bar, raw_max)
    
    y_min = 0
    y_range = y_max if y_max > 0 else 1
    
    text = f'Kruskal-Wallis p = {p_val:.3f}' if p_val >= 0.001 else 'Kruskal-Wallis p < 0.001'
    if p_val < 0.05:
        color = '#D32F2F'
        if p_val < 0.001: text += ' ***'
        elif p_val < 0.01: text += ' **'
        else: text += ' *'
    else:
        color = 'black'
        text += ' (ns)'
        
    y_bracket = y_max + y_range * 0.1
    ax.text(1, y_bracket, text, ha='center', va='bottom', color=color, fontweight='bold', fontsize=14)
    
    sns.despine(top=True, right=True)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.tick_params(width=2, length=6)
    ax.set_ylim(0, y_bracket + y_range * 0.2)
    
    plt.xlabel("Group")
    plt.ylabel(f"{item.capitalize()} Subtest Score")
    plt.title(f"{item.capitalize()} Score by 3 Groups", pad=20)
    
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
        
    plt.tight_layout()
    plt.savefig(os.path.join(out_3g_dir, f"{item}_3groups.png"), dpi=600, bbox_inches='tight')
    plt.close()

print("Successfully finished drawing all plots.")



