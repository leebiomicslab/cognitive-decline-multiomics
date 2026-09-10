import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.outliers_influence import variance_inflation_factor

warnings.filterwarnings('ignore')

# High-impact journal aesthetics
sns.set_theme(style="ticks", context="paper")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# Paths
data_path = str(CLINICAL_CSV)
out_dir = FIGURES_DIR / "phenotyping"
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv(data_path)

# Map gender 2->0 if needed (1: male, 2: female -> 0)
if 'gender' in df.columns:
    df['gender'] = df['gender'].map({1: 1, 2: 0, '1': 1, '2': 0})

models = {
    "Model B (Proteomics)": ['age', 'gender', 'edu', 'bmi', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'Statin']
}

for model_name, features in models.items():
    print(f"Processing {model_name}...")
    df_vif = df[features].dropna()
    df_vif = df_vif.apply(pd.to_numeric, errors='coerce').dropna()
    
    if len(df_vif) > 10:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(df_vif)
        X_vif = pd.DataFrame(X_scaled, columns=df_vif.columns, index=df_vif.index)
        X_vif.insert(0, 'const', 1)

        vif_data = pd.DataFrame()
        vif_data['Feature'] = X_vif.columns
        vif_data['VIF'] = [variance_inflation_factor(X_vif.values, i) for i in range(X_vif.shape[1])]
        
        vif_data = vif_data[vif_data['Feature'] != 'const']
        vif_data = vif_data.sort_values(by='VIF', ascending=False)

        # Plot
        fig, ax = plt.subplots(figsize=(8, max(5, len(features)*0.35)))
        sns.barplot(data=vif_data, x='VIF', y='Feature', color='#2b6a99', ax=ax, width=0.7, edgecolor='none')
        ax.axvline(x=5, color='#d9534f', linestyle='--', linewidth=1.5, label='Threshold (VIF=5)')
        
        ax.set_title(f'Variance Inflation Factor (Model C)', fontweight='bold', fontsize=14, pad=15)
        ax.set_xlabel('VIF Score', fontweight='bold', fontsize=12)
        ax.set_ylabel('Clinical Variable', fontweight='bold', fontsize=12)
        
        # High impact adjustments
        sns.despine(trim=False)
        ax.tick_params(axis='both', which='major', labelsize=11)
        ax.xaxis.grid(True, linestyle=':', alpha=0.6)
        ax.yaxis.grid(False)
        ax.legend(frameon=False, loc='lower right')
        
        plt.tight_layout()
        save_name = 'model_c'
        plt.savefig(os.path.join(out_dir, f'clinical_vif_{save_name}.png'), dpi=600, bbox_inches='tight')
        plt.close()
    else:
        print(f"Not enough data for {model_name}")

print("Done!")


