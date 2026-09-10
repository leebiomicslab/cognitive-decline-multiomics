"""
Purpose:          Pearson and Spearman correlation matrices for continuous and
                  categorical clinical covariates; VIF calculation for all
                  candidate model variables.
Manuscript:       Methods — Covariate collinearity assessment
Figure/Table:     Supplementary Fig. 4 (VIF barplot)
Input:            vascular.csv, vascular_with_medication.xlsx, common_samples.csv
Output:           results/figures/phenotyping/clinical_continuous_correlation.png
                  results/figures/phenotyping/clinical_categorical_overlap.png
                  results/figures/phenotyping/clinical_vif_barplot.png
                  results/tables/clinical_vif_values.csv
Main dependencies: pandas, numpy, matplotlib, seaborn, statsmodels, scikit-learn
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.outliers_influence import variance_inflation_factor

warnings.filterwarnings('ignore')
sns.set_theme(style="ticks", context="paper", font_scale=1.3)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# Paths (from config/paths.py)
out_dir = FIGURES_DIR / "phenotyping"
os.makedirs(out_dir, exist_ok=True)

print("Loading datasets...")
vascular = pd.read_csv(VASCULAR_CSV)
med = pd.read_excel(MEDICATION_XLSX)
common = pd.read_csv(COMMON_SAMPLES_CSV)

def normalize_id(df):
    if 'ID' in df.columns:
        df.rename(columns={'ID': 'id'}, inplace=True)
    if 'id' in df.columns:
        df['id'] = df['id'].astype(str).str.strip()
    return df

vascular = normalize_id(vascular)
med = normalize_id(med)
common = normalize_id(common)

def get_actual_cols(df, expected_cols):
    actual = {}
    for exp in expected_cols:
        found = False
        if exp in df.columns:
            actual[exp] = exp
            found = True
        else:
            for col in df.columns:
                if str(col).lower().strip() == exp.lower().strip():
                    actual[exp] = col
                    found = True
                    break
    return actual

expected_vas = ['age', 'gender', 'edu', 'bmi', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'CRP', 'egfr', 'moca_1']
actual_vas = get_actual_cols(vascular, expected_vas)

expected_med = ['Statin', 'dm_med', 'htn_med', 'chol_med', 'PPI', 'Antibiotics']
actual_med = get_actual_cols(med, expected_med)

df = pd.merge(vascular[['id'] + list(actual_vas.values())], med[['id'] + list(actual_med.values())], on='id', how='inner')
rename_dict = {v: k for k, v in actual_vas.items()}
rename_dict.update({v: k for k, v in actual_med.items()})
df.rename(columns=rename_dict, inplace=True)

if 'moca_1' in df.columns:
    df.rename(columns={'moca_1': 'moca'}, inplace=True)

df = df[df['id'].isin(common['id'])]

num_cols = ['age', 'edu', 'bmi', 'CRP', 'egfr', 'moca']
cat_cols = ['gender', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'Statin', 'dm_med', 'htn_med', 'chol_med', 'PPI', 'Antibiotics']

print(f"Total samples matched: {len(df)}")

for col in num_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

for col in cat_cols:
    s = pd.to_numeric(df[col], errors='coerce')
    df[col] = np.where(s.notna() & (s == s // 1), s.astype('Int64').astype(str), df[col])
    df[col] = df[col].replace(['nan', 'NaN', 'None', 'NA', '<NA>'], np.nan)

df_num = df.copy()
# Gender mapping back to zero indexing
df_num['gender'] = df_num['gender'].map({'1': 1, '2': 0, 1:1, 2:0})
for col in cat_cols:
    if col != 'gender':
        df_num[col] = pd.to_numeric(df_num[col], errors='coerce')

# 1. Continuous Variables Correlation Matrix
print("Generating Continuous Correlation Matrix...")
plt.figure(figsize=(8, 6))
corr_num = df_num[num_cols].corr(method='pearson')
mask_num = np.triu(np.ones_like(corr_num, dtype=bool))
sns.heatmap(corr_num, annot=True, cmap='coolwarm', vmin=-1, vmax=1, fmt='.2f', square=True, mask=mask_num,
            cbar_kws={'shrink': 0.8}, linewidths=0.5, linecolor='white')
plt.title('Continuous clinical variable correlation matrix', fontweight='bold', pad=15, fontsize=14)
plt.figtext(0.5, -0.05, "MoCA was included for descriptive assessment of cognitive relevance\nand was not included as a covariate in omics association models.",
            ha="center", fontsize=11, style='italic')
plt.savefig(str(out_dir / 'clinical_continuous_correlation.png'), dpi=600, bbox_inches='tight')
plt.close()

# 2. Categorical/Binary Variables Overlap Check
print("Generating Categorical Overlap Matrix...")
plt.figure(figsize=(10, 8))
corr_cat = df_num[cat_cols].corr(method='spearman')
mask_cat = np.triu(np.ones_like(corr_cat, dtype=bool))
sns.heatmap(corr_cat, annot=True, cmap='viridis', vmin=-1, vmax=1, fmt='.2f', square=True, mask=mask_cat,
            cbar_kws={'shrink': 0.8}, linewidths=0.5, linecolor='white')
plt.title('Binary covariate association matrix', fontweight='bold', pad=15, fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.savefig(str(out_dir / 'clinical_categorical_overlap.png'), dpi=600, bbox_inches='tight')
plt.close()

# 3. VIF Calculation
print("Calculating VIF...")
features_vif = num_cols + cat_cols
df_vif_data = df_num[features_vif].dropna()
print(f"Samples remaining for VIF after dropping NAs: {len(df_vif_data)}")

if len(df_vif_data) > 10:
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df_vif_data)
    X_vif = pd.DataFrame(X_scaled, columns=df_vif_data.columns, index=df_vif_data.index)
    X_vif.insert(0, 'const', 1)

    vif_data = pd.DataFrame()
    vif_data['Feature'] = X_vif.columns
    vif_data['VIF'] = [variance_inflation_factor(X_vif.values, i) for i in range(X_vif.shape[1])]

    vif_data = vif_data[vif_data['Feature'] != 'const']
    vif_data = vif_data.sort_values(by='VIF', ascending=False)

    plt.figure(figsize=(10, 6))
    sns.barplot(data=vif_data, x='VIF', y='Feature', palette='magma')
    plt.axvline(x=5, color='r', linestyle='--', label='Threshold (VIF=5)')
    plt.title('Variance Inflation Factor (VIF)', fontweight='bold', pad=15)
    plt.xlabel('VIF Score', fontweight='bold')
    plt.ylabel('Clinical Variable', fontweight='bold')
    plt.tight_layout()
    plt.savefig(str(out_dir / 'clinical_vif_barplot.png'), dpi=300, bbox_inches='tight')
    plt.close()

    os.makedirs(TABLES_DIR, exist_ok=True)
    vif_data.to_csv(str(TABLES_DIR / 'clinical_vif_values.csv'), index=False)
else:
    print("Not enough samples to calculate VIF reliably.")

print("All collinearity checks finished successfully!")
