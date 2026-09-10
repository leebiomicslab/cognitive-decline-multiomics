"""
Purpose:          Generate Table 1 — baseline characteristics for all participants
                  stratified by cognitive phenotype (C1/C2/C3).
Manuscript:       Methods — Participant characteristics
Figure/Table:     Supplementary Table 1
Input:            vascular.csv, vascular_with_medication.xlsx,
                  k-means_3group.csv, common_samples.csv
Output:           results/tables/Table1_Baseline_Characteristics.csv
Main dependencies: pandas, numpy, scipy
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import pandas as pd
import numpy as np
from scipy.stats import kruskal, chi2_contingency

# Paths (from config/paths.py)
out_dir = TABLES_DIR
os.makedirs(out_dir, exist_ok=True)

print("Loading datasets...")
vascular = pd.read_csv(VASCULAR_CSV)
med = pd.read_excel(MEDICATION_XLSX)
group = pd.read_csv(GROUP_CSV)
common = pd.read_csv(COMMON_SAMPLES_CSV)

def normalize_id(df):
    if 'ID' in df.columns:
        df.rename(columns={'ID': 'id'}, inplace=True)
    if 'id' in df.columns:
        df['id'] = df['id'].astype(str).str.strip()
    return df

vascular = normalize_id(vascular)
med = normalize_id(med)
group = normalize_id(group)
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
        if not found:
            print(f"WARNING: column '{exp}' not found in dataframe!")
    return actual

expected_vas = ['age', 'gender', 'edu', 'bmi', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'CRP', 'egfr', 'moca_1']
actual_vas = get_actual_cols(vascular, expected_vas)

expected_med = ['Statin', 'dm_med', 'htn_med', 'chol_med', 'PPI', 'Antibiotics']
actual_med = get_actual_cols(med, expected_med)

# Merge mappings
df = pd.merge(group[['id', 'group']], vascular[['id'] + list(actual_vas.values())], on='id', how='inner')
df = pd.merge(df, med[['id'] + list(actual_med.values())], on='id', how='inner')

# Rename back to standard keys
rename_dict = {v: k for k, v in actual_vas.items()}
rename_dict.update({v: k for k, v in actual_med.items()})
df.rename(columns=rename_dict, inplace=True)

# For Cognitive profile section, rename moca_1 to moca
if 'moca_1' in df.columns:
    df.rename(columns={'moca_1': 'moca'}, inplace=True)

print(f"Rows before specific sample intersect: {len(df)}")
df = df[df['id'].isin(common['id'])]
print(f"Total target samples for Table 1: {len(df)}")

num_cols = ['age', 'edu', 'bmi', 'CRP', 'egfr', 'moca']
cat_cols = ['gender', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'Statin', 'dm_med', 'htn_med', 'chol_med', 'PPI', 'Antibiotics']
binary_cols = ['smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'Statin', 'dm_med', 'htn_med', 'chol_med', 'PPI', 'Antibiotics']

for col in cat_cols:
    if col in df.columns:
        s = pd.to_numeric(df[col], errors='coerce')
        df[col] = np.where(s.notna() & (s == s // 1), s.astype('Int64').astype(str), df[col])
        df[col] = df[col].replace(['nan', 'NaN', 'None', 'NA', '<NA>'], np.nan)

# Apply value labels for gender
if 'gender' in df.columns:
    df['gender'] = df['gender'].map({'1': 'Male', '2': 'Female'}).fillna(df['gender'])

groups = ['C1', 'C2', 'C3']
results = []

def format_num(series):
    return f"{series.mean():.2f} ± {series.std():.2f}"

def format_cat(series, val):
    count = (series == val).sum()
    pct = count / len(series.dropna()) * 100 if len(series.dropna()) > 0 else 0
    return f"{count} ({pct:.1f}%)"

sections = {
    'Section 1: Demographics': ['age', 'gender', 'edu', 'bmi', 'smoking', 'alcohol'],
    'Section 2: Comorbidity': ['HTN', 'DM', 'CVD'],
    'Section 3: Systemic state': ['CRP', 'egfr'],
    'Section 4: Medication (summary)': ['Statin', 'dm_med', 'htn_med', 'chol_med', 'PPI', 'Antibiotics'],
    'Section 5: Cognitive profile': ['moca']
}

for sec_name, vars_list in sections.items():
    # Append section header
    results.append({
        'Variable': f"[{sec_name}]",
        'Overall': '', 'C1': '', 'C2': '', 'C3': '', 'P-value': ''
    })
    
    for var in vars_list:
        if var not in df.columns:
            continue
            
        row = {'Variable': var}
        valid_data = df.dropna(subset=[var, 'group'])
        if len(valid_data) == 0:
            continue

        # Numeric handling (Kruskal-Wallis)
        if var in num_cols:
            row['Overall'] = format_num(valid_data[var])
            group_data = []
            for g in groups:
                g_data = valid_data[valid_data['group'] == g][var]
                row[g] = format_num(g_data)
                group_data.append(g_data.values)
            
            if all(len(g) > 0 for g in group_data):
                try:
                    stat, p = kruskal(*group_data)
                    row['P-value'] = f"{p:.4e}" if p < 0.001 else f"{p:.3f}"
                except:
                    row['P-value'] = "NA"
            else:
                row['P-value'] = "NA"
            results.append(row)
            
        # Categorical handling (Chi-Square)
        elif var in cat_cols:
            unique_vals = valid_data[var].dropna().unique()
            unique_vals = sorted([v for v in unique_vals if str(v).lower() not in ['nan', 'none', 'na', '<na>']])
            
            crosstab = pd.crosstab(valid_data[var], valid_data['group'])
            try:
                chi2, p, dof, expected = chi2_contingency(crosstab)
                pval_str = f"{p:.4e}" if p < 0.001 else f"{p:.3f}"
            except:
                pval_str = "NA"
                
            if var in binary_cols:
                yes_val = next((v for v in unique_vals if str(v) in ['1', 'yes', 'Yes']), None)
                row['Variable'] = f"{var} (Yes)"
                if yes_val is not None:
                    row['Overall'] = format_cat(valid_data[var], yes_val)
                    for g in groups:
                        g_data = valid_data[valid_data['group'] == g][var]
                        row[g] = format_cat(g_data, yes_val)
                else:
                    row['Overall'] = "0 (0.0%)"
                    for g in groups:
                        row[g] = "0 (0.0%)"
                row['P-value'] = pval_str
                results.append(row)
            else:
                row['Overall'] = f"Missing: {len(df) - len(valid_data)}"
                for g in groups:
                    row[g] = ""
                row['P-value'] = pval_str
                results.append(row)
                
                for val in unique_vals:
                    sub_row = {'Variable': f"  {val}"}
                    sub_row['Overall'] = format_cat(valid_data[var], val)
                    for g in groups:
                        g_data = valid_data[valid_data['group'] == g][var]
                        sub_row[g] = format_cat(g_data, val)
                    sub_row['P-value'] = ""
                    results.append(sub_row)

res_df = pd.DataFrame(results, columns=['Variable', 'Overall', 'C1', 'C2', 'C3', 'P-value'])
out_file = str(out_dir / "Table1_Baseline_Characteristics.csv")
res_df.to_csv(out_file, index=False, encoding='utf-8-sig')
print(f"Table 1 generated successfully at {out_file}")
