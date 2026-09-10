import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Table 1 after propensity score matching
Manuscript:       Methods - PSM sensitivity
Figure/Table:     Supplementary Table S2
"""

import pandas as pd
import numpy as np
from scipy import stats
import os

psm_file = str(PSM_MATCHED_CSV)
cog_file = str(COGNITIVE_XLSX)
out_file = str(RESULTS_DIR)

os.makedirs(os.path.dirname(out_file), exist_ok=True)

df_psm = pd.read_csv(psm_file)
df_cog = pd.read_excel(cog_file)

if 'ID' in df_cog.columns:
    df_cog = df_cog.rename(columns={'ID': 'id'})

# Merge moca_1 as moca
df_psm = pd.merge(df_psm, df_cog[['id', 'moca_1']], on='id', how='left')
df_psm = df_psm.rename(columns={'moca_1': 'moca'})

def format_continuous(series):
    return f"{series.mean():.2f} 簣 {series.std():.2f}"

def format_categorical(series, val):
    count = (series == val).sum()
    total = len(series.dropna())
    pct = (count / total) * 100 if total > 0 else 0
    return f"{count} ({pct:.1f}%)"

def p_value_continuous(df, col):
    g1 = df[df['group'] == 'C1'][col].dropna()
    g2 = df[df['group'] == 'C2'][col].dropna()
    g3 = df[df['group'] == 'C3'][col].dropna()
    if len(g1)>0 and len(g2)>0 and len(g3)>0:
        # Use ANOVA
        _, p = stats.f_oneway(g1, g2, g3)
        return p
    return np.nan

def p_value_categorical(df, col):
    ct = pd.crosstab(df['group'], df[col])
    if ct.shape[0] > 1 and ct.shape[1] > 1:
        _, p, _, _ = stats.chi2_contingency(ct)
        return p
    return np.nan

def format_p(p):
    if pd.isna(p): return ""
    if p < 0.001: return "<0.001"
    return f"{p:.3f}"

rows = []
rows.append(["Variable", "Overall", "C1", "C2", "C3", "P-value"])

groups = ['C1', 'C2', 'C3']

# Helper to add section
def add_section(name):
    rows.append([f"[{name}]", "", "", "", "", ""])

# Helper to add continuous
def add_continuous(col, label):
    overall = format_continuous(df_psm[col].dropna())
    c1 = format_continuous(df_psm[df_psm['group'] == 'C1'][col].dropna())
    c2 = format_continuous(df_psm[df_psm['group'] == 'C2'][col].dropna())
    c3 = format_continuous(df_psm[df_psm['group'] == 'C3'][col].dropna())
    p = p_value_continuous(df_psm, col)
    rows.append([label, overall, c1, c2, c3, format_p(p)])

# Helper to add categorical
def add_categorical(col, label, val_to_show=1):
    overall = format_categorical(df_psm[col], val_to_show)
    c1 = format_categorical(df_psm[df_psm['group'] == 'C1'][col], val_to_show)
    c2 = format_categorical(df_psm[df_psm['group'] == 'C2'][col], val_to_show)
    c3 = format_categorical(df_psm[df_psm['group'] == 'C3'][col], val_to_show)
    p = p_value_categorical(df_psm, col)
    rows.append([f"{label} (Yes)", overall, c1, c2, c3, format_p(p)])

# Add Demographics
add_section("Section 1: Demographics")
add_continuous("age", "age")

# Gender requires special handling for Female/Male
p_gender = p_value_categorical(df_psm, "gender")
rows.append(["gender", "Missing: 0", "", "", "", format_p(p_gender)])
female_overall = format_categorical(df_psm["gender"], 2)
female_c1 = format_categorical(df_psm[df_psm['group'] == 'C1']["gender"], 2)
female_c2 = format_categorical(df_psm[df_psm['group'] == 'C2']["gender"], 2)
female_c3 = format_categorical(df_psm[df_psm['group'] == 'C3']["gender"], 2)
rows.append(["  Female", female_overall, female_c1, female_c2, female_c3, ""])

male_overall = format_categorical(df_psm["gender"], 1)
male_c1 = format_categorical(df_psm[df_psm['group'] == 'C1']["gender"], 1)
male_c2 = format_categorical(df_psm[df_psm['group'] == 'C2']["gender"], 1)
male_c3 = format_categorical(df_psm[df_psm['group'] == 'C3']["gender"], 1)
rows.append(["  Male", male_overall, male_c1, male_c2, male_c3, ""])

add_continuous("edu", "edu")
add_continuous("bmi", "bmi")
add_categorical("smoking", "smoking")
add_categorical("alcohol", "alcohol")

# Add Comorbidity
add_section("Section 2: Comorbidity")
add_categorical("HTN", "HTN")
add_categorical("DM", "DM")
add_categorical("CVD", "CVD")

# Add Systemic state
add_section("Section 3: Systemic state")
add_continuous("CRP", "CRP")
add_continuous("egfr", "egfr")

# Add Medication
add_section("Section 4: Medication (summary)")
add_categorical("Statin", "Statin")
add_categorical("dm_med", "dm_med")
add_categorical("htn_med", "htn_med")
add_categorical("chol_med", "chol_med")
add_categorical("PPI", "PPI")
add_categorical("Antibiotics", "Antibiotics")

# Add Cognitive profile
add_section("Section 5: Cognitive profile")
if 'moca' in df_psm.columns:
    add_continuous("moca", "moca")

# Write to CSV
import csv
with open(out_file, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerows(rows)

print(f"Table 1 generated successfully at: {out_file}")



