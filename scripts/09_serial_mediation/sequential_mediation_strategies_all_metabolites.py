import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from scipy.stats import norm
import warnings

warnings.filterwarnings('ignore')

# 1. PATHS
RATIO_CSV = CISS_RATIO_CSV
MET_CSV   = METABOLOMICS_CSV
# VASCULAR_CSV loaded from config/paths.py
CLIN_MERGED_CSV = CLINICAL_CSV
OUT_DIR   = str(MODELING_DIR)

os.makedirs(OUT_DIR, exist_ok=True)

# 2. LOAD DATA
ratio_df = pd.read_csv(RATIO_CSV, index_col=0)

# Load moca_1 from vascular
vasc_df = pd.read_csv(VASCULAR_CSV)
vasc_df = vasc_df[['id', 'moca_1']].set_index('id')

# Load rest of covariates from clinical_merged_371
clin_df = pd.read_csv(CLIN_MERGED_CSV)
req_clin_cols = ['id', 'age', 'gender', 'edu', 
                 'bmi', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 
                 'CRP', 'egfr', 'Statin', 'PPI', 'Antibiotics', 
                 'dm_med', 'htn_med', 'chol_med']
clin_df = clin_df[req_clin_cols].set_index('id')

# Combine them
clin_df = clin_df.join(vasc_df, how='inner')

# Load RAW metabolite matrix (rows=features, cols=samples)
met_raw = pd.read_csv(MET_CSV, index_col=0)
met_T_raw = met_raw.T.copy() # samples x features

# Define names
IAA    = "3-indoleacetic acid (Indole-3-acetic acid)"
IPA    = "3-indolepropionic acid (Indole-3-propionic acid)_2"
IS1    = "3-indoxyl sulfate (3-indoxylsulfuric acid)_1"
IS2    = "3-indoxyl sulfate (3-indoxylsulfuric acid)_2"
I3CA   = next((c for c in met_T_raw.columns if "carboxaldehyde" in c.lower()), "Indole-3-carboxaldehyde")
KYN    = "Kynurenine"
TRP    = "Tryptophan"
PAGLN  = "Phenylacetylglutamine (PAGln)"
PAGL1  = "Phenylacetylglycine (PAGly)_1"
PAGL2  = "Phenylacetylglycine (PAGly)_2"
ILA    = "Indolelactic acid"
GB1    = [c for c in met_T_raw.columns if "butyrobetaine_1" in c][0]
GB2    = [c for c in met_T_raw.columns if "butyrobetaine_2" in c][0]

# Compute averages on RAW data
met_T_raw["IS_avg"]       = met_T_raw[[IS1, IS2]].mean(1)
met_T_raw["PAGly_avg"]    = met_T_raw[[PAGL1, PAGL2]].mean(1)
met_T_raw["ILA_avg"]      = met_T_raw[ILA]
met_T_raw["TMAO_pre_avg"] = met_T_raw[[GB1, GB2]].mean(1)

PANELS = [
    (IAA,           "IAA"),
    (IPA,           "IPA"),
    (IS1,           "IS_1"),
    (IS2,           "IS_2"),
    ("IS_avg",      "IS_avg"),
    (I3CA,          "I3CA"),
    (KYN,           "KYN"),
    (TRP,           "Trp"),
    (PAGLN,         "PAGln"),
    (PAGL1,         "PAGly_1"),
    (PAGL2,         "PAGly_2"),
    ("PAGly_avg",   "PAGly_avg"),
    (ILA,           "ILA"),
    ("ILA_avg",     "ILA_avg"),
    (GB1,           "TMAO_pre_1"),
    (GB2,           "TMAO_pre_2"),
    ("TMAO_pre_avg","TMAO_pre_avg"),
]

met_features = [c[0] for c in PANELS]
met_labels = [c[1] for c in PANELS]

# Combine raw
met_selected_raw = met_T_raw[met_features].copy()
met_selected_raw.columns = met_labels

base_df = ratio_df.join(clin_df, how='inner').join(met_selected_raw, how='inner')
base_df.rename(columns={
    "Micro_Ratio": "microbiome_ratio",
    "Prot_Ratio": "protein_ratio",
    "moca_1": "MoCA",
    "edu": "education"
}, inplace=True)

# Treat categorical covariates
for cat_col in ['gender', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'Statin', 'PPI', 'Antibiotics', 'dm_med', 'htn_med', 'chol_med']:
    base_df[cat_col] = base_df[cat_col].astype(str)

# 3. MEDIATION FUNCTION
def serial_mediation_fast(df, x, m1, m2, y, covars=None, n_boot=5000, seed=42):
    if covars is None: covars = []
    np.random.seed(seed)
    
    df_used = pd.get_dummies(df[[x, m1, m2, y] + covars], drop_first=True, dtype=float)
    cov_cols = [c for c in df_used.columns if c not in [x, m1, m2, y]]
    
    X1 = sm.add_constant(df_used[[x] + cov_cols]).values
    Y1 = df_used[m1].values
    X2 = sm.add_constant(df_used[[x, m1] + cov_cols]).values
    Y2 = df_used[m2].values
    X3 = sm.add_constant(df_used[[x, m1, m2] + cov_cols]).values
    Y3 = df_used[y].values
    
    def fit_coef(X_mat, Y_vec, var_idx):
        return np.linalg.lstsq(X_mat, Y_vec, rcond=None)[0][var_idx]
        
    a1 = fit_coef(X1, Y1, 1)
    d  = fit_coef(X2, Y2, 1)
    a2 = fit_coef(X2, Y2, 2)
    c_prime = fit_coef(X3, Y3, 1)
    f  = fit_coef(X3, Y3, 2)
    b  = fit_coef(X3, Y3, 3)
    
    indirect_serial = a1 * a2 * b
    n = len(Y1)
    boot_serial = np.zeros(n_boot)
    
    for i in range(n_boot):
        idx = np.random.randint(0, n, n)
        ba1 = fit_coef(X1[idx], Y1[idx], 1)
        ba2 = fit_coef(X2[idx], Y2[idx], 2)
        bb  = fit_coef(X3[idx], Y3[idx], 3)
        boot_serial[i] = ba1 * ba2 * bb
        
    ci_low, ci_high = np.percentile(boot_serial, [2.5, 97.5])
    p = 2 * min((boot_serial <= 0).mean(), (boot_serial >= 0).mean())
    
    return {
        "a1_X_to_M1": a1,
        "a2_M1_to_M2": a2,
        "b_M2_to_Y": b,
        "c_prime_direct_X_to_Y": c_prime,
        "indirect_serial": indirect_serial,
        "CI95_lower": ci_low,
        "CI95_upper": ci_high,
        "p_value": p
    }

# 4. RUN ALL METABOLITES ACROSS 3 STRATEGIES & 4 COVARIATE MODELS
cov_models = {
    "Model_1_Core": ['age', 'gender', 'education', 'bmi', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD'],
    "Model_2_Sens_Inflam_Renal": ['age', 'gender', 'education', 'bmi', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'CRP', 'egfr'],
    "Model_3_Sens_Micro_Drug": ['age', 'gender', 'education', 'bmi', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'Statin', 'PPI', 'Antibiotics'],
    "Model_4_Sens_Chronic_Med": ['age', 'gender', 'education', 'bmi', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'dm_med', 'htn_med', 'chol_med']
}

# We will collect results for each strategy
results_raw = []
results_log2 = []
results_rint = []

print("Running mediation for all 17 metabolites across 3 strategies and 4 covariate models...")

for met_name in met_labels:
    # We will standardize continuous predictors (ratio, moca, age, edu, bmi, crp, egfr, met_col) for each model execution to keep scale comparable
    print(f"Processing {met_name}")
    
    for mod_name, cov_list in cov_models.items():
        # Prepare Base Data for this specific covariate model (so dropna only drops what's needed)
        cols_to_keep = ['microbiome_ratio', 'protein_ratio', 'MoCA', met_name] + cov_list
        df_base = base_df[cols_to_keep].dropna().copy()
        
        # Transformations
        s_raw = df_base[met_name]
        s_log2 = np.log2(s_raw + 1)
        
        ranks = s_raw.rank()
        p_vals = (ranks - 3/8) / (len(ranks) + 1/4)
        s_rint = norm.ppf(p_vals)
        
        df_base[met_name + "_raw"] = s_raw
        df_base[met_name + "_log2"] = s_log2
        df_base[met_name + "_rint"] = s_rint
        
        # Standardize continuous variables
        continuous_vars = ['microbiome_ratio', 'protein_ratio', 'MoCA', 'age', 'education', 'bmi']
        if 'CRP' in cov_list: continuous_vars.append('CRP')
        if 'egfr' in cov_list: continuous_vars.append('egfr')
        
        for strategy, col_suffix in [("raw", "_raw"), ("log2", "_log2"), ("rint", "_rint")]:
            df_strat = df_base.copy()
            met_col = met_name + col_suffix
            
            cols_to_std = continuous_vars + [met_col]
            scaler = StandardScaler()
            df_strat[cols_to_std] = scaler.fit_transform(df_strat[cols_to_std])
            
            # FORWARD
            res_fw = serial_mediation_fast(
                df_strat, x="microbiome_ratio", m1=met_col, m2="protein_ratio", y="MoCA", covars=cov_list
            )
            res_fw["Metabolite"] = met_name
            res_fw["Direction"] = "Forward (Micro->Met->Prot->MoCA)"
            res_fw["Covariate_Model"] = mod_name
            res_fw["N"] = len(df_strat)
            
            # REVERSE
            res_rv = serial_mediation_fast(
                df_strat, x="protein_ratio", m1=met_col, m2="microbiome_ratio", y="MoCA", covars=cov_list
            )
            res_rv["Metabolite"] = met_name
            res_rv["Direction"] = "Reverse (Prot->Met->Micro->MoCA)"
            res_rv["Covariate_Model"] = mod_name
            res_rv["N"] = len(df_strat)
            
            if strategy == "raw":
                results_raw.extend([res_fw, res_rv])
            elif strategy == "log2":
                results_log2.extend([res_fw, res_rv])
            elif strategy == "rint":
                results_rint.extend([res_fw, res_rv])

# 5. SAVE FILES
out_cols = ["Covariate_Model", "Metabolite", "Direction", "N", "a1_X_to_M1", "a2_M1_to_M2", "b_M2_to_Y", "c_prime_direct_X_to_Y", "indirect_serial", "CI95_lower", "CI95_upper", "p_value"]

df_out_raw = pd.DataFrame(results_raw)[out_cols]
csv_raw = os.path.join(OUT_DIR, "sequential_mediation_all_metabolites_raw.csv")
df_out_raw.to_csv(csv_raw, index=False)

df_out_log2 = pd.DataFrame(results_log2)[out_cols]
csv_log2 = os.path.join(OUT_DIR, "sequential_mediation_all_metabolites_log2.csv")
df_out_log2.to_csv(csv_log2, index=False)

df_out_rint = pd.DataFrame(results_rint)[out_cols]
csv_rint = os.path.join(OUT_DIR, "sequential_mediation_all_metabolites_rint.csv")
df_out_rint.to_csv(csv_rint, index=False)

print("\nSaved 3 strategy files:")
print(f"1. {csv_raw}")
print(f"2. {csv_log2}")
print(f"3. {csv_rint}")

