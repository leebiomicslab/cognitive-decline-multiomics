import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
import warnings

warnings.filterwarnings('ignore')
np.random.seed(42)

# 1. LOAD DATA
RATIO_CSV = CISS_RATIO_CSV
MET_CSV   = METABOLOMICS_CSV
CLIN_CSV  = VASCULAR_CSV
OUT_DIR   = str(MODELING_DIR)

os.makedirs(OUT_DIR, exist_ok=True)

# Load ratio
ratio_df = pd.read_csv(RATIO_CSV, index_col=0) # Index is 'id'

# Load clinical
clin_df = pd.read_csv(CLIN_CSV)
clin_df = clin_df[['id', 'age', 'gender', 'edu', 'moca_1']].copy()
clin_df.set_index('id', inplace=True)

# Load metabolite and transform
met_raw = pd.read_csv(MET_CSV, index_col=0)
met_T = np.log2(met_raw + 1).T

# Prepare specific metabolites
IAA    = "3-indoleacetic acid (Indole-3-acetic acid)"
IPA    = "3-indolepropionic acid (Indole-3-propionic acid)_2"
IS1    = "3-indoxyl sulfate (3-indoxylsulfuric acid)_1"
IS2    = "3-indoxyl sulfate (3-indoxylsulfuric acid)_2"
I3CA   = "Indole-3-carboxaldehyde"
KYN    = "Kynurenine"
TRP    = "Tryptophan"
PAGLN  = "Phenylacetylglutamine (PAGln)"
PAGL1  = "Phenylacetylglycine (PAGly)_1"
PAGL2  = "Phenylacetylglycine (PAGly)_2"
ILA    = "Indolelactic acid"
GB1    = [c for c in met_T.columns if "butyrobetaine_1" in c][0]
GB2    = [c for c in met_T.columns if "butyrobetaine_2" in c][0]

met_T["IS_avg"]       = met_T[[IS1, IS2]].mean(1)
met_T["PAGly_avg"]    = met_T[[PAGL1, PAGL2]].mean(1)
met_T["ILA_avg"]      = met_T[ILA]
met_T["TMAO_pre_avg"] = met_T[[GB1, GB2]].mean(1)

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

met_selected = met_T[met_features].copy()
met_selected.columns = met_labels

# Base dataframe
base_df = ratio_df.join(clin_df, how='inner').join(met_selected, how='inner')
base_df.rename(columns={
    "Micro_Ratio": "microbiome_ratio",
    "Prot_Ratio": "protein_ratio",
    "moca_1": "MoCA",
    "edu": "education"
}, inplace=True)
base_df['gender'] = base_df['gender'].astype(str)

# 3. SERIAL MEDIATION FUNCTION (FAST NUMPY LSTSQ)
def serial_mediation_fast(df, x, m1, m2, y, covars=None, n_boot=5000, seed=42):
    if covars is None: covars = []
    np.random.seed(seed)
    
    # One-hot encode covariances if they are objects (like gender)
    df_used = pd.get_dummies(df[[x, m1, m2, y] + covars], drop_first=True, dtype=float)
    x_col = x
    m1_col = m1
    m2_col = m2
    y_col = y
    
    cov_cols = [c for c in df_used.columns if c not in [x, m1, m2, y]]
    
    X1 = sm.add_constant(df_used[[x_col] + cov_cols]).values
    Y1 = df_used[m1_col].values
    
    X2 = sm.add_constant(df_used[[x_col, m1_col] + cov_cols]).values
    Y2 = df_used[m2_col].values
    
    X3 = sm.add_constant(df_used[[x_col, m1_col, m2_col] + cov_cols]).values
    Y3 = df_used[y_col].values
    
    def fit_coef(X_mat, Y_vec, var_idx):
        # lstsq returns (x, res, rank, s)
        return np.linalg.lstsq(X_mat, Y_vec, rcond=None)[0][var_idx]
        
    a1 = fit_coef(X1, Y1, 1) # x is index 1
    d  = fit_coef(X2, Y2, 1) # x is index 1
    a2 = fit_coef(X2, Y2, 2) # m1 is index 2
    c_prime = fit_coef(X3, Y3, 1) # x is index 1
    f  = fit_coef(X3, Y3, 2) # m1 is index 2
    b  = fit_coef(X3, Y3, 3) # m2 is index 3
    
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

covariates = ["age", "gender", "education"]
results_list = []

print(f"Running FAST mediation analysis for 17 metabolites with N_boot=5000...")
for met_name in met_labels:
    cols_to_keep = ['microbiome_ratio', 'protein_ratio', 'MoCA', 'age', 'gender', 'education', met_name]
    df = base_df[cols_to_keep].dropna()
    
    cols_to_std = ['microbiome_ratio', 'protein_ratio', 'MoCA', 'age', 'education', met_name]
    scaler = StandardScaler()
    df[cols_to_std] = scaler.fit_transform(df[cols_to_std])
    
    print(f" > Processing {met_name} (N={len(df)})")
    
    # FORWARD
    res_fw = serial_mediation_fast(
        df, x="microbiome_ratio", m1=met_name, m2="protein_ratio", y="MoCA", covars=covariates
    )
    res_fw["Metabolite"] = met_name
    res_fw["Model"] = "Forward (Micro->Met->Prot->MoCA)"
    res_fw["N"] = len(df)
    results_list.append(res_fw)
    
    # REVERSE
    res_rv = serial_mediation_fast(
        df, x="protein_ratio", m1=met_name, m2="microbiome_ratio", y="MoCA", covars=covariates
    )
    res_rv["Metabolite"] = met_name
    res_rv["Model"] = "Reverse (Prot->Met->Micro->MoCA)"
    res_rv["N"] = len(df)
    results_list.append(res_rv)

out_df = pd.DataFrame(results_list)
out_cols = ["Metabolite", "Model", "N", "a1_X_to_M1", "a2_M1_to_M2", "b_M2_to_Y", "c_prime_direct_X_to_Y", "indirect_serial", "CI95_lower", "CI95_upper", "p_value"]
out_df = out_df[out_cols]

csv_path = os.path.join(OUT_DIR, "sequential_mediation_all_metabolites.csv")
out_df.to_csv(csv_path, index=False)
print(f"Finished! Results saved to {csv_path}")

sig_df = out_df[out_df['p_value'] < 0.05]
if len(sig_df) > 0:
    print("\n--- SIGNIFICANT MEDIATIONS ---")
    for _, row in sig_df.iterrows():
        print(f"{row['Metabolite']} ({row['Model']}): p={row['p_value']:.4f}, CI=[{row['CI95_lower']:.4f}, {row['CI95_upper']:.4f}]")
else:
    print("\n--- NO SIGNIFICANT MEDIATIONS FOUND (p < 0.05) ---")

