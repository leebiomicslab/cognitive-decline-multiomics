"""
Purpose:          1:1 greedy propensity score matching (PSM) of C2 and C3 to C1,
                  matching on: age, gender, edu, bmi, smoking, alcohol, HTN, DM, CVD.
                  Uses logistic regression to estimate propensity scores.
Manuscript:       Methods — Propensity score matching (sensitivity analysis)
Figure/Table:     Supplementary Figs. 7, 8
Input:            clinical_merged_371.csv
Output:           data/clinical_matched_psm.csv, results/psm/psm_balance_summary.txt
Main dependencies: pandas, numpy, scikit-learn, scipy
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy.stats import ttest_ind, chi2_contingency

out_dir = PSM_DIR
os.makedirs(out_dir, exist_ok=True)

print("Loading clinical data...")
df = pd.read_csv(CLINICAL_CSV)

# We want to match on: age, gender, edu, bmi, smoking, alcohol, HTN, DM, CVD
# Continuous: age, edu, bmi
# Categorical: gender, smoking, alcohol, HTN, DM, CVD
covariates = ['age', 'gender', 'edu', 'bmi', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD']
# ensure gender is dummy coded or treated properly. 
# in this dataset gender is 1 and 2, let's keep it as is, or one-hot encode.
# Logistic regression handles binary categorical natively if we just pass it in.

# Drop missing
df_clean = df.dropna(subset=['group'] + covariates).copy()

# C1 is reference
c1_df = df_clean[df_clean['group'] == 'C1'].copy()
c2_df = df_clean[df_clean['group'] == 'C2'].copy()
c3_df = df_clean[df_clean['group'] == 'C3'].copy()

print(f"Original group sizes - C1: {len(c1_df)}, C2: {len(c2_df)}, C3: {len(c3_df)}")

# Function to perform 1:1 matching without replacement
def match_groups(ref_df, target_df, covs):
    # Combine data
    combined = pd.concat([ref_df, target_df], ignore_index=True)
    # Target = 1 if in target_df, 0 if in ref_df
    y = np.concatenate([np.zeros(len(ref_df)), np.ones(len(target_df))])
    
    # Scale continuous and pass categorical as is (since they are binary/ordinal numbers)
    X = combined[covs].copy()
    
    # Standardize continuous variables for better logistic regression convergence
    scaler = StandardScaler()
    cont_vars = ['age', 'edu', 'bmi']
    X[cont_vars] = scaler.fit_transform(X[cont_vars])
    
    # Fit Logistic Regression
    lr = LogisticRegression(penalty=None, solver='lbfgs', max_iter=1000)
    lr.fit(X, y)
    
    # Calculate propensity scores
    ps = lr.predict_proba(X)[:, 1]
    
    # Logit of PS is better for matching
    eps = 1e-6
    ps = np.clip(ps, eps, 1 - eps)
    logit_ps = np.log(ps / (1 - ps))
    
    # Split back
    ref_logits = logit_ps[:len(ref_df)].reshape(-1, 1)
    target_logits = logit_ps[len(ref_df):].reshape(-1, 1)
    
    # Find nearest neighbors
    nn = NearestNeighbors(n_neighbors=1, algorithm='brute')
    nn.fit(target_logits)
    
    # Match greedily without replacement
    matched_target_indices = []
    distances = []
    
    # To do greedy matching properly: calculate pairwise distances
    # For simplicity, we just iterate through C1 and find the closest available in target
    available_target_indices = set(range(len(target_df)))
    
    import scipy.spatial.distance as distance
    dist_matrix = distance.cdist(ref_logits, target_logits, metric='euclidean')
    
    # Flatten and sort to match closest pairs first
    n_ref, n_target = dist_matrix.shape
    flat_indices = np.argsort(dist_matrix, axis=None)
    
    matched_ref = set()
    matched_tgt = set()
    
    ref_to_tgt_map = {}
    
    for idx in flat_indices:
        r = idx // n_target
        t = idx % n_target
        
        if r not in matched_ref and t not in matched_tgt:
            matched_ref.add(r)
            matched_tgt.add(t)
            ref_to_tgt_map[r] = t
            
        if len(matched_ref) == n_ref:
            break
            
    # Reconstruct matched target df
    target_matched_indices = [ref_to_tgt_map[i] for i in range(len(ref_df))]
    target_matched_df = target_df.iloc[target_matched_indices].copy()
    
    return target_matched_df

print("Matching C2 to C1...")
c2_matched = match_groups(c1_df, c2_df, covariates)

print("Matching C3 to C1...")
c3_matched = match_groups(c1_df, c3_df, covariates)

print(f"Matched group sizes - C1: {len(c1_df)}, C2: {len(c2_matched)}, C3: {len(c3_matched)}")

# Combine all matched
matched_cohort = pd.concat([c1_df, c2_matched, c3_matched], ignore_index=True)
output_path = PSM_MATCHED_CSV
matched_cohort.to_csv(str(output_path), index=False)
print(f"Saved matched cohort to {output_path}")

# Evaluate balance
def calculate_smd(g1, g2, is_continuous):
    if is_continuous:
        mean1, mean2 = g1.mean(), g2.mean()
        var1, var2 = g1.var(), g2.var()
        smd = (mean1 - mean2) / np.sqrt((var1 + var2) / 2)
    else:
        # For categorical, calculate Cohen's w or just use proportion differences
        p1, p2 = g1.mean(), g2.mean()
        # Ensure proportions are between 0 and 1
        # if categorical like gender is 1/2, convert to 0/1
        if set(g1.unique()).issubset({1, 2}):
            p1 = (g1 == 2).mean()
            p2 = (g2 == 2).mean()
        smd = (p1 - p2) / np.sqrt((p1*(1-p1) + p2*(1-p2))/2 + 1e-6)
    return abs(smd)

print("\n--- Balance Check (Standardized Mean Differences) ---")
print("Target is SMD < 0.1 for good balance, < 0.2 for acceptable balance.")

def print_balance(df_before, df_after, ref_grp='C1', cmp_grp='C2'):
    b_ref = df_before[df_before['group'] == ref_grp]
    b_cmp = df_before[df_before['group'] == cmp_grp]
    
    a_ref = df_after[df_after['group'] == ref_grp]
    a_cmp = df_after[df_after['group'] == cmp_grp]
    
    print(f"\n{ref_grp} vs {cmp_grp}:")
    print(f"{'Covariate':<10} | {'SMD Before':<15} | {'SMD After':<15}")
    print("-" * 45)
    for cov in covariates:
        is_cont = cov in ['age', 'edu', 'bmi']
        smd_b = calculate_smd(b_ref[cov], b_cmp[cov], is_cont)
        smd_a = calculate_smd(a_ref[cov], a_cmp[cov], is_cont)
        print(f"{cov:<10} | {smd_b:<15.4f} | {smd_a:<15.4f}")

print_balance(df_clean, matched_cohort, 'C1', 'C2')
print_balance(df_clean, matched_cohort, 'C1', 'C3')

# Write a basic Table 1 summary to text
summary_path = os.path.join(out_dir, "psm_balance_summary.txt")
with open(summary_path, 'w') as f:
    f.write("Matched Cohort Size: C1={}, C2={}, C3={}\n\n".format(len(c1_df), len(c2_matched), len(c3_matched)))
    for cov in covariates:
        f.write(f"--- {cov} ---\n")
        f.write(f"C1 mean (std): {matched_cohort[matched_cohort['group'] == 'C1'][cov].mean():.2f} ({matched_cohort[matched_cohort['group'] == 'C1'][cov].std():.2f})\n")
        f.write(f"C2 mean (std): {matched_cohort[matched_cohort['group'] == 'C2'][cov].mean():.2f} ({matched_cohort[matched_cohort['group'] == 'C2'][cov].std():.2f})\n")
        f.write(f"C3 mean (std): {matched_cohort[matched_cohort['group'] == 'C3'][cov].mean():.2f} ({matched_cohort[matched_cohort['group'] == 'C3'][cov].std():.2f})\n\n")

print(f"\nSaved balance summary to {summary_path}")
