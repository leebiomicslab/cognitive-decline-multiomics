"""
Purpose:          Merge clinical vascular data, medication data, and cognitive
                  phenotype assignments into a single metadata file used as
                  primary covariate input for microbiome differential abundance.
Manuscript:       Methods — Clinical covariates
Figure/Table:     Prerequisite step (generates clinical_merged_371.csv)
Input:            vascular.csv, vascular_with_medication.xlsx,
                  k-means_3group.csv, common_samples.csv
Output:           data/clinical_merged_371.csv
Main dependencies: pandas
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import pandas as pd

vasc = pd.read_csv(VASCULAR_CSV)
med = pd.read_excel(MEDICATION_XLSX)
group = pd.read_csv(GROUP_CSV)
common = pd.read_csv(COMMON_SAMPLES_CSV)

for df in [vasc, med, group, common]:
    if 'ID' in df.columns: df.rename(columns={'ID': 'id'}, inplace=True)
    if 'id' in df.columns: df['id'] = df['id'].astype(str).str.strip()

target_vasc = ['id', 'age', 'gender', 'edu', 'bmi', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'CRP', 'egfr']
target_med = ['id', 'Statin', 'PPI', 'Antibiotics', 'dm_med', 'htn_med', 'chol_med']

vasc_cols = {col.lower().strip(): col for col in vasc.columns}
cv_vasc = [vasc_cols[x.lower().strip()] for x in target_vasc if x.lower().strip() in vasc_cols]

med_cols = {col.lower().strip(): col for col in med.columns}
cv_med = [med_cols[x.lower().strip()] for x in target_med if x.lower().strip() in med_cols]

clin = pd.merge(group[['id', 'group']], vasc[cv_vasc], on='id', how='inner')
clin = pd.merge(clin, med[cv_med], on='id', how='inner')

rename_mapping = {vc: target_vasc[i] for i, vc in enumerate(cv_vasc)}
rename_mapping.update({mc: target_med[i] for i, mc in enumerate(cv_med)})
clin.rename(columns=rename_mapping, inplace=True)

clin = clin[clin['id'].isin(common['id'])]

cat_cols = ['gender', 'smoking', 'alcohol', 'HTN', 'DM', 'CVD', 'Statin', 'PPI', 'Antibiotics', 'dm_med', 'htn_med', 'chol_med']
num_cols = ['age', 'edu', 'bmi', 'CRP', 'egfr']

for c in cat_cols + num_cols:
    clin[c] = pd.to_numeric(clin[c], errors='coerce')
    
# Flatten entirely without dropping NAs to allow R sequential models to purge selectively
clin.to_csv(str(CLINICAL_CSV), index=False)
print("Consolidated metadata saved to:", CLINICAL_CSV)
