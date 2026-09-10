"""
Purpose:          Pairwise logistic regression classification (C1 vs C2, C1 vs C3)
                  using domain-specific CISS ratio features derived from concordant
                  stable features across both comparisons.
                  Also constructs the Cognitive Impairment Signature Score (CISS)
                  for Clinical, Microbiome, and Proteomic domains.
Manuscript:       Methods — Classification and CISS construction
Figure/Table:     Figure 5a–f; Supplementary Table 5
Input:            clinical_merged_371.csv, asv_genus_table.csv (CLR-transformed),
                  data_log2_mapping_deduplicated.csv,
                  coef_stability.csv from upstream clinica/ and binary_advanced/ runs
Output:           results/modeling/classification/feature/logistic/ratio_modeling_cont_clin/
                    C1_vs_C{2,3}/{config}/summary_metrics.csv  (Supplementary Table 5)
                    sample_ratio_features.csv                   (Figure 5d, 5e, 5f data)
                    used_features.csv                           (Figure 5c data)
                    performance_bar_C1_vs_C{2,3}.png            (Figure 5a–b)
                    roc_auc_C1_vs_C{2,3}.png                   (Figure 5a–b)
Main dependencies: scikit-learn, pandas, numpy, matplotlib, seaborn

Cross-validation design:
  Outer: 5-fold x 20 repeats (RepeatedStratifiedKFold, RANDOM_STATE=42)
  Inner: 4-fold (StratifiedKFold, for hyperparameter tuning)
  Hyperparameter: C in {0.001, 0.01, 0.1, 1, 10, 100}, scored by AUC
  Threshold: F1-maximizing, selected on training-fold predictions only

Status: FINAL IMPLEMENTATION — Figure 5 / Supplementary Table 5 / CISS
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, roc_auc_score, confusion_matrix, roc_curve, auc

import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

# --- Paths (from config/paths.py) ---
OUT_BASE = MODELING_DIR / "classification" / "feature" / "logistic" / "ratio_modeling_cont_clin"
OUT_BASE.mkdir(parents=True, exist_ok=True)

CLINICAL  = CLINICAL_CSV
MICRO_ASV = MICRO_GENUS_COUNTS
PROT_EXPR = PROTEOMICS_LOG2

# Pre-computed feature stability dirs (produced by upstream classification runs)
CLIN_DIR = MODELING_DIR / "classification" / "feature" / "logistic" / "clinica"
ADV_DIR  = MODELING_DIR / "classification" / "feature" / "logistic" / "binary_advanced"

OUTER_SPLITS  = 5
OUTER_REPEATS = 20
INNER_SPLITS  = 4
RANDOM_STATE  = 42
CLR_EPS       = 1e-6

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'axes.labelweight': 'bold', 'axes.titleweight': 'bold'})
GROUP_COLORS = {'micro': '#1F77B4', 'clin': '#FF7F0E', 'prot': '#2CA02C', 
                'micro_clin': '#9467BD', 'micro_prot': '#8C564B', 'clin_prot': '#E377C2',
                'combined': '#D62728'}

# --- 1. Identify Stable Features ---
print("Identifying overlapping stable features from scatter plot logic...")
df1 = pd.read_csv(CLIN_DIR / "C1_vs_C2" / "clin_combined" / "coef_stability.csv").rename(columns={"Unnamed: 0": "Feature", "mean": "mean_c12"})
df2 = pd.read_csv(CLIN_DIR / "C1_vs_C3" / "clin_combined" / "coef_stability.csv").rename(columns={"Unnamed: 0": "Feature", "mean": "mean_c13"})
if "Feature" not in df1.columns: df1.rename(columns={df1.columns[0]: "Feature"}, inplace=True)
if "Feature" not in df2.columns: df2.rename(columns={df2.columns[0]: "Feature"}, inplace=True)
merged_df = pd.merge(df1, df2, on="Feature", suffixes=("_c12", "_c13"))

pos_features_all = merged_df[(merged_df["mean_c12"] > 0) & (merged_df["mean_c13"] > 0)]["Feature"].tolist()
neg_features_all = merged_df[(merged_df["mean_c12"] < 0) & (merged_df["mean_c13"] < 0)]["Feature"].tolist()

# Classify features
df_clin1 = pd.read_csv(CLIN_DIR / "C1_vs_C2" / "clinical" / "coef_stability.csv")
df_clin2 = pd.read_csv(CLIN_DIR / "C1_vs_C3" / "clinical" / "coef_stability.csv")
clinical_features = set(df_clin1.iloc[:, 0].unique()) | set(df_clin2.iloc[:, 0].unique())

df_micro1 = pd.read_csv(ADV_DIR / "C1_vs_C2" / "micro" / "coef_stability.csv")
df_micro2 = pd.read_csv(ADV_DIR / "C1_vs_C3" / "micro" / "coef_stability.csv")
microbiota_features = set(df_micro1.iloc[:, 0].unique()) | set(df_micro2.iloc[:, 0].unique())

pos_clin = [f for f in pos_features_all if f in clinical_features]
neg_clin = [f for f in neg_features_all if f in clinical_features]

CONTINUOUS_CLIN = ['age', 'edu', 'egfr', 'bmi']
pos_clin = [f for f in pos_clin if f in CONTINUOUS_CLIN]
neg_clin = [f for f in neg_clin if f in CONTINUOUS_CLIN]
pos_micro = [f for f in pos_features_all if f in microbiota_features]
neg_micro = [f for f in neg_features_all if f in microbiota_features]
pos_prot = [f for f in pos_features_all if f not in clinical_features and f not in microbiota_features]
neg_prot = [f for f in neg_features_all if f not in clinical_features and f not in microbiota_features]

print(f"Clin  Pos: {len(pos_clin)}, Neg: {len(neg_clin)}")
print(f"Micro Pos: {len(pos_micro)}, Neg: {len(neg_micro)}")
print(f"Prot  Pos: {len(pos_prot)}, Neg: {len(neg_prot)}")

# Output what features were used
feature_out = []
for f in pos_clin: feature_out.append({"Category": "Clin", "Sign": "Positive", "Feature": f})
for f in neg_clin: feature_out.append({"Category": "Clin", "Sign": "Negative", "Feature": f})
for f in pos_micro: feature_out.append({"Category": "Micro", "Sign": "Positive", "Feature": f})
for f in neg_micro: feature_out.append({"Category": "Micro", "Sign": "Negative", "Feature": f})
for f in pos_prot: feature_out.append({"Category": "Prot", "Sign": "Positive", "Feature": f})
for f in neg_prot: feature_out.append({"Category": "Prot", "Sign": "Negative", "Feature": f})
pd.DataFrame(feature_out).to_csv(OUT_BASE / "used_features.csv", index=False)

# --- 2. Load Raw Data Array ---
def get_clr(df, eps=CLR_EPS):
    df_imp = df.replace(0, eps).astype(float)
    log_df = np.log(df_imp)
    return log_df.sub(log_df.mean(axis=1), axis=0)

print("Loading raw data...")
clin = pd.read_csv(CLINICAL).rename(columns={'id': 'SampleID', 'group': 'Group'})
clin['SampleID'] = clin['SampleID'].astype(str)
y_all = clin[['SampleID', 'Group']].dropna().set_index('SampleID')

asv_raw = pd.read_csv(MICRO_ASV)
asv_raw.set_index(asv_raw.columns[0], inplace=True)
asv_raw = asv_raw.loc[:, ~asv_raw.columns.duplicated()].T
asv_raw.index = asv_raw.index.astype(str)
asv_clr = get_clr(asv_raw)

prot_raw = pd.read_csv(PROT_EXPR)
prot_raw = prot_raw.rename(columns={prot_raw.columns[0]: 'SampleID'})
prot_raw['SampleID'] = prot_raw['SampleID'].astype(str)
prot_raw.set_index('SampleID', inplace=True)

common_all = y_all.index.intersection(asv_clr.index).intersection(prot_raw.index)

# Merge all valid features
valid_features = set(pos_features_all + neg_features_all)
clin_cols = [c for c in clin.columns if c in valid_features]
X_clin = clin.set_index("SampleID")[clin_cols].loc[common_all]
X_micro = asv_clr[[c for c in asv_clr.columns if c in valid_features]].loc[common_all]
X_prot = prot_raw[[c for c in prot_raw.columns if c in valid_features]].loc[common_all]

X_master = pd.concat([X_clin, X_micro, X_prot], axis=1)
# Handle NAs for numeric data (fill with median)
X_master = X_master.fillna(X_master.median())

# Ensure all tracked features are in columns (create 0-cols if missing)
for f in pos_features_all + neg_features_all:
    if f not in X_master.columns:
        X_master[f] = 0.0

# --- Transformers ---
class DfStandardScaler(BaseEstimator, TransformerMixin):
    def fit(self, X, y=None):
        self.means = X.mean(axis=0)
        self.stds = X.std(axis=0).replace(0, 1)
        return self
    def transform(self, X):
        return (X - self.means) / self.stds

class RatioTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, config):
        self.config = config
    def fit(self, X, y=None):
        return self
    def transform(self, X):
        # X is DataFrame from DfStandardScaler
        out = pd.DataFrame(index=X.index)
        if self.config in ['micro', 'micro_clin', 'micro_prot', 'combined']:
             out['Micro_Ratio'] = X[pos_micro].mean(axis=1).fillna(0) - X[neg_micro].mean(axis=1).fillna(0)
        if self.config in ['clin', 'micro_clin', 'clin_prot', 'combined']:
             out['Clin_Ratio'] = X[pos_clin].mean(axis=1).fillna(0) - X[neg_clin].mean(axis=1).fillna(0)
        if self.config in ['prot', 'micro_prot', 'clin_prot', 'combined']:
             out['Prot_Ratio'] = X[pos_prot].mean(axis=1).fillna(0) - X[neg_prot].mean(axis=1).fillna(0)
        return out.astype(float).values

# Compute and output global sample ratios
scaler_global = DfStandardScaler().fit(X_master)
X_scaled_global = scaler_global.transform(X_master)
transformer_global = RatioTransformer('combined')
ratios_global = pd.DataFrame(transformer_global.transform(X_scaled_global), index=X_master.index, columns=['Micro_Ratio', 'Clin_Ratio', 'Prot_Ratio'])
ratios_global['Group'] = y_all.loc[ratios_global.index, 'Group'].values
ratios_global.to_csv(OUT_BASE / "sample_ratio_features.csv")

# Metrics logic
def binary_metrics_at_threshold(y_true, y_prob, threshold=0.5):
    y_pred = (y_prob >= threshold).astype(int)
    cm = confusion_matrix(y_true, y_pred)
    if cm.shape != (2, 2): return dict(ACC=np.nan, SEN=np.nan, SPE=np.nan, F1=np.nan, MCC=np.nan)
    tn, fp, fn, tp = cm.ravel()
    return {
        'ACC': accuracy_score(y_true, y_pred),
        'SEN': tp / (tp + fn) if (tp + fn) > 0 else 0.0,
        'SPE': tn / (tn + fp) if (tn + fp) > 0 else 0.0,
        'F1' : f1_score(y_true, y_pred, zero_division=0),
        'MCC': matthews_corrcoef(y_true, y_pred)
    }

def find_best_threshold(y_true, y_prob, n_thresholds=200):
    thresholds = np.linspace(0, 1, n_thresholds)
    best_th, best_f1 = 0.5, 0.0
    for th in thresholds:
        y_pred = (y_prob >= th).astype(int)
        f = f1_score(y_true, y_pred, zero_division=0)
        if f > best_f1: best_f1, best_th = f, th
    return best_th, best_f1

PARAM_GRID = {'clf__C': [0.001, 0.01, 0.1, 1, 10, 100]}

def run_pipeline(pos_class, neg_class, config):
    print(f"\n[Running] {pos_class} vs {neg_class} | {config}")
    
    mask = y_all.loc[common_all, 'Group'].isin([pos_class, neg_class])
    samples = np.array(common_all[mask])
    y_sub = y_all.loc[samples, 'Group'].map({pos_class: 1, neg_class: 0}).values
    X_sub = X_master.loc[samples]

    outer_cv = RepeatedStratifiedKFold(n_splits=OUTER_SPLITS, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE)
    inner_cv = StratifiedKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    fold_rows = []
    y_prob_aggregated = np.zeros(len(samples))
    y_prob_counts = np.zeros(len(samples))

    for fold_idx, (tr_idx, te_idx) in enumerate(outer_cv.split(samples, y_sub)):
        X_tr, X_te = X_sub.iloc[tr_idx], X_sub.iloc[te_idx]
        y_tr, y_te = y_sub[tr_idx], y_sub[te_idx]
        
        pipe = Pipeline([
            ('scaler', DfStandardScaler()),
            ('ratio', RatioTransformer(config=config)),
            ('clf', LogisticRegression(penalty='l2', solver='lbfgs', class_weight='balanced', max_iter=2000, random_state=RANDOM_STATE))
        ])
        
        grid = GridSearchCV(pipe, PARAM_GRID, cv=inner_cv, scoring='roc_auc', n_jobs=-1)
        grid.fit(X_tr, y_tr)
        best = grid.best_estimator_
        
        y_prob_tr = best.predict_proba(X_tr)[:, 1]
        best_th, _ = find_best_threshold(y_tr, y_prob_tr)
        
        y_prob_te = best.predict_proba(X_te)[:, 1]
        
        # Accumulate predictions for average ROC
        y_prob_aggregated[te_idx] += y_prob_te
        y_prob_counts[te_idx] += 1
        
        m = binary_metrics_at_threshold(y_te, y_prob_te, best_th)
        m['AUC'] = roc_auc_score(y_te, y_prob_te) if len(np.unique(y_te)) > 1 else np.nan
        m['fold'] = fold_idx
        fold_rows.append(m)
        
    y_prob_avg = y_prob_aggregated / y_prob_counts
    final_roc_auc = roc_auc_score(y_sub, y_prob_avg) if len(np.unique(y_sub)) > 1 else np.nan
    
    # Also calculate false positive rates and true positive rates for plotting
    fpr, tpr, _ = roc_curve(y_sub, y_prob_avg)
    
    fold_df = pd.DataFrame(fold_rows)
    metrics_cols = ['ACC', 'SEN', 'SPE', 'F1', 'MCC', 'AUC']
    summary = fold_df[metrics_cols].agg(['mean', 'std']).T
    summary.columns = ['Mean', 'Std']
    
    return summary, fpr, tpr, final_roc_auc

# Master Loop
comparisons = [('C1', 'C2'), ('C1', 'C3')]
configs = ['micro', 'clin', 'prot', 'micro_clin', 'micro_prot', 'clin_prot', 'combined']

all_summaries = {}
all_rocs = {}

for pos, neg in comparisons:
    all_rocs[f"{pos}_vs_{neg}"] = {}
    for cfg in configs:
        key = f"{pos}vs{neg}__{cfg}"
        sum_df, fpr, tpr, auc_val = run_pipeline(pos, neg, cfg)
        all_summaries[key] = sum_df
        all_rocs[f"{pos}_vs_{neg}"][cfg] = (fpr, tpr, auc_val)
        
        # Save summary CSV
        cfg_out = OUT_BASE / f"{pos}_vs_{neg}" / cfg
        cfg_out.mkdir(parents=True, exist_ok=True)
        sum_df.to_csv(cfg_out / "summary_metrics.csv")

# ─── 3. Visualization ──────────────────────────────────────────────────────────
print("\nCreating Visualization...")

for pos, neg in comparisons:
    # --- Performance Bar Plot ---
    metrics = ['AUC', 'ACC', 'SEN', 'SPE', 'F1', 'MCC']
    plot_data = []
    for cfg in configs:
        k = f"{pos}vs{neg}__{cfg}"
        df = all_summaries[k]
        for m in metrics:
            plot_data.append({'Config': cfg, 'Metric': m, 'Mean': df.loc[m, 'Mean'], 'Std': df.loc[m, 'Std']})
    
    pdf = pd.DataFrame(plot_data)
    
    plt.figure(figsize=(14, 6))
    sns.barplot(data=pdf, x='Metric', y='Mean', hue='Config', palette=GROUP_COLORS, capsize=0.1, errorbar=None, edgecolor='k')
    
    # Custom error bars overlay
    ax = plt.gca()
    num_metrics = len(metrics)
    num_hues = len(configs)
    bar_width = 0.8 / num_hues
    
    for i, p in enumerate(ax.patches):
        cfg_idx = i // num_metrics
        met_idx = i % num_metrics
        if cfg_idx < num_hues:
            val = p.get_height()
            x = p.get_x() + p.get_width() / 2
            # Retrieve std
            cfg_name = configs[cfg_idx]
            met_name = metrics[met_idx]
            std_val = pdf[(pdf['Config']==cfg_name) & (pdf['Metric']==met_name)]['Std'].values[0]
            ax.errorbar(x, val, yerr=std_val, color='black', capsize=3, elinewidth=1, linestyle='none')
            
    plt.title(f"Model Performance using Log-Ratio Features ({pos} vs {neg})", fontsize=15, pad=15, fontweight='bold')
    plt.ylabel("Score", fontweight='bold')
    plt.xlabel("")
    plt.ylim(0, 1.05)
    plt.legend(title='Configuration', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(OUT_BASE / f"performance_bar_{pos}_vs_{neg}.png", dpi=600, bbox_inches='tight')
    plt.close()

    # --- ROC AUC Overlay Plot ---
    plt.figure(figsize=(7, 7))
    for cfg in configs:
        fpr, tpr, roc_auc = all_rocs[f"{pos}_vs_{neg}"][cfg]
        plt.plot(fpr, tpr, lw=2, color=GROUP_COLORS[cfg], label=f"{cfg.upper()} (AUC = {roc_auc:.3f})")
        
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([-0.02, 1.0])
    plt.ylim([0.0, 1.02])
    plt.xlabel('False Positive Rate', fontweight='bold', fontsize=12)
    plt.ylabel('True Positive Rate', fontweight='bold', fontsize=12)
    plt.title(f'Log-Ratio Model ROC Curves\n({pos} vs {neg})', fontweight='bold', fontsize=14)
    plt.legend(loc="lower right", frameon=True, fontsize=10)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_BASE / f"roc_auc_{pos}_vs_{neg}.png", dpi=600, bbox_inches='tight')
    plt.close()

print(f"\n[DONE] Finished writing results and plots to {OUT_BASE}")
