"""
Purpose:          Nested CV Clinical vs Omics Comparison
                  Compares:
                  1. Clinical-only
                  2. Multi-omics Only (micro + prot)
                  3. Integrated (Clinical + Multi-omics)
                  Across C1 vs C2 and C1 vs C3 comparisons.
Manuscript:       Supplementary Table 5 — Classification Performance Metrics Across Feature Configurations
Figure/Table:     Supplementary Table 5
Input:            clinical_merged_371.csv, asv_genus_table.csv, data_log2_mapping_deduplicated.csv,
                  microbiota_da_results_final.csv, nominal_significant_proteins.csv
Output:           results/modeling/classification/feature/logistic/clinica/master_comparison.csv
Status:           EXACT PROVENANCE MATCH FOR SUPPLEMENTARY TABLE 5
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.metrics import (
    accuracy_score, f1_score, matthews_corrcoef,
    roc_auc_score, confusion_matrix
)
from pathlib import Path

# ─── Config ────────────────────────────────────────────────────────────────────
OUT_BASE = OUT_CLASSIFICATION_DIR / "feature" / "logistic" / "clinica"

CLINICAL     = CLINICAL_MERGED_FILE
MICRO_ASV    = ASV_GENUS_FILE
PROT_EXPR    = PROTEOMICS_MAPPING_FILE
MICRO_DA_RES = RESULTS_DIR / "microbiota_covariate_analysis" / "Model_A_Main" / "microbiota_da_results_final.csv"
PROT_DA_RES  = RESULTS_DIR / "proteomics_covariate_analysis" / "Main_Model" / "nominal_significant_proteins.csv"

# Clinical features to include
CLINICAL_COLS = [
    'age', 'gender', 'edu', 'bmi', 'smoking', 'alcohol', 
    'HTN', 'DM', 'CVD', 'CRP', 'egfr', 'Statin', 'PPI', 
    'Antibiotics', 'dm_med', 'htn_med', 'chol_med'
]

OUTER_SPLITS  = 5
OUTER_REPEATS = 20
INNER_SPLITS  = 4
RANDOM_STATE  = 42
CLR_EPS       = 1e-6

# ─── SignConsistencySelector ───────────────────────────────────────────────────
class SignConsistencySelector(BaseEstimator, TransformerMixin):
    def __init__(self, threshold=0.75, n_subsamples=20, random_state=42):
        self.threshold    = threshold
        self.n_subsamples = n_subsamples
        self.random_state = random_state

    def fit(self, X, y):
        rng   = np.random.default_rng(self.random_state)
        n     = X.shape[0]
        coefs = []
        for _ in range(self.n_subsamples):
            idx = rng.choice(n, size=int(n * 0.8), replace=False)
            if len(np.unique(y[idx])) < 2:
                continue
            lr = LogisticRegression(
                C=0.1, penalty='l2', solver='lbfgs',
                class_weight='balanced', max_iter=2000,
                random_state=self.random_state
            )
            lr.fit(X[idx], y[idx])
            coefs.append(lr.coef_.ravel())

        if len(coefs) == 0:
            self.selected_mask_ = np.ones(X.shape[1], dtype=bool)
            self.sign_consistency_ = np.ones(X.shape[1])
            return self

        coef_arr = np.array(coefs)
        mean_sign = np.sign(coef_arr.mean(axis=0))
        consistency = (np.sign(coef_arr) == mean_sign).mean(axis=0)
        self.sign_consistency_ = consistency
        self.selected_mask_    = consistency >= self.threshold
        if self.selected_mask_.sum() == 0:
            self.selected_mask_[np.argmax(consistency)] = True
        return self

    def transform(self, X):
        return X[:, self.selected_mask_]

    def get_selected_names(self, feature_names):
        return np.array(feature_names)[self.selected_mask_]

PARAM_GRID = {
    'selector__threshold': [0, 0.55, 0.65, 0.75],
    'clf__C':              [0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 1, 5, 10],
}

INNER_SCORING = 'roc_auc'

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial'],
                     'axes.labelweight': 'bold',
                     'axes.titleweight': 'bold'})
GROUP_COLORS = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'}

# ─── Helpers ───────────────────────────────────────────────────────────────────
def get_clr(df, eps=CLR_EPS):
    df_imp = df.replace(0, eps).astype(float)
    log_df = np.log(df_imp)
    return log_df.sub(log_df.mean(axis=1), axis=0)

def compute_features(clin_df, asv_tr, asv_te, prot_tr, prot_te, m_items, p_items, feature_set):
    c_tr = clin_df.loc[asv_tr.index, CLINICAL_COLS]
    c_te = clin_df.loc[asv_te.index, CLINICAL_COLS]

    asv_tr_clr = get_clr(asv_tr)
    asv_te_clr = get_clr(asv_te)
    m_valid = [x for x in m_items if x in asv_tr_clr.columns]
    m_tr, m_te = asv_tr_clr[m_valid], asv_te_clr[m_valid]

    p_valid = [x for x in p_items if x in prot_tr.columns]
    p_tr, p_te = prot_tr[p_valid], prot_te[p_valid]

    if feature_set == 'clinical':
        return c_tr, c_te
    elif feature_set == 'micro':
        return m_tr, m_te
    elif feature_set == 'prot':
        return p_tr, p_te
    elif feature_set == 'combined':
        return pd.concat([m_tr, p_tr], axis=1), pd.concat([m_te, p_te], axis=1)
    elif feature_set == 'clin_micro':
        return pd.concat([c_tr, m_tr], axis=1), pd.concat([c_te, m_te], axis=1)
    elif feature_set == 'clin_prot':
        return pd.concat([c_tr, p_tr], axis=1), pd.concat([c_te, p_te], axis=1)
    elif feature_set == 'clin_combined':
        return pd.concat([c_tr, m_tr, p_tr], axis=1), pd.concat([c_te, m_te, p_te], axis=1)
    else:
        raise ValueError(f"Unknown feature set: {feature_set}")

def binary_metrics_at_threshold(y_true, y_prob, threshold=0.5):
    y_pred = (y_prob >= threshold).astype(int)
    cm = confusion_matrix(y_true, y_pred)
    if cm.shape != (2, 2):
        return dict(ACC=np.nan, SEN=np.nan, SPE=np.nan, F1=np.nan, MCC=np.nan)
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
        f1_val = f1_score(y_true, y_pred, zero_division=0)
        if f1_val > best_f1:
            best_f1 = f1_val
            best_th = th
    return best_th, best_f1

# ─── Load Raw Data ─────────────────────────────────────────────────────────────
def main():
    print("=== Loading Data for Supplementary Table 5 Nested CV ===")
    clin_raw = pd.read_csv(CLINICAL).rename(columns={'id': 'SampleID', 'group': 'Group'})
    clin_raw['SampleID'] = clin_raw['SampleID'].astype(str)
    clin_features = clin_raw.set_index('SampleID')
    y_all = clin_raw[['SampleID', 'Group']].dropna().set_index('SampleID')

    asv_raw = pd.read_csv(MICRO_ASV)
    asv_raw.set_index(asv_raw.columns[0], inplace=True)
    asv_raw = asv_raw.loc[:, ~asv_raw.columns.duplicated()].T
    asv_raw.index = asv_raw.index.astype(str)

    prot_raw = pd.read_csv(PROT_EXPR)
    prot_raw = prot_raw.rename(columns={prot_raw.columns[0]: 'SampleID'})
    prot_raw['SampleID'] = prot_raw['SampleID'].astype(str)
    prot_raw.set_index('SampleID', inplace=True)

    common_all = y_all.index.intersection(asv_raw.index).intersection(prot_raw.index)
    micro_da = pd.read_csv(MICRO_DA_RES)
    m_items  = micro_da.loc[micro_da['Union_Sig'] == True, 'taxon'].unique()

    prot_da  = pd.read_csv(PROT_DA_RES)
    p_items  = prot_da['Protein'].unique()

    def run_nested_cv(pos_class, neg_class, feature_set):
        label    = f"{pos_class}_vs_{neg_class}"
        run_name = f"{label}__{feature_set}"
        out_dir  = OUT_BASE / label / feature_set
        out_dir.mkdir(parents=True, exist_ok=True)

        mask    = y_all.loc[common_all, 'Group'].isin([pos_class, neg_class])
        samples = np.array(common_all[mask])
        y_sub   = y_all.loc[samples, 'Group'].map({pos_class: 1, neg_class: 0}).values

        outer_cv = RepeatedStratifiedKFold(n_splits=OUTER_SPLITS, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE)
        inner_cv = StratifiedKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=RANDOM_STATE)

        fold_rows, coef_rows = [], []
        sel_counts = None
        total = OUTER_SPLITS * OUTER_REPEATS

        for fold_idx, (tr_idx, te_idx) in enumerate(outer_cv.split(samples, y_sub)):
            if (fold_idx + 1) % 20 == 0 or fold_idx == 0:
                print(f"    [{run_name}] fold {fold_idx+1}/{total}", flush=True)

            tr_ids, te_ids = samples[tr_idx].tolist(), samples[te_idx].tolist()
            y_tr, y_te = y_sub[tr_idx], y_sub[te_idx]

            X_tr_df, X_te_df = compute_features(
                clin_features, asv_raw.loc[tr_ids], asv_raw.loc[te_ids],
                prot_raw.loc[tr_ids], prot_raw.loc[te_ids],
                m_items, p_items, feature_set
            )
            X_te_df = X_te_df.reindex(columns=X_tr_df.columns)
            feature_names = list(X_tr_df.columns)
            X_tr, X_te = X_tr_df.values.astype(float), X_te_df.values.astype(float)

            if sel_counts is None: sel_counts = pd.Series(0, index=feature_names)

            pipe = Pipeline([
                ('scaler',   StandardScaler()),
                ('selector', SignConsistencySelector(random_state=RANDOM_STATE)),
                ('clf',      LogisticRegression(penalty='l2', solver='lbfgs', class_weight='balanced', max_iter=10000, random_state=RANDOM_STATE))
            ])
            grid = GridSearchCV(pipe, PARAM_GRID, cv=inner_cv, scoring=INNER_SCORING, n_jobs=-1)
            grid.fit(X_tr, y_tr)
            best = grid.best_estimator_

            selector = best.named_steps['selector']
            selected_names = selector.get_selected_names(feature_names)
            sel_counts[selected_names] += 1

            y_prob_tr = best.predict_proba(X_tr)[:, 1]
            best_th, _ = find_best_threshold(y_tr, y_prob_tr)
            y_prob_te = best.predict_proba(X_te)[:, 1]
            
            auc_val = roc_auc_score(y_te, y_prob_te) if len(np.unique(y_te)) > 1 else np.nan
            m = binary_metrics_at_threshold(y_te, y_prob_te, best_th)
            m['AUC'] = auc_val
            m['train_AUC'] = roc_auc_score(y_tr, y_prob_tr) if len(np.unique(y_tr)) > 1 else np.nan
            tr_m = binary_metrics_at_threshold(y_tr, y_prob_tr, best_th)
            m['train_ACC'], m['train_SEN'], m['train_SPE'], m['train_F1'], m['train_MCC'] = tr_m['ACC'], tr_m['SEN'], tr_m['SPE'], tr_m['F1'], tr_m['MCC']
            m['inner_best_cv_score'], m['best_C'], m['best_sign_threshold'] = grid.best_score_, grid.best_params_.get('clf__C'), grid.best_params_.get('selector__threshold')
            m['n_features_selected'], m['best_threshold'], m['fold'] = int(selector.selected_mask_.sum()), best_th, fold_idx
            fold_rows.append(m)

            clf_coef = best.named_steps['clf'].coef_.ravel()
            full_coef = np.zeros(len(feature_names))
            full_coef[selector.selected_mask_] = clf_coef
            coef_rows.append(pd.Series(full_coef, index=feature_names, name=f'fold_{fold_idx}'))

        fold_df = pd.DataFrame(fold_rows)
        fold_df.to_csv(out_dir / "fold_metrics.csv", index=False)
        metrics_cols = ['ACC', 'SEN', 'SPE', 'F1', 'MCC', 'AUC']
        summary = fold_df[metrics_cols].agg(['mean', 'std', 'median']).T
        summary.columns = ['Mean', 'Std', 'Med']
        summary.to_csv(out_dir / "summary_metrics.csv")
        
        sel_freq = (sel_counts / total).sort_values(ascending=False)
        sel_freq.to_csv(out_dir / "feature_selection_frequency.csv", header=True)

        coef_df = pd.DataFrame(coef_rows)
        ever_selected = sel_counts[sel_counts > 0].index
        coef_df_sel = coef_df[ever_selected]
        coef_stat = pd.DataFrame({
            'mean': coef_df_sel.mean(),
            'std': coef_df_sel.std(),
            'sign_consistency': (np.sign(coef_df_sel) == np.sign(coef_df_sel.mean())).mean(),
            'abs_mean': coef_df_sel.abs().mean(),
            'selection_frequency': sel_freq[ever_selected],
        }).sort_values('abs_mean', ascending=False)
        coef_stat.to_csv(out_dir / "coef_stability.csv")

        top_n = min(40, len(coef_stat))
        plot_df = coef_stat.head(top_n)
        fig, ax = plt.subplots(figsize=(8, max(4, len(plot_df) * 0.38)))
        y_pos = np.arange(len(plot_df))
        colors = ['#C0392B' if v > 0 else '#2980B9' for v in plot_df['mean']]
        ax.barh(y_pos, plot_df['mean'], xerr=plot_df['std'], color=colors, alpha=0.75, error_kw=dict(ecolor='black', lw=1, capsize=3))
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"{feat[:40]} (sc={sc:.0%}, sel={sf:.0%})" for feat, sc, sf in zip(plot_df.index, plot_df['sign_consistency'], plot_df['selection_frequency'])], fontsize=7)
        ax.axvline(0, color='black', lw=0.8, ls='--')
        ax.set_xlabel('Coefficient')
        ax.set_title(f"{run_name}\nCoef Stability", fontsize=10, fontweight='bold')
        plt.tight_layout()
        plt.savefig(out_dir / "coef_stability.png", dpi=600, bbox_inches='tight')
        plt.close()

        return summary

    comparisons  = [('C1', 'C2'), ('C1', 'C3')]
    feature_sets = ['clinical', 'combined', 'clin_combined']
    all_summaries = {}

    for pos, neg in comparisons:
        for fset in feature_sets:
            key = f"{pos}vs{neg}__{fset}"
            print(f"\nRunning: {key}")
            summary = run_nested_cv(pos, neg, fset)
            all_summaries[key] = summary

    rows = []
    for key, s in all_summaries.items():
        comp, fset = key.split('__')
        row = {'Comparison': comp, 'FeatureSet': fset}
        for met in ['ACC', 'SEN', 'SPE', 'F1', 'MCC', 'AUC']:
            row[f"{met}_Mean"], row[f"{met}_Std"] = s.loc[met, 'Mean'], s.loc[met, 'Std']
        rows.append(row)

    master_df = pd.DataFrame(rows)
    master_df.to_csv(OUT_BASE / "master_comparison.csv", index=False)
    print("\n[DONE] Supplementary Table 5 results saved to", OUT_BASE / "master_comparison.csv")

if __name__ == "__main__":
    main()
