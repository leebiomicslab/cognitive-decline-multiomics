"""
Status:           LEGACY / EXPLORATORY ANALYSIS
                  Not used for final Figure 5 or final CISS results.
                  This script performed classification using module eigengene features.
                  The module eigengene approach was explored as an alternative to the
                  ratio-feature approach. Final results use ratio features constructed
                  by feature_ratio_modeling_cont_clin.py.
                  OUTER_SPLITS = 4 in this script (differs from final script's 5-fold);
                  this discrepancy was previously flagged as D1 but is now RESOLVED
                  because the final implementation correctly uses 5-fold outer CV.
"""
"""
Purpose:          Binary classification of cognitive phenotype groups (C1 vs C2,
                  C1 vs C3) using L2-regularised logistic regression. Input features
                  combine module eigengenes (Leiden clusters, PCA-derived) and
                  raw DA-feature expression levels from both proteomics and
                  microbiome. Sign-consistency feature selection applied in each
                  outer fold. Computes CISS (Cognitive Impairment Signature Score)
                  via per-sample classification probabilities.
Manuscript:       Methods ??Multi-omics classification and CISS computation
Figure/Table:     Fig. 6b, 6c; Supplementary Table 5
Input:            CLINICAL_CSV, MICRO_GENUS_COUNTS, MICRO_MODULE assignments,
                  PROTEOMICS_LOG2, PROT_MODULE assignments
Output:           results/modeling/classification/module_feature/<comparison>/
Main dependencies: pandas, numpy, scikit-learn, matplotlib, seaborn

# NOTE ??DISCREPANCY (D1):
#   Code: OUTER_SPLITS = 4 (nested CV, 4 outer folds x 20 repeats = 80 folds)
#   Manuscript describes: "5-fold x 20 repeats" ??awaiting author confirmation.
#   DO NOT change OUTER_SPLITS without explicit author sign-off.
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
from sklearn.decomposition import PCA
from pathlib import Path

# ??? Config ?????????????????????????????????????????????????????????????????????????????????
# Paths (from config/paths.py)
OUT_BASE = MODELING_DIR / "classification" / "module_feature"

CLINICAL     = CLINICAL_CSV
MICRO_ASV    = MICRO_GENUS_COUNTS
MICRO_MODULE = FASTSPAR_MICRO_DIR / "network_topology" / "leiden_clusters_24_taxa_G_all.csv"
PROT_EXPR    = PROTEOMICS_LOG2
PROT_MODULE  = FASTSPAR_PROT_DIR / "pearson" / "module_protein_assignment.csv"

OUTER_SPLITS  = 4
OUTER_REPEATS = 20
INNER_SPLITS  = 4
RANDOM_STATE  = 42
CLR_EPS       = 1e-6

# ??? SignConsistencySelector ???????????????????????????????????????????????????
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
                C=0.1, penalty='elasticnet', solver='saga', l1_ratio=0.5,
                class_weight='balanced', max_iter=10000,
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


# ??? Param Grid ?????????
PARAM_GRID = {
    'selector__threshold': [0.85, 0.90, 0.95],
    'clf__C':              [1e-4, 5e-4, 1e-3, 1e-2, 1e-1, 1]
}

INNER_SCORING = 'roc_auc'

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial'],
                     'axes.labelweight': 'bold',
                     'axes.titleweight': 'bold'})
GROUP_COLORS = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'}

# ??? Helpers ???????????????????????????????????????????????????????????????????
def get_clr(df, eps=CLR_EPS):
    df_imp = df.replace(0, eps).astype(float)
    log_df = np.log(df_imp)
    return log_df.sub(log_df.mean(axis=1), axis=0)

def fit_pca_me(module_df, tr_mat, te_mat, module_col, id_col, prefix):
    me_tr, me_te = {}, {}
    for cid in sorted(module_df[module_col].unique()):
        name = f"{prefix}{cid+1}"
        items = module_df.loc[module_df[module_col] == cid, id_col].tolist()
        valid = [x for x in items if x in tr_mat.columns]
        if not valid: continue
        X_tr = tr_mat.loc[:, valid]
        X_te = te_mat.loc[:, valid]
        pca = PCA(n_components=1).fit(X_tr)
        if pca.components_[0].mean() < 0:
            pca.components_ *= -1
        me_tr[name] = pd.Series(pca.transform(X_tr).ravel(), index=tr_mat.index)
        me_te[name] = pd.Series(pca.transform(X_te).ravel(), index=te_mat.index)
    if not me_tr:
        return pd.DataFrame(index=tr_mat.index), pd.DataFrame(index=te_mat.index)
    return pd.DataFrame(me_tr), pd.DataFrame(me_te)

def compute_features(asv_tr, asv_te, prot_tr, prot_te, micro_m, prot_m, feature_set):
    """Extract Raw features + Module Eigengenes."""
    asv_tr_clr = get_clr(asv_tr)
    asv_te_clr = get_clr(asv_te)

    # 1. Raw features
    m_items = micro_m['Taxon'].unique()
    m_valid = [x for x in m_items if x in asv_tr_clr.columns]
    p_items = prot_m['Protein'].unique()
    p_valid = [x for x in p_items if x in prot_tr.columns]

    m_tr_raw, m_te_raw = asv_tr_clr[m_valid], asv_te_clr[m_valid]
    p_tr_raw, p_te_raw = prot_tr[p_valid], prot_te[p_valid]

    # 2. Module Eigengenes
    m_tr_me, m_te_me = fit_pca_me(micro_m, asv_tr_clr, asv_te_clr, 'Leiden_Cluster', 'Taxon', 'M')
    p_tr_me, p_te_me = fit_pca_me(prot_m, prot_tr, prot_te, 'Module', 'Protein', 'P')

    # Drop M1, M2, P1, P4 (as per module analysis script)
    drop_cols = ['M1', 'M2', 'P1', 'P4']
    m_tr_me = m_tr_me.drop(columns=[c for c in drop_cols if c in m_tr_me.columns], errors='ignore')
    m_te_me = m_te_me.drop(columns=[c for c in drop_cols if c in m_te_me.columns], errors='ignore')
    p_tr_me = p_tr_me.drop(columns=[c for c in drop_cols if c in p_tr_me.columns], errors='ignore')
    p_te_me = p_te_me.drop(columns=[c for c in drop_cols if c in p_te_me.columns], errors='ignore')

    tr_parts, te_parts = [], []
    if feature_set in ['micro', 'combined']:
        tr_parts.extend([m_tr_me, m_tr_raw])
        te_parts.extend([m_te_me, m_te_raw])
    if feature_set in ['prot', 'combined']:
        tr_parts.extend([p_tr_me, p_tr_raw])
        te_parts.extend([p_te_me, p_te_raw])

    return pd.concat(tr_parts, axis=1), pd.concat(te_parts, axis=1)


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

# ??? Load Raw Data Once ????????????????????????????????????????????????????????
print("=== Loading Raw Data ===")
clin = pd.read_csv(CLINICAL).rename(columns={'id': 'SampleID', 'group': 'Group'})
clin['SampleID'] = clin['SampleID'].astype(str)
y_all = clin[['SampleID', 'Group']].dropna().set_index('SampleID')

asv_raw = pd.read_csv(MICRO_ASV)
asv_raw.set_index(asv_raw.columns[0], inplace=True)
asv_raw = asv_raw.loc[:, ~asv_raw.columns.duplicated()].T
asv_raw.index = asv_raw.index.astype(str)

prot_raw = pd.read_csv(PROT_EXPR)
prot_raw = prot_raw.rename(columns={prot_raw.columns[0]: 'SampleID'})
prot_raw['SampleID'] = prot_raw['SampleID'].astype(str)
prot_raw.set_index('SampleID', inplace=True)

common_all = y_all.index.intersection(asv_raw.index).intersection(prot_raw.index)
micro_m = pd.read_csv(MICRO_MODULE)
prot_m  = pd.read_csv(PROT_MODULE)

# ??? Nested CV Runner ??????????????????????????????????????????????????????????
def run_nested_cv(pos_class, neg_class, feature_set):
    label    = f"{pos_class}_vs_{neg_class}"
    run_name = f"{label}__{feature_set}"
    out_dir  = OUT_BASE / label / feature_set
    out_dir.mkdir(parents=True, exist_ok=True)

    mask    = y_all.loc[common_all, 'Group'].isin([pos_class, neg_class])
    samples = np.array(common_all[mask])
    y_sub   = y_all.loc[samples, 'Group'].map({pos_class: 1, neg_class: 0}).values

    outer_cv = RepeatedStratifiedKFold(
        n_splits=OUTER_SPLITS, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE)
    inner_cv = StratifiedKFold(
        n_splits=INNER_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    fold_rows  = []
    coef_rows  = []
    sel_counts = None
    total = OUTER_SPLITS * OUTER_REPEATS

    for fold_idx, (tr_idx, te_idx) in enumerate(outer_cv.split(samples, y_sub)):
        if (fold_idx + 1) % 20 == 0 or fold_idx == 0:
            print(f"    [{run_name}] fold {fold_idx+1}/{total}", flush=True)

        tr_ids = samples[tr_idx].tolist()
        te_ids = samples[te_idx].tolist()
        y_tr, y_te = y_sub[tr_idx], y_sub[te_idx]

        asv_tr  = asv_raw.loc[tr_ids]
        asv_te  = asv_raw.loc[te_ids]
        prot_tr = prot_raw.loc[tr_ids]
        prot_te = prot_raw.loc[te_ids]

        X_tr_df, X_te_df = compute_features(
            asv_tr, asv_te, prot_tr, prot_te,
            micro_m, prot_m, feature_set
        )
        X_te_df = X_te_df.reindex(columns=X_tr_df.columns)
        feature_names = list(X_tr_df.columns)
        X_tr = X_tr_df.values.astype(float)
        X_te = X_te_df.values.astype(float)

        if sel_counts is None:
            sel_counts = pd.Series(0, index=feature_names)

        pipe = Pipeline([
            ('scaler',   StandardScaler()),
            ('selector', SignConsistencySelector(random_state=RANDOM_STATE)),
            ('clf',      LogisticRegression(
                penalty='l2',
                solver='lbfgs',
                class_weight='balanced',
                max_iter=10000,
                random_state=RANDOM_STATE
            ))
        ])
        grid = GridSearchCV(pipe, PARAM_GRID, cv=inner_cv,
                            scoring=INNER_SCORING, n_jobs=-1)
        grid.fit(X_tr, y_tr)
        best = grid.best_estimator_

        selector: SignConsistencySelector = best.named_steps['selector']
        selected_names = selector.get_selected_names(feature_names)
        sel_counts[selected_names] += 1

        y_prob_tr = best.predict_proba(X_tr)[:, 1]
        best_th, _ = find_best_threshold(y_tr, y_prob_tr)

        y_prob_te = best.predict_proba(X_te)[:, 1]
        auc_val   = roc_auc_score(y_te, y_prob_te) if len(np.unique(y_te)) > 1 else np.nan
        m = binary_metrics_at_threshold(y_te, y_prob_te, best_th)
        m['AUC'] = auc_val

        m['train_AUC'] = roc_auc_score(y_tr, y_prob_tr) if len(np.unique(y_tr)) > 1 else np.nan
        tr_m = binary_metrics_at_threshold(y_tr, y_prob_tr, best_th)
        m['train_ACC'] = tr_m['ACC']
        m['train_SEN'] = tr_m['SEN']
        m['train_SPE'] = tr_m['SPE']
        m['train_F1']  = tr_m['F1']
        m['train_MCC'] = tr_m['MCC']

        m['inner_best_cv_score']   = grid.best_score_
        m['best_C']                = grid.best_params_.get('clf__C')
        m['best_sign_threshold']   = grid.best_params_.get('selector__threshold')
        m['n_features_selected']   = int(selector.selected_mask_.sum())
        m['best_threshold']        = best_th
        m['fold']                  = fold_idx
        fold_rows.append(m)

        clf_coef = best.named_steps['clf'].coef_.ravel()
        full_coef = np.zeros(len(feature_names))
        full_coef[selector.selected_mask_] = clf_coef
        coef_row = pd.Series(full_coef, index=feature_names, name=f'fold_{fold_idx}')
        coef_rows.append(coef_row)

    fold_df = pd.DataFrame(fold_rows)
    fold_df.to_csv(out_dir / "fold_metrics.csv", index=False)

    metrics_cols = ['ACC', 'SEN', 'SPE', 'F1', 'MCC', 'AUC']
    sub = fold_df[metrics_cols]
    summary = sub.agg(['mean', 'std', 'median']).T
    summary.columns = ['Mean', 'Std', 'Med']
    summary.to_csv(out_dir / "summary_metrics.csv")
    print(f"\n  Summary [{run_name}]:")
    print(summary[['Mean', 'Std']].to_string())

    param_freq = (
        fold_df[['best_C', 'best_sign_threshold']]
        .value_counts()
        .reset_index(name='Count')
        .sort_values('Count', ascending=False)
    )
    param_freq.to_csv(out_dir / "best_param_frequency.csv", index=False)

    th_mean = fold_df['best_threshold'].mean()
    th_std  = fold_df['best_threshold'].std()
    n_feat_mean = fold_df['n_features_selected'].mean()
    n_feat_std  = fold_df['n_features_selected'].std()
    print(f"    Optimal threshold: {th_mean:.3f} 簣 {th_std:.3f}")
    print(f"    Features selected: {n_feat_mean:.1f} 簣 {n_feat_std:.1f}")

    sel_freq = (sel_counts / total).sort_values(ascending=False)
    sel_freq.index.name = 'Feature'
    sel_freq.name = 'selection_frequency'
    sel_freq.to_csv(out_dir / "feature_selection_frequency.csv", header=True)

    coef_df  = pd.DataFrame(coef_rows)
    ever_selected = sel_counts[sel_counts > 0].index
    coef_df_sel = coef_df[ever_selected]

    coef_stat = pd.DataFrame({
        'mean':                coef_df_sel.mean(),
        'std':                 coef_df_sel.std(),
        'sign_consistency':    (np.sign(coef_df_sel) == np.sign(coef_df_sel.mean())).mean(),
        'abs_mean':            coef_df_sel.abs().mean(),
        'selection_frequency': sel_freq[ever_selected],
    }).sort_values('abs_mean', ascending=False)
    coef_stat.index.name = 'Feature'
    coef_stat.to_csv(out_dir / "coef_stability.csv")

    top_n   = min(40, len(coef_stat))
    plot_df = coef_stat.head(top_n)

    fig, ax = plt.subplots(figsize=(8, max(4, len(plot_df) * 0.38)))
    y_pos   = np.arange(len(plot_df))
    colors  = ['#C0392B' if v > 0 else '#2980B9' for v in plot_df['mean']]
    ax.barh(y_pos, plot_df['mean'], xerr=plot_df['std'],
            color=colors, alpha=0.75, error_kw=dict(ecolor='black', lw=1, capsize=3))
    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [str(feat)[:40] + f"  (sc={sc:.0%}, sel={sf:.0%})" for feat, sc, sf in
         zip(plot_df.index, plot_df['sign_consistency'], plot_df['selection_frequency'])],
        fontsize=7)
    ax.axvline(0, color='black', lw=0.8, ls='--')
    ax.set_xlabel('Coefficient (mean 簣 std across folds)')
    ax.set_title(
        f"{run_name}\nCoefficient Stability - Top {top_n} Features\n"
        f"(sc=sign_consistency, sel=selection_freq)",
        fontsize=10, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_dir / "coef_stability.png", dpi=600, bbox_inches='tight')
    plt.close()

    col = GROUP_COLORS[pos_class]

    fig, axes = plt.subplots(1, len(metrics_cols), figsize=(18, 4), sharey=False)
    for ax, met in zip(axes, metrics_cols):
        vals = fold_df[met].dropna()
        ax.violinplot(vals, positions=[0], showmedians=True, widths=0.7)
        ax.scatter([0]*len(vals), vals, s=8, alpha=0.25, color=col, zorder=3)
        ax.axhline(vals.mean(), color='black', lw=1.5, ls='--')
        ax.set_title(
            f"{met}\n{vals.mean():.3f}簣{vals.std():.3f}",
            fontsize=9, fontweight='bold')
        ax.set_xticks([])
    fig.suptitle(f"{run_name.replace('__','  |  ')} ??Distribution (n={total} folds)",
                 fontsize=12, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(out_dir / "violin_distribution.png", dpi=600, bbox_inches='tight')
    plt.close()

    melt = fold_df[metrics_cols].melt(var_name='Metric', value_name='Score')
    plt.figure(figsize=(10, 5))
    sns.boxplot(data=melt, x='Metric', y='Score',
                order=metrics_cols, color=col,
                linewidth=1.5, fliersize=2)
    sns.stripplot(data=melt, x='Metric', y='Score',
                  order=metrics_cols, color='black',
                  size=2, alpha=0.2, jitter=True)
    plt.axhline(0, color='black', lw=0.8, ls='--')
    plt.title(f"{run_name} (n={total} outer folds)", fontsize=12, fontweight='bold')
    plt.ylabel("Score"); plt.grid(axis='y', alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_dir / "boxplot.png", dpi=600, bbox_inches='tight')
    plt.close()

    return summary

# ??? Master Loop ???????????????????????????????????????????????????????????????
if __name__ == "__main__":
    comparisons  = [('C1', 'C2'), ('C1', 'C3')]
    feature_sets = ['combined']  # Because the goal is to combine modules + features. 
    # If the user wants micro/prot separately, we can leave all 3. Let's do all 3.
    feature_sets = ['micro', 'prot', 'combined']
    all_summaries = {}

    for pos, neg in comparisons:
        for fset in feature_sets:
            key = f"{pos}vs{neg}__{fset}"
            print(f"\n{'#'*65}")
            print(f"  Running: {key}")
            print(f"{'#'*65}")
            summary = run_nested_cv(pos, neg, fset)
            all_summaries[key] = summary

    rows = []
    for key, s in all_summaries.items():
        comp, fset = key.split('__')
        row = {'Comparison': comp, 'FeatureSet': fset}
        for met in ['ACC', 'SEN', 'SPE', 'F1', 'MCC', 'AUC']:
            row[f"{met}_Mean"] = s.loc[met, 'Mean']
            row[f"{met}_Std"]  = s.loc[met, 'Std']
        rows.append(row)

    master_df = pd.DataFrame(rows)
    master_df.to_csv(OUT_BASE / "master_comparison.csv", index=False)

    print("\n\n" + "="*70)
    print("MASTER COMPARISON TABLE")
    print("="*70)
    print(master_df.to_string(index=False))
    print(f"\n[DONE] All results saved to {OUT_BASE}")

