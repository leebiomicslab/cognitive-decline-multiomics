"""
Status:           LEGACY / EXPLORATORY ANALYSIS
                  Not used for final Figure 5 or final CISS results.
                  This script explored CISS construction using DA effect sizes
                  (e.g., ANCOM-BC2 log-fold change) as the basis for feature weighting,
                  rather than logistic regression coefficient concordance.
                  The final manuscript uses logistic regression coefficient concordance,
                  implemented in feature_ratio_modeling_cont_clin.py.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Compute CISS ratio features from classification probabilities
Manuscript:       Methods - CISS construction
Figure/Table:     Fig. 6a
"""

from pathlib import Path
"""
DA-Based Ratio Modeling
========================
敺?DA ??蝯?? C2/C3 ?孵?銝?港??賡＊???孵噩嚗?  - ?: ANCOMBC2 p_groupC2 < 0.05 OR p_groupC3 < 0.05嚗fc ???  - ?: ANCOVA P_Value < 0.05 銝??喳?銝??頛? 95% CI 銝楊?塚?Beta ???
???亦?銝??log-ratio嚗ean(Pos ?孵噩) - Mean(Neg ?孵噩)
??拙?ratio (Micro_Ratio, Prot_Ratio) 撱箸芋

蝯?: 7 蝔桅?蝵?? C1vsC2/C1vsC3 ??Performance Bar + AUC ??
      boxplot by group + MoCA correlation
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from scipy import stats
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.metrics import (accuracy_score, f1_score, matthews_corrcoef,
                             roc_auc_score, confusion_matrix, roc_curve)
import warnings
warnings.filterwarnings('ignore')

# ??? Paths ????????????????????????????????????????????????????????????????????
OUT_BASE   = MODELING_DIR / "classification" / "feature" / "logistic" / "ratio_modeling_da"
OUT_BASE.mkdir(parents=True, exist_ok=True)

MICRO_DA   = MICRO_ASSOC_DIR / "Model_A_Main" / "microbiota_da_results_final.csv"
PROT_DA    = PROT_ASSOC_DIR / "Main_Model" / "proteomics_ancova_detailed_results.csv"
CLINICAL   = CLINICAL_CSV
MICRO_ASV  = MICRO_GENUS_COUNTS
PROT_EXPR  = PROTEOMICS_LOG2
VASCULAR   = BASE_DIR / "data" / "vascular.csv"

OUTER_SPLITS  = 5
OUTER_REPEATS = 20
INNER_SPLITS  = 4
RANDOM_STATE  = 42
CLR_EPS       = 1e-6

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial'],
                     'axes.labelweight': 'bold', 'axes.titleweight': 'bold',
                     'axes.spines.top': False, 'axes.spines.right': False})

GROUP_ORDER  = ['C1', 'C2', 'C3']
GROUP_COLORS_GRP = {'C1': '#ee7a5b', 'C2': '#eca362', 'C3': '#6fcbbc'}
CFG_COLORS   = {'micro': '#1F77B4', 'prot': '#2CA02C',
                'micro_prot': '#9467BD', 'combined': '#D62728'}
CONFIGS      = ['micro', 'prot', 'micro_prot']   # combined == micro_prot here

# ??? 1. Feature Selection from DA Results ?????????????????????????????????????
print("=" * 60)
print("Selecting features from DA results...")

## Microbiome: either C2 or C3 significant (p<0.05), same LFC direction
micro_da = pd.read_csv(MICRO_DA)
micro_da['sig_C2'] = micro_da['p_groupC2'] < 0.05
micro_da['sig_C3'] = micro_da['p_groupC3'] < 0.05
micro_da['either_sig'] = micro_da['sig_C2'] | micro_da['sig_C3']
micro_da['same_dir']   = (micro_da['lfc_groupC2'].apply(lambda x: x > 0) ==
                          micro_da['lfc_groupC3'].apply(lambda x: x > 0))
micro_sig = micro_da[micro_da['either_sig'] & micro_da['same_dir']].copy()
micro_sig['avg_lfc'] = (micro_sig['lfc_groupC2'] + micro_sig['lfc_groupC3']) / 2

pos_micro = micro_sig[micro_sig['avg_lfc'] > 0]['taxon'].tolist()
neg_micro = micro_sig[micro_sig['avg_lfc'] < 0]['taxon'].tolist()
print(f"\nMicro  Pos ({len(pos_micro)}): {pos_micro}")
print(f"Micro  Neg ({len(neg_micro)}): {neg_micro}")

## Proteomics: P_Value < 0.05嚗??喳?銝??頛? 95% CI 銝楊?塚?Beta ?孵?銝??prot_da = pd.read_csv(PROT_DA)
prot_da['sig_C2_ci'] = prot_da['CI_C2_Low'] * prot_da['CI_C2_High'] > 0   # CI doesn't cross 0
prot_da['sig_C3_ci'] = prot_da['CI_C3_Low'] * prot_da['CI_C3_High'] > 0
prot_da['either_ci'] = prot_da['sig_C2_ci'] | prot_da['sig_C3_ci']         # at least one CI significant
prot_da['p_sig']     = prot_da['P_Value'] < 0.05                            # overall ANCOVA p < 0.05
prot_da['same_dir']  = (prot_da['Beta_C2'].apply(lambda x: x > 0) ==
                        prot_da['Beta_C3'].apply(lambda x: x > 0))
prot_sig = prot_da[prot_da['p_sig'] & prot_da['either_ci'] & prot_da['same_dir']].copy()
prot_sig['avg_beta'] = (prot_sig['Beta_C2'] + prot_sig['Beta_C3']) / 2

# avg_beta > 0 ??higher in C2/C3 vs C1 ??Negative marker for C1
# avg_beta < 0 ??lower  in C2/C3 vs C1 ??Positive marker for C1
pos_prot = prot_sig[prot_sig['avg_beta'] < 0]['Protein'].tolist()
neg_prot = prot_sig[prot_sig['avg_beta'] > 0]['Protein'].tolist()

print(f"\nProt   Pos ({len(pos_prot)}): {pos_prot}")
print(f"Prot   Neg ({len(neg_prot)}): {neg_prot}")

# ??? Export feature list ???????????????????????????????????????????????????????
feat_rows = []
for f in pos_micro: feat_rows.append({'Category': 'Micro', 'Sign': 'Positive', 'Feature': f})
for f in neg_micro: feat_rows.append({'Category': 'Micro', 'Sign': 'Negative', 'Feature': f})
for f in pos_prot:  feat_rows.append({'Category': 'Prot',  'Sign': 'Positive', 'Feature': f})
for f in neg_prot:  feat_rows.append({'Category': 'Prot',  'Sign': 'Negative', 'Feature': f})
pd.DataFrame(feat_rows).to_csv(OUT_BASE / "da_selected_features.csv", index=False)
print(f"\nSaved: da_selected_features.csv")

# ??? 2. Load Raw Expression Data ??????????????????????????????????????????????
def get_clr(df, eps=CLR_EPS):
    df_imp = df.replace(0, eps).astype(float)
    log_df = np.log(df_imp)
    return log_df.sub(log_df.mean(axis=1), axis=0)

print("\nLoading raw data...")
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
print(f"Common samples: {len(common_all)}")

# Build master feature table (all needed columns)
all_features = set(pos_micro + neg_micro + pos_prot + neg_prot)
X_micro_df = asv_clr[[c for c in asv_clr.columns if c in all_features]].loc[common_all]
X_prot_df  = prot_raw[[c for c in prot_raw.columns if c in all_features]].loc[common_all]
X_master   = pd.concat([X_micro_df, X_prot_df], axis=1)
X_master   = X_master.fillna(X_master.median())

# Fill any still-missing feature columns with 0
for f in all_features:
    if f not in X_master.columns:
        print(f"  [WARN] Feature not found in data: {f}")
        X_master[f] = 0.0

# ??? 3. Transformers ??????????????????????????????????????????????????????????
class DfStandardScaler(BaseEstimator, TransformerMixin):
    def fit(self, X, y=None):
        self.means_ = X.mean(axis=0)
        self.stds_  = X.std(axis=0).replace(0, 1)
        return self
    def transform(self, X):
        return (X - self.means_) / self.stds_

class RatioTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, config='micro_prot'):
        self.config = config
    def fit(self, X, y=None):
        return self
    def transform(self, X):
        out = pd.DataFrame(index=X.index)
        pm = [f for f in pos_micro if f in X.columns]
        nm = [f for f in neg_micro if f in X.columns]
        pp = [f for f in pos_prot  if f in X.columns]
        np_ = [f for f in neg_prot  if f in X.columns]
        if self.config in ['micro', 'micro_prot']:
            out['Micro_Ratio'] = X[pm].mean(axis=1).fillna(0) - X[nm].mean(axis=1).fillna(0)
        if self.config in ['prot', 'micro_prot']:
            out['Prot_Ratio']  = X[pp].mean(axis=1).fillna(0) - X[np_].mean(axis=1).fillna(0)
        return out.astype(float).values

# Compute global sample ratios (all samples, for boxplot/correlation)
scaler_g = DfStandardScaler().fit(X_master)
X_scaled_g = scaler_g.transform(X_master)
rt = RatioTransformer('micro_prot')
ratio_vals = rt.transform(X_scaled_g)
df_ratio = pd.DataFrame(ratio_vals, index=X_master.index, columns=['Micro_Ratio', 'Prot_Ratio'])
df_ratio['Group'] = y_all.loc[df_ratio.index, 'Group'].values
df_ratio.to_csv(OUT_BASE / "sample_ratio_features.csv")
print("Saved: sample_ratio_features.csv")

# ??? 4. Metrics & CV ??????????????????????????????????????????????????????????
def binary_metrics_at_threshold(y_true, y_prob, threshold=0.5):
    y_pred = (y_prob >= threshold).astype(int)
    cm = confusion_matrix(y_true, y_pred)
    if cm.shape != (2, 2):
        return dict(ACC=np.nan, SEN=np.nan, SPE=np.nan, F1=np.nan, MCC=np.nan)
    tn, fp, fn, tp = cm.ravel()
    return {'ACC': accuracy_score(y_true, y_pred),
            'SEN': tp / (tp + fn) if (tp + fn) > 0 else 0.0,
            'SPE': tn / (tn + fp) if (tn + fp) > 0 else 0.0,
            'F1' : f1_score(y_true, y_pred, zero_division=0),
            'MCC': matthews_corrcoef(y_true, y_pred)}

def find_best_threshold(y_true, y_prob, n=200):
    best_th, best_f1 = 0.5, 0.0
    for th in np.linspace(0, 1, n):
        y_pred = (y_prob >= th).astype(int)
        f = f1_score(y_true, y_pred, zero_division=0)
        if f > best_f1:
            best_f1, best_th = f, th
    return best_th

PARAM_GRID = {'clf__C': [0.001, 0.01, 0.1, 1, 10, 100]}

def run_pipeline(pos_class, neg_class, config):
    print(f"  [{pos_class} vs {neg_class} | {config}]", flush=True)
    mask    = y_all.loc[common_all, 'Group'].isin([pos_class, neg_class])
    samples = np.array(common_all[mask])
    y_sub   = y_all.loc[samples, 'Group'].map({pos_class: 1, neg_class: 0}).values
    X_sub   = X_master.loc[samples]

    outer_cv = RepeatedStratifiedKFold(n_splits=OUTER_SPLITS, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE)
    inner_cv = StratifiedKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    fold_rows       = []
    y_prob_agg      = np.zeros(len(samples))
    y_prob_cnt      = np.zeros(len(samples))

    for fold_idx, (tr_idx, te_idx) in enumerate(outer_cv.split(samples, y_sub)):
        X_tr, X_te = X_sub.iloc[tr_idx], X_sub.iloc[te_idx]
        y_tr, y_te = y_sub[tr_idx], y_sub[te_idx]

        pipe = Pipeline([
            ('scaler', DfStandardScaler()),
            ('ratio',  RatioTransformer(config=config)),
            ('clf',    LogisticRegression(penalty='l2', solver='lbfgs',
                                          class_weight='balanced', max_iter=2000,
                                          random_state=RANDOM_STATE))
        ])
        grid = GridSearchCV(pipe, PARAM_GRID, cv=inner_cv, scoring='roc_auc', n_jobs=-1)
        grid.fit(X_tr, y_tr)
        best = grid.best_estimator_

        y_prob_tr = best.predict_proba(X_tr)[:, 1]
        best_th   = find_best_threshold(y_tr, y_prob_tr)
        y_prob_te = best.predict_proba(X_te)[:, 1]

        y_prob_agg[te_idx] += y_prob_te
        y_prob_cnt[te_idx] += 1

        m = binary_metrics_at_threshold(y_te, y_prob_te, best_th)
        m['AUC']  = roc_auc_score(y_te, y_prob_te) if len(np.unique(y_te)) > 1 else np.nan
        m['fold'] = fold_idx
        fold_rows.append(m)

    y_prob_avg  = y_prob_agg / y_prob_cnt
    fpr, tpr, _ = roc_curve(y_sub, y_prob_avg)
    final_auc   = roc_auc_score(y_sub, y_prob_avg) if len(np.unique(y_sub)) > 1 else np.nan

    fold_df  = pd.DataFrame(fold_rows)
    metrics_cols = ['ACC', 'SEN', 'SPE', 'F1', 'MCC', 'AUC']
    summary  = fold_df[metrics_cols].agg(['mean', 'std']).T
    summary.columns = ['Mean', 'Std']

    out_dir = OUT_BASE / f"{pos_class}_vs_{neg_class}" / config
    out_dir.mkdir(parents=True, exist_ok=True)
    fold_df.to_csv(out_dir / "fold_metrics.csv", index=False)
    summary.to_csv(out_dir / "summary_metrics.csv")

    return summary, fpr, tpr, final_auc

# ??? 5. Master Loop ???????????????????????????????????????????????????????????
print("\nRunning pipelines...")
comparisons  = [('C1', 'C2'), ('C1', 'C3')]
all_summaries  = {}
all_rocs       = {}

for pos, neg in comparisons:
    cmp_key = f"{pos}_vs_{neg}"
    all_rocs[cmp_key] = {}
    for cfg in CONFIGS:
        key = f"{pos}vs{neg}__{cfg}"
        s, fpr, tpr, auc_val = run_pipeline(pos, neg, cfg)
        all_summaries[key] = s
        all_rocs[cmp_key][cfg] = (fpr, tpr, auc_val)

# ??? 6. Performance Bar Plots ?????????????????????????????????????????????????
print("\nCreating plots...")
metrics_order = ['AUC', 'ACC', 'SEN', 'SPE', 'F1', 'MCC']

for pos, neg in comparisons:
    plot_data = []
    for cfg in CONFIGS:
        k  = f"{pos}vs{neg}__{cfg}"
        df = all_summaries[k]
        for m in metrics_order:
            plot_data.append({'Config': cfg, 'Metric': m,
                              'Mean': df.loc[m, 'Mean'], 'Std': df.loc[m, 'Std']})
    pdf = pd.DataFrame(plot_data)

    fig, ax = plt.subplots(figsize=(12, 5))
    sns.barplot(data=pdf, x='Metric', y='Mean', hue='Config',
                palette={k: CFG_COLORS[k] for k in CONFIGS},
                hue_order=CONFIGS,
                capsize=0.0, errorbar=None, edgecolor='k', ax=ax)

    # Manual error bars
    bars = [p for p in ax.patches if p.get_width() > 0]
    n_hues = len(CONFIGS)
    n_mets = len(metrics_order)
    for i, p in enumerate(bars):
        cfg_i = i // n_mets
        met_i = i % n_mets
        if cfg_i >= n_hues: continue
        cfg_name = CONFIGS[cfg_i]
        met_name = metrics_order[met_i]
        row = pdf[(pdf['Config'] == cfg_name) & (pdf['Metric'] == met_name)]
        if row.empty: continue
        std_val = row['Std'].values[0]
        x = p.get_x() + p.get_width() / 2
        ax.errorbar(x, p.get_height(), yerr=std_val,
                    color='black', capsize=3, elinewidth=1, linestyle='none')

    ax.set_title(f"DA-Based Ratio Model Performance ({pos} vs {neg})", fontsize=14, pad=12)
    ax.set_ylabel("Score"); ax.set_xlabel("")
    ax.set_ylim(0, 1.08)
    ax.legend(title='Config', bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(OUT_BASE / f"performance_bar_{pos}_vs_{neg}.png", dpi=600, bbox_inches='tight')
    plt.close()

# ??? 7. ROC AUC Overlay ???????????????????????????????????????????????????????
for pos, neg in comparisons:
    fig, ax = plt.subplots(figsize=(7, 7))
    for cfg in CONFIGS:
        fpr, tpr, auc_val = all_rocs[f"{pos}_vs_{neg}"][cfg]
        ax.plot(fpr, tpr, lw=2, color=CFG_COLORS[cfg],
                label=f"{cfg.upper()}  (AUC = {auc_val:.3f})")
    ax.plot([0, 1], [0, 1], 'navy', lw=1.5, ls='--')
    ax.set_xlim([-0.02, 1.0]); ax.set_ylim([0.0, 1.02])
    ax.set_xlabel('False Positive Rate'); ax.set_ylabel('True Positive Rate')
    ax.set_title(f'DA-Based Ratio ROC Curves\n({pos} vs {neg})', fontsize=14)
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_BASE / f"roc_auc_{pos}_vs_{neg}.png", dpi=600, bbox_inches='tight')
    plt.close()

# ??? 8. Group Boxplot ?????????????????????????????????????????????????????????
ratio_cols   = ['Micro_Ratio', 'Prot_Ratio']
ratio_labels = {'Micro_Ratio': 'Microbiome Ratio (DA)', 'Prot_Ratio': 'Proteomics Ratio (DA)'}

fig, axes = plt.subplots(1, 2, figsize=(10, 6))
for ax, col in zip(axes, ratio_cols):
    data_by_grp = [df_ratio.loc[df_ratio['Group'] == g, col].dropna().values for g in GROUP_ORDER]
    bp = ax.boxplot(data_by_grp, patch_artist=True, notch=False,
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(linewidth=1.5),
                    capprops=dict(linewidth=1.5),
                    flierprops=dict(marker='o', markersize=3, alpha=0.4, linestyle='none'))
    for patch, grp in zip(bp['boxes'], GROUP_ORDER):
        patch.set_facecolor(GROUP_COLORS_GRP[grp]); patch.set_alpha(0.8)

    for i, (grp, vals) in enumerate(zip(GROUP_ORDER, data_by_grp)):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                   color=GROUP_COLORS_GRP[grp], edgecolor='k', linewidth=0.4,
                   s=18, alpha=0.5, zorder=3)

    pairs = [('C1', 'C2', 1, 2), ('C1', 'C3', 1, 3), ('C2', 'C3', 2, 3)]
    y_all_vals = np.concatenate(data_by_grp)
    y_top  = y_all_vals.max()
    y_bot  = y_all_vals.min()
    tick   = (y_top - y_bot) * 0.07
    y_sig  = y_top + tick * 0.5

    for g1, g2, x1, x2 in pairs:
        v1 = df_ratio.loc[df_ratio['Group'] == g1, col].dropna().values
        v2 = df_ratio.loc[df_ratio['Group'] == g2, col].dropna().values
        _, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
        ax.plot([x1, x1, x2, x2], [y_sig, y_sig + tick*0.3, y_sig + tick*0.3, y_sig], color='black', lw=1.2)
        ax.text((x1 + x2) / 2, y_sig + tick*0.35, sig, ha='center', va='bottom', fontsize=10)
        y_sig += tick * 1.6

    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(GROUP_ORDER, fontsize=11)
    ax.set_title(ratio_labels[col], fontsize=12, pad=10)
    ax.set_ylabel("Log-Ratio Score"); ax.axhline(0, color='grey', lw=0.8, ls='--', alpha=0.6)

legend_patches = [mpatches.Patch(facecolor=GROUP_COLORS_GRP[g], edgecolor='k', label=g) for g in GROUP_ORDER]
fig.legend(handles=legend_patches, title='Group', loc='upper right', fontsize=10)
fig.suptitle("DA-Derived Ratio Distribution by Cognitive Group", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(OUT_BASE / "ratio_boxplot_groups.png", dpi=600, bbox_inches='tight')
plt.close()
print("Saved: ratio_boxplot_groups.png")

# ??? 9. MoCA Correlation ??????????????????????????????????????????????????????
df_vasc = pd.read_csv(VASCULAR).rename(columns={'id': 'SampleID'})
df_vasc['SampleID'] = df_vasc['SampleID'].astype(str)
df_vasc = df_vasc[['SampleID', 'moca_1']].dropna()

df_ratio_reset = df_ratio.reset_index().rename(columns={'index': 'SampleID'})
df_corr = pd.merge(df_ratio_reset, df_vasc, on='SampleID').dropna(subset=['moca_1'])
print(f"Samples with MoCA: {len(df_corr)}")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
corr_rows = []
for ax, col in zip(axes, ratio_cols):
    x = df_corr[col].values
    y = df_corr['moca_1'].values
    rho, pval = stats.spearmanr(x, y)
    pstr = 'p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'
    corr_rows.append({'Ratio': col, 'Spearman_rho': rho, 'p_value': pval})

    for grp in GROUP_ORDER:
        mask = df_corr['Group'] == grp
        ax.scatter(x[mask], y[mask], c=GROUP_COLORS_GRP[grp], edgecolor='k',
                   linewidth=0.3, s=30, alpha=0.7, label=grp, zorder=3)

    m, b = np.polyfit(x, y, 1)
    x_line = np.linspace(x.min(), x.max(), 200)
    ax.plot(x_line, m * x_line + b, color='black', lw=1.8, ls='--', alpha=0.8)

    ax.set_xlabel(ratio_labels[col], fontsize=10)
    ax.set_ylabel("MoCA Score", fontsize=10)
    ax.set_title(f"{ratio_labels[col]}\n? = {rho:.3f},  {pstr}", fontsize=11, pad=8)
    ax.grid(alpha=0.25)

handles = [mpatches.Patch(facecolor=GROUP_COLORS_GRP[g], edgecolor='k', label=g) for g in GROUP_ORDER]
fig.legend(handles=handles, title='Group', loc='upper right', fontsize=10)
fig.suptitle("DA-Derived Ratios vs. MoCA Score (Spearman)", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(OUT_BASE / "ratio_moca_correlation.png", dpi=600, bbox_inches='tight')
plt.close()
print("Saved: ratio_moca_correlation.png")

corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(OUT_BASE / "ratio_moca_spearman.csv", index=False)
print("\nSpearman Correlation:")
print(corr_df.to_string(index=False))
print(f"\n[DONE] All outputs saved to {OUT_BASE}")



