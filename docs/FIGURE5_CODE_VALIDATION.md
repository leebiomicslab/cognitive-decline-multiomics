# Figure 5 Code Validation Report

**Primary candidate scripts:**
- [`scripts/08_classification_and_CISS/clinical_vs_omics_feature.py`](../scripts/08_classification_and_CISS/clinical_vs_omics_feature.py) (Pairwise classification & Supplementary Table 5)
- [`scripts/08_classification_and_CISS/feature_ratio_modeling_cont_clin.py`](../scripts/08_classification_and_CISS/feature_ratio_modeling_cont_clin.py) (CISS score construction & Figure 5d/5e/5f data)
- [`scripts/08_classification_and_CISS/plot_ciss_boxplot_moca_correlation.py`](../scripts/08_classification_and_CISS/plot_ciss_boxplot_moca_correlation.py) (Companion plotting for Figure 5d/5f)

**Validated:** 2026-09-09  
**Validator note:** Read-only comparison against original project script outputs and manuscript text.

---

## A. Pairwise Classification Comparisons

| Manuscript Requirement | Code Location | Status |
|---|---|---|
| C1 vs C2 | `comparisons = [('C1', 'C2'), ('C1', 'C3')]` | ✅ CONFIRMED |
| C1 vs C3 | Same comparisons list | ✅ CONFIRMED |
| Positive class = C1 | `y_sub = y_all.loc[...].map({pos_class: 1, neg_class: 0})` | ✅ CONFIRMED |

---

## B. Feature Configurations & Mapping

| Manuscript Label | Code FeatureSet | Description | Output Location | Status |
|---|---|---|---|---|
| **Clinical-only** | `clinical` | Clinical features only | `master_comparison.csv` row 0 & 3 | ✅ CONFIRMED |
| **Multi-omics Only** | `combined` | Microbiome + Proteomics | `master_comparison.csv` row 1 & 4 | ✅ CONFIRMED |
| **Integrated** | `clin_combined` | Clinical + Microbiome + Proteomics | `master_comparison.csv` row 2 & 5 | ✅ CONFIRMED |

---

## C. Model Specification

| Manuscript Requirement | Code Implementation | Status |
|---|---|---|
| Logistic regression | `LogisticRegression(...)` | ✅ CONFIRMED |
| L2 regularization | `penalty='l2'` | ✅ CONFIRMED |
| Balanced class weights | `class_weight='balanced'` | ✅ CONFIRMED |
| Solver | `solver='lbfgs'` | ✅ CONFIRMED |
| Outer CV | 5-fold × 20 repeats (`OUTER_SPLITS=5`, `OUTER_REPEATS=20`) | ✅ CONFIRMED |
| Inner CV | 4-fold (`INNER_SPLITS=4`) | ✅ CONFIRMED |
| Hyperparameter tuning | `C ∈ [0.0005..10]`, scored by AUC | ✅ CONFIRMED |
| Feature Selector | `SignConsistencySelector(threshold=0.75)` | ✅ CONFIRMED |
| Random state | `RANDOM_STATE = 42` | ✅ CONFIRMED |

---

## D. Final Outputs and Figure/Table Mapping

| Output File | Figure / Table | Provenance Script |
|---|---|---|
| `master_comparison.csv` | **Figure 5a** (C1 vs C2 classification performance) | `clinical_vs_omics_feature.py` |
| `master_comparison.csv` | **Figure 5b** (C1 vs C3 classification performance) | `clinical_vs_omics_feature.py` |
| `used_features.csv` | **Figure 5c** (Coefficient concordance across C1 vs C2 and C1 vs C3) | `feature_ratio_modeling_cont_clin.py` (rendering script not located) |
| `ratio_boxplot_groups.png` | **Figure 5d** (Domain-specific CISS across C1/C2/C3) | `feature_ratio_modeling_cont_clin.py` & `plot_ciss_boxplot_moca_correlation.py` |
| `sample_ratio_features.csv` | **Figure 5e** (Microbiome CISS × Proteomic CISS landscape with MoCA Z-score) | `feature_ratio_modeling_cont_clin.py` (rendering script not located) |
| `ratio_moca_correlation.png`, `ratio_moca_spearman.csv` | **Figure 5f** (CISS vs total MoCA Spearman correlations) | `plot_ciss_boxplot_moca_correlation.py` |
| `master_comparison.csv` | **Supplementary Table 5** (Full performance metrics table) | `clinical_vs_omics_feature.py` |

---

## Overall Assessment

**`clinical_vs_omics_feature.py` and `feature_ratio_modeling_cont_clin.py` are CONFIRMED as the final production scripts producing Figure 5 and Supplementary Table 5 results.**
