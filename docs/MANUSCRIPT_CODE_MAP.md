# Manuscript–Code Map

> Formal mapping of manuscript figures and supplementary tables to analysis scripts.
> Last updated: 2026-09-09

---

## Main Figures

### Figure 1 — Study Design / Analysis Workflow

Conceptual overview figure. No corresponding analysis code.

---

### Figure 2 — Cognitive Phenotyping + Microbiome Diversity

| Panel | Description | Script | Output |
|-------|-------------|--------|--------|
| **2a** | PCA of cognitive features | `scripts/10_figure_generation/plot_pca_cognitive.py` | `results/figures/` |
| **2b** | Total MoCA score boxplot | `scripts/10_figure_generation/plot_moca_boxplot.py` | `results/figures/phenotyping/moca_boxplot.png` |
| **2c** | MoCA subtests pairwise comparisons | `scripts/10_figure_generation/plot_moca_group_comparisons.py` | `results/figures/phenotyping/moca_subtests_*.png` |
| **2d** | Phylum alpha diversity | `scripts/10_figure_generation/plot_alpha_diversity.py` | `results/figures/` |
| **2e** | Genus alpha diversity | `scripts/10_figure_generation/plot_alpha_diversity.py` | `results/figures/` |
| **2f** | Full-cohort Bray-Curtis PCoA | `scripts/10_figure_generation/plot_bray_curtis_pcoa.py` | `results/figures/` |
| **2g** | PSM phylum alpha diversity | `scripts/02_psm_and_microbiome_diversity/plot_bray_curtis_pcoa_psm.py` | `results/figures/` |
| **2h** | PSM genus alpha diversity | `scripts/02_psm_and_microbiome_diversity/plot_bray_curtis_pcoa_psm.py` | `results/figures/` |
| **2i** | PSM Bray-Curtis PCoA | `scripts/02_psm_and_microbiome_diversity/plot_bray_curtis_pcoa_psm.py` | `results/figures/` |

---

### Figure 3 — Proteomic + Microbiome Association and Sensitivity

| Panel | Description | Script | Output |
|-------|-------------|--------|--------|
| **3a** | Proteomics forest plot | `scripts/03_proteomics_association/run_covariate_models.py` | `results/proteomics_covariate_analysis/` |
| **3b** | Proteomics model overlap | `scripts/10_figure_generation/plot_venn_ABCD_labeled.py` | `results/figures/` |
| **3c** | Proteomics P-value trajectories | `scripts/03_proteomics_association/run_covariate_models.py` | `results/proteomics_covariate_analysis/` |
| **3d** | Microbiome forest plot | `scripts/04_microbiome_differential_abundance/run_microbiome_covariate_models.R` | `results/microbiota_covariate_analysis/` |
| **3e** | Microbiome model overlap | `scripts/10_figure_generation/plot_venn_ABCD_labeled.py` | `results/figures/` |
| **3f** | C2 vs C1 microbial P-value trajectories | `scripts/04_microbiome_differential_abundance/run_microbiome_covariate_models.R` | `results/microbiota_covariate_analysis/` |
| **3g** | C3 vs C1 microbial P-value trajectories | `scripts/04_microbiome_differential_abundance/run_microbiome_covariate_models.R` | `results/microbiota_covariate_analysis/` |

---

### Figure 4 — Protein / Microbial Networks, Module Eigengenes, Cross-omics Correlations

| Panel | Description | Script | Output |
|-------|-------------|--------|--------|
| **4a** | Protein networks | `scripts/05_protein_network/run_proteomics_network_analysis.py` | `results/proteomics_network/` |
| **4b** | Microbial networks | `scripts/06_microbial_network/analyze_network_topology.py` | `results/microbiota_network/` |
| **4c** | Protein module eigengene scores | `scripts/10_figure_generation/plot_proteomics_modules.py` | `results/figures/` |
| **4d** | Microbial module eigengene scores | `scripts/10_figure_generation/plot_leiden_24_taxa.py` | `results/figures/` |
| **4e** | Cross-omics ME×ME correlations | UNRESOLVED (rendering script not located) | `results/fastspar/microbiome_x_proteomics_ME_corr.png` |

---

### Figure 5 — Pairwise Classification + CISS

| Panel | Description | Script | Output File |
|-------|-------------|--------|-------------|
| **5a** | C1 vs C2 classification performance | `scripts/08_classification_and_CISS/clinical_vs_omics_feature.py` | `results/modeling/classification/feature/logistic/clinica/master_comparison.csv` |
| **5b** | C1 vs C3 classification performance | `scripts/08_classification_and_CISS/clinical_vs_omics_feature.py` | `results/modeling/classification/feature/logistic/clinica/master_comparison.csv` |
| **5c** | Coefficient concordance across C1 vs C2 and C1 vs C3 | `scripts/08_classification_and_CISS/feature_ratio_modeling_cont_clin.py` | `used_features.csv` |
| **5d** | Domain-specific CISS across C1/C2/C3 | `scripts/08_classification_and_CISS/feature_ratio_modeling_cont_clin.py` & `plot_ciss_boxplot_moca_correlation.py` | `ratio_boxplot_groups.png` |
| **5e** | Microbiome CISS × Proteomic CISS landscape with MoCA Z-score | `scripts/08_classification_and_CISS/feature_ratio_modeling_cont_clin.py` | `sample_ratio_features.csv` |
| **5f** | CISS vs total MoCA Spearman correlations | `scripts/08_classification_and_CISS/plot_ciss_boxplot_moca_correlation.py` | `ratio_moca_correlation.png`, `ratio_moca_spearman.csv` |

---

### Figure 6 — GBB / Metabolite Correlations + Serial Mediation

| Panel | Description | Script | Output |
|-------|-------------|--------|--------|
| **6a** | CISS × targeted metabolite heatmap | `scripts/10_figure_generation/plot_metabolite_violins_heatmap_v2.py` | `results/figures/metabolic/` |
| **6b** | Mediation robustness Models 1–4 | `scripts/10_figure_generation/plot_mediation_sensitivity_forest.py` | `results/mediation/` |
| **6c** | Forward mediation | `scripts/09_serial_mediation/sequential_mediation.py` | `results/mediation/` |
| **6d** | Reverse mediation | `scripts/09_serial_mediation/sequential_mediation.py` | `results/mediation/` |

---

## Supplementary Tables

### Supplementary Table 1 — Full-Cohort Baseline Characteristics
- **Script:** `scripts/01_cognitive_phenotyping/generate_table1.py`
- **Output:** `results/tables/Table1_Baseline_Characteristics.csv`

### Supplementary Table 2 — PSM Baseline Characteristics
- **Script:** `scripts/02_psm_and_microbiome_diversity/generate_psm_table1.py`
- **Output:** `results/tables/Table1_PSM_Baseline_Characteristics.csv`

### Supplementary Table 3 — Proteomic Network Topology Metrics
- **Script:** `scripts/03_proteomics_association/analyze_proteomics_topology.py`
- **Output:** `results/proteomics_network/network_metrics_*.csv`

### Supplementary Table 4 — Microbial Network Topology Metrics
- **Script:** `scripts/06_microbial_network/analyze_network_topology.py`
- **Output:** `results/microbiota_network/network_metrics_*.csv`

### Supplementary Table 5 — Classification Performance Metrics
- **Script:** [`scripts/08_classification_and_CISS/clinical_vs_omics_feature.py`](../scripts/08_classification_and_CISS/clinical_vs_omics_feature.py)
- **Output:** `results/modeling/classification/feature/logistic/clinica/master_comparison.csv`

**100% Exact Match to Manuscript (`SUPPLEMENTARY.docx` Table 5):**

| Comparison | Feature Configuration | Script FeatureSet | ACC | SEN | SPE | F1 | MCC | AUC |
|------------|-----------------------|-------------------|-----|-----|-----|----|-----|-----|
| **C1 vs C2** | Clinical-only | `clinical` | 0.533 ± 0.088 | 0.668 ± 0.204 | 0.483 ± 0.148 | 0.434 ± 0.097 | 0.142 ± 0.154 | 0.617 ± 0.104 |
| | Multi-omics Only | `combined` | 0.698 ± 0.079 | 0.448 ± 0.172 | 0.793 ± 0.123 | 0.441 ± 0.126 | 0.252 ± 0.162 | 0.652 ± 0.091 |
| | Integrated | `clin_combined` | 0.683 ± 0.080 | 0.484 ± 0.167 | 0.759 ± 0.121 | 0.451 ± 0.113 | 0.246 ± 0.158 | 0.668 ± 0.085 |
| **C1 vs C3** | Clinical-only | `clinical` | 0.723 ± 0.071 | 0.530 ± 0.160 | 0.772 ± 0.091 | 0.435 ± 0.107 | 0.273 ± 0.143 | 0.707 ± 0.079 |
| | Multi-omics Only | `combined` | 0.712 ± 0.069 | 0.462 ± 0.178 | 0.776 ± 0.106 | 0.386 ± 0.107 | 0.222 ± 0.132 | 0.702 ± 0.074 |
| | Integrated | `clin_combined` | 0.802 ± 0.054 | 0.483 ± 0.149 | 0.883 ± 0.067 | 0.494 ± 0.124 | 0.384 ± 0.153 | 0.779 ± 0.071 |
