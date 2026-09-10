# Code Inventory

Generated: 2026-09-09 | Repository: `nature_aging_code_release`
Source scan: Actual Public Release Repository Tree

---

## Legend

| Column | Meaning |
|--------|---------|
| **Script Path** | Relative repository path of the script |
| **Lang** | Programming language (Python / R / Batch) |
| **Purpose** | Key scientific / computational function |
| **Manuscript Section** | Relevant section in the manuscript |
| **Figure / Table** | Output figure or table in manuscript |
| **Input Files** | Required input data files (noting restricted upstream inputs) |
| **Output Files** | Expected / generated output files |
| **Restricted Input?** | Requires participant-level data or precomputed restricted inputs |
| **Status** | Implementation status: **Final**, **Legacy**, **Helper** |

---

## 1. Upstream Microbiome Preprocessing (`scripts/00_microbiome_preprocessing/`)

| Script Path | Lang | Purpose | Manuscript Section | Figure / Table | Input Files | Output Files | Restricted Input? | Status |
|-------------|------|---------|--------------------|----------------|-------------|--------------|-------------------|--------|
| `scripts/00_microbiome_preprocessing/README.md` | Doc | Upstream 16S rRNA preprocessing specification (nf-core/ampliseq v2.11.0, DADA2, SILVA 138.1) | Methods — Microbiome preprocessing | — | Raw 16S FASTQ files (restricted) | `asv_genus_table_reads.csv` (expected output) | **YES** | Upstream Spec |

---

## 2. Cognitive Phenotyping (`scripts/01_cognitive_phenotyping/`)

| Script Path | Lang | Purpose | Manuscript Section | Figure / Table | Input Files | Output Files | Restricted Input? | Status |
|-------------|------|---------|--------------------|----------------|-------------|--------------|-------------------|--------|
| `scripts/01_cognitive_phenotyping/generate_table1.py` | Python | Baseline demographic and clinical characteristics Table 1 (Kruskal-Wallis, Chi-square) | Methods — Clinical characterization | Supp. Table 1 | `vascular.csv`, `medication.xlsx`, `k-means_3group.csv` (restricted upstream input, not included) | `Table1_Baseline_Characteristics.csv` (expected output) | **YES** | **Final** |
| `scripts/01_cognitive_phenotyping/clinical_collinearity_check.py` | Python | Continuous and categorical covariate correlation matrix & VIF audit | Methods — Covariate collinearity | Supp. (VIF) | `clinical_merged_371.csv` | `clinical_vif_values.csv`, heatmap PNGs (expected outputs) | **YES** | **Helper** |
| `scripts/01_cognitive_phenotyping/clinical_vif_4models.py` | Python | Multicollinearity VIF audit across all 4 covariate adjustment models (A–D) | Methods — Covariate collinearity | Supp. Fig. 4 | `clinical_merged_371.csv` | VIF plot PNGs per model (expected outputs) | **YES** | **Helper** |
| `scripts/01_cognitive_phenotyping/clinical_vif_proteomics_model_b.py` | Python | Multicollinearity VIF audit specifically for proteomics Model B covariates | Methods — Covariate collinearity | Supp. | `clinical_merged_371.csv` | Model B VIF plot PNG (expected output) | **YES** | **Helper** |

---

## 3. PSM and Microbiome Diversity (`scripts/02_psm_and_microbiome_diversity/`)

| Script Path | Lang | Purpose | Manuscript Section | Figure / Table | Input Files | Output Files | Restricted Input? | Status |
|-------------|------|---------|--------------------|----------------|-------------|--------------|-------------------|--------|
| `scripts/02_psm_and_microbiome_diversity/prepare_microbiome_metadata.py` | Python | Merges clinical data, K-means cluster assignments (`k-means_3group.csv`), and medication data | Methods — Data preprocessing | — | `vascular.csv`, `medication.xlsx`, `k-means_3group.csv` (restricted upstream input, not included) | `clinical_merged_371.csv` (expected output) | **YES** | **Helper / Preprocessing** |
| `scripts/02_psm_and_microbiome_diversity/run_3group_psm.py` | Python | 1:1 Propensity Score Matching across C1/C2/C3 (n=49 per phenotype) | Methods — Propensity score matching | Supp. Table 2, Supp. Fig. 2 | `clinical_merged_371.csv` | `psm_matched_cohort.csv` (expected output) | **YES** | **Final** |
| `scripts/02_psm_and_microbiome_diversity/generate_psm_table1.py` | Python | Baseline clinical characteristics Table 1 for PSM-matched cohort | Methods — Propensity score matching | Supp. Table 2 | `psm_matched_cohort.csv` | `PSM_Table1_Baseline_Characteristics.csv` (expected output) | **YES** | **Final** |
| `scripts/02_psm_and_microbiome_diversity/plot_bray_curtis_pcoa_psm.py` | Python | Bray-Curtis PCoA beta diversity visualization for PSM-matched cohort | Methods — Microbiome diversity | Supp. Fig. 2 | `asv_genus_table.csv`, `psm_matched_cohort.csv` | `pcoa_bray_curtis_psm.png` (expected output) | **YES** | **Final** |

---

## 4. Proteomics Association (`scripts/03_proteomics_association/`)

| Script Path | Lang | Purpose | Manuscript Section | Figure / Table | Input Files | Output Files | Restricted Input? | Status |
|-------------|------|---------|--------------------|----------------|-------------|--------------|-------------------|--------|
| `scripts/03_proteomics_association/run_covariate_models.py` | Python | **Main proteomics ANCOVA** — OLS with HC3 SE, omnibus F-test, BH-FDR, adjusted marginal means across Models A–D | Methods — Proteomics association | Fig. 3a–c, Supp. Fig. 5, Supp. Table 3 | `data_log2_mapping_deduplicated.csv`, `clinical_merged_371.csv` | `proteomics_ancova_detailed_results.csv` (expected output), forest plots, rank plots | **YES** | **Final** |
| `scripts/03_proteomics_association/analyze_proteomics_topology.py` | Python | Network topology metrics (degree, betweenness centrality, modularity) for protein networks | Methods — Proteomics association | Supp. Table 3 | Protein network GraphML files | `protein_network_metrics.csv` (expected output) | **YES** | **Final** |

---

## 5. Microbiome Differential Abundance (`scripts/04_microbiome_differential_abundance/`)

| Script Path | Lang | Purpose | Manuscript Section | Figure / Table | Input Files | Output Files | Restricted Input? | Status |
|-------------|------|---------|--------------------|----------------|-------------|--------------|-------------------|--------|
| `scripts/04_microbiome_differential_abundance/run_microbiome_covariate_models.R` | R | **Main microbiome DA** — ANCOM-BC2 + ALDEx2 (128 MC), Models A–D across taxonomic levels | Methods — Microbiome differential abundance | Fig. 3d–g, Supp. Table 4 | `asv_genus_table_reads.csv`, `clinical_merged_371.csv` | `microbiota_da_results_final.csv` (expected output), Venn plots, forest plots | **YES** | **Final** |
| `scripts/04_microbiome_differential_abundance/run_models.bat` | Batch | Orchestration wrapper script executing `run_microbiome_covariate_models.R` | Methods — Orchestration | — | R environment | Model execution logs | **NO** | **Helper / Orchestration** |

---

## 6. Protein Network (`scripts/05_protein_network/`)

| Script Path | Lang | Purpose | Manuscript Section | Figure / Table | Input Files | Output Files | Restricted Input? | Status |
|-------------|------|---------|--------------------|----------------|-------------|--------------|-------------------|--------|
| `scripts/05_protein_network/run_proteomics_network_analysis.py` | Python | **Protein co-expression network** — Pearson r, Leiden clustering (resolution=1.0, seed=42), module eigengene calculation | Methods — Protein network | Fig. 4a, 4c | `data_log2_mapping_deduplicated.csv`, `clinical_merged_371.csv` | `ME_proteomics_ALL.csv` (expected output), GraphML files | **YES** | **Final** |

---

## 7. Microbial Network (`scripts/06_microbial_network/`)

| Script Path | Lang | Purpose | Manuscript Section | Figure / Table | Input Files | Output Files | Restricted Input? | Status |
|-------------|------|---------|--------------------|----------------|-------------|--------------|-------------------|--------|
| `scripts/06_microbial_network/analyze_network_topology.py` | Python | **Microbial network topology** — FastSpar matrices, Leiden community detection (resolution=1.1, seed=42), module eigengenes | Methods — Microbial network | Fig. 4b, 4d, Supp. Table 4 | FastSpar matrices (`median_correlation.tsv`), `asv_genus_table.csv`, `clinical_merged_371.csv` | `microbial_network_topology.csv` (expected output), GraphML files | **YES** | **Final** |

---

## 8. Classification and CISS (`scripts/08_classification_and_CISS/`)

| Script Path | Lang | Purpose | Manuscript Section | Figure / Table | Input Files | Output Files | Restricted Input? | Status |
|-------------|------|---------|--------------------|----------------|-------------|--------------|-------------------|--------|
| `scripts/08_classification_and_CISS/clinical_vs_omics_feature.py` | Python | **Primary pairwise classification** — L2 logistic regression under 5-fold CV × 20 repeats (Clinical-only, Multi-omics Only, Integrated) | Methods — Classification | Fig. 5a, 5b, Supp. Table 5 | `data_log2_mapping_deduplicated.csv`, `asv_genus_table.csv`, `clinical_merged_371.csv` | `master_comparison.csv` (expected output), ROC curve PNGs | **YES** | **Final** |
| `scripts/08_classification_and_CISS/feature_ratio_modeling_cont_clin.py` | Python | **Primary CISS construction** — Selects concordant features across pairwise models, builds domain CISS scores, exports ratio features | Methods — CISS construction | Fig. 5c, 5d, 5e, 5f | `data_log2_mapping_deduplicated.csv`, `asv_genus_table.csv`, `clinical_merged_371.csv` | `sample_ratio_features.csv` (expected output), `ratio_moca_spearman.csv` (expected output) | **YES** | **Final** |
| `scripts/08_classification_and_CISS/plot_ciss_boxplot_moca_correlation.py` | Python | **Companion CISS plotting** — Renders domain-specific CISS boxplots by phenotype and CISS vs MoCA Spearman correlations | Methods — CISS visualization | Fig. 5d, 5f | `sample_ratio_features.csv`, `ratio_moca_spearman.csv` | `ciss_boxplots.png`, `ciss_moca_correlations.png` (expected outputs) | **YES** | **Final / Companion Plotting** |
| `scripts/08_classification_and_CISS/binary_classification_feature.py` | Python | Precursor 5-fold CV feature-level L2 logistic regression classifier | Methods — Sensitivity analysis | — | `clinical_merged_371.csv`, omics matrices | `feature_classification_metrics.csv` (expected output) | **YES** | **Legacy / Sensitivity** |
| `scripts/08_classification_and_CISS/binary_classification_module_feature.py` | Python | Precursor module-level & combined feature classification experiment | Methods — Sensitivity analysis | — | `clinical_merged_371.csv`, module eigengene matrices | `module_feature_classification_metrics.csv` (expected output) | **YES** | **Legacy / Sensitivity** |
| `scripts/08_classification_and_CISS/ratio_modeling_da.py` | Python | Precursor DA-feature ratio modeling experiment | Methods — CISS construction | — | `clinical_merged_371.csv`, DA features | `da_ratio_modeling_results.csv` (expected output) | **YES** | **Legacy / Sensitivity** |

---

## 9. Serial Mediation (`scripts/09_serial_mediation/`)

| Script Path | Lang | Purpose | Manuscript Section | Figure / Table | Input Files | Output Files | Restricted Input? | Status |
|-------------|------|---------|--------------------|----------------|-------------|--------------|-------------------|--------|
| `scripts/09_serial_mediation/sequential_mediation.py` | Python | **Primary serial mediation** — Forward & reverse models (Group → Microbiome CISS → GBB → Protein CISS → MoCA), 5000 bootstrap iterations | Methods — Serial mediation | Fig. 6a–d, Supp. Table 6 | `sample_ratio_features.csv`, `clinical_merged_371.csv`, metabolite data | `mediation_results_summary.csv` (expected output), path diagram PNGs | **YES** | **Final** |
| `scripts/09_serial_mediation/sequential_mediation_strategies_all_metabolites.py` | Python | **Sensitivity serial mediation** — Evaluates all candidate metabolite mediators and strategy routing configurations | Methods — Mediation sensitivity | Fig. 6b, Supp. Fig. 11 | `sample_ratio_features.csv`, `clinical_merged_371.csv`, metabolite matrix | `mediation_sensitivity_all_metabolites.csv` (expected output) | **YES** | **Sensitivity / Helper** |

---

## 10. Figure Generation (`scripts/10_figure_generation/`)

| Script Path | Lang | Purpose | Manuscript Section | Figure / Table | Input Files | Output Files | Restricted Input? | Status |
|-------------|------|---------|--------------------|----------------|-------------|--------------|-------------------|--------|
| `scripts/10_figure_generation/plot_alpha_diversity.py` | Python | Phylum & genus alpha diversity (Observed, Shannon, Chao1, Pielou) + rarefaction | Methods — Diversity | Fig. 2d, 2e | Microbiome abundance tables | Alpha diversity PNG figures | **YES** | **Final** |
| `scripts/10_figure_generation/plot_bray_curtis_pcoa.py` | Python | Full-cohort Bray-Curtis PCoA beta diversity & PERMANOVA | Methods — Diversity | Fig. 2f | `asv_genus_table.csv`, `clinical_merged_371.csv` | `pcoa_bray_curtis.png` | **YES** | **Final** |
| `scripts/10_figure_generation/plot_cognitive_scores.py` | Python | MoCA subtest score distributions across C1/C2/C3 | Methods — Cognitive phenotyping | Fig. 2c | `clinical_merged_371.csv` | Cognitive score distribution PNGs | **YES** | **Final** |
| `scripts/10_figure_generation/plot_combined_network_C123.py` | Python | Combined multi-omics network topology across C1/C2/C3 | Methods — Network analysis | Fig. 4a, 4b | Multi-omics GraphML files | Combined network PNG figure | **YES** | **Final** |
| `scripts/10_figure_generation/plot_da_heatmap.py` | Python | Heatmap of differential microbial taxa by group | Methods — Differential abundance | Fig. 3d | `microbiota_da_results_final.csv` | Microbiome DA heatmap PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_da_heatmap_proteomics.py` | Python | Heatmap of differential plasma proteins by group | Methods — Proteomics association | Fig. 3a | `proteomics_ancova_detailed_results.csv` | Proteomics DA heatmap PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_leiden_24_taxa.py` | Python | Abundance profiles of 24 key taxa in microbial Leiden modules | Methods — Microbial network | Fig. 4d | Microbial module assignments | 24-taxa profile PNG figure | **YES** | **Final** |
| `scripts/10_figure_generation/plot_leiden_me_loadings.py` | Python | Leiden module eigengene PCA loadings | Methods — Network analysis | Supp. Fig. 9 | Module eigengenes | Module loadings PNG plot | **YES** | **Final** |
| `scripts/10_figure_generation/plot_me_moca_correlation.py` | Python | Module eigengene × MoCA score correlation analysis | Methods — Cross-omics integration | Supplementary analysis — exact figure/panel to be confirmed | Module eigengenes, MoCA scores | ME-MoCA correlation PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_mediation_sensitivity_forest.py` | Python | Forest plot of serial mediation direct/indirect effects across sensitivity models | Methods — Serial mediation | Fig. 6b | Mediation results CSV | Mediation sensitivity forest PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_metabolite_violins_heatmap_v2.py` | Python | Targeted metabolite violin plots & CISS correlation heatmap | Methods — Targeted metabolomics | Fig. 6a | Targeted metabolite matrix, CISS scores | Metabolite violins & heatmap PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_moca_boxplot.py` | Python | Total MoCA score boxplots across C1/C2/C3 | Methods — Cognitive phenotyping | Fig. 2b | `clinical_merged_371.csv` | MoCA total boxplot PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_moca_group_comparisons.py` | Python | Pairwise MoCA subtest group comparisons | Methods — Cognitive phenotyping | Fig. 2c detail | `clinical_merged_371.csv` | Subtest comparison PNGs | **YES** | **Final** |
| `scripts/10_figure_generation/plot_module_dim_reduction.py` | Python | UMAP/PCA visualization of multi-omics module eigengenes | Methods — Network analysis | Fig. 4c, 4d | Module eigengenes | UMAP/PCA module PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_network_comparisons.py` | Python | Network topology metric comparisons across groups | Methods — Network topology | Supp. Fig. 8 | Topology CSV files | Network comparison PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_path_diagram_4models.py` | Python | Serial mediation path diagrams across Models 1–4 | Methods — Serial mediation | Fig. 6c, 6d | Mediation path coefficients | Mediation path diagram PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_pca_cognitive.py` | Python | PCA of MoCA subtest scores colored by phenotype | Methods — Cognitive phenotyping | Fig. 2a | `clinical_merged_371.csv` | Cognitive PCA PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_proteomics_modules.py` | Python | Protein module eigengenes across C1/C2/C3 | Methods — Protein network | Fig. 4c | `ME_proteomics_ALL.csv` | Protein module eigengene PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_proteomics_pca.py` | Python | Full-cohort plasma proteomics PCA plot | Methods — Proteomics association | Supp. Fig. 3 | Proteomics intensity matrix | Proteomics PCA PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_psm_networks.py` | Python | PSM-matched microbial network plots | Methods — Propensity score matching | Supp. Fig. 10 | PSM FastSpar matrices | PSM microbial network PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_psm_networks_proteomics.py` | Python | PSM-matched protein network plots | Methods — Propensity score matching | Supp. Fig. 10 | PSM protein network files | PSM protein network PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_venn_ABCD_labeled.py` | Python | Venn diagram of candidate feature overlaps across Models A–D | Methods — Differential abundance | Fig. 3b, 3e | Model A–D results | Venn diagram PNG | **YES** | **Final** |
| `scripts/10_figure_generation/plot_venn_intersection.py` | Python | Cohort multi-omics sample intersection Venn diagram | Methods — Cohort summary | Methods | Cohort sample IDs | Sample intersection Venn PNG | **YES** | **Final** |
