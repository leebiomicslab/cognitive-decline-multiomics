# Missing Code Report

This document identifies analysis steps described in the manuscript for which
**custom code was not located** in the repository at the time of code release preparation.

Do NOT create placeholder or synthetic scripts for these steps.
They require the actual source code or confirmation of upstream input status.

---

## UNRESOLVED CUSTOM CODE

### MISSING_001 — K-means Cognitive Phenotyping Script

**Manuscript section**: Methods — Cognitive phenotype discovery

**What is missing**:
The K-means clustering script that performs:
- MoCA subtest preprocessing and standardization (`StandardScaler`)
- K-means clustering with k = 2–6 evaluation
- Silhouette score, Calinski–Harabasz index, Davies–Bouldin index computation
- Selection of final k = 3
- PCA visualization of clusters
- Generation of `k-means_3group.csv` (the cluster assignment file)

**Current Status**:
`k-means_3group.csv` is a required restricted upstream input and is not included in the public repository. Downstream characterization scripts (`plot_pca_cognitive.py`, `plot_cognitive_scores.py`, `plot_moca_boxplot.py`) use these precomputed cluster assignments. The original K-means cluster-generation script was not retained.

---

### MISSING_002 — Cross-omics ME×ME Correlation Heatmap Script (Figure 4e)

**Manuscript section**: Methods — Cross-omics module analysis

**What is missing**:
A dedicated script to compute the protein module eigengene × microbiome module eigengene
Pearson correlations stratified by C1/C2/C3 and render the Figure 4e heatmap.

**Current Status**:
Output files exist in original analysis outputs (`microbiome_x_proteomics_ME_r_C*.csv`, `microbiome_x_proteomics_ME_pval_C*.csv`), but the dedicated rendering script for Figure 4e has not been located.

---

### MISSING_003 — Figure 5e Multi-Omics Landscape Rendering Script

**Manuscript section**: Methods — CISS landscape analysis

**What is missing**:
A custom rendering script that plots the 3D / 2D landscape of Microbiome CISS × Proteomic CISS with MoCA Z-score overlay for Figure 5e.

**Current Status**:
The underlying data file `sample_ratio_features.csv` is an expected/generated output exported by `feature_ratio_modeling_cont_clin.py`, but the exact rendering script for panel 5e has not been located.

---

### MISSING_004 — PSM Network Reconstruction Code

**Manuscript section**: Supplementary — PSM network sensitivity

**What is missing**:
A standalone script that rebuilds protein and microbial networks using only the PSM-matched cohort (n = 49 per phenotype) and computes network topology.

**Current Status**:
Downstream plotting scripts (`plot_psm_networks.py` and `plot_psm_networks_proteomics.py`) exist for visualization, but the network reconstruction pipeline for PSM samples is not confirmed.

---

### MISSING_005 — Microbiome Rarefaction Depth Confirmation

**Manuscript section**: Methods — Alpha diversity

**What is missing**:
`plot_alpha_diversity.py` implements rarefaction dynamically using the minimum observed read depth across samples. The exact fixed rarefaction depth integer (if manually set) is unconfirmed.

---

## RESOLVED ITEMS ✅

### RESOLVED_001 — CISS Construction Script

**Confirmed Final Script**: [`scripts/08_classification_and_CISS/feature_ratio_modeling_cont_clin.py`](scripts/08_classification_and_CISS/feature_ratio_modeling_cont_clin.py)

**Status**: **RESOLVED**  
All 10 manuscript CISS construction steps (concordant coefficient selection across pairwise C1vsC2 and C1vsC3, domain split into Clinical / Microbiome / Proteomics, standardization, positive mean minus negative mean) are fully implemented and verified in `feature_ratio_modeling_cont_clin.py`.

---

### RESOLVED_002 — FastSpar Execution

**Status**: **EXTERNAL SERVER EXECUTION — NOT DISTRIBUTED**  
FastSpar correlation estimation was executed on an institutional computing cluster. Server-specific execution wrapper scripts are not distributed as custom repository code. The repository provides downstream topology analysis scripts ([`scripts/06_microbial_network/analyze_network_topology.py`](scripts/06_microbial_network/analyze_network_topology.py)), expected FastSpar input/output specifications, and methodology documentation in [`scripts/06_microbial_network/README.md`](scripts/06_microbial_network/README.md).

---

## Summary Table

| ID | Analysis Step | Status / Action Owner |
|----|--------------|----------------------|
| **MISSING_001** | K-means cognitive phenotyping script | 🔴 UNRESOLVED — `k-means_3group.csv` is a required restricted upstream input and is not included in the public repository |
| **MISSING_002** | Cross-omics ME×ME heatmap (Figure 4e) | 🟡 UNRESOLVED — Outputs exist; rendering script not located |
| **MISSING_003** | Figure 5e CISS landscape rendering script | 🟡 UNRESOLVED — `sample_ratio_features.csv` is an expected/generated output |
| **MISSING_004** | PSM network reconstruction pipeline | 🟢 UNRESOLVED — Plotting scripts provided |
| **MISSING_005** | Rarefaction depth parameter | 🟢 UNRESOLVED — Dynamic minimum depth implemented |
| **RESOLVED_001**| CISS construction script | ✅ RESOLVED — `feature_ratio_modeling_cont_clin.py` |
| **RESOLVED_002**| FastSpar execution script | ✅ RESOLVED — Documented as external server execution |
