# Microbial Network Analysis

## Overview

This directory contains downstream analysis code for microbial co-occurrence network characterisation.
Network construction (FastSpar correlation estimation and bootstrap P-value calculation) was performed
on an institutional computing server and is **not** distributed in this repository.

---

## FastSpar Execution (External Server — Not Distributed)

FastSpar correlation estimation and bootstrap-based empirical P-value calculation were executed on an
institutional computing server. Server-specific execution scripts and scheduler scripts (e.g., SLURM
batch files) are **not** distributed in this repository.

The repository provides:
1. Downstream microbial network analysis code (topology, module detection, visualization)
2. Input specification for FastSpar
3. Expected FastSpar output structure
4. Methods-consistent description of the procedure
5. Software citation and version information (where confirmed)

**Do not assume that the raw-to-final FastSpar workflow can be reproduced from this repository alone.**

---

## FastSpar: Input Specification

FastSpar operates on a **count table** (samples × genera) in TSV/CSV format.

- **Input file:** ASV genus-level count table (pre-rarefied or filtered)
  - Expected path: `data/MICROBIOTA/asv_genus_table_reads.csv`
  - Format: rows = genera, columns = samples (or transposed; verify with your FastSpar version)
- **Zero handling:** FastSpar uses a compositional framework internally; zero imputation is not required
  but low-abundance filtering is recommended (minimum prevalence threshold applied in preprocessing)

---

## FastSpar: Parameters

The following parameters were used. Where a parameter could not be confirmed from a server run log,
it is marked **TO BE CONFIRMED FROM SERVER RUN LOG**.

| Parameter | Value | Confirmation Status |
|-----------|-------|---------------------|
| Software | FastSpar | ✅ Confirmed (downstream code reads FastSpar output format) |
| Version | TO BE CONFIRMED FROM SERVER RUN LOG | ⚠️ |
| Number of iterations | TO BE CONFIRMED FROM SERVER RUN LOG | ⚠️ |
| Number of bootstrap replicates | TO BE CONFIRMED FROM SERVER RUN LOG | ⚠️ |
| Minimum sample count filter | TO BE CONFIRMED FROM SERVER RUN LOG | ⚠️ |
| Empirical P-value threshold | < 0.05 (confirmed from downstream analysis code) | ✅ |
| Random seed | TO BE CONFIRMED FROM SERVER RUN LOG | ⚠️ |

---

## FastSpar: Expected Output Structure

The downstream analysis scripts in this directory expect the following FastSpar output files:

```
results/
└── fastspar/
    └── microbiota/
        ├── correlations.tsv      # Spearman-like correlation matrix (genera × genera)
        └── pvalues.tsv           # Empirical P-values from bootstrap resampling
```

---

## FastSpar: Citation

> Watts SC, Ritchie SC, Inouye M, Holt KE. FastSpar: rapid and scalable correlation estimation
> for compositional data. *Bioinformatics*. 2019;35(6):1064–1066.
> https://doi.org/10.1093/bioinformatics/bty734

---

## Downstream Analysis Scripts (Provided in This Repository)

| Script | Description | Output |
|--------|-------------|--------|
| `analyze_network_topology.py` | Network topology metrics (degree, betweenness, modularity, Leiden community detection) using FastSpar correlation + P-value matrices | `results/microbiota_network/` |
| `plot_da_heatmap.py` → `scripts/10_figure_generation/` | Heatmap of DA taxa by group | `results/figures/` |
| `plot_leiden_24_taxa.py` → `scripts/10_figure_generation/` | Leiden module visualisation for 24 DA taxa | `results/figures/` |
| `plot_leiden_me_loadings.py` → `scripts/10_figure_generation/` | Module eigengene PCA loadings | `results/figures/` |

---

## Preprocessing Notes

The microbiome preprocessing pipeline (QIIME2 → DADA2 → taxonomic assignment) was performed with
proprietary institutional bioinformatics infrastructure and is not included in this repository.
The downstream analysis assumes the following preprocessed input files are available:

- `data/MICROBIOTA/asv_genus_table.csv` — CLR-ready count table (samples × genera)
- `data/MICROBIOTA/asv_genus_table_reads.csv` — raw read counts for FastSpar input
