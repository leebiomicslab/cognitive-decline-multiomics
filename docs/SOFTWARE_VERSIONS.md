# Software Environment

## Python

| Component | Version used | Notes |
|-----------|-------------|-------|
| Python    | 3.10        | Conda-managed |
| numpy     | ≥ 1.24      | |
| pandas    | ≥ 2.0       | |
| scipy     | ≥ 1.10      | |
| scikit-learn | ≥ 1.3   | |
| statsmodels | ≥ 0.14   | OLS, HC3 robust SE, marginal means |
| matplotlib | ≥ 3.7      | |
| seaborn   | ≥ 0.12      | |
| networkx  | ≥ 3.1       | |
| python-igraph | ≥ 0.10  | igraph Python bindings |
| leidenalg | ≥ 0.10      | Leiden community detection (RBConfigurationVertexPartition) |
| openpyxl  | ≥ 3.1       | |
| adjustText | ≥ 0.8     | Label placement for proteomics volcano plots |

> **Note**: Exact patch versions used during analysis are to be confirmed from the original
> execution environment (conda `conda env export`). The versions listed above are
> the minimum confirmed to work.

---

## R

The following R packages are required for microbiome differential abundance (Step 04):

| Package   | Version used        | Notes |
|-----------|---------------------|-------|
| R         | ≥ 4.3.0             | |
| tidyverse | version to be confirmed | |
| phyloseq  | version to be confirmed | |
| ANCOMBC   | version to be confirmed | ANCOM-BC2 (`ancombc2()`) |
| ALDEx2    | version to be confirmed | 128 Monte Carlo instances |
| VennDiagram | version to be confirmed | |
| emmeans   | version to be confirmed | Adjusted marginal means for CLR boxplots |

### Installing R packages

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("phyloseq", "ANCOMBC", "ALDEx2"))
install.packages(c("tidyverse", "VennDiagram", "emmeans"))
```

---

## FastSpar (Microbial Co-abundance Networks)

FastSpar is used for estimating sparse correlations in microbial data via bootstrapping.

- **Version**: to be confirmed from the execution log
- **Repository**: https://github.com/scwatts/fastspar
- **License**: GNU GPL v3

### Installation (Conda — recommended)

```bash
conda install -c bioconda fastspar
```

### Usage summary

The FastSpar runs that produced the network inputs (`results/fastspar/microbiota/`) were
executed as follows (commands to be confirmed from execution logs):

```bash
# Example — actual parameters to be verified
fastspar --otu_table <clr_genus_table.tsv> \
         --correlation <correlation_matrix.tsv> \
         --covariance <covariance_matrix.tsv> \
         --bootstrap_number <N> \
         --threads <N>
```

> The GraphML files in `results/fastspar/` are **intermediate outputs** from FastSpar
> that serve as input to downstream network analysis scripts.
> They are **not** participant-level raw data.

---

## Upstream Preprocessing Software (Proprietary / External Pipelines)

These tools produced the omics matrices used as inputs to this repository.
Their source code is **not** included here.

### Microbiome 16S rRNA

| Software | Version | Reference |
|----------|---------|-----------|
| nf-core/ampliseq | 2.11.0 | https://nf-co.re/ampliseq |
| Cutadapt | version to be confirmed | Adapter trimming |
| DADA2 | version to be confirmed | ASV calling |
| SILVA reference | v138.1 | Taxonomic classification |
| Docker | version to be confirmed | Container runtime |

**Command** (to be confirmed from execution log):
```bash
nextflow run nf-core/ampliseq \
  -r 2.11.0 \
  --input samplesheet.csv \
  --FW_primer <forward_primer> \
  --RV_primer <reverse_primer> \
  --dada_ref_taxonomy silva=138.1 \
  --outdir results/ \
  -profile docker
```

### Proteomics (MS/MS)

| Software | Version |
|----------|---------|
| PEAKS Studio | 11.5 |

Downstream analysis begins from the log2-transformed protein intensity matrix output by PEAKS.

### Targeted Metabolomics

| Software | Version |
|----------|---------|
| Progenesis QI | version to be confirmed |
| statTarget | version to be confirmed |

---

## Operating System

Analyses were performed on Windows 11 with Anaconda Python 3.10.
R analyses were executed via the R terminal.
FastSpar was executed in WSL2 (Ubuntu) or equivalent Linux environment.
