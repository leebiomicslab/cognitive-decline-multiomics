# Environment Setup

## Python Environment

Install using Conda (recommended):

```bash
conda env create -f environment.yml
conda activate nature_aging_multiomics
```

Or using pip:

```bash
pip install -r requirements.txt
```

## R Environment

Install Bioconductor packages:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("phyloseq", "ANCOMBC", "ALDEx2"))
install.packages(c("tidyverse", "VennDiagram", "emmeans"))
```

Verify installation:

```r
library(ANCOMBC)
library(ALDEx2)
library(phyloseq)
sessionInfo()
```

## FastSpar

FastSpar is an external command-line tool for microbial correlation estimation.

```bash
# Install via Conda
conda install -c bioconda fastspar

# Verify
fastspar --help
```

See [docs/SOFTWARE_VERSIONS.md](../docs/SOFTWARE_VERSIONS.md) for version details.
See [scripts/00_microbiome_preprocessing/README.md](../scripts/00_microbiome_preprocessing/README.md)
for nf-core/ampliseq setup.
