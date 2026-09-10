# paths.example.py
# ============================================================
# Template for configuring all data paths used by this pipeline.
# Copy this file to paths.py (never commit paths.py) and update
# BASE_DIR to your local data directory.
# ============================================================
from pathlib import Path

# ── Project layout ────────────────────────────────────────────
# This repository (code only)
REPO_DIR = Path(__file__).resolve().parent.parent

# Base directory where restricted participant-level data reside.
# This must be set by each user; it is NOT included in the repository.
# Example:  BASE_DIR = Path("/data/tvgh_cognitive")
BASE_DIR = Path("REPLACE_WITH_YOUR_DATA_PATH")

# ── Input files ───────────────────────────────────────────────
# Clinical data (merged, de-identified integer participant IDs)
CLINICAL_CSV          = BASE_DIR / "clinical_merged_371.csv"

# Medication covariates
MEDICATION_XLSX       = BASE_DIR / "vascular_with_medication.xlsx"

# Cognitive phenotype cluster assignments from K-means step
GROUP_CSV             = BASE_DIR / "k-means_3group.csv"

# Samples common to all three omics layers
COMMON_SAMPLES_CSV    = BASE_DIR / "common_samples.csv"

# Proteomics: log2-transformed participant × protein matrix
PROTEOMICS_LOG2       = BASE_DIR / "PROTEOMICS" / "data_log2_mapping_deduplicated.csv"

# Microbiome: genus-level raw read counts (genus × sample)
MICRO_GENUS_COUNTS    = BASE_DIR / "MICROBIOTA" / "asv_genus_table.csv"
MICRO_GENUS_READS     = BASE_DIR / "MICROBIOTA" / "asv_genus_table_reads.csv"

# Cognitive sub-test scores (MoCA domains)
COGNITIVE_XLSX        = BASE_DIR / "Cognitive.xlsx"

# Metabolomics feature matrix (log2 values)
METABOLOMICS_CSV      = BASE_DIR / "METABOLITE" / "clean_metabolite_matrix_all_features_final_v4.csv"

# ── FastSpar / network input (intermediate results) ───────────
FASTSPAR_MICRO_DIR    = REPO_DIR / "results" / "fastspar" / "microbiota"
FASTSPAR_PROT_DIR     = REPO_DIR / "results" / "fastspar" / "proteomics"

# ── Output directories ────────────────────────────────────────
RESULTS_DIR           = REPO_DIR / "results"
FIGURES_DIR           = RESULTS_DIR / "figures"
TABLES_DIR            = RESULTS_DIR / "tables"

for p in [RESULTS_DIR, FIGURES_DIR, TABLES_DIR]:
    p.mkdir(parents=True, exist_ok=True)
