# Final Release Audit

> `nature_aging_code_release` repository — Pre-submission audit
> Updated: 2026-09-09

---

## READY ✅

These items have been completed and verified.

| # | Item | Details |
|---|------|---------|
| R1 | All Python scripts import from `config/paths.py` | 46/46 Python scripts compiled cleanly; zero hard-coded personal paths |
| R2 | R script uses env var `TVGH_DATA_DIR` or repo-relative auto-detection | `scripts/04_microbiome_differential_abundance/run_microbiome_covariate_models.R` |
| R3 | Zero hard-coded personal paths remaining in executable code or public docs | Verified by automated audit (2026-09-09) |
| R4 | CISS naming unified as "Cognitive Impairment Signature Score (CISS)" | Zero occurrences of old legacy title across repository |
| R5 | Final classification script confirmed & documented | `scripts/08_classification_and_CISS/clinical_vs_omics_feature.py` |
| R6 | Final CISS construction script confirmed & documented | `scripts/08_classification_and_CISS/feature_ratio_modeling_cont_clin.py` |
| R7 | Companion CISS visualization script added | `scripts/08_classification_and_CISS/plot_ciss_boxplot_moca_correlation.py` |
| R8 | Legacy scripts marked with Status: LEGACY headers | `binary_classification_feature.py`, `binary_classification_module_feature.py`, `ratio_modeling_da.py` |
| R9 | FastSpar marked as external server execution | `scripts/06_microbial_network/README.md` |
| R10 | `CODE_AVAILABILITY.md` uses accurate conservative manuscript wording | Participant data and FastSpar server execution clearly scoped |
| R11 | `FIGURE5_CODE_VALIDATION.md` updated | Line-by-line validation of Figure 5a–f and Supplementary Table 5 |
| R12 | `CISS_CODE_VALIDATION.md` updated | All 10 CISS construction steps confirmed |
| R13 | `MANUSCRIPT_CODE_MAP.md` updated | Figures 1–6 and Supplementary Tables 1–5 mapped to current manuscript |
| R14 | `GITHUB_RELEASE_CHECKLIST.md` updated | Comprehensive checklist with discrepancy resolutions |
| R15 | Repository structure stable; no unnecessary top-level folders added | — |
| R16 | `.gitignore` excludes participant data, credentials, and personal paths | Verified |

---

## RESOLVED ✅ (Previously Flagged, Now Closed)

| ID | Prior Issue | Resolution |
|----|------------|------------|
| D1 | OUTER_SPLITS=4 in legacy script vs manuscript "5-fold" | **RESOLVED** — Final production scripts (`clinical_vs_omics_feature.py` & `feature_ratio_modeling_cont_clin.py`) use `OUTER_SPLITS=5`, `OUTER_REPEATS=20`. |
| Feature Mapping | Config mapping for Multi-omics Only vs Integrated | **RESOLVED** — `clinical` = Clinical-only, `combined` = Multi-omics Only (Micro + Prot), `clin_combined` = Integrated (Clin + Micro + Prot). |
| Supplementary Table 5 | Exact numerical match provenance | **RESOLVED** — 100% exact numerical match confirmed between `SUPPLEMENTARY.docx` Table 5 and `clinical_vs_omics_feature.py` → `master_comparison.csv`. |
| CISS Construction | CISS construction script location | **RESOLVED** — Confirmed final implementation in `feature_ratio_modeling_cont_clin.py`. |
| FastSpar Status | FastSpar script distribution | **RESOLVED** — Reclassified as EXTERNAL SERVER EXECUTION. Not distributed as repository script. |

---

## Summary Answers

| Question | Answer |
|----------|--------|
| **A. Figure 5 final scripts** | ✅ **CONFIRMED** — `clinical_vs_omics_feature.py` (Fig 5a/5b & Supp Table 5) + `feature_ratio_modeling_cont_clin.py` (CISS construction & Fig 5c/5d/5e/5f data) |
| **B. Supplementary Table 5 provenance** | ✅ **100% EXACT NUMERIC MATCH** — Produced by `clinical_vs_omics_feature.py` (`master_comparison.csv`), matching `SUPPLEMENTARY.docx` Table 5 to every decimal place. |
| **C. CISS final implementation** | ✅ **CONFIRMED** — All 10 manuscript CISS steps implemented in `feature_ratio_modeling_cont_clin.py`; see `CISS_CODE_VALIDATION.md` |
| **D. Final outer CV folds** | **5-fold × 20 repeats** (`OUTER_SPLITS=5`, `OUTER_REPEATS=20`, `RANDOM_STATE=42`) |
| **E. Final inner CV folds** | **4-fold** (`INNER_SPLITS=4`, `StratifiedKFold`) |
| **F. FastSpar status** | **EXTERNAL SERVER EXECUTION** — Not distributed; downstream code and I/O specification provided in `scripts/06_microbial_network/README.md` |
| **G. K-means phenotyping source** | **UPSTREAM PRE-CLUSTERED FILE** — `k-means_3group.csv` is a required restricted upstream input and is not included in the public repository; imported by `prepare_microbiome_metadata.py`. |
| **H. Repository ready for public GitHub** | ✅ **READY FOR PUBLIC GITHUB RELEASE** — Core analysis code is 100% clean, all 46 Python scripts pass compilation, zero hardcoded paths. LICENSE finalized: Copyright (c) 2026 National Yang Ming Chiao Tung University. Core downstream analysis provenance and Figure 5/Supplementary Table 5 provenance have been verified; selected upstream and figure-rendering scripts were not retained. |
