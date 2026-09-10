# GitHub Release Checklist
> Multi-omics cognitive aging manuscript — Code Availability package
> Last updated: 2026-09-09

---

## ✅ Completed Items

| # | Task | Status |
|---|------|--------|
| 1 | Repository skeleton created (`scripts/`, `config/`, `data/`, `results/`, `docs/`) | ✅ |
| 2 | `config/paths.py` — single-file path configuration, no hard-coded personal paths | ✅ |
| 3 | `scripts/_path_helper.py` — portable import bootstrap for all scripts | ✅ |
| 4 | All 40+ Python scripts refactored to import from `config/paths.py` | ✅ |
| 5 | R script (`run_microbiome_covariate_models.R`) refactored; now reads from env var `TVGH_DATA_DIR` or auto-detects repo structure | ✅ |
| 6 | Automated audit: **zero** hard-coded personal absolute paths remain in executable code | ✅ |
| 7 | CISS terminology unified as **"Cognitive Impairment Signature Score (CISS)"** throughout all scripts and documentation | ✅ |
| 8 | Standardised docstring header added to all core analysis scripts (Purpose, Manuscript section, Figure/Table, Input, Output, Dependencies) | ✅ |
| 9 | `CODE_INVENTORY.md` — full table mapping each script to manuscript figure/table | ✅ |
| 10 | `MISSING_CODE_REPORT.md` — lists unresolved custom scripts | ✅ |
| 11 | `data/README.md` — describes expected data structure; no raw data committed | ✅ |
| 12 | FastSpar documented as external server execution in `scripts/06_microbial_network/README.md` | ✅ |
| 13 | Final classification (`clinical_vs_omics_feature.py`) and CISS construction (`feature_ratio_modeling_cont_clin.py`) scripts confirmed | ✅ |

---

## ⚠️ Known Discrepancies & Status

### D1 — Nested CV Outer Splits (RESOLVED)
| Field | Value |
|-------|-------|
| **File** | `scripts/08_classification_and_CISS/clinical_vs_omics_feature.py` & `feature_ratio_modeling_cont_clin.py` |
| **Status** | **RESOLVED** — Final production scripts use `OUTER_SPLITS = 5` (5-fold × 20 repeats = 100 outer splits), consistent with manuscript. |

### D2 — Cognitive phenotype K-means k=3 source
| Field | Value |
|-------|-------|
| **File** | `scripts/01_cognitive_phenotyping/` |
| **Issue** | K-means clustering script (`k-means_3group.csv` generator) is **unresolved** |
| **Impact** | `k-means_3group.csv` is a required restricted upstream input and is not included in the public repository |
| **Action** | Document as upstream input file requirement |

### D3 — FastSpar Server Execution Parameters
| Field | Value |
|-------|-------|
| **File** | `scripts/06_microbial_network/README.md` |
| **Issue** | FastSpar executed externally on computing server; exact version/bootstrap parameters marked `TO BE CONFIRMED FROM SERVER RUN LOG` |
| **Action** | Confirm server run log parameters before publication |

---

## 📋 TODO Before Final Public Release

- [ ] **License**: Add `LICENSE` file. Recommended: MIT or Apache 2.0.
- [ ] **README.md**: Update top-level README with final bioRxiv / journal DOI once assigned.
- [ ] **Data Availability**: Add dataset accession numbers / institutional access instructions.

---

## 🔒 Privacy & Data Safety Audit (Passed)

Run date: 2026-09-09

```
Zero personal absolute paths in executable code ✅
Zero participant-level raw data committed ✅
Zero credentials or credentials patterns found ✅
```

---

## 📁 Repository Structure Summary

```
nature_aging_code_release/
├── config/
│   ├── paths.py            ← single path config (set BASE_DIR here)
│   └── paths.example.py   ← template for new users
├── data/
│   └── README.md           ← data structure documentation (no raw data)
├── scripts/
│   ├── _path_helper.py    ← portable import bootstrap
│   ├── 00_microbiome_preprocessing/
│   ├── 01_cognitive_phenotyping/
│   ├── 02_psm_and_microbiome_diversity/
│   ├── 03_proteomics_association/
│   ├── 04_microbiome_differential_abundance/
│   ├── 05_protein_network/
│   ├── 06_microbial_network/
│   ├── 07_cross_omics_module/
│   ├── 08_classification_and_CISS/
│   ├── 09_serial_mediation/
│   └── 10_figure_generation/
├── results/                ← generated at runtime (gitignored)
├── docs/
│   ├── MANUSCRIPT_CODE_MAP.md
│   ├── FIGURE5_CODE_VALIDATION.md
│   └── CISS_CODE_VALIDATION.md
├── CODE_AVAILABILITY.md
├── CODE_INVENTORY.md
├── MISSING_CODE_REPORT.md
├── GITHUB_RELEASE_CHECKLIST.md  ← this file
├── FINAL_RELEASE_AUDIT.md
└── README.md
```
