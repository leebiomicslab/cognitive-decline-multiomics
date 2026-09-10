# Data Directory

> **No participant-level data are distributed with this repository.**

This directory is a placeholder. All omics and clinical data used in the study are
subject to institutional data governance and ethics requirements and cannot be made
publicly available without explicit approval from the relevant ethics board and
data access committee.

---

## What data files are required to run the analysis?

The scripts in this repository expect the following input files to be placed
in an accessible data directory (configured via `config/paths.example.py`).

### 1. Clinical / Demographic Data

| File | Format | Description |
|------|--------|-------------|
| `clinical_merged_371.csv` | CSV | Merged clinical covariates for n = 371 participants with multi-omics overlap. Columns: `id`, `age`, `gender`, `edu`, `bmi`, `smoking`, `alcohol`, `HTN`, `DM`, `CVD`, `CRP`, `egfr`, `group`, `Statin`, `dm_med`, `htn_med`, `chol_med`, `PPI`, `Antibiotics` |
| `vascular_with_medication.xlsx` | Excel | Source medication data (subset of columns) |
| `k-means_3group.csv` | CSV | K-means cluster assignments. Columns: `id`, `group` (values: C1, C2, C3) |
| `common_samples.csv` | CSV | Participant IDs present in both proteomics and microbiome datasets. Column: `id` |
| `Cognitive.xlsx` | Excel | MoCA subtest raw scores. Columns: `id`, `trail`, `cube`, `clock_shape`, `clock_number`, `clock_time`, `naming`, `con1`, `con2`, `con3`, `lang1`, `lang2`, `abstract`, `delay`, `moca_1` |

### 2. Proteomics

| File | Format | Description |
|------|--------|-------------|
| `PROTEOMICS/data_log2_mapping_deduplicated.csv` | CSV | Log2-transformed protein intensities. Rows = participant IDs, Columns = protein identifiers. Produced by PEAKS Studio 11.5 and de-duplicated. |

### 3. Microbiome

| File | Format | Description |
|------|--------|-------------|
| `MICROBIOTA/asv_genus_table.csv` | CSV | Genus-level relative abundance table (genus × sample). Produced by nf-core/ampliseq v2.11.0 with SILVA v138.1. |
| `MICROBIOTA/asv_genus_table_reads.csv` | CSV | Genus-level raw read count table (genus × sample). |

### 4. Metabolomics

| File | Format | Description |
|------|--------|-------------|
| `METABOLITE/clean_metabolite_matrix_all_features_final_v4.csv` | CSV | Targeted metabolite feature matrix. Rows = metabolite features, Columns = participant IDs. Produced by Progenesis QI / statTarget preprocessing. |

---

## Schema for Dummy / Example Data

For code testing purposes, researchers may supply synthetic data matching these schemas.
No real participant data should be used without ethics approval.

### `clinical_merged_371.csv` example rows

```csv
id,group,age,gender,edu,bmi,smoking,alcohol,HTN,DM,CVD,CRP,egfr,Statin,dm_med,htn_med,chol_med,PPI,Antibiotics
E0001,C3,68,1,12,24.5,0,0,1,0,0,1.2,72.3,0,0,1,0,0,0
E0002,C2,72,2,9,26.1,0,1,1,1,0,2.1,65.8,1,1,1,0,0,0
```

- `id`: participant identifier (de-identified string, e.g. "E0001")
- `group`: C1, C2, or C3 (cognitive phenotype)
- `gender`: 1 = male, 2 = female
- `smoking`, `alcohol`, `HTN`, `DM`, `CVD`, `Statin`, `dm_med`, `htn_med`, `chol_med`, `PPI`, `Antibiotics`: binary (0/1)

---

## Data Availability

Participant-level clinical data, proteomics matrices, and microbiome tables
from this study **are not deposited in a public repository** at this time,
due to IRB restrictions and patient privacy requirements.

Researchers interested in accessing the data for replication should contact the
corresponding author. Data may be shared under a data use agreement following
institutional review.

Aggregate summary statistics sufficient to reproduce key figures may be
provided upon reasonable request.
