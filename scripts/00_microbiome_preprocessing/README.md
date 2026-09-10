# 00 — Microbiome Upstream Preprocessing

This step covers the processing of raw 16S rRNA amplicon sequencing data
(FASTQ files) into genus-level feature tables used by all downstream analyses.

**This step is NOT run from this repository.**
The pipeline was executed externally using nf-core/ampliseq v2.11.0.

---

## Software

| Tool | Version | Purpose |
|------|---------|---------|
| nf-core/ampliseq | 2.11.0 | Workflow manager |
| Nextflow | ≥ 22.0 | Pipeline executor |
| Docker | version to be confirmed | Container runtime |
| Cutadapt | (via ampliseq) | Primer trimming |
| DADA2 | (via ampliseq) | ASV denoising and chimera removal |
| SILVA | v138.1 | Taxonomic classification reference |

## Input

- Raw paired-end FASTQ files from 16S rRNA V3-V4 amplicon sequencing
- A sample sheet (`samplesheet.csv`) mapping sample IDs to FASTQ file pairs
- Primer sequences (to be confirmed from lab protocol)

## Example command (to be confirmed from execution log)

```bash
nextflow run nf-core/ampliseq \
    -r 2.11.0 \
    --input samplesheet.csv \
    --FW_primer <forward_primer_sequence> \
    --RV_primer <reverse_primer_sequence> \
    --dada_ref_taxonomy silva=138.1 \
    --outdir results/ampliseq/ \
    --trunclenf <to_be_confirmed> \
    --trunclenr <to_be_confirmed> \
    -profile docker \
    -resume
```

> ⚠ The exact command including primer sequences, truncation lengths, and additional
> parameters must be confirmed from the execution log before public release.

## Output (relevant for downstream analysis)

The following files were produced and are used as inputs to the analysis scripts:

| File | Location in this repo | Description |
|------|-----------------------|-------------|
| `asv_genus_table.csv` | `data/MICROBIOTA/` | Genus-level relative abundance (genus × sample) |
| `asv_genus_table_reads.csv` | `data/MICROBIOTA/` | Genus-level raw read counts (genus × sample) |
| `ASV_tax.user.tsv` | `data/MICROBIOTA/` | Full taxonomy table (not used directly in analysis) |

> These files contain participant-level microbiome data and are **not included** in this
> public repository. See [data/README.md](../../data/README.md).

## Reference

nf-core/ampliseq:
- Ewels et al. (2020). The nf-core framework for community-curated bioinformatics pipelines.
  *Nature Biotechnology*, 38, 276–278. https://doi.org/10.1038/s41587-020-0439-x
- SILVA: Quast et al. (2013). The SILVA ribosomal RNA gene database project.
  *Nucleic Acids Research*, 41, D590–D596. https://doi.org/10.1093/nar/gks1219
