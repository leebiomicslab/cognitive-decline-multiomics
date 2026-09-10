# Purpose:           Microbial differential abundance analysis using ANCOM-BC2 + ALDEx2
#                    across four covariate model configurations (Models A-D).
#                    Candidate taxa = union of ANCOM-BC2 and ALDEx2 p < 0.05 hits.
# Manuscript:        Methods — Microbiome differential abundance analysis
# Figure/Table:      Fig. 5a-c; Supplementary Fig. 6; Supplementary Table 4
# Input:             clinical_merged_371.csv, asv_genus_table_reads.csv
# Output:            results/microbiota_covariate_analysis/<Model>/*
# Main dependencies: phyloseq, ANCOMBC (>= 2.0), ALDEx2, emmeans, VennDiagram, tidyverse
#
# CONFIGURATION:
#   Set BASE_DIR below to the data root on your system.
#   All input CSVs must be placed as specified under data/ and
#   results/ are written relative to REPO_DIR (auto-detected).
#
# NOTE — DISCREPANCY (D5): ALDEx2 run with mc.samples=128 (see line 73).
#   Manuscript section does not specify; confirm for reproducibility.

options(warn=-1)
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
  library(ALDEx2)
  library(VennDiagram)
  library(grid)
  library(emmeans)
})

# ── Path configuration ────────────────────────────────────────────────────────
# Set BASE_DIR to the directory containing your data files.
# Alternatively, set the env var TVGH_DATA_DIR before running:
#   Rscript run_microbiome_covariate_models.R
if (Sys.getenv("TVGH_DATA_DIR") != "") {
  base_dir <- Sys.getenv("TVGH_DATA_DIR")
} else {
  # Default: repo results + data/ are siblings of this script's grandparent
  base_dir <- file.path(dirname(dirname(dirname(normalizePath(sys.frame(1)$ofile,
               mustWork = FALSE)))), "data")
  # Fallback: use working directory
  if (!dir.exists(base_dir)) base_dir <- file.path(getwd(), "data")
}

repo_dir    <- dirname(dirname(dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))))
if (!dir.exists(repo_dir)) repo_dir <- getwd()

data_dir        <- base_dir
micro_dir       <- file.path(data_dir, "MICROBIOTA")
out_dir_base    <- file.path(repo_dir, "results", "microbiota_covariate_analysis")
if (!dir.exists(out_dir_base)) dir.create(out_dir_base, recursive = TRUE)

clin    <- read.csv(file.path(data_dir, "clinical_merged_371.csv"),
                    colClasses = c("id" = "character"))
reads_raw <- read.csv(file.path(micro_dir, "asv_genus_table_reads.csv"),
                      check.names = FALSE)
rownames(reads_raw) <- reads_raw$Genus
reads_mat <- reads_raw %>% dplyr::select(-Genus) %>% as.matrix()

for(cc in c("gender","smoking","alcohol","HTN","DM","CVD","Statin","PPI","Antibiotics","dm_med","htn_med","chol_med")) {
  clin[[cc]] <- as.factor(clin[[cc]])
}
clin$group <- factor(clin$group, levels = c("C1","C2","C3"))

run_model <- function(formula_rhs, model_name, base_out) {
  cat(sprintf("\n====  %s  ====\n", model_name))
  out_dir <- file.path(base_out, model_name)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

  vars <- all.vars(as.formula(paste("~", formula_rhs)))
  meta_df <- clin %>% drop_na(all_of(vars)) %>% tibble::column_to_rownames("id")
  valid <- intersect(colnames(reads_mat), rownames(meta_df))
  m_reads <- reads_mat[, valid]
  m_meta <- meta_df[valid, ]
  cat(sprintf("  Samples=%d\n", length(valid)))

  tse <- phyloseq(otu_table(m_reads, taxa_are_rows=TRUE), sample_data(m_meta))
  keep <- taxa_sums(tse) > 0 & rowSums(otu_table(tse) > 0) >= 0.1*nsamples(tse)
  tse_f <- prune_taxa(keep, tse)
  cat(sprintf("  Taxa=%d\n", ntaxa(tse_f)))

  # ── ANCOM-BC2 ────────────────────────────────────────────────────
  cat("  ANCOM-BC2...\n")
  ancom_res <- NULL
  tryCatch({
    out <- ancombc2(data=tse_f, fix_formula=formula_rhs, p_adj_method="BH",
                    pseudo_sens=TRUE, prv_cut=0, lib_cut=0,
                    alpha=0.05, n_cl=1, verbose=FALSE)
    ancom_res <- out$res
    cat(sprintf("  ANCOM-BC2 OK: %d taxa\n", nrow(ancom_res)))
  }, error=function(e) cat(sprintf("  ANCOM-BC2 FAIL: %s\n", e$message)))

  if(is.null(ancom_res)) {
    cat("  Skipping — no ANCOM result\n")
    return(invisible(NULL))
  }

  p_c2 <- grep("p_.*groupC2|p_val.*groupC2", colnames(ancom_res), value=TRUE)[1]
  p_c3 <- grep("p_.*groupC3|p_val.*groupC3", colnames(ancom_res), value=TRUE)[1]
  lfc_c2 <- grep("lfc.*groupC2|beta.*groupC2", colnames(ancom_res), value=TRUE)[1]
  lfc_c3 <- grep("lfc.*groupC3|beta.*groupC3", colnames(ancom_res), value=TRUE)[1]
  se_c2 <- grep("^se_groupC2", colnames(ancom_res), value=TRUE)[1]
  se_c3 <- grep("^se_groupC3", colnames(ancom_res), value=TRUE)[1]

  # ── ALDEx2 ────────────────────────────────────────────────────────
  cat("  ALDEx2...\n")
  aldex_hits <- c()
  tryCatch({
    dmat_local <- model.matrix(as.formula(paste("~", formula_rhs)), data=m_meta)
    xc_obj <- aldex.clr(as.matrix(otu_table(tse_f)), dmat_local, mc.samples=128, denom="all", verbose=FALSE)
    xg_obj <- aldex.glm(xc_obj, dmat_local)
    ap2 <- grep("groupC2.*pval", colnames(xg_obj), value=TRUE)[1]
    ap3 <- grep("groupC3.*pval", colnames(xg_obj), value=TRUE)[1]
    if(!is.na(ap2) && !is.na(ap3)) {
      aldex_hits <- unique(c(rownames(xg_obj)[!is.na(xg_obj[,ap2]) & xg_obj[,ap2]<0.05],
                             rownames(xg_obj)[!is.na(xg_obj[,ap3]) & xg_obj[,ap3]<0.05]))
    }
    cat(sprintf("  ALDEx2 OK: %d hits\n", length(aldex_hits)))
    
    aldex_res <- as.data.frame(xg_obj)
    aldex_res$taxon <- rownames(aldex_res)
    write.csv(aldex_res, file.path(out_dir, "ALDEx2_full_results.csv"), row.names = FALSE)
    if(!is.na(ap2) && !is.na(ap3)) {
      aldex_res$padj_C2 <- p.adjust(aldex_res[[ap2]], method = "BH")
      aldex_res$padj_C3 <- p.adjust(aldex_res[[ap3]], method = "BH")
      write.csv(aldex_res, file.path(out_dir, "ALDEx2_full_results_with_BH.csv"), row.names = FALSE)
    }
  }, error=function(e) cat(sprintf("  ALDEx2 FAIL: %s\n", e$message)))

  ancom_hits <- ancom_res$taxon[!is.na(ancom_res[[p_c2]]) & (ancom_res[[p_c2]]<0.05 | ancom_res[[p_c3]]<0.05)]
  union_set <- unique(c(ancom_hits, aldex_hits))
  cat(sprintf("  Union=%d (ANCOM=%d, ALDEx2=%d)\n", length(union_set), length(ancom_hits), length(aldex_hits)))

  cat("  CSV...\n")
  res_out <- ancom_res
  res_out$ALDEx2_Sig <- res_out$taxon %in% aldex_hits
  res_out$Union_Sig <- res_out$taxon %in% union_set
  write.csv(res_out, file.path(out_dir, "microbiota_da_results_final.csv"), row.names=FALSE)
  cat("  CSV saved\n")

  cat("  Venn...\n")
  tryCatch({
    png(file.path(out_dir, "Figure1_VennDiagram.png"), width=8, height=6, units="in", res=600)
    vd <- venn.diagram(x=list("ANCOM-BC2\n(P<0.05)"=unique(ancom_hits), "ALDEx2\n(P<0.05)"=unique(aldex_hits)),
                       filename=NULL, fill=c("#70a1d7","#f4a261"), alpha=0.5,
                       cex=1.5, cat.cex=1.3, cat.fontface="bold",
                       cat.pos=c(-20,20), disable.logging=TRUE)
    grid.draw(vd)
    dev.off()
    cat("  Venn OK\n")
  }, error=function(e) { try(dev.off(), silent=TRUE); cat(sprintf("  Venn FAIL: %s\n", e$message)) })

  if(length(union_set) == 0) {
    cat("  No union genera — skip Fig2 & Fig3\n")
    return(invisible(NULL))
  }

  cat("  Forest...\n")
  tryCatch({
    bs <- res_out[res_out$taxon %in% union_set, ]
    fd <- rbind(data.frame(taxon=bs$taxon, Comparison="C1 vs C2", LFC=bs[[lfc_c2]], AncomSE=bs[[se_c2]], P=bs[[p_c2]], stringsAsFactors=FALSE),
                data.frame(taxon=bs$taxon, Comparison="C1 vs C3", LFC=bs[[lfc_c3]], AncomSE=bs[[se_c3]], P=bs[[p_c3]], stringsAsFactors=FALSE))
    fd$CI_lo <- fd$LFC - 1.96 * fd$AncomSE
    fd$CI_hi <- fd$LFC + 1.96 * fd$AncomSE
    fd$CI_Cross <- factor(ifelse(fd$CI_lo <= 0 & fd$CI_hi >= 0, "95% CI crosses 0", "95% CI does not cross 0"),
                          levels=c("95% CI does not cross 0", "95% CI crosses 0"))
    fd <- fd[order(fd$Comparison, -abs(fd$LFC)), ]

    dodge <- position_dodge(width = 0.6)

    pf <- ggplot(fd, aes(x=LFC, y=reorder(taxon, LFC), color=Comparison, group=Comparison)) +
      geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=1) +
      geom_errorbar(aes(xmin=CI_lo, xmax=CI_hi), width=0.3, alpha=0.9, linewidth=1.2, position=dodge) +
      geom_point(aes(fill=interaction(Comparison, CI_Cross)), shape=21, size=4.5, stroke=1.2, position=dodge) +
      scale_color_manual(values=c("C1 vs C2"="#eca362", "C1 vs C3"="#6fcbbc")) +
      scale_fill_manual(values=c("C1 vs C2.95% CI does not cross 0"="#eca362", 
                                 "C1 vs C3.95% CI does not cross 0"="#6fcbbc",
                                 "C1 vs C2.95% CI crosses 0"="white",
                                 "C1 vs C3.95% CI crosses 0"="white"), guide="none") +
      theme_classic(base_size=18) +
      theme(
        axis.text.y = element_text(size=20, face="bold", color="black"),
        axis.text.x = element_text(size=18, color="black"),
        axis.title.x = element_text(size=20, face="bold", color="black"),
        axis.title.y = element_blank(),
        title = element_text(size=22, face="bold"),
        panel.grid.major.x = element_line(color="grey80", linetype="dotted", linewidth=1),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=16)
      ) +
      labs(title="Effect Sizes of Significant Taxa",
           x="Adjusted log2 Fold Change (Beta)") +
      geom_point(aes(x=0, y=0, shape=CI_Cross), alpha=0) +
      scale_shape_manual(values=c("95% CI does not cross 0"=19, "95% CI crosses 0"=21), drop=FALSE, na.translate=FALSE) +
      guides(
        color = guide_legend(order=1, override.aes=list(shape=15, size=6)),
        shape = guide_legend(order=2, override.aes=list(alpha=1, color="gray40", fill=c("gray40", "white"), size=5))
      )
      
    ggsave(file.path(out_dir, "Figure2_ForestPlot.png"), pf, width=14, height=max(8, length(union_set)*0.6), dpi=600)
    cat("  Forest OK\n")
  }, error=function(e) cat(sprintf("  Forest FAIL: %s\n", e$message)))

  cat("  Boxplot...\n")
  tryCatch({
    otu_m <- as.matrix(otu_table(tse_f))
    sample_mat <- t(otu_m)                # samples x taxa
    
    clr_sample <- t(apply(sample_mat + 1, 1, function(x) {
      log(x) - mean(log(x))
    }))
    
    clr_df <- as.data.frame(clr_sample)
    avail <- intersect(union_set, colnames(clr_df))

    long_df <- clr_df[, avail, drop=FALSE]
    long_df$id <- rownames(long_df)
    long_df <- long_df %>%
      left_join(m_meta %>% rownames_to_column("id"), by="id") %>%
      pivot_longer(cols=all_of(avail), names_to="Genus", values_to="CLR")

    adj_list <- lapply(avail, function(g) {
      tryCatch({
        ds <- long_df[long_df$Genus==g, ]
        fit <- lm(as.formula(paste("CLR ~", formula_rhs)), data=ds)
        em <- as.data.frame(emmeans(fit, ~group))
        names(em)[names(em)=="emmean"] <- "AdjMean"
        names(em)[names(em)=="SE"] <- "AdjSE"
        em$Genus <- g
        em
      }, error=function(e) NULL)
    })
    adj_means <- do.call(rbind, Filter(Negate(is.null), adj_list))

    if(!is.null(adj_means) && nrow(adj_means)>0) {
      gc <- c("C1"="#ee7a5b","C2"="#eca362","C3"="#6fcbbc")
      pb <- ggplot() +
        geom_jitter(data=long_df, aes(x=group, y=CLR, color=group), width=0.2, alpha=0.3, size=1) +
        geom_violin(data=long_df, aes(x=group, y=CLR, fill=group), alpha=0.1, color=NA) +
        geom_errorbar(data=adj_means, aes(x=group, ymin=AdjMean-AdjSE, ymax=AdjMean+AdjSE), width=0.2, linewidth=0.8, color="grey30") +
        geom_point(data=adj_means, aes(x=group, y=AdjMean, color=group), size=4, shape=18) +
        geom_line(data=adj_means, aes(x=group, y=AdjMean, group=Genus), color="black", alpha=0.6, linetype="dashed", linewidth=0.5) +
        facet_wrap(~Genus, scales="free_y", ncol=5) +
        scale_color_manual(values=gc) + scale_fill_manual(values=gc) +
        theme_bw(base_size=11) +
        theme(panel.grid.minor=element_blank(), legend.position="bottom", strip.text=element_text(size=8)) +
        labs(title="Differential Genera — Union Set (Covariate-Adjusted Marginal Means)",
             subtitle=sprintf("%s | N=%d genera", model_name, length(avail)),
             y="CLR Abundance", x="Group")
      rows <- ceiling(length(avail)/5)
      ggsave(file.path(out_dir, "Figure3_AdjustedBoxplot.png"), pb, width=18, height=max(8, 3*rows), dpi=600)
      cat("  Boxplot OK\n")
    } else {
      cat("  Boxplot: no adj_means\n")
    }

  }, error=function(e) cat(sprintf("  Boxplot FAIL: %s\n", e$message)))

  cat(sprintf("  [%s DONE]\n", model_name))
}

models <- list(
  "Model_A_Main" = "group + age + gender + edu + bmi + smoking + alcohol + HTN + DM + CVD",
  "Model_B_Sens1_Systemic" = "group + age + gender + edu + bmi + smoking + alcohol + HTN + DM + CVD + CRP + egfr",
  "Model_C_Sens2_MicrobiomeDrugs" = "group + age + gender + edu + bmi + smoking + alcohol + HTN + DM + CVD + Statin + PPI + Antibiotics",
  "Model_D_Sens3_TxSummary" = "group + age + gender + edu + bmi + smoking + alcohol + HTN + DM + CVD + dm_med + htn_med + chol_med"
)
for(nm in names(models)) {
  run_model(models[[nm]], nm, out_dir_base)
}
cat("\n[ALL MODELS DONE]\n")