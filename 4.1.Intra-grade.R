##############################
# Intra-Stratum (Intra-Grade) Pathway & Gene Analysis
##############################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
})

utils::globalVariables(c("rho", "p", "FDR", "feature", "abs_rho", "type"))

# -----------------------------
# 0. Create directories
# -----------------------------
dir.create("exports/intra-grade", showWarnings = FALSE)
dir.create("exports/intra-grade/plots", showWarnings = FALSE)
plot_dir <- "exports/intra-grade/plots"

# -----------------------------
# 1. Load data
# -----------------------------
stratified_results <- readRDS("exports/embeddings/gsva_xcell/stratified_gsva_xcell.rds")
progression_list <- readRDS("exports/progression/progression_scores.rds")
cell_fraction_results <- readRDS("exports/embeddings/xcell/stratified_xcell_only.rds")
tpm <- read.csv("exports/final-supplementary-data/tpm.csv", row.names = 1, check.names = FALSE)
tpm_log <- log2(as.matrix(tpm) + 1)

# -----------------------------
# 2. Feature–Progression Association Function
# -----------------------------
run_feature_assoc <- function(feature_mat, prog_df) {
  common_samples <- intersect(prog_df$rna_barcode, colnames(feature_mat))
  prog <- prog_df$progression_score[match(common_samples, prog_df$rna_barcode)]
  feature_sub <- feature_mat[, common_samples, drop = FALSE]
  
  res <- t(apply(feature_sub, 1, function(f) {
    ct <- suppressWarnings(cor.test(f, prog, method = "spearman"))
    c(rho = as.numeric(ct$estimate), p = as.numeric(ct$p.value))
  }))
  
  as.data.frame(res) %>%
    tibble::rownames_to_column("feature") %>%
    mutate(
      rho = as.numeric(rho),
      p = as.numeric(p),
      FDR = p.adjust(p, method = "fdr"),
      abs_rho = abs(rho)
    ) %>%
    arrange(FDR)
}

# -----------------------------
# 3. Visualization Functions
# -----------------------------
plot_volcano <- function(df, title, top_n = 15) {
  top_features <- df %>% slice_min(FDR, n = top_n)
  
  ggplot(df, aes(x = rho, y = -log10(FDR))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_point(data = top_features, color = "red", size = 2) +
    geom_text_repel(
      data = top_features, 
      aes(label = feature), 
      size = 2.5, 
      max.overlaps = 30
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal() +
    labs(
      title = title,
      x = "Spearman ρ (Progression Association)",
      y = "-log10(FDR)"
    )
}

plot_top_features <- function(df, title, top_n = 20) {
  top <- df %>% 
    slice_min(FDR, n = top_n) %>%
    mutate(direction = ifelse(rho > 0, "Increases", "Decreases"))
  
  ggplot(top, aes(x = reorder(feature, rho), y = rho, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Increases" = "firebrick", "Decreases" = "steelblue")) +
    theme_minimal() +
    labs(
      title = title,
      x = "Feature (Pathway / Cell Type)",
      y = "Spearman ρ",
      fill = "With Progression"
    )
}

# -----------------------------
# 4. Run intra-stratum analysis
# -----------------------------
for (stratum in names(progression_list)) {
  
  message("Processing stratum: ", stratum)
  
  prog_df <- progression_list[[stratum]]
  feature_mat <- stratified_results[[stratum]]
  
  # Pathway/Cell-type association
  result <- run_feature_assoc(feature_mat, prog_df)
  
  # Annotate type (GSVA vs xCell)
  xcell_features <- rownames(cell_fraction_results[[stratum]])
  result$type <- ifelse(result$feature %in% xcell_features, "xCell", "GSVA")
  
  # Save results
  write.csv(
    result, 
    paste0("exports/intra-grade/pathway_progression_", stratum, ".csv"),
    row.names = FALSE
  )
  
  # Volcano plot
  pdf(file.path(plot_dir, paste0("pathway_progression_volcano_", stratum, ".pdf")), width = 13, height = 9)
  print(plot_volcano(result, paste0("Pathway Association – ", stratum)))
  dev.off()
  
  # Top features bar plot
  pdf(file.path(plot_dir, paste0("pathway_progression_top_", stratum, ".pdf")), width = 13, height = 12)
  print(plot_top_features(result, paste0("Top 20 Features – ", stratum)))
  dev.off()
  
  # Save top 50 summary
  top_summary <- result %>%
    slice_min(FDR, n = 50) %>%
    select(feature, type, rho, p, FDR)
  
  write.csv(
    top_summary,
    paste0("exports/intra-grade/top50_", stratum, ".csv"),
    row.names = FALSE
  )
  
  # -----------------------------
  # Gene-level association (secondary validation)
  # -----------------------------
  result_genes <- run_feature_assoc(tpm_log, prog_df)
  
  write.csv(
    result_genes,
    paste0("exports/intra-grade/genes_", stratum, ".csv"),
    row.names = FALSE
  )
  
  pdf(file.path(plot_dir, paste0("genes_progression_volcano_", stratum, ".pdf")), width = 13, height = 9)
  print(plot_volcano(result_genes, paste0("Gene-Level Validation – ", stratum), top_n = 20))
  dev.off()
}

message("✓ Intra-grade analysis complete!")
message("Results saved in: exports/intra-grade/")
message("Plots saved in: ", plot_dir)
