##############################
# Pathway/Cell-Type–Progression Association (Corrected)
##############################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

utils::globalVariables(c("rho", "p", "FDR", "feature", "abs_rho", "type"))

dir.create("exports/gene-level", showWarnings = FALSE)
dir.create("exports/gene-level/plots", showWarnings = FALSE)
plot_dir <- "exports/gene-level/plots"

# -----------------------------
# 1. Load GSVA+xCell scores and progression
# -----------------------------
stratified_results <- readRDS("exports/embeddings/gsva_xcell/stratified_gsva_xcell.rds")
progression_list <- readRDS("exports/progression/progression_scores.rds")

# Optional: separate GSVA from xCell for interpretation
cell_fraction_results <- readRDS("exports/embeddings/xcell/stratified_xcell_only.rds")

# -----------------------------
# 2. Function: Feature–Progression Association
# -----------------------------
run_feature_assoc <- function(feature_mat, prog_df) {
  # Align samples
  common_samples <- intersect(prog_df$rna_barcode, colnames(feature_mat))
  prog <- prog_df$progression_score[match(common_samples, prog_df$rna_barcode)]
  feature_sub <- feature_mat[, common_samples, drop = FALSE]
  
  # Spearman correlation per feature (pathway or cell type)
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
# 3. Visualization: Volcano + Ranked Bar
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
# 4. Run per stratum
# -----------------------------
for (stratum in names(progression_list)) {
  
  message("Processing stratum: ", stratum)
  
  prog_df <- progression_list[[stratum]]
  feature_mat <- stratified_results[[stratum]]
  
  # Run association
  result <- run_feature_assoc(feature_mat, prog_df)
  
  # Annotate feature type (GSVA pathway vs xCell cell type)
  xcell_features <- rownames(cell_fraction_results[[stratum]])
  result$type <- ifelse(result$feature %in% xcell_features, "xCell", "GSVA")
  
  # Save results
  write.csv(
    result, 
    paste0("exports/gene-level/pathway_progression_", stratum, ".csv"),
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
  
  # Save top features summary
  top_summary <- result %>%
    slice_min(FDR, n = 50) %>%
    select(feature, type, rho, p, FDR)
  
  write.csv(
    top_summary,
    paste0("exports/gene-level/top50_", stratum, ".csv"),
    row.names = FALSE
  )
}

message("✓ Pathway-progression associations complete.")
message("  Results: exports/gene-level/")
message("  Plots: plots/")


# -----------------------------
# 5. OPTIONAL: Gene-level analysis (as secondary validation)
# -----------------------------
message("\n--- Running secondary gene-level analysis ---")

tpm <- read.csv("exports/final-supplementary-data/tpm.csv", row.names = 1, check.names = FALSE)
tpm_log <- log2(as.matrix(tpm) + 1)

for (stratum in names(progression_list)) {
  
  prog_df <- progression_list[[stratum]]
  
  # Only use samples from this stratum
  result_genes <- run_feature_assoc(tpm_log, prog_df)
  
  write.csv(
    result_genes,
    paste0("exports/gene-level/genes_", stratum, ".csv"),
    row.names = FALSE
  )
  
  # Top genes volcano
  pdf(file.path(plot_dir, paste0("genes_progression_volcano_", stratum, ".pdf")), width = 13, height = 9)
  print(plot_volcano(
    result_genes, 
    paste0("Gene-Level Validation – ", stratum),
    top_n = 20
  ))
  dev.off()
}

message("✓ Gene-level validation complete (secondary analysis).")


# -----------------------------
# 6. Cross-stratum comparison (identify shared vs unique drivers)
# -----------------------------
message("\n--- Comparing progression drivers across strata ---")

all_results <- lapply(names(progression_list), function(s) {
  read.csv(paste0("exports/gene-level/pathway_progression_", s, ".csv")) %>%
    mutate(stratum = s) %>%
    filter(FDR < 0.05, abs(rho) > 0.3)  # Significant + strong associations
})

combined <- bind_rows(all_results)

# Features appearing in multiple strata
shared_features <- combined %>%
  group_by(feature) %>%
  summarise(
    n_strata = n_distinct(stratum),
    mean_rho = mean(rho),
    min_FDR = min(FDR),
    strata_list = paste(unique(stratum), collapse = "; ")
  ) %>%
  arrange(desc(n_strata), min_FDR)

write.csv(
  shared_features,
  "exports/gene-level/shared_features_across_strata.csv",
  row.names = FALSE
)

# Plot heatmap of top shared features
top_shared <- head(shared_features$feature, 20)

heatmap_data <- combined %>%
  filter(feature %in% top_shared) %>%
  select(feature, stratum, rho) %>%
  tidyr::pivot_wider(names_from = stratum, values_from = rho, values_fill = 0)

write.csv(
  heatmap_data,
  "exports/gene-level/shared_features_heatmap.csv",
  row.names = FALSE
)

message("✓ Cross-stratum comparison complete.")
message("\n=== ANALYSIS COMPLETE ===")
message("Key outputs:")
message("  1. Pathway associations per stratum")
message("  2. Gene-level validation (secondary)")
message("  3. Shared features across strata")
message("  4. Volcano plots and bar charts in plots/")