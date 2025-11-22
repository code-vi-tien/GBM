# Load libraries
library(GSVA)
library(GSEABase)
library(phateR)
library(xCell)
library(pracma)

# =====================================
# Step 0: Load TPM and metadata
# =====================================
tpm_mat <- read.csv("exports/final/tpm.csv", row.names=1, check.names=FALSE)
meta_df <- read.csv("exports/final/metadata.csv", check.names=FALSE)

# Check column names
if(!all(colnames(tpm_mat) %in% meta_df$rna_barcode)){
  stop("Some TPM columns do not match metadata rna_barcode")
}

# Log2-transform TPM
tpm_log <- log2(tpm_mat + 1)
# =====================================
# DIAGNOSTIC: Check metadata columns
# =====================================
cat("=== METADATA DIAGNOSTICS ===\n")
cat("Total samples in metadata:", nrow(meta_df), "\n")
cat("\nColumn names in metadata:\n")
print(colnames(meta_df))

cat("\n--- Checking 'grade' column ---\n")
if("grade" %in% colnames(meta_df)){
  cat("Unique grades:\n")
  print(table(meta_df$grade, useNA = "ifany"))
} else {
  cat("⚠ 'grade' column not found!\n")
  cat("Available columns:", paste(colnames(meta_df), collapse=", "), "\n")
}

cat("\n--- Checking 'subtype' column ---\n")
if("subtype" %in% colnames(meta_df)){
  cat("Unique molecular subtypes:\n")
  print(table(meta_df$subtype, useNA = "ifany"))
} else {
  cat("⚠ 'subtype' column not found!\n")
  cat("Available columns:", paste(colnames(meta_df), collapse=", "), "\n")
}

# =====================================
# Step 1: Stratify samples (FIXED)
# =====================================
cat("\n=== CREATING STRATA ===\n")

# Check if required columns exist
required_cols <- c("grade", "subtype")
missing_cols <- setdiff(required_cols, colnames(meta_df))

if(length(missing_cols) > 0){
  stop("Missing required columns: ", paste(missing_cols, collapse=", "))
}

# Remove samples with NA in stratification variables
meta_df_clean <- meta_df[!is.na(meta_df$grade) & !is.na(meta_df$subtype), ]

cat("Samples before NA removal:", nrow(meta_df), "\n")
cat("Samples after NA removal:", nrow(meta_df_clean), "\n")

if(nrow(meta_df_clean) == 0){
  stop("No samples remaining after removing NAs from grade/subtype!")
}

# Create stratum column
meta_df_clean$stratum <- paste(meta_df_clean$grade, meta_df_clean$subtype, sep="_")
strata <- unique(meta_df_clean$stratum)

cat("\nStrata created:\n")
print(table(meta_df_clean$stratum))
cat("\nTotal strata:", length(strata), "\n\n")

if(length(strata) == 0){
  stop("No strata were created! Check your metadata.")
}

# =====================================
# Step 2: Load MSigDB GMT files
# =====================================
cat("=== LOADING GENE SETS ===\n")
g_h <- getGmt("supplementary-data/h.all.v2025.1.Hs.symbols.gmt")
g_c2 <- getGmt("supplementary-data/c2.all.v2025.1.Hs.symbols.gmt")

# Convert GeneSetCollection to list of gene vectors
gsc_all <- GeneSetCollection(c(g_h, g_c2))

# Convert to named list of gene vectors
genesets_list <- lapply(gsc_all, geneIds)
names(genesets_list) <- sapply(gsc_all, setName)

# =====================================
# DIAGNOSTIC: Check gene ID overlap
# =====================================
cat("\nChecking gene ID format...\n")
cat("First 5 genes in expression data:\n")
print(head(rownames(tpm_log), 5))

cat("\nFirst 5 genes in first gene set:\n")
print(head(genesets_list[[1]], 5))

# Get all unique genes in gene sets
all_geneset_genes <- unique(unlist(genesets_list))
cat("\nTotal unique genes in gene sets:", length(all_geneset_genes), "\n")
cat("Total genes in expression data:", nrow(tpm_log), "\n")

# Check overlap
overlap <- intersect(rownames(tpm_log), all_geneset_genes)
cat("Genes overlapping:", length(overlap), "\n")

if(length(overlap) < 100){
  cat("\n⚠ WARNING: Very few overlapping genes!\n")
  cat("Gene ID format might be mismatched.\n")
}

# Filter gene sets to only include genes present in expression data
cat("\nFiltering gene sets to match expression data...\n")
genesets_filtered <- lapply(genesets_list, function(gs) {
  intersect(gs, rownames(tpm_log))
})

# Remove empty gene sets (names automatically preserved)
keep_sets <- sapply(genesets_filtered, length) >= 5
genesets_filtered <- genesets_filtered[keep_sets]

cat("Gene sets after filtering:", length(genesets_filtered), "\n")
cat("First 3 gene set names:\n")
print(names(genesets_filtered)[1:3])
cat("Gene sets have names? ", !is.null(names(genesets_filtered)), "\n")

# Check size distribution
set_sizes <- sapply(genesets_filtered, length)
if(length(set_sizes) > 0){
  cat("Gene set size range:", min(set_sizes), "-", max(set_sizes), "\n")
} else {
  cat("⚠ No gene sets remaining after filtering!\n")
}
cat("Gene sets with ≥5 genes:", sum(set_sizes >= 5), "\n\n")

# =====================================
# Step 3: Stratified GSVA + xCell
# =====================================
cat("=== RUNNING STRATIFIED ANALYSIS ===\n")
stratified_results <- list()
cell_fraction_results <- list()

for(s in strata){
  
  cat("\n--- Processing stratum:", s, "---\n")
  
  # Subset TPM using cleaned metadata
  samples_in_stratum <- meta_df_clean$rna_barcode[meta_df_clean$stratum == s]
  
  # Ensure samples exist in TPM matrix
  samples_in_stratum <- intersect(samples_in_stratum, colnames(tpm_log))
  
  if(length(samples_in_stratum) == 0){
    cat("  ⚠ No samples found for this stratum. Skipping.\n")
    next
  }
  
  tpm_sub <- tpm_log[, samples_in_stratum, drop=FALSE]
  
  cat("  Samples:", length(samples_in_stratum), "\n")
  
  # --- GSVA with GSVA 2.x syntax ---
  # Create parameter object with filtered gene sets
  cat("  Checking row names after subsetting...\n")
  cat("    First 5 row names of tpm_sub:", head(rownames(tpm_sub), 5), "\n")
  cat("    Row names exist?", !is.null(rownames(tpm_sub)), "\n")
  cat("    Testing overlap with gene sets:", 
      length(intersect(rownames(tpm_sub), genesets_filtered[[1]])), "\n")

  # Convert to matrix explicitly with row names preserved
  tpm_sub_mat <- as.matrix(tpm_sub)
  cat("    Row names after as.matrix():", head(rownames(tpm_sub_mat), 5), "\n")

  cat("  Converting gene sets to GeneSetCollection...\n")
  gsc <- GeneSetCollection(lapply(names(genesets_filtered), function(nm) {
    GeneSet(genesets_filtered[[nm]], setName = nm)
  }))
  cat("  GeneSetCollection created with", length(gsc), "gene sets\n")

  gsva_par <- gsvaParam(
    exprData = as.matrix(tpm_sub),
    geneSets = genesets_filtered,
    kcdf = "Gaussian",
    minSize = 5,
    maxSize = 5000
  )
  
  # Run GSVA
  cat("  Running GSVA...\n")
  gsva_res <- gsva(gsva_par, verbose = FALSE)
  cat("  GSVA complete. Gene sets scored:", nrow(gsva_res), "\n")
  
  # --- xCell enrichment ---
  cat("  Running xCell...\n")
  xcell_res <- xCellAnalysis(as.matrix(tpm_sub))
  xcell_res <- as.matrix(xcell_res)
  cat("  xCell complete. Cell types scored:", nrow(xcell_res), "\n")
  
  # --- Combine GSVA + xCell ---
  combined <- rbind(gsva_res, xcell_res)
  
  # --- Save per stratum ---
  stratified_results[[s]] <- combined
  cell_fraction_results[[s]] <- xcell_res
  
  cat("  Total features:", nrow(combined), "\n")
}

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Strata processed:", length(stratified_results), "\n")
cat("Stratum names:", paste(names(stratified_results), collapse=", "), "\n")

# =====================================
# Step 4: PHATE embedding per stratum
# =====================================
cat("\n=== CREATING PHATE EMBEDDINGS ===\n")
phate_embeddings <- list()

for(s in names(stratified_results)){
  cat("\nPHATE embedding for stratum:", s, "\n")
  
  # PHATE expects samples x features
  phate_obj <- phate(t(stratified_results[[s]]))
  phate_embeddings[[s]] <- phate_obj$embedding
  
  # Export embeddings
  write.csv(
    phate_obj$embedding,
    paste0("exports/phate_embedding_", s, ".csv"),
    row.names = TRUE
  )
  
  cat("  Saved:", paste0("exports/phate_embedding_", s, ".csv"), "\n")
}

# =====================================
# Save combined results
# =====================================
cat("\n=== SAVING RESULTS ===\n")
saveRDS(stratified_results, "exports/stratified_gsva_xcell.rds")
saveRDS(cell_fraction_results, "exports/stratified_xcell_only.rds")
saveRDS(phate_embeddings, "exports/phate_embeddings.rds")

cat("\n✓ Analysis complete!\n")
cat("  - Stratified GSVA+xCell results saved\n")
cat("  - PHATE embeddings exported per stratum\n")
cat("\nFinal check:\n")
cat("  stratified_results length:", length(stratified_results), "\n")
cat("  cell_fraction_results length:", length(cell_fraction_results), "\n")
cat("  phate_embeddings length:", length(phate_embeddings), "\n")