# Load libraries
library(GSVA)
library(GSEABase)
library(phateR)
library(xCell)

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
# Step 1: Stratify samples
# =====================================
meta_df$stratum <- paste(meta_df$grade, meta_df$molecular_subtype, sep="_")
strata <- unique(meta_df$stratum)

# =====================================
# Step 2: Load MSigDB GMT files
# =====================================
g_h <- getGmt("supplementary-data/h.all.v2025.1.Hs.symbols.gmt")
g_c2 <- getGmt("supplementary-data/c2.all.v2025.1.Hs.symbols.gmt")

# Convert GeneSetCollection to list of gene vectors
genesets_list <- lapply(c(g_h, g_c2), geneIds)

# =====================================
# DIAGNOSTIC: Check gene ID overlap
# =====================================
cat("Checking gene ID format...\n")
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
  cat("Common issues:\n")
  cat("- Expression data uses Ensembl IDs (ENSG...) but gene sets use symbols\n")
  cat("- Version differences in gene symbols\n")
  cat("- Case sensitivity issues\n\n")
  
  # Attempt to detect Ensembl IDs
  if(any(grepl("^ENSG", rownames(tpm_log)))){
    stop("Expression data appears to use Ensembl IDs, but gene sets use gene symbols. Please convert your TPM matrix row names to gene symbols first.")
  }
}

# Filter gene sets to only include genes present in expression data
cat("\nFiltering gene sets to match expression data...\n")
genesets_filtered <- filterGeneSets(genesets_list, rownames(tpm_log))
cat("Gene sets after filtering:", length(genesets_filtered), "\n")

# Check size distribution
set_sizes <- sapply(genesets_filtered, length)
cat("Gene set size range:", min(set_sizes), "-", max(set_sizes), "\n")
cat("Gene sets with ≥5 genes:", sum(set_sizes >= 5), "\n\n")

# =====================================
# Step 3: Stratified GSVA + xCell
# =====================================
stratified_results <- list()
cell_fraction_results <- list()

for(s in strata){
  
  cat("Processing stratum:", s, "\n")
  
  # Subset TPM
  samples_in_stratum <- meta_df$rna_barcode[meta_df$stratum == s]
  tpm_sub <- tpm_log[, samples_in_stratum, drop=FALSE]
  
  cat("  Samples:", length(samples_in_stratum), "\n")
  
  # --- GSVA with GSVA 2.x syntax ---
  # Create parameter object with filtered gene sets
  gsva_par <- gsvaParam(
    exprData = as.matrix(tpm_sub),
    geneSets = genesets_filtered,
    kcdf = "Gaussian",
    minSize = 5,  # Increased from 1 to avoid tiny gene sets
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
  
  cat("  Total features:", nrow(combined), "\n\n")
}

# =====================================
# Step 4: PHATE embedding per stratum
# =====================================
phate_embeddings <- list()

for(s in names(stratified_results)){
  cat("PHATE embedding for stratum:", s, "\n")
  
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
saveRDS(stratified_results, "exports/stratified_gsva_xcell.rds")
saveRDS(cell_fraction_results, "exports/stratified_xcell_only.rds")
saveRDS(phate_embeddings, "exports/phate_embeddings.rds")

cat("\n✓ Analysis complete!\n")
cat("  - Stratified GSVA+xCell results saved\n")
cat("  - PHATE embeddings exported per stratum\n")