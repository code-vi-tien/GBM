# =====================================
# BLAST Embedding QC Script + Visualizations
# =====================================
library(cluster)      # silhouette
library(RANN)         # kNN
library(phateR)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

# -------------------------------
# Load data
# -------------------------------
phate_embeddings <- readRDS("exports/embeddings/phate_embeddings.rds")
stratified_results <- readRDS("exports/embeddings/gsva_xcell/stratified_gsva_xcell.rds")
meta_df <- read.csv("exports/final-supplementary-data/metadata.csv", check.names=FALSE)

# -------------------------------
# Create stratum
# -------------------------------
meta_df_clean <- meta_df[!is.na(meta_df$grade) & !is.na(meta_df$molecular_subtype), ]
meta_df_clean$stratum <- paste(meta_df_clean$grade, meta_df_clean$molecular_subtype, sep="_")

# -------------------------------
# Parameters
# -------------------------------
n_subsample <- 0.8
k_nn <- 5
set.seed(123)

# Create folder for QC plots
dir.create("exports/embeddings/qc", recursive = TRUE, showWarnings = FALSE)

# Initialize storage
qc_metrics <- list()

# ============================================================
# Helper Functions
# ============================================================
# ---- 1. Silhouette score -----------------------------------
compute_silhouette <- function(emb, labels) {
  if(length(unique(labels)) <= 1) return(NA)

  dmat <- dist(emb)
  sil <- silhouette(as.numeric(as.factor(labels)), dmat)
  mean(sil[, 3])
}

# ---- 2. Subsampling stability -------------------------------
compute_subsample_corr <- function(emb, frac = 0.8) {
  n <- nrow(emb)
  if(n < 10) return(NA)

  n_sub <- ceiling(n * frac)
  idx <- sample(1:n, n_sub)
  emb_sub <- emb[idx, , drop = FALSE]

  ph <- phate(emb_sub)

  dims <- min(ncol(emb_sub), ncol(ph$embedding))
  cors <- sapply(1:dims, function(i)
    cor(emb_sub[, i], ph$embedding[, i], method = "spearman")
  )

  mean(cors)
}

# ---- 3. Local continuity (kNN overlap) -----------------------
compute_knn_overlap <- function(emb, k = 10) {
  if(nrow(emb) <= k + 1) return(NA)

  nn_orig <- nn2(emb, k = k)$nn.idx

  idx <- sample(1:nrow(emb), ceiling(nrow(emb) * 0.8))
  emb_sub <- emb[idx, , drop = FALSE]
  nn_sub <- nn2(emb_sub, k = k)$nn.idx

  min_rows <- min(nrow(nn_orig), nrow(nn_sub))

  overlaps <- sapply(1:min_rows, function(i) {
    length(intersect(nn_orig[i, ], nn_sub[i, ])) / k
  })

  mean(overlaps)
}

# ---- 4. Embedding–feature correlation ------------------------
compute_embedding_feature_corr <- function(emb, feats) {
  feats <- as.matrix(feats)
  if(nrow(feats) == 1) {
    return(mean(cor(emb, feats[1,], method = "spearman")))
  }
  cor_vals <- sapply(1:ncol(emb), function(i) {
    mean(apply(feats, 1, function(f) cor(emb[,i], f, method="spearman")))
  })
  mean(cor_vals)
}

# ---- 5. Grade correlation ------------------------------------
compute_grade_corr <- function(emb, grade_vec) {
  if(length(unique(grade_vec)) <= 1) return(NA)
  grade_vec <- as.numeric(factor(grade_vec))
  mean(apply(emb, 2, function(x) cor(x, grade_vec, method="spearman")))
}

# ---- 6. PHATE disconnected-component check ------------------
detect_phate_disconnect <- function(emb) {
  # compute kNN graph connectivity
  kn <- nn2(emb, k = 10)
  adj <- matrix(0, nrow=nrow(emb), ncol=nrow(emb))
  for(i in 1:nrow(emb)) adj[i, kn$nn.idx[i, ]] <- 1
  g <- igraph::graph_from_adjacency_matrix(adj, mode="undirected")
  comps <- igraph::components(g)$no
  comps
}

# ---- 7. Safe plot assembly -----------------------------------
save_qc_plots <- function(p_list, outfile) {
  p_list <- p_list[!sapply(p_list, is.null)]
  if(length(p_list) == 0) return(NULL)

  pdf(outfile, width=10, height=8)
  do.call(gridExtra::grid.arrange, c(p_list, ncol=1))
  dev.off()
}

# ============================================================
#           Main QC Loop
# ============================================================

qc_metrics <- list()

adj_full <- as.matrix(read.csv("exports/final-supplementary-data/adjacency_matrix.csv", row.names = 1))
adj_full <- 1 * ((adj_full + t(adj_full)) > 0)   # ensure global symmetry

for(s in names(phate_embeddings)) {
  cat("\n=== QC for stratum:", s, "===\n")

  emb <- phate_embeddings[[s]]
  if(nrow(emb) < 5) {
    cat(" ⚠ Too few samples, skipping\n")
    next
  }

  samples <- rownames(emb)

  meta_sub <- meta_df_clean[meta_df_clean$stratum == s, ]
  rownames(meta_sub) <- meta_sub$rna_barcode
  meta_sub <- meta_sub[samples, , drop=FALSE]

  feats <- stratified_results[[s]]

  # -------------------------
  # Compute metrics
  # -------------------------
  sil <- compute_silhouette(emb, meta_sub$molecular_subtype)
  subs <- compute_subsample_corr(emb)
  knn  <- compute_knn_overlap(emb, k=10)
  efc  <- compute_embedding_feature_corr(emb, feats)
  gcor <- compute_grade_corr(emb, meta_sub$grade)
  comps <- detect_phate_disconnect(emb)

  cat(" Silhouette:", sil, "\n")
  cat(" Subsample correlation:", subs, "\n")
  cat(" kNN continuity:", knn, "\n")
  cat(" Embedding-feature corr:", efc, "\n")
  cat(" Grade corr:", gcor, "\n")
  cat(" PHATE graph components:", comps, "\n")

  # -------------------------
  # Construct plots
  # -------------------------
  p_sil <- if(!is.na(sil)) {
    ggplot() + 
      geom_bar(aes(x="Silhouette", y=sil), stat="identity") +
      ylim(0,1) +
      ggtitle(paste(s, "- Silhouette"))
  } else NULL

  p_sub <- ggplot() +
    geom_bar(aes(x="Subsample", y=subs), stat="identity") +
    ggtitle(paste(s, "- Subsample Stability"))

  p_knn <- ggplot() +
    geom_bar(aes(x="kNN overlap", y=knn), stat="identity") +
    ggtitle(paste(s, "- Local Continuity"))

  plot_file <- paste0("exports/embeddings/qc/qc_", s, ".pdf")
  save_qc_plots(list(p_sil, p_sub, p_knn), plot_file)

  # -------------------------
  # Store metrics
  # -------------------------
  qc_metrics[[s]] <- list(
    silhouette = sil,
    subsampling = subs,
    knn_continuity = knn,
    emb_feat_cor = efc,
    grade_cor = gcor,
    phate_components = comps
  )
}

# Save metrics
qc_df <- do.call(rbind, lapply(names(qc_metrics), function(nm) {
  c(stratum = nm, unlist(qc_metrics[[nm]]))
}))
write.csv(qc_df, "exports/embeddings/qc/embeddings_qc_summary.csv", row.names=FALSE)
