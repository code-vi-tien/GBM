# Load libraries
library(igraph)

# =====================================
# Step 0: Load TPM and metadata
# =====================================
phate_embeddings <- readRDS("exports/embeddings/phate_embeddings.rds")
meta_df <- read.csv("exports/final-supplementary-data/metadata.csv", check.names = FALSE)
adj_mat <- read.csv("exports/final-supplementary-data/adjacency_matrix.csv", 
                    row.names = 1, check.names = FALSE)
adj_mat <- as.matrix(adj_mat)


# =====================================
# Step 1: Build patient groups using adjacency matrix
# =====================================
message("Constructing patient trajectories from adjacency matrix...")

g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
comp <- components(g)

patient_ids <- paste0("patient_", comp$membership)
names(patient_ids) <- names(comp$membership)


# =====================================
# Step 2: Process PHATE embeddings per stratum
# =====================================
progression_raw <- list()

for(s in names(phate_embeddings)) {
  
  message("Processing stratum: ", s)
  
  emb <- phate_embeddings[[s]]
  samps <- rownames(emb)  # rna_barcodes
  
  # match rna_barcode to metadata
  meta_s <- meta_df[match(samps, meta_df$rna_barcode), ]
  
  # append patient_id
  meta_s$patient_id <- patient_ids[meta_s$sample_barcode]

  # ensure surgery_number comes from metadata
  surgery_number <- meta_s$surgery_number

  # build raw dataframe
  df <- data.frame(
    rna_barcode     = samps,
    sample_barcode  = meta_s$sample_barcode,
    patient_id      = meta_s$patient_id,
    surgery_number  = surgery_number,   # from metadata
    phate1          = emb[,1],
    phate2          = emb[,2],
    stringsAsFactors = FALSE
  )
  
  # raw progression axis = PHATE1
  df$raw_axis <- df$phate1
  
  progression_raw[[s]] <- df
}


# =====================================
# Step 3: EXPORT the raw axis before anchoring
# =====================================
saveRDS(progression_raw, "exports/progression_raw_axis.rds")

for(s in names(progression_raw)){
  write.csv(
    progression_raw[[s]],
    paste0("exports/progression_raw_axis_", s, ".csv"),
    row.names = FALSE
  )
}

message("✓ Raw progression axis exported.")


# =====================================
# Step 4: Anchor direction using paired surgeries
# =====================================
progression_anchored <- progression_raw

for(s in names(progression_anchored)) {
  
  message("Anchoring direction for stratum: ", s)
  
  df <- progression_anchored[[s]]
  
  deltas <- c()
  
  # For each patient trajectory (connected component)
  for(pid in unique(df$patient_id)) {
    
    sub <- df[df$patient_id == pid, ]
    
    if(nrow(sub) >= 2) {
      sub <- sub[order(sub$surgery_number), ]
      
      # change from earliest to latest surgery
      delta <- tail(sub$raw_axis, 1) - head(sub$raw_axis, 1)
      deltas <- c(deltas, delta)
    }
  }
  
  # If most deltas are negative → flip direction
  if(median(deltas, na.rm = TRUE) < 0) {
    df$raw_axis <- -df$raw_axis
  }
  
  progression_anchored[[s]] <- df
}

# =====================================
# Normalize axis to [0,1]
# =====================================
progression_norm <- progression_anchored

for(s in names(progression_norm)) {
  df <- progression_norm[[s]]
  
  min_val <- min(df$raw_axis)
  max_val <- max(df$raw_axis)
  
  df$progression_score <- (df$raw_axis - min_val) / (max_val - min_val)
  
  progression_norm[[s]] <- df
}

# =====================================
# Exports
# =====================================
saveRDS(progression_norm, "exports/progression/progression_scores.rds")

for(s in names(progression_norm)){
  write.csv(
    progression_norm[[s]],
    paste0("exports/progression/progression_scores_", s, ".csv"),
    row.names = FALSE
  )
}

message("✓ Progression direction anchored and normalized.")
message("✓ All results saved successfully.")
