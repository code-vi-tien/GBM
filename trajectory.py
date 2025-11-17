import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import StandardScaler

# ----------------------------
# Utility functions
# ----------------------------
def normalize_sample_barcode(s: str) -> str:
    """Normalize sample barcode to match metadata."""
    if pd.isna(s):
        return s
    s2 = s.replace('.', '-')
    parts = s2.split('-')
    if len(parts) >= 4:
        return '-'.join(parts[:4])
    return s2

def load_tpm(tpm_path: Path) -> pd.DataFrame:
    df = pd.read_csv(tpm_path, sep='\t', header=0, index_col=0)
    df.index.name = 'Gene_symbol'
    return df

def load_metadata(meta_path: Path) -> pd.DataFrame:
    df = pd.read_csv(meta_path, sep='\t', dtype=str)
    df['sample_barcode'] = df['sample_barcode'].str.strip().map(normalize_sample_barcode)
    df['case_barcode'] = df['case_barcode'].str.strip()
    return df.drop_duplicates(subset='sample_barcode')

def load_pairs(pairs_path: Path) -> pd.DataFrame:
    df = pd.read_csv(pairs_path, sep='\t', dtype=str)
    df['rna_a'] = df['rna_a'].map(normalize_sample_barcode)
    df['rna_b'] = df['rna_b'].map(normalize_sample_barcode)
    if 'surgical_interval_mo' in df.columns:
        df['surgical_interval_mo'] = pd.to_numeric(df['surgical_interval_mo'], errors='coerce')
    else:
        df['surgical_interval_mo'] = np.nan
    return df

def identify_static_samples(meta_df: pd.DataFrame, pairs_df: pd.DataFrame, tpm_cols_normalized: set) -> set:
    """Samples in metadata + TPM but not in any longitudinal pair."""
    pairs_samples = set(pairs_df['rna_a'].dropna()).union(set(pairs_df['rna_b'].dropna()))
    meta_samples = set(meta_df['sample_barcode'].unique())
    return (meta_samples & tpm_cols_normalized) - pairs_samples

# ----------------------------
# Main DPT pipeline
# ----------------------------
def run_dpt_pipeline(tpm_path, meta_path, pairs_path, outdir: Path, n_pcs=30, random_state=42):
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading TPM...")
    tpm = load_tpm(Path(tpm_path))
    tpm.columns = [normalize_sample_barcode(c) for c in tpm.columns]

    print("Loading metadata...")
    meta = load_metadata(Path(meta_path))

    print("Loading longitudinal pairs...")
    pairs = load_pairs(Path(pairs_path))

    print("Identifying static samples...")
    static_samples = identify_static_samples(meta, pairs, set(tpm.columns))
    print(f"Found {len(static_samples)} static samples.")

    # ----------------------------
    # Build the manifold (PCA + diffusion map)
    # ----------------------------
    print("Preparing data for Scanpy...")
    # Only use static + longitudinal samples
    # Ensure all sample barcodes are strings and drop NaNs
    meta_samples = meta['sample_barcode'].dropna().astype(str)
    all_samples = sorted(list(static_samples) + list(set(meta_samples) - set(static_samples)))

    tpm = tpm.loc[:, tpm.columns.intersection(all_samples)]
    # log-transform
    adata = sc.AnnData(tpm.T)
    adata.var_names_make_unique()
    adata.obs = meta.set_index('sample_barcode').loc[adata.obs_names]

    print("Scaling data...")
    sc.pp.log1p(adata)  # log1p(TPM)
    sc.pp.scale(adata)

    print(f"Computing PCA ({n_pcs} components)...")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack', random_state=random_state)

    print("Computing neighbors...")
    sc.pp.neighbors(adata, n_pcs=n_pcs)

    print("Computing diffusion map...")
    sc.tl.diffmap(adata)

    # ----------------------------
    # Diffusion pseudotime
    # ----------------------------
    print("Selecting root sample for DPT (earliest static or longitudinal)...")
    # Choose a root: first static sample if exists, else first sample overall
    root_sample = list(static_samples)[0] if len(static_samples) > 0 else adata.obs_names[0]
    # Ensure root sample exists
    if root_sample not in adata.obs_names:
        root_sample = adata.obs_names[0]

    # Get integer index of root sample
    root_locs = np.where(adata.obs_names == root_sample)[0]
    if len(root_locs) == 0:
        raise ValueError(f"Root sample {root_sample} not found in adata.obs_names")
    root_idx = int(root_locs[0])  # take first occurrence

    # Set the root in .uns BEFORE calling dpt
    adata.uns['iroot'] = root_idx

    # Compute DPT pseudotime (legacy API)
    sc.tl.dpt(adata, n_dcs=10, n_branchings=0)


    # ----------------------------
    # Save outputs
    # ----------------------------
    print("Saving embedding and pseudotime...")
    embedding = pd.DataFrame(
        adata.obsm['X_diffmap'],
        index=adata.obs_names,
        columns=[f'DiffMap{i+1}' for i in range(adata.obsm['X_diffmap'].shape[1])]
    )
    embedding['pseudotime'] = adata.obs['dpt_pseudotime']

    # Merge metadata
    embedding = embedding.merge(meta, left_index=True, right_on='sample_barcode', how='left')
    embedding.set_index('sample_barcode', inplace=True)
    embedding.to_csv(outdir / 'dpt_trajectory.tsv', sep='\t')

    print(f"DPT trajectory saved: {outdir / 'dpt_trajectory.tsv'}")

    return adata, embedding

# ----------------------------
# Example usage
# ----------------------------
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="DPT-based trajectory")
    parser.add_argument('--tpm', default='exports/final/cleaned_tpm.tsv')
    parser.add_argument('--meta', default='exports/final/subtyped_full_metadata.tsv')
    parser.add_argument('--pairs', default='exports/final/relation_matrix.tsv')
    parser.add_argument('--outdir', default='results/dpt')
    parser.add_argument('--n_pcs', type=int, default=30)
    args = parser.parse_args()

    run_dpt_pipeline(args.tpm, args.meta, args.pairs, Path(args.outdir), n_pcs=args.n_pcs)
