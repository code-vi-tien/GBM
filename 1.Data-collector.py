import pandas as pd
import numpy as np


# ================================================================
# Utility functions
# ================================================================
# Normalization of strings
def norm(x):
    return x.strip().lower() if isinstance(x, str) else x

# Map values
map_idh = {
    "idhmut": "mut", "idh-mutant": "mut", "mut": "mut",
    "idhwt": "wt", "idh-wildtype": "wt", "wt": "wt"
}
map_codel = {
    "codel": "codel", "codeleted": "codel",
    "1p/19q-codeleted": "codel",
    "noncodel": "intact"
}
# Fallback: WHO classification parsing
def parse_who(x):
    if not isinstance(x, str):
        return np.nan

    s = x.lower()
    if "oligodendro" in s:
        return "Oligodendroglioma"
    if "astro" in s:
        return "Astrocytoma"
    if "glioblastoma" in s or "gbm" in s:
        return "Primary GBM"
    return np.nan

# Convert tumor_barcode → sample_barcode Ex: GLSS-19-0267-TP-RNA-XXXX  → GLSS-19-0267-TP
def tumor_to_sample(tumor_bc):
    return "-".join(tumor_bc.split("-")[:4])

def barcode_to_tpm_format(bc):
    return bc.replace("-", ".")

def annotate_molecular_subtype(df, allowed_grades, allowed_idh, parse_who_fn):
    df = df.copy()
    
    # Initialize molecular_subtype column
    df['molecular_subtype'] = pd.NA
    
    # Only consider rows with allowed grade and IDH
    mask = df['grade'].isin(allowed_grades) & df['idh_status'].isin(allowed_idh)
    
    # Assign molecular_subtype based on idh_codel_subtype
    df.loc[mask & (df['idh_codel_subtype'] == 'IDHmut-codel'), 'molecular_subtype'] = 'Oligodendroglioma'
    df.loc[mask & (df['idh_codel_subtype'] == 'IDHmut-noncodel'), 'molecular_subtype'] = 'Astrocytoma'
    df.loc[mask & (df['idh_codel_subtype'] == 'IDHwt'), 'molecular_subtype'] = 'Primary GBM'
    
    # For unresolved rows, use WHO classification
    unresolved = df['molecular_subtype'].isna()
    df.loc[unresolved, 'molecular_subtype'] = df.loc[unresolved, 'who_classification'].apply(parse_who_fn)
    
    return df

def map_subtype(row):
    if row['idh_codel_subtype'] == 'IDHmut-codel':
        return 'Oligo'
    elif row['idh_codel_subtype'] == 'IDHmut-noncodel':
        return 'Astro'
    elif row['idh_codel_subtype'] == 'IDHwt':
        return 'GBM'
    else:
        return None
    

# ================================================================
# Load input files
# ================================================================
pairs = pd.read_csv("og-supplementary-data/pairs.csv")
clinical = pd.read_csv("og-supplementary-data/clinical_metadata.csv")
analyte = pd.read_csv("og-supplementary-data/analysis_analyte_set.csv")
aliquot_batch = pd.read_csv("og-supplementary-data/biospecimen_aliquots.csv")
gene_tpm = pd.read_csv('og-supplementary-data/gene_tpm_all.tsv', sep='\t')
# Normalize tpm
new_columns = [
    col.replace('.', '-') if col != gene_tpm.columns[0] else col 
    for col in gene_tpm.columns
]
gene_tpm.columns = new_columns

# ================================================================
# Extract paired samples and annotate subtypes
# ================================================================
# 1. Add rna_barcode and batch data
analyte_unique = analyte.drop_duplicates(subset=['sample_barcode'])
clinical_with_rna = clinical.merge(analyte[['sample_barcode', 'rna_barcode']], on='sample_barcode', how='left')

# 2. Add aliquot_batch
aliquot_batch_clean = aliquot_batch.drop(
    columns=["aliquot_barcode", "aliquot_uuid_short", "aliquot_analyte_type", 
             "aliquot_analysis_type", "aliquot_portion"]
).drop_duplicates(subset=['sample_barcode'])

clinical_with_rna = clinical_with_rna.merge(aliquot_batch_clean, on='sample_barcode', how='left')

# 3. Extract tumor → sample barcodes
all_tumor = pd.concat([
    pairs["tumor_barcode_a"],
    pairs["tumor_barcode_b"]
], ignore_index=True).dropna().astype(str).str.strip().unique()

# 4. Subset the clinical metadata
paired_meta = clinical_with_rna[
    clinical_with_rna["rna_barcode"].notna() & clinical_with_rna["rna_barcode"].isin(all_tumor)
].copy()
cols = list(paired_meta.columns)
cols.insert(3, cols.pop(cols.index('rna_barcode')))
paired_meta = paired_meta[cols]
# Drop duplicated rna_barcodes in paired
paired_meta = paired_meta.drop_duplicates(subset='rna_barcode', keep='first')

nonpaired_meta = clinical_with_rna[
    clinical_with_rna["rna_barcode"].isna() | ~clinical_with_rna["rna_barcode"].isin(all_tumor)
].copy()
cols = list(nonpaired_meta.columns)
cols.insert(3, cols.pop(cols.index('rna_barcode')))
nonpaired_meta = nonpaired_meta[cols]

print(f"✔ Paired samples: {paired_meta.shape[0]}")
print(f"✔ Non-paired samples: {nonpaired_meta.shape[0]}")
#nonpaired_meta.to_csv("raw_nonpaired_metadata.csv", index=False)
#paired_meta.to_csv("raw_paired_metadata.csv", index=False)

# ================================================================
# Annotate molecular subtype for paired and non-paired and Filter
# ================================================================
allowed_grades = ["II", "III", "IV"]
allowed_idh = ["IDHwt", "IDHmut"]

rna_barcodes = gene_tpm.columns[1:]  # assuming first column is gene names
rna_barcodes = [bc.replace('.', '-') for bc in rna_barcodes]
rna_barcodes_set = set(rna_barcodes)

paired = annotate_molecular_subtype(
    paired_meta,
    allowed_grades,
    allowed_idh,
    parse_who
)
paired = paired[paired['rna_barcode'].isin(rna_barcodes_set)].copy()
print(f'Paired after filtered with tpm: {paired.shape}')
#paired.to_csv("subtyped_paired_metadata.csv", index=False)

nonpaired = annotate_molecular_subtype(
    nonpaired_meta,
    allowed_grades,
    allowed_idh,
    parse_who
)
nonpaired = nonpaired[nonpaired['rna_barcode'].isin(rna_barcodes_set)].copy()
print(f'Nonpaired after filtered with tpm: {nonpaired.shape}')
#nonpaired.to_csv("subtyped_nonpaired_metadata.csv", index=False)

required_cols = [
    'case_barcode','sample_barcode','rna_barcode','molecular_subtype'
]

paired = paired.dropna(subset=['rna_barcode', "case_barcode", 'sample_barcode'])
print(f'after drop na {paired.shape}')
nonpaired = nonpaired.dropna(subset=['rna_barcode', "case_barcode", 'sample_barcode'])
print(f'after drop na {nonpaired.shape}')

# Sampling according to target
target = {
    ('II','Oligo'): 12,
    ('II','Astro'): 20,
    ('II','GBM'): 3,
    ('III','Oligo'): 15,
    ('III','Astro'): 14,
    ('III','GBM'): 9,
    ('IV','Oligo'): 0,
    ('IV','Astro'): 7,
    ('IV','GBM'): 52
}

nonpaired = nonpaired.sort_values(
    by=['case_barcode', 'surgery_number'], # Use the actual chronological column here
    ascending=True
)
nonpaired = nonpaired.groupby('case_barcode', as_index=False).first()

nonpaired['subtype'] = nonpaired.apply(map_subtype, axis=1)

samples_list = []

for (grade, subtype), n_samples in target.items():
    if n_samples == 0:
        continue

    # Subset the compressed DataFrame
    subset = nonpaired[
        (nonpaired['grade'] == grade) &
        (nonpaired['subtype'] == subtype)
    ]
    
    # Identify unique cases (already unique in this DataFrame)
    unique_cases = subset['case_barcode'].unique()
    
    # Determine how many cases to select
    n_select = min(len(unique_cases), n_samples)
    
    # Sample the cases randomly
    sampled_cases = np.random.default_rng(seed=42).choice(
        unique_cases,
        size=n_select,
        replace=False
    )
    
    # Select the corresponding rows (which are the full samples)
    selected_rows = subset[subset['case_barcode'].isin(sampled_cases)]
    samples_list.append(selected_rows)

nonpaired = pd.concat(samples_list).reset_index(drop=True)
print(f'after selected nonpaired: {nonpaired.shape}')

# Final export of the result containing only the selected first samples with all original columns
#nonpaired.to_csv('selected_nonpaired_metadata.csv', index=False)


# ================================================================
# Combine subtype metadata
# ================================================================
subtyped_full = pd.concat([paired, nonpaired]).reset_index(drop=True)

# Export
#subtyped_full.to_csv("subtyped_full_metadata.csv", index=False)
print("✔ Exported: subtyped_full_metadata.csv")
print(f'full: {subtyped_full.shape}')


# ================================================================
# Prepare TPM matrix (sample-matched)
# ================================================================
required_barcodes = subtyped_full['rna_barcode'].tolist()
print(len(required_barcodes))

columns_to_keep = [gene_tpm.columns[0]] + required_barcodes
print(len(columns_to_keep))

available_barcodes = gene_tpm.columns.intersection(columns_to_keep)
print(available_barcodes.shape)

gene_tpm = gene_tpm[available_barcodes]
print(gene_tpm.shape)

#gene_tpm.to_csv("filtered_tpm.csv",index=False)
print("✔ Exported: filtered_tpm.csv")


# ================================================================
# Gene filtering
# ================================================================
expr = gene_tpm.set_index("Gene_symbol")

mask_half = (expr > 0).sum(axis=1) >= expr.shape[1] * 0.5
mask_mean = expr.mean(axis=1) >= 1

expr_filtered = expr[mask_half & mask_mean]
print(expr_filtered.shape)
#expr_filtered.to_csv("cleaned_tpm.csv")
print("✔ Exported: cleaned_tpm.csv")


# ================================================================
# Build relation matrix (for DM / OT)
# ================================================================
# Get sorted list of unique barcodes
all_rna_barcodes = sorted(subtyped_full['rna_barcode'].unique())
N = len(all_rna_barcodes)

# Map barcode → index
barcode_to_index = {barcode: i for i, barcode in enumerate(all_rna_barcodes)}

# Initialize adjacency matrix
relation_matrix_np = np.zeros((N, N), dtype=int)

# Filter pairs to valid barcodes
valid_pairs = pairs[
    (pairs['tumor_barcode_a'].isin(all_rna_barcodes)) &
    (pairs['tumor_barcode_b'].isin(all_rna_barcodes))
]

# Vectorized approach: get indices
idx_a = valid_pairs['tumor_barcode_a'].map(barcode_to_index).to_numpy()
idx_b = valid_pairs['tumor_barcode_b'].map(barcode_to_index).to_numpy()

# Populate adjacency matrix
relation_matrix_np[idx_a, idx_b] = 1
relation_matrix_np[idx_b, idx_a] = 1  # ensure undirected

# Convert to DataFrame
adjacency_matrix = pd.DataFrame(
    relation_matrix_np,
    index=all_rna_barcodes,
    columns=all_rna_barcodes
)

print(f"Adjacency matrix shape: {adjacency_matrix.shape}")

# Export
adjacency_matrix.to_csv("adjacency_matrix.csv")
print("✔ Exported: adjacency_matrix.csv")