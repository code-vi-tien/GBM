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
pairs = pd.read_csv("supplementary-data/pairs.csv")
clinical = pd.read_csv("supplementary-data/clinical_metadata.csv")
analyte = pd.read_csv("supplementary-data/analysis_analyte_set.csv")
aliquot_batch = pd.read_csv("supplementary-data/biospecimen_aliquots.csv")
gene_tpm = pd.read_csv('supplementary-data/gene_tpm_all.tsv', sep='\t')


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
# Annotate molecular subtype for paired and non-paired
# ================================================================
allowed_grades = ["II", "III", "IV"]
allowed_idh = ["IDHwt", "IDHmut"]

paired = annotate_molecular_subtype(
    paired_meta,
    allowed_grades,
    allowed_idh,
    parse_who
)
#paired.to_csv("subtyped_paired_metadata.csv", index=False)
nonpaired = annotate_molecular_subtype(
    nonpaired_meta,
    allowed_grades,
    allowed_idh,
    parse_who
)
#nonpaired.to_csv("subtyped_nonpaired_metadata.csv", index=False)

required_cols = [
    'case_barcode','sample_barcode','rna_barcode','molecular_subtype'
]

paired = paired.dropna(subset=['rna_barcode', "case_barcode", 'sample_barcode'])
nonpaired = nonpaired.dropna(subset=['rna_barcode', "case_barcode", 'sample_barcode'])

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

nonpaired['subtype'] = nonpaired.apply(map_subtype, axis=1)

samples_list = []

for (grade, subtype), n_samples in target.items():
    if n_samples == 0:
        continue
    subset = nonpaired[(nonpaired['grade'] == grade) & (nonpaired['subtype'] == subtype)]
    
    # If there are fewer rows than n_samples, take all rows
    n_select = min(len(subset), n_samples)
    
    selected = subset.sample(n=n_select, random_state=42)
    samples_list.append(selected)

nonpaired = pd.concat(samples_list)

# Reset index
nonpaired = nonpaired.reset_index(drop=True)


# ================================================================
# Combine subtype metadata
# ================================================================
subtyped_full = pd.concat([paired, nonpaired]).reset_index(drop=True)

# Export
subtyped_full.to_csv("subtyped_full_metadata.csv", index=False)
print("✔ Exported: subtyped_full_metadata.csv")
print(subtyped_full.isna().sum())

'''
# ================================================================
# Prepare TPM matrix (sample-matched)
# ================================================================
subtyped_full["rna_barcode_tpm"] = subtyped_full["rna_barcode"].apply(barcode_to_tpm_format)

# pull columns
cols = ["Gene_symbol"] + [c for c in subtyped_full["rna_barcode_tpm"] if c in tpm.columns]
tpm_samples = tpm[cols]

tpm_samples.to_csv("filtered_tpm.csv",index=False)
print("✔ Exported: filtered_tpm.csv")


# ================================================================
# Gene filtering
# ================================================================
expr = tpm_samples.set_index("Gene_symbol")

mask_half = (expr > 0).sum(axis=1) >= expr.shape[1] * 0.5
mask_mean = expr.mean(axis=1) >= 1

expr_filtered = expr[mask_half & mask_mean]
expr_filtered.to_csv("cleaned_tpm.csv")
print("✔ Exported: cleaned_tpm.csv")
'''
'''
# ================================================================
# Build relation matrix (for DM / OT)
# ================================================================
pairs["sample_a"] = pairs["tumor_barcode_a"].apply(tumor_to_sample)
pairs["sample_b"] = pairs["tumor_barcode_b"].apply(tumor_to_sample)

# join original hyphenated RNA barcodes
rmat = pairs.merge(
    meta[["sample_barcode","rna_barcode"]],
    left_on="sample_a", right_on="sample_barcode",
    how="left"
).rename(columns={"rna_barcode":"rna_a"})

rmat = rmat.merge(
    meta[["sample_barcode","rna_barcode"]],
    left_on="sample_b", right_on="sample_barcode",
    how="left"
).rename(columns={"rna_barcode":"rna_b"})

relation = rmat[[
    "case_barcode","rna_a","rna_b","surgical_interval_mo"
]].dropna()

relation.to_csv("relation_matrix.csv", index=False)
print("✔ Exported: relation_matrix.csv")
'''