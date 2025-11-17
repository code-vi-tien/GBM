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


# ================================================================
# Load input files
# ================================================================
pairs = pd.read_csv("supplementary-data/pairs.csv")
clinical = pd.read_csv("supplementary-data/clinical_metadata.csv")
analyte = pd.read_csv("supplementary-data/analysis_analyte_set.csv")
tpm = pd.read_csv("supplementary-data/gene_tpm_all.tsv", sep="\t")


# ================================================================
# Extract paired samples and annotate subtypes
# ================================================================
# 1. Extract tumor → sample barcodes
all_tumor = pd.concat([
    pairs["tumor_barcode_a"],
    pairs["tumor_barcode_b"]
], ignore_index=True).dropna()

sample_barcodes = list({tumor_to_sample(t) for t in all_tumor})

# 2. Subset the clinical metadata
paired_meta = clinical[clinical["sample_barcode"].isin(sample_barcodes)].copy()
nonpaired_meta = clinical[~clinical["sample_barcode"].isin(sample_barcodes)].copy()

print(f"✔ Paired samples: {paired_meta.shape[0]}")
print(f"✔ Non-paired samples: {nonpaired_meta.shape[0]}")

# 3. Normalize paired metadata
paired = paired_meta.copy()
paired["idh_norm"] = paired["idh_status"].apply(norm)
paired["codel_norm"] = paired["codel_status"].apply(norm)
paired["idh_clean"] = paired["idh_norm"].map(map_idh)
paired["codel_clean"] = paired["codel_norm"].map(map_codel)

paired["molecular_subtype"] = np.nan
paired.loc[(paired["idh_clean"]=="mut") & (paired["codel_clean"]=="codel"),  "molecular_subtype"] = "Oligodendroglioma"
paired.loc[(paired["idh_clean"]=="mut") & (paired["codel_clean"]=="intact"), "molecular_subtype"] = "Astrocytoma"
paired.loc[(paired["idh_clean"]=="wt"), "molecular_subtype"] = "Primary GBM"
paired.loc[paired["molecular_subtype"].isna(), "molecular_subtype"] = paired["who_classification"].apply(parse_who)


# ================================================================
# Subtype non-paired samples
# ================================================================
allowed_grades = ["II", "III", "IV"]
allowed_idh = ["IDHwt", "IDHmut"]

filtered_non = nonpaired_meta.loc[
    nonpaired_meta["grade"].isin(allowed_grades) &
    nonpaired_meta["idh_status"].isin(allowed_idh)
].copy()

filtered_non["idh_norm"]   = filtered_non["idh_status"].apply(norm)
filtered_non["codel_norm"] = filtered_non["codel_status"].apply(norm)
filtered_non["idh_clean"]  = filtered_non["idh_norm"].map(map_idh)
filtered_non["codel_clean"]= filtered_non["codel_norm"].map(map_codel)

filtered_non["molecular_subtype"] = np.nan
filtered_non.loc[(filtered_non["idh_clean"]=="mut") & (filtered_non["codel_clean"]=="codel"),"molecular_subtype"] = "Oligodendroglioma"
filtered_non.loc[(filtered_non["idh_clean"]=="mut") & (filtered_non["codel_clean"]=="intact"),"molecular_subtype"] = "Astrocytoma"
filtered_non.loc[(filtered_non["idh_clean"]=="wt"), "molecular_subtype"] = "Primary GBM"
filtered_non.loc[filtered_non["molecular_subtype"].isna(),"molecular_subtype"] = filtered_non["who_classification"].apply(parse_who)


# ================================================================
# Combine subtype metadata
# ================================================================
subtyped_full = pd.concat([paired, filtered_non]).reset_index(drop=True)

cols = [
    "case_barcode","sample_barcode","histology",
    "grade","idh_status","codel_status",
    "who_classification","mgmt_methylation",
    "molecular_subtype"
]

subtyped_full = subtyped_full[cols]

subtyped_full.to_csv("subtyped_full_metadata.tsv", sep="\t", index=False)
print("✔ Exported: subtyped_full_metadata.tsv")


# ================================================================
# Attach RNA barcodes from analyte set
# ================================================================
meta = subtyped_full.merge(
    analyte[["sample_barcode","rna_barcode"]],
    on="sample_barcode", how="right"
)

missing = meta["rna_barcode"].isna().sum()
print(f"✔ RNA barcode missing: {missing}")
if missing > 0:
    print(meta[meta["rna_barcode"].isna()])
    raise ValueError("Some samples missing RNA barcode.")


# ================================================================
# Prepare TPM matrix (sample-matched)
# ================================================================
meta["rna_barcode_tpm"] = meta["rna_barcode"].apply(barcode_to_tpm_format)

# pull columns
cols = ["Gene_symbol"] + [c for c in meta["rna_barcode_tpm"] if c in tpm.columns]
tpm_samples = tpm[cols]

tpm_samples.to_csv("filtered_tpm.tsv", sep="\t", index=False)
print("✔ Exported: filtered_tpm.tsv")


# ================================================================
# Gene filtering
# ================================================================
expr = tpm_samples.set_index("Gene_symbol")

mask_half = (expr > 0).sum(axis=1) >= expr.shape[1] * 0.5
mask_mean = expr.mean(axis=1) >= 1

expr_filtered = expr[mask_half & mask_mean]
expr_filtered.to_csv("cleaned_tpm.tsv", sep="\t")
print("✔ Exported: cleaned_tpm.tsv")


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

relation.to_csv("relation_matrix.tsv", sep="\t", index=False)
print("✔ Exported: relation_matrix.tsv")
