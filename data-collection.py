import pandas as pd
import numpy as np
import re

# --------------------------------------------------------
# Load data
# --------------------------------------------------------
pairs = pd.read_csv("supplementary-data/pairs.csv")
clinical = pd.read_csv("supplementary-data/clinical_metadata.csv")

# --------------------------------------------------------
# 1. Extract all unique tumor barcodes
# --------------------------------------------------------
all_tumor = pd.concat([
    pairs["tumor_barcode_a"],
    pairs["tumor_barcode_b"]
], ignore_index=True).dropna()

all_tumor = all_tumor.unique()

# --------------------------------------------------------
# 2. Convert tumor_barcode → sample_barcode
# --------------------------------------------------------
def tumor_to_sample(tumor):
    parts = tumor.split("-")
    return "-".join(parts[:4])  # KEEP FIRST 4 PARTS ONLY

sample_barcodes = list({tumor_to_sample(t) for t in all_tumor})

# --------------------------------------------------------
# 3. Inner join to clinical metadata → metadata for paired samples
# --------------------------------------------------------
paired_metadata = clinical[clinical["sample_barcode"].isin(sample_barcodes)].copy()

# --------------------------------------------------------
# 4. Extract the non-paired sample metadata
# --------------------------------------------------------
nonpaired_metadata = clinical[~clinical["sample_barcode"].isin(sample_barcodes)].copy()

print("Paired samples:", paired_metadata.shape)
print("Non-paired samples:", nonpaired_metadata.shape)
#paired_metadata.to_csv("paired_metadata.csv", index=False)
#nonpaired_metadata.to_csv("nonpaired_metadata.csv", index=False)

# --------------------------------------------------------
# 6. Sanity checking
# --------------------------------------------------------
bad_tumor_format = [t for t in all_tumor if len(t.split("-")) < 4]
if bad_tumor_format:
    print("❌ Bad tumor barcode format:", bad_tumor_format)
else:
    print("✓ All tumor barcodes have ≥ 4 hyphen-separated parts")

# ---- Sanity Check 2: uniqueness mapping ----
tumor_to_sample_map = {t: tumor_to_sample(t) for t in all_tumor}
if len(tumor_to_sample_map.values()) != len(set(tumor_to_sample_map.values())):
    print("⚠ Warning: multiple tumor barcodes map to the SAME sample barcode (this is normal for replicates).")
else:
    print("✓ 1:1 tumor→sample mapping")

# ---- Sanity Check 3: sample exists in clinical metadata ----
missing_samples = [s for s in sample_barcodes if s not in clinical["sample_barcode"].values]
if missing_samples:
    print("❌ Some sample_barcodes not found in clinical metadata:", missing_samples)
else:
    print("✓ All sample_barcodes found in clinical_metadata")

# ---- Sanity Check 4: paired_meta should contain all samples ----
if len(paired_metadata) != len(sample_barcodes):
    print(f"⚠ Mismatch: {len(sample_barcodes)} extracted samples vs {len(paired_metadata)} clinical matches")
else:
    print("✓ paired_metadata perfectly matches extracted sample_barcodes")

print("Sanity checks complete.")

# --------------------------------------------------------
# 5. Filter non-paired metadata
# --------------------------------------------------------
allowed_grades = ["II", "III", "IV"]
allowed_idh = ["IDHwt", "IDHmut"]

filtered_nonpaired = (
    nonpaired_metadata
    .loc[
        nonpaired_metadata["grade"].isin(allowed_grades) &
        nonpaired_metadata["idh_status"].isin(allowed_idh)
    ][[
        "case_barcode",
        "sample_barcode",
        "histology",
        "grade",
        "idh_status",
        "codel_status",
        "who_classification",
        "mgmt_methylation"
    ]]
)

#filtered_nonpaired.to_csv("nonpaired_filtered.csv", index=False)

# --------------------------------------------------------
# Subtyping and classifications of tumors 🚨🚨🚨
# --------------------------------------------------------
df = filtered_nonpaired.copy()
# --- Normalization ---
def norm(x):
    if isinstance(x, str):
        return x.strip().lower()
    return x

df["idh_norm"] = df["idh_status"].apply(norm)
df["codel_norm"] = df["codel_status"].apply(norm)

map_idh = {
    "idhmut": "mut",
    "idh-mutant": "mut",
    "mut": "mut",
    "idhwt": "wt",
    "idh-wildtype": "wt",
    "wt": "wt"
}

map_codel = {
    "codel": "codel",
    "codeleted": "codel",
    "1p/19q-codeleted": "codel",
    "noncodel": "intact",
    "intact": "intact"
}

df["idh_clean"] = df["idh_norm"].map(map_idh)
df["codel_clean"] = df["codel_norm"].map(map_codel)

# --- Build subtype from IDH + CODEL ---
df["molecular_subtype"] = np.nan

df.loc[(df["idh_clean"] == "mut") & (df["codel_clean"] == "codel"),
       "molecular_subtype"] = "Oligodendroglioma"

df.loc[(df["idh_clean"] == "mut") & (df["codel_clean"] == "intact"),
       "molecular_subtype"] = "Astrocytoma"

df.loc[(df["idh_clean"] == "wt"),
       "molecular_subtype"] = "Primary GBM"

# --- Fallback: WHO classification parsing ---
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

df.loc[df["molecular_subtype"].isna(), "molecular_subtype"] = \
    df["who_classification"].apply(parse_who)

print(df["molecular_subtype"].value_counts(dropna=False))
print(pd.crosstab(df["grade"], df["molecular_subtype"], dropna=False))