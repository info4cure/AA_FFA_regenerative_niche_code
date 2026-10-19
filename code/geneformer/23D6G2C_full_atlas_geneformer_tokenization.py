#!/usr/bin/env python3

# ==============================================================================
# 23D6G2C — FULL-ATLAS GENEFORMER V2 TOKENIZATION
#
# Purpose
# -------
# Tokenize the frozen 23D6G1C full-atlas H5AD using the same Geneformer V2
# specification validated in 23D4C, while:
#
#   1. preserving exact cell identity explicitly;
#   2. preserving biological / execution metadata;
#   3. retaining uncropped token lengths for exact truncation adjudication;
#   4. performing structural validation only.
#
# NO Geneformer model inference.
# NO embeddings.
# NO perturbation.
# ==============================================================================

from pathlib import Path
import shutil
import sys

import anndata as ad
import numpy as np
import pandas as pd
from datasets import load_from_disk
from geneformer import TranscriptomeTokenizer


# ==============================================================================
# 1. PATHS
# ==============================================================================

ROOT = Path.cwd()

BASE = (
    ROOT
    / "PIPELINE_V2/results/disease_biology"
    / "23_ADVANCED_ANALYSIS_ROADMAP"
    / "23D_AI_PERTURBATION"
    / "23D6_PERTURBATION"
    / "23D6G_FULL_ATLAS_REASSESSMENT"
)

SRC = (
    BASE
    / "23D6G1C_FULL_ATLAS_H5AD"
    / "objects"
    / "23D6G1C_geneformer_full_atlas.h5ad"
)

OUT = BASE / "23D6G2C_FULL_ATLAS_TOKENIZATION"
INPUT = OUT / "input"
DATA = OUT / "dataset"
TAB = OUT / "tables"

for p in [OUT, INPUT, DATA, TAB]:
    p.mkdir(parents=True, exist_ok=True)

WORK_H5AD = INPUT / "23D6G2C_geneformer_full_atlas.h5ad"

DATASET_PATH = DATA / "23D6G2C_geneformer_full_atlas.dataset"


print("=" * 80)
print("23D6G2C — FULL-ATLAS GENEFORMER V2 TOKENIZATION")
print("STRUCTURAL VALIDATION + EXACT TRUNCATION AUDIT")
print("NO EMBEDDINGS — NO MODEL INFERENCE — NO PERTURBATION")
print("=" * 80)


# ==============================================================================
# 2. INPUT VALIDATION
# ==============================================================================

print("\n===== 1. INPUT VALIDATION =====")

assert SRC.exists(), f"Missing source H5AD: {SRC}"

adata = ad.read_h5ad(SRC)

print("source:", SRC)
print("shape:", adata.shape)
print("cells:", adata.n_obs)
print("genes:", adata.n_vars)

assert adata.n_obs == 38699, (
    f"Unexpected number of cells: {adata.n_obs}; expected 38699"
)

required_obs = [
    "patient_id",
    "disease",
    "module",
    "L2",
    "execution_universe_id",
    "full_atlas_order_index",
    "compartment",
    "stratum",
    "n_counts",
    "filter_pass",
]

missing_obs = [x for x in required_obs if x not in adata.obs.columns]

print("required obs missing:", missing_obs)

assert not missing_obs

assert "ensembl_id" in adata.var.columns

print("INPUT_VALIDATION=PASS")


# ==============================================================================
# 3. PRESERVE EXACT CELL ID
# ==============================================================================

print("\n===== 2. CELL-ID PRESERVATION =====")

assert adata.obs_names.is_unique, "H5AD obs_names are not unique."

adata.obs["cell_id"] = adata.obs_names.astype(str)

assert adata.obs["cell_id"].is_unique
assert adata.obs["cell_id"].notna().all()

print("unique cell IDs:", adata.obs["cell_id"].nunique())
print("CELL_ID_PRESERVATION=PASS")


# ==============================================================================
# 4. FILTER-PASS AUDIT
# ==============================================================================

print("\n===== 3. FILTER_PASS AUDIT =====")

filter_counts = adata.obs["filter_pass"].value_counts(dropna=False)

print(filter_counts.to_string())

n_filter_pass = int((adata.obs["filter_pass"] == 1).sum())

print("filter_pass == 1:", n_filter_pass)
print("total cells:", adata.n_obs)

assert n_filter_pass == adata.n_obs, (
    "Not all 23D6G1C cells have filter_pass == 1."
)

print("FILTER_PASS_AUDIT=PASS")


# ==============================================================================
# 5. ORDER INDEX AUDIT
# ==============================================================================

print("\n===== 4. ORDER INDEX AUDIT =====")

order_index = pd.to_numeric(
    adata.obs["full_atlas_order_index"],
    errors="raise",
).astype(int)

assert order_index.is_unique
assert order_index.min() == 1
assert order_index.max() == adata.n_obs
assert set(order_index.tolist()) == set(range(1, adata.n_obs + 1))

print("min:", order_index.min())
print("max:", order_index.max())
print("unique:", order_index.nunique())

print("ORDER_INDEX_AUDIT=PASS")


# ==============================================================================
# 6. WRITE DEDICATED TOKENIZATION INPUT
# ==============================================================================

print("\n===== 5. WRITE TOKENIZATION INPUT =====")

# Never alter the frozen 23D6G1C source.
# The only addition is explicit cell_id metadata.

if WORK_H5AD.exists():
    WORK_H5AD.unlink()

adata.write_h5ad(WORK_H5AD)

print("written:", WORK_H5AD)
print("size bytes:", WORK_H5AD.stat().st_size)

del adata

print("TOKENIZATION_INPUT=PASS")


# ==============================================================================
# 7. TOKENIZER CONFIGURATION
# ==============================================================================

print("\n===== 6. TOKENIZER CONFIGURATION =====")

custom_attrs = {
    "cell_id": "cell_id",
    "patient_id": "patient_id",
    "disease": "disease",
    "module": "module",
    "L2": "L2",
    "execution_universe_id": "execution_universe_id",
    "full_atlas_order_index": "full_atlas_order_index",
    "compartment": "compartment",
    "stratum": "stratum",
}

tk = TranscriptomeTokenizer(
    custom_attr_name_dict=custom_attrs,
    nproc=1,
    chunk_size=256,
    model_input_size=4096,
    special_token=True,
    collapse_gene_ids=True,
    use_h5ad_index=False,
    keep_counts=False,
    model_version="V2",
)

print("model_version:", tk.model_version)
print("model_input_size:", tk.model_input_size)
print("special_token:", tk.special_token)
print("collapse_gene_ids:", tk.collapse_gene_ids)
print("use_h5ad_index:", tk.use_h5ad_index)
print("keep_counts:", tk.keep_counts)
print("chunk_size:", tk.chunk_size)

assert tk.model_version == "V2"
assert tk.model_input_size == 4096
assert tk.special_token is True
assert tk.collapse_gene_ids is True
assert tk.use_h5ad_index is False
assert tk.keep_counts is False

print("TOKENIZER_CONFIGURATION=PASS")


# ==============================================================================
# 8. TOKENIZE FULL ATLAS
# ==============================================================================

print("\n===== 7. TOKENIZE FULL ATLAS =====")

# tokenize_files performs the same tokenization used by tokenize_data().
# create_dataset is called explicitly so that the pre-crop token sequence
# can be retained for exact truncation auditing.

tokenized_cells, cell_metadata, tokenized_counts = tk.tokenize_files(
    INPUT,
    file_format="h5ad",
    input_identifier="23D6G2C_geneformer_full_atlas",
)

print("raw tokenized cells:", len(tokenized_cells))

assert len(tokenized_cells) > 0

ds = tk.create_dataset(
    tokenized_cells,
    cell_metadata,
    tokenized_counts,
    use_generator=False,
    keep_uncropped_input_ids=True,
)

print("dataset rows:", len(ds))
print("dataset columns:", ds.column_names)

required_dataset_cols = [
    "input_ids",
    "input_ids_uncropped",
    "length",
    "length_uncropped",
    "cell_id",
    "patient_id",
    "disease",
    "module",
    "L2",
    "execution_universe_id",
    "full_atlas_order_index",
    "compartment",
    "stratum",
]

missing_ds = [x for x in required_dataset_cols if x not in ds.column_names]

print("required dataset columns missing:", missing_ds)

assert not missing_ds

print("TOKENIZATION=PASS")


# ==============================================================================
# 9. SAVE DATASET
# ==============================================================================

print("\n===== 8. SAVE TOKENIZED DATASET =====")

if DATASET_PATH.exists():
    shutil.rmtree(DATASET_PATH)

ds.save_to_disk(str(DATASET_PATH))

print("saved:", DATASET_PATH)

# Reload from disk so all downstream validation addresses the persisted object.

del ds
del tokenized_cells
del tokenized_counts
del cell_metadata

ds = load_from_disk(str(DATASET_PATH))

print("reloaded rows:", len(ds))
print("reloaded columns:", ds.column_names)

print("DATASET_SAVE_RELOAD=PASS")


# ==============================================================================
# 10. CELL-IDENTITY VALIDATION
# ==============================================================================

print("\n===== 9. CELL-IDENTITY VALIDATION =====")

tok_cell_ids = np.asarray(ds["cell_id"], dtype=str)

assert len(tok_cell_ids) == len(ds)
assert len(set(tok_cell_ids.tolist())) == len(tok_cell_ids)

source = ad.read_h5ad(WORK_H5AD, backed="r")

source_ids = source.obs_names.astype(str).to_numpy()

source.file.close()

source_set = set(source_ids.tolist())
token_set = set(tok_cell_ids.tolist())

missing_after_tokenization = sorted(source_set - token_set)
unexpected_after_tokenization = sorted(token_set - source_set)

print("source cells:", len(source_ids))
print("tokenized cells:", len(tok_cell_ids))
print("missing after tokenization:", len(missing_after_tokenization))
print("unexpected after tokenization:", len(unexpected_after_tokenization))

assert len(missing_after_tokenization) == 0, (
    f"{len(missing_after_tokenization)} source cells absent after tokenization."
)

assert len(unexpected_after_tokenization) == 0

assert len(tok_cell_ids) == 38699

print("CELL_IDENTITY_VALIDATION=PASS")


# ==============================================================================
# 11. ORDER PRESERVATION
# ==============================================================================

print("\n===== 10. ORDER PRESERVATION =====")

tok_order = np.asarray(
    ds["full_atlas_order_index"],
    dtype=int,
)

expected_order = np.arange(1, len(ds) + 1, dtype=int)

order_preserved = bool(np.array_equal(tok_order, expected_order))

print("order preserved:", order_preserved)
print("first order index:", tok_order[0])
print("last order index:", tok_order[-1])

assert order_preserved

print("ORDER_PRESERVATION=PASS")


# ==============================================================================
# 12. LENGTH CONSISTENCY
# ==============================================================================

print("\n===== 11. TOKEN LENGTH CONSISTENCY =====")

input_ids = ds["input_ids"]
input_ids_uncropped = ds["input_ids_uncropped"]

stored_lengths = np.asarray(ds["length"], dtype=int)
stored_uncropped = np.asarray(ds["length_uncropped"], dtype=int)

observed_lengths = np.asarray(
    [len(x) for x in input_ids],
    dtype=int,
)

observed_uncropped = np.asarray(
    [len(x) for x in input_ids_uncropped],
    dtype=int,
)

length_match = bool(
    np.array_equal(stored_lengths, observed_lengths)
)

uncropped_match = bool(
    np.array_equal(stored_uncropped, observed_uncropped)
)

print("cropped length field matches input_ids:", length_match)
print(
    "uncropped length field matches input_ids_uncropped:",
    uncropped_match,
)

assert length_match
assert uncropped_match

print("TOKEN_LENGTH_CONSISTENCY=PASS")


# ==============================================================================
# 13. MODEL-LIMIT AUDIT
# ==============================================================================

print("\n===== 12. MODEL-LIMIT AUDIT =====")

n = len(ds)

n_lt = int(np.sum(observed_lengths < 4096))
n_eq = int(np.sum(observed_lengths == 4096))
n_gt = int(np.sum(observed_lengths > 4096))

print("n:", n)
print("min:", int(observed_lengths.min()))
print("median:", float(np.median(observed_lengths)))
print("mean:", float(np.mean(observed_lengths)))
print("max:", int(observed_lengths.max()))
print("<4096:", n_lt)
print("==4096:", n_eq)
print(">4096:", n_gt)
print("fraction ==4096:", n_eq / n)

assert n_gt == 0
assert observed_lengths.max() <= 4096

print("MODEL_LIMIT_AUDIT=PASS")


# ==============================================================================
# 14. EXACT TRUNCATION AUDIT
# ==============================================================================

print("\n===== 13. EXACT TRUNCATION AUDIT =====")

# input_ids_uncropped contains gene tokens BEFORE CLS/EOS insertion and
# BEFORE cropping. With special_token=True, at most 4094 gene tokens can
# be retained in the final sequence:
#
#   CLS + 4094 genes + EOS = 4096 tokens
#
# Therefore:
#
#   length_uncropped > 4094  -> true gene-token truncation
#   length_uncropped ==4094  -> exactly fills model capacity, no truncation
#   length_uncropped < 4094  -> no truncation

true_truncated = observed_uncropped > 4094
exact_capacity = observed_uncropped == 4094

n_true_truncated = int(true_truncated.sum())
n_exact_capacity = int(exact_capacity.sum())

print("uncropped min:", int(observed_uncropped.min()))
print("uncropped median:", float(np.median(observed_uncropped)))
print("uncropped mean:", float(np.mean(observed_uncropped)))
print("uncropped max:", int(observed_uncropped.max()))

print("true truncated (>4094 gene tokens):", n_true_truncated)
print("exact capacity (=4094 gene tokens):", n_exact_capacity)
print("fraction true truncated:", n_true_truncated / n)

# Every genuinely truncated sequence must finish at model length 4096.

if n_true_truncated > 0:
    assert np.all(observed_lengths[true_truncated] == 4096)

# Every final sequence must equal:
# min(raw gene-token count, 4094) + CLS + EOS.

expected_final_lengths = (
    np.minimum(observed_uncropped, 4094) + 2
)

assert np.array_equal(
    observed_lengths,
    expected_final_lengths,
)

print("EXACT_TRUNCATION_AUDIT=PASS")


# ==============================================================================
# 15. TOKEN BOUNDARY STRUCTURE
# ==============================================================================

print("\n===== 14. TOKEN BOUNDARY STRUCTURE =====")

first_tokens = np.asarray(
    [x[0] for x in input_ids],
    dtype=int,
)

last_tokens = np.asarray(
    [x[-1] for x in input_ids],
    dtype=int,
)

unique_first = sorted(set(first_tokens.tolist()))
unique_last = sorted(set(last_tokens.tolist()))

print("unique first tokens:", unique_first)
print("unique last tokens:", unique_last)

first_token_2 = bool(np.all(first_tokens == 2))
last_token_3 = bool(np.all(last_tokens == 3))

print("ALL_FIRST_TOKEN_2 =", first_token_2)
print("ALL_LAST_TOKEN_3 =", last_token_3)

assert first_token_2
assert last_token_3

print("TOKEN_BOUNDARY_STRUCTURE=PASS")


# ==============================================================================
# 16. CELL-LEVEL TABLE
# ==============================================================================

print("\n===== 15. CELL-LEVEL TOKENIZATION TABLE =====")

cell_table = pd.DataFrame({
    "cell_id": tok_cell_ids,
    "full_atlas_order_index": tok_order,
    "patient_id": ds["patient_id"],
    "disease": ds["disease"],
    "module": ds["module"],
    "L2": ds["L2"],
    "execution_universe_id": ds["execution_universe_id"],
    "compartment": ds["compartment"],
    "stratum": ds["stratum"],
    "length_uncropped": observed_uncropped,
    "length": observed_lengths,
    "at_model_limit": observed_lengths == 4096,
    "true_truncated": true_truncated,
    "exact_gene_capacity": exact_capacity,
})

cell_table = cell_table.sort_values(
    "full_atlas_order_index"
).reset_index(drop=True)

cell_table.to_csv(
    TAB / "23D6G2C_cell_level_tokenization.csv",
    index=False,
)

print("rows:", len(cell_table))
print("true truncated:", int(cell_table["true_truncated"].sum()))

print("CELL_LEVEL_TABLE=PASS")


# ==============================================================================
# 17. SUMMARY FUNCTION
# ==============================================================================

def length_summary(df, group_cols):
    rows = []

    grouped = df.groupby(
        group_cols,
        dropna=False,
        observed=True,
    )

    for key, x in grouped:
        if not isinstance(key, tuple):
            key = (key,)

        lengths = x["length"].to_numpy()
        uncropped = x["length_uncropped"].to_numpy()

        row = dict(zip(group_cols, key))

        row.update({
            "n": len(x),
            "min": int(lengths.min()),
            "median": float(np.median(lengths)),
            "mean": float(np.mean(lengths)),
            "q95": float(np.quantile(lengths, 0.95)),
            "max": int(lengths.max()),
            "n_eq_4096": int(np.sum(lengths == 4096)),
            "fraction_eq_4096": float(np.mean(lengths == 4096)),
            "n_true_truncated": int(np.sum(uncropped > 4094)),
            "fraction_true_truncated": float(
                np.mean(uncropped > 4094)
            ),
            "max_length_uncropped": int(uncropped.max()),
        })

        rows.append(row)

    return pd.DataFrame(rows)


# ==============================================================================
# 18. LENGTH BY DISEASE
# ==============================================================================

print("\n===== 18. LENGTH BY DISEASE =====")

by_disease = length_summary(
    cell_table,
    ["disease"],
).sort_values("disease")

by_disease.to_csv(
    TAB / "23D6G2C_length_by_disease.csv",
    index=False,
)

print(by_disease.to_string(index=False))


# ==============================================================================
# 19. LENGTH BY COMPARTMENT
# ==============================================================================

print("\n===== 19. LENGTH BY COMPARTMENT =====")

by_compartment = length_summary(
    cell_table,
    ["compartment"],
).sort_values("compartment")

by_compartment.to_csv(
    TAB / "23D6G2C_length_by_compartment.csv",
    index=False,
)

print(by_compartment.to_string(index=False))


# ==============================================================================
# 20. LENGTH BY L2
# ==============================================================================

print("\n===== 20. LENGTH BY L2 =====")

by_L2 = length_summary(
    cell_table,
    ["compartment", "L2"],
).sort_values(
    ["fraction_true_truncated", "n"],
    ascending=[False, False],
)

by_L2.to_csv(
    TAB / "23D6G2C_length_by_L2.csv",
    index=False,
)

print(by_L2.to_string(index=False))


# ==============================================================================
# 21. CELLS AT MODEL LIMIT / TRUE TRUNCATION
# ==============================================================================

print("\n===== 21. CELLS AT MODEL LIMIT =====")

at_limit = cell_table.loc[
    cell_table["length"] == 4096
].copy()

at_limit = at_limit.sort_values(
    "full_atlas_order_index"
)

at_limit.to_csv(
    TAB / "23D6G2C_cells_at_model_limit.csv",
    index=False,
)

print("cells at 4096:", len(at_limit))
print(
    "true truncated:",
    int(at_limit["true_truncated"].sum()),
)

if len(at_limit) > 0:
    cols = [
        "cell_id",
        "patient_id",
        "disease",
        "compartment",
        "L2",
        "length",
        "length_uncropped",
        "true_truncated",
    ]
    print(
        at_limit[cols]
        .head(50)
        .to_string(index=False)
    )


# ==============================================================================
# 22. FORMAL CHECKS
# ==============================================================================

print("\n===== 22. FORMAL CHECKS =====")

n_cells = len(cell_table)
n_eq = int(np.sum(cell_table["length"] == 4096))
n_gt = int(np.sum(cell_table["length"] > 4096))
n_true_truncated = int(
    np.sum(cell_table["length_uncropped"] > 4094)
)

formal_checks = pd.DataFrame([
    {
        "check": "tokenized_cells_match_source",
        "status": (
            "PASS"
            if len(cell_table) == len(source_ids)
            else "FAIL"
        ),
    },
    {
        "check": "cell_ids_unique",
        "status": (
            "PASS"
            if cell_table["cell_id"].is_unique
            else "FAIL"
        ),
    },
    {
        "check": "cell_ids_match_source",
        "status": (
            "PASS"
            if set(cell_table["cell_id"]) == set(source_ids)
            else "FAIL"
        ),
    },
    {
        "check": "full_atlas_order_preserved",
        "status": (
            "PASS"
            if (
                cell_table["full_atlas_order_index"].is_monotonic_increasing
                and cell_table["cell_id"].tolist() == list(source_ids)
            )
            else "FAIL"
        ),
    },
    {
        "check": "length_field_matches_input_ids",
        "status": (
            "PASS"
            if length_match
            else "FAIL"
        ),
    },
    {
        "check": "max_length_le_4096",
        "status": (
            "PASS"
            if observed_lengths.max() <= 4096
            else "FAIL"
        ),
    },
    {
        "check": "no_sequence_gt_4096",
        "status": (
            "PASS"
            if n_gt == 0
            else "FAIL"
        ),
    },
    {
        "check": "uncropped_length_available",
        "status": (
            "PASS"
            if "length_uncropped" in ds.column_names
            else "FAIL"
        ),
    },
    {
        "check": "true_truncation_quantified",
        "status": "PASS",
    },
    {
        "check": "all_first_token_2",
        "status": (
            "PASS"
            if first_token_2
            else "FAIL"
        ),
    },
    {
        "check": "all_last_token_3",
        "status": (
            "PASS"
            if last_token_3
            else "FAIL"
        ),
    },
    {
        "check": "no_geneformer_model_execution",
        "status": "PASS",
    },
    {
        "check": "no_perturbation",
        "status": "PASS",
    },
])

formal_checks.to_csv(
    TAB / "23D6G2C_formal_checks.csv",
    index=False,
)

print(formal_checks.to_string(index=False))

assert (formal_checks["status"] == "PASS").all()

print("FORMAL_CHECKS=PASS")


# ==============================================================================
# 23. AUDIT
# ==============================================================================

print("\n===== 23. AUDIT =====")

audit_lines = [
    "MODULE=23D6G2C",
    "TITLE=Full-atlas Geneformer V2 tokenization and truncation audit",
    "SOURCE_MODULE=23D6G1C",
    "SOURCE_H5AD=23D6G1C_geneformer_full_atlas.h5ad",
    "MODEL_VERSION=Geneformer_V2",
    f"SOURCE_CELLS={len(source_ids)}",
    f"TOKENIZED_CELLS={len(cell_table)}",
    "MODEL_INPUT_SIZE=4096",
    f"MIN_LENGTH={int(observed_lengths.min())}",
    f"MEDIAN_LENGTH={float(np.median(observed_lengths))}",
    f"MEAN_LENGTH={float(np.mean(observed_lengths))}",
    f"MAX_LENGTH={int(observed_lengths.max())}",
    f"N_EQ_4096={n_eq}",
    f"N_GT_4096={n_gt}",
    f"FRACTION_EQ_4096={n_eq / len(observed_lengths)}",
    f"N_TRUE_TRUNCATED={n_true_truncated}",
    f"FRACTION_TRUE_TRUNCATED={n_true_truncated / len(cell_table)}",
    f"MAX_LENGTH_UNCROPPED={int(cell_table['length_uncropped'].max())}",
    f"UNIQUE_FIRST_TOKENS={unique_first}",
    f"UNIQUE_LAST_TOKENS={unique_last}",
    "LENGTH_FIELD_MATCHES_INPUT_IDS=TRUE",
    "NO_SEQUENCE_EXCEEDS_MODEL_LIMIT=TRUE",
    "UNCROPPED_LENGTH_PRESERVED=TRUE",
    "TRUE_TRUNCATION_DEFINED_AS_LENGTH_UNCROPPED_GT_4094=TRUE",
    "ALL_FIRST_TOKEN_2=TRUE",
    "ALL_LAST_TOKEN_3=TRUE",
    "TOKEN_BOUNDARY_STRUCTURE=PASS",
    "CELL_IDS_PRESERVED=TRUE",
    "FULL_ATLAS_ORDER_PRESERVED=TRUE",
    "GENEFORMER_MODEL_EXECUTION_PERFORMED=FALSE",
    "PERTURBATION_PERFORMED=FALSE",
    "STATUS=PASS",
]

audit = "\n".join(audit_lines) + "\n"

audit_path = OUT / "23D6G2C_audit.txt"
audit_path.write_text(audit)

print(audit)


# ==============================================================================
# 24. FINAL INTERPRETATION
# ==============================================================================

print("""
================================================================================
23D6G2C COMPLETE

INTERPRETATION:
- Full-atlas Geneformer V2 tokenization completed and passed structural
  validation.
- Cell identity and full-atlas ordering were explicitly preserved through
  tokenization.
- Stored token lengths match the observed input_ids lengths.
- No token sequence exceeds the 4096-token model input limit.
- Uncropped token lengths were retained before Geneformer's 4094-gene-token
  crop plus CLS/EOS insertion.
- Therefore, true boundary truncation is directly quantified as
  length_uncropped > 4094 rather than inferred solely from final length=4096.
- Cells with final length 4096 are explicitly identified.
- All tokenized sequences reproduce the expected Geneformer V2 boundary
  structure: first token=2 (CLS) and last token=3 (EOS).
- Cell-level disease, compartment, L2, stratum, patient and ordering metadata
  are preserved for downstream stratified analyses.
- No Geneformer model inference was performed in this validation step.
- No perturbation analysis was performed.

OUTPUTS:
- dataset/23D6G2C_geneformer_full_atlas.dataset
- tables/23D6G2C_cell_level_tokenization.csv
- tables/23D6G2C_length_by_disease.csv
- tables/23D6G2C_length_by_compartment.csv
- tables/23D6G2C_length_by_L2.csv
- tables/23D6G2C_cells_at_model_limit.csv
- tables/23D6G2C_formal_checks.csv
- 23D6G2C_audit.txt

STATUS=PASS
================================================================================
""")
