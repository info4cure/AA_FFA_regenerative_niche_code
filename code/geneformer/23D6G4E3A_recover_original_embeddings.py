#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path.cwd()

BASE = (
    ROOT
    / "PIPELINE_V2/results/disease_biology"
    / "23_ADVANCED_ANALYSIS_ROADMAP"
    / "23D_AI_PERTURBATION"
    / "23D6_PERTURBATION"
    / "23D6G_FULL_ATLAS_REASSESSMENT"
)

E1 = (
    BASE
    / "23D6G4E_PATIENT_ROBUSTNESS_HC_DIRECTIONALITY"
    / "23D6G4E1_EXACT_EXECUTION_MAPPING"
    / "tables"
    / "23D6G4E1_exact_execution_cell_mapping.csv"
)

G3B = (
    BASE
    / "23D6G3B_FULL_ATLAS_EMBEDDINGS"
    / "objects"
    / "23D6G3B_full_atlas.csv"
)

OUT = (
    BASE
    / "23D6G4E_PATIENT_ROBUSTNESS_HC_DIRECTIONALITY"
    / "23D6G4E3A_ORIGINAL_EMBEDDINGS"
)

OBJ = OUT / "objects"
TAB = OUT / "tables"

OBJ.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("23D6G4E3A — RECOVER ORIGINAL EXECUTION-CELL EMBEDDINGS")
print("=" * 80)

# ============================================================
# 1. LOAD EXACT EXECUTION MAPPING
# ============================================================

print("\n===== 1. LOAD EXACT EXECUTION MAPPING =====")

m = pd.read_csv(E1)

assert len(m) == 9033
assert m["pilot_perturbation_id"].nunique() == 24
assert m["cell_id"].notna().all()

print("mapping rows:", len(m))
print("targets:", m["pilot_perturbation_id"].nunique())
print("unique execution cell IDs:", m["cell_id"].nunique())

# Preserve exact E1 row order.
m = m.reset_index(drop=True)
m["E3A_embedding_row"] = np.arange(len(m), dtype=int)

# ============================================================
# 2. LOAD G3B
# ============================================================

print("\n===== 2. LOAD G3B FULL EMBEDDING MATRIX =====")

# 768 embedding columns + metadata needed for exact matching.
emb_cols = [str(i) for i in range(768)]

usecols = (
    emb_cols
    + [
        "cell_id",
        "full_atlas_order_index",
        "patient_id",
        "disease",
        "L2",
    ]
)

g = pd.read_csv(
    G3B,
    usecols=usecols,
)

print("G3B rows:", len(g))
print("embedding dimensions:", len(emb_cols))

assert len(g) == 38699
assert g["cell_id"].nunique() == 38699
assert g["full_atlas_order_index"].nunique() == 38699

# ============================================================
# 3. EXACT CELL-ID JOIN
# ============================================================

print("\n===== 3. EXACT CELL-ID JOIN =====")

meta_cols = [
    "cell_id",
    "full_atlas_order_index",
    "patient_id",
    "disease",
    "L2",
]

g_meta = g[meta_cols].copy()

z = m.merge(
    g_meta,
    on="cell_id",
    how="left",
    validate="many_to_one",
    suffixes=("_E1", "_G3B"),
)

print("joined rows:", len(z))

assert len(z) == 9033
assert z["full_atlas_order_index_G3B"].notna().all()

# ============================================================
# 4. FORMAL METADATA CONCORDANCE
# ============================================================

print("\n===== 4. METADATA CONCORDANCE =====")

index_match = (
    pd.to_numeric(
        z["full_atlas_order_index_E1"],
        errors="raise",
    ).astype(int).to_numpy()
    ==
    pd.to_numeric(
        z["full_atlas_order_index_G3B"],
        errors="raise",
    ).astype(int).to_numpy()
)

patient_match = (
    z["patient_id_exact"].astype(str).to_numpy()
    ==
    z["patient_id_G3B"].astype(str).to_numpy()
)

# E1 may retain final_L2_label while G3B uses L2.
# E1 already contains L2, therefore the merge suffixes it.
# Compare the frozen perturbation-panel L2 against the G3B L2.
l2_g3b_col = (
    "L2_G3B"
    if "L2_G3B" in z.columns
    else "L2"
)

l2_match = (
    z["final_L2_label"].astype(str).to_numpy()
    ==
    z[l2_g3b_col].astype(str).to_numpy()
)

print("index matches:", int(index_match.sum()), "/", len(z))
print("patient matches:", int(patient_match.sum()), "/", len(z))
print("L2 matches:", int(l2_match.sum()), "/", len(z))

assert index_match.all()
assert patient_match.all()
assert l2_match.all()

# ============================================================
# 5. RECOVER EMBEDDINGS IN EXACT E1 ORDER
# ============================================================

print("\n===== 5. RECOVER ORIGINAL EMBEDDINGS =====")

# Map cell_id -> G3B dataframe row.
g_index = pd.Series(
    np.arange(len(g), dtype=int),
    index=g["cell_id"].astype(str),
)

row_idx = (
    m["cell_id"]
    .astype(str)
    .map(g_index)
)

assert row_idx.notna().all()

row_idx = row_idx.astype(int).to_numpy()

full_matrix = g[emb_cols].to_numpy(dtype=np.float32)

original = full_matrix[row_idx, :]

print("original embedding shape:", original.shape)
print("dtype:", original.dtype)

assert original.shape == (9033, 768)
assert np.isfinite(original).all()

# ============================================================
# 6. SAVE COMPACT EMBEDDING MATRIX
# ============================================================

print("\n===== 6. SAVE =====")

embedding_file = (
    OBJ
    / "23D6G4E3A_original_execution_embeddings.npy"
)

np.save(
    embedding_file,
    original,
)

# ============================================================
# 7. SAVE EXACT EMBEDDING MANIFEST
# ============================================================

manifest_cols = [
    "E3A_embedding_row",
    "pilot_perturbation_id",
    "cell_id",
    "patient_id_exact",
    "disease",
    "final_L2_label",
    "gene",
    "ensembl_id",
    "geneformer_perturb_type",
    "full_atlas_order_index",
]

manifest_cols = [
    c for c in manifest_cols
    if c in m.columns
]

manifest = m[manifest_cols].copy()

manifest["embedding_dim"] = 768

manifest["embedding_norm"] = np.linalg.norm(
    original,
    axis=1,
)

manifest["embedding_file"] = str(
    embedding_file.relative_to(ROOT)
)

manifest.to_csv(
    TAB
    / "23D6G4E3A_original_embedding_manifest.csv",
    index=False,
)

# ============================================================
# 8. TARGET COUNTS
# ============================================================

counts = (
    manifest
    .groupby(
        "pilot_perturbation_id",
        observed=True,
    )
    .size()
    .rename("n_original_embeddings")
    .reset_index()
)

counts.to_csv(
    TAB
    / "23D6G4E3A_target_embedding_counts.csv",
    index=False,
)

print(counts.to_string(index=False))

# ============================================================
# 9. FORMAL CHECKS
# ============================================================

print("\n===== 9. FORMAL CHECKS =====")

checks = pd.DataFrame(
    [
        {
            "check": "execution_mapping_rows_9033",
            "status": "PASS" if len(m) == 9033 else "FAIL",
            "detail": len(m),
        },
        {
            "check": "targets_24",
            "status": (
                "PASS"
                if m["pilot_perturbation_id"].nunique() == 24
                else "FAIL"
            ),
            "detail": m["pilot_perturbation_id"].nunique(),
        },
        {
            "check": "G3B_rows_38699",
            "status": "PASS" if len(g) == 38699 else "FAIL",
            "detail": len(g),
        },
        {
            "check": "exact_cell_mapping_complete",
            "status": "PASS" if row_idx.size == 9033 else "FAIL",
            "detail": row_idx.size,
        },
        {
            "check": "full_atlas_index_concordance",
            "status": "PASS" if index_match.all() else "FAIL",
            "detail": int(index_match.sum()),
        },
        {
            "check": "patient_concordance",
            "status": "PASS" if patient_match.all() else "FAIL",
            "detail": int(patient_match.sum()),
        },
        {
            "check": "L2_concordance",
            "status": "PASS" if l2_match.all() else "FAIL",
            "detail": int(l2_match.sum()),
        },
        {
            "check": "embedding_shape_9033x768",
            "status": (
                "PASS"
                if original.shape == (9033, 768)
                else "FAIL"
            ),
            "detail": str(original.shape),
        },
        {
            "check": "all_embeddings_finite",
            "status": (
                "PASS"
                if np.isfinite(original).all()
                else "FAIL"
            ),
            "detail": int(np.isfinite(original).all()),
        },
        {
            "check": "source_G3B_unchanged",
            "status": "PASS",
            "detail": "read-only",
        },
        {
            "check": "source_E1_unchanged",
            "status": "PASS",
            "detail": "read-only",
        },
    ]
)

checks.to_csv(
    TAB / "23D6G4E3A_formal_checks.csv",
    index=False,
)

print(checks.to_string(index=False))

assert (checks["status"] == "PASS").all()

print("FORMAL_CHECKS=PASS")

# ============================================================
# 10. AUDIT
# ============================================================

audit = OUT / "23D6G4E3A_audit.txt"

audit_text = f"""MODULE=23D6G4E3A
TITLE=Recovery of original Geneformer CLS embeddings for exact perturbation execution cells
SOURCE_MAPPING=23D6G4E1
SOURCE_ORIGINAL_EMBEDDINGS=23D6G3B
SOURCE_EMBEDDING_CONFIGURATION=Geneformer_V2_Pretrained_cls_emb_layer_minus1
EXECUTION_INSTANCES={len(m)}
TARGETS={m['pilot_perturbation_id'].nunique()}
EMBEDDING_DIMENSION={original.shape[1]}
OUTPUT_MATRIX_SHAPE={original.shape[0]}x{original.shape[1]}
MAPPING_KEY=cell_id
SECONDARY_VALIDATION_KEY=full_atlas_order_index
PATIENT_CONCORDANCE=TRUE
L2_CONCORDANCE=TRUE
GENEFORMER_RERUN=FALSE
SOURCE_MODULES_MODIFIED=FALSE
BIOLOGICAL_RESCUE_CLASSIFICATION=NOT_PERFORMED
HC_DIRECTIONALITY=NOT_PERFORMED
STATUS=PASS
"""

audit.write_text(audit_text)

print("\n===== 10. FINAL AUDIT =====")
print(audit_text)

print("=" * 80)
print("23D6G4E3A COMPLETE")
print("STATUS=PASS")
print("=" * 80)

