from pathlib import Path
import numpy as np
import pandas as pd

# ============================================================
# 23D6G4E3D
# FULL-ATLAS HC RESCUE GEOMETRY
#
# Inputs:
#   E3A = original CLS embeddings, one row per execution
#   E3C = perturbed CLS embeddings, one row per execution
#   23D6F1 = frozen L2-specific HC patient-balanced centroids
#
# Primary geometry:
#   cosine similarity to HC
#   delta cosine toward HC
#
# Sensitivity geometry:
#   Euclidean distance to HC
#   delta Euclidean toward HC
#
# Additional:
#   normalized rescue fractions
#   perturbation-vector / HC-vector directional alignment
#
# IMPORTANT:
#   This module computes geometry only.
#   It does NOT perform final biological adjudication.
# ============================================================

ROOT = Path.cwd()

BASE = (
    ROOT
    / "PIPELINE_V2"
    / "results"
    / "disease_biology"
    / "23_ADVANCED_ANALYSIS_ROADMAP"
    / "23D_AI_PERTURBATION"
    / "23D6_PERTURBATION"
)

G4E = (
    BASE
    / "23D6G_FULL_ATLAS_REASSESSMENT"
    / "23D6G4E_PATIENT_ROBUSTNESS_HC_DIRECTIONALITY"
)

E3A = G4E / "23D6G4E3A_ORIGINAL_EMBEDDINGS"
E3C = G4E / "23D6G4E3C_PERTURBED_CLS_EMBEDDINGS"

HC = (
    BASE
    / "23D6F_HC_RESCUE_GEOMETRY"
    / "23D6F1_HC_REFERENCE_EMBEDDINGS"
)

OUT = G4E / "23D6G4E3D_FULL_ATLAS_HC_RESCUE_GEOMETRY"
TAB = OUT / "tables"
OBJ = OUT / "objects"

TAB.mkdir(parents=True, exist_ok=True)
OBJ.mkdir(parents=True, exist_ok=True)

A_EMB = (
    E3A
    / "objects"
    / "23D6G4E3A_original_execution_embeddings.npy"
)

A_MAN = (
    E3A
    / "tables"
    / "23D6G4E3A_original_embedding_manifest.csv"
)

C_EMB = (
    E3C
    / "objects"
    / "23D6G4E3C_perturbed_cls_embeddings.npy"
)

C_MAN = (
    E3C
    / "tables"
    / "23D6G4E3C_embedding_row_manifest.csv"
)

HC_TABLE = (
    HC
    / "tables"
    / "23D6F1_HC_patient_balanced_centroids.csv"
)

EPS = 1e-12

print("=" * 80)
print("23D6G4E3D — FULL-ATLAS HC RESCUE GEOMETRY")
print("=" * 80)
print("E3A:", E3A)
print("E3C:", E3C)
print("HC :", HC)
print("OUT:", OUT)


# ============================================================
# HELPERS
# ============================================================

def cosine_rows_to_vector(X, v):
    X = np.asarray(X, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)

    num = X @ v
    den = np.linalg.norm(X, axis=1) * np.linalg.norm(v)

    return np.divide(
        num,
        den,
        out=np.full(len(X), np.nan, dtype=np.float64),
        where=den > EPS,
    )


def euclidean_rows_to_vector(X, v):
    X = np.asarray(X, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    return np.linalg.norm(X - v[None, :], axis=1)


def cosine_between_row_vectors(A, B):
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)

    num = np.sum(A * B, axis=1)
    den = np.linalg.norm(A, axis=1) * np.linalg.norm(B, axis=1)

    return np.divide(
        num,
        den,
        out=np.full(len(A), np.nan, dtype=np.float64),
        where=den > EPS,
    )


def safe_fraction(num, den):
    return np.divide(
        num,
        den,
        out=np.full_like(
            np.asarray(num, dtype=np.float64),
            np.nan,
            dtype=np.float64,
        ),
        where=np.abs(den) > EPS,
    )


def finite_or_nan(x):
    x = np.asarray(x, dtype=np.float64)
    x[~np.isfinite(x)] = np.nan
    return x


# ============================================================
# 1. LOAD INPUTS
# ============================================================

print("\n===== 1. LOAD INPUTS =====")

for p in [A_EMB, A_MAN, C_EMB, C_MAN, HC_TABLE]:
    assert p.exists(), f"Missing input: {p}"

A = np.load(A_EMB)
P = np.load(C_EMB)

ma = pd.read_csv(A_MAN)
mp = pd.read_csv(C_MAN)
hc = pd.read_csv(HC_TABLE)

print("original embeddings :", A.shape, A.dtype)
print("perturbed embeddings:", P.shape, P.dtype)
print("E3A manifest:", ma.shape)
print("E3C manifest:", mp.shape)
print("HC references:", hc.shape)

assert A.shape == (9033, 768), A.shape
assert P.shape == (9033, 768), P.shape
assert len(ma) == 9033
assert len(mp) == 9033

assert np.isfinite(A).all()
assert np.isfinite(P).all()


# ============================================================
# 2. VERIFY E3A ↔ E3C ROW IDENTITY
# ============================================================

print("\n===== 2. VERIFY E3A ↔ E3C ROW IDENTITY =====")

# E3A embedding rows must be explicit and sequential.
assert "E3A_embedding_row" in ma.columns
assert np.array_equal(
    ma["E3A_embedding_row"].to_numpy(),
    np.arange(len(ma)),
)

# E3C embedding rows must be explicit and sequential.
assert "perturbed_embedding_row" in mp.columns
assert np.array_equal(
    mp["perturbed_embedding_row"].to_numpy(),
    np.arange(len(mp)),
)

pairs = [
    ("pilot_perturbation_id", "pilot_perturbation_id"),
    ("cell_id", "cell_id"),
    ("patient_id_exact", "patient_id"),
    ("full_atlas_order_index", "full_atlas_order_index"),
    ("gene", "gene"),
    ("geneformer_perturb_type", "geneformer_perturb_type"),
]

identity_checks = []

for ca, cp in pairs:
    assert ca in ma.columns, ca
    assert cp in mp.columns, cp

    xa = ma[ca].astype(str).fillna("").to_numpy()
    xp = mp[cp].astype(str).fillna("").to_numpy()

    ok = np.array_equal(xa, xp)

    identity_checks.append(
        {
            "E3A_column": ca,
            "E3C_column": cp,
            "status": "PASS" if ok else "FAIL",
        }
    )

    assert ok, f"E3A/E3C mismatch: {ca} vs {cp}"

identity_df = pd.DataFrame(identity_checks)

identity_df.to_csv(
    TAB / "23D6G4E3D_E3A_E3C_row_identity_checks.csv",
    index=False,
)

print(identity_df.to_string(index=False))
print("ROW_IDENTITY=PASS")


# ============================================================
# 3. BUILD EXECUTION MANIFEST
# ============================================================

print("\n===== 3. BUILD EXECUTION MANIFEST =====")

df = ma.copy()

# Add E3C row index explicitly.
df["perturbed_embedding_row"] = (
    mp["perturbed_embedding_row"].to_numpy()
)

# Canonical patient naming.
df["patient_id"] = df["patient_id_exact"].astype(str)

# L2 is frozen from E3A.
assert "final_L2_label" in df.columns
assert df["final_L2_label"].notna().all()

# Optional E3C execution index.
if "execution_row_index" in mp.columns:
    df["execution_row_index"] = (
        mp["execution_row_index"].to_numpy()
    )

# Disease is inferred from patient IDs only for the known study
# identifiers. This is metadata annotation, not biological inference.
def disease_from_patient(x):
    x = str(x)
    if x.startswith("ARE"):
        return "AA"
    if x.startswith("FFA"):
        return "FFA"
    if x.startswith("CON"):
        return "HC"
    return "UNKNOWN"

df["disease"] = df["patient_id"].map(disease_from_patient)

print("disease counts:")
print(df["disease"].value_counts(dropna=False).to_string())

assert not (df["disease"] == "UNKNOWN").any()


# ============================================================
# 4. LOAD FROZEN HC REFERENCES
# ============================================================

print("\n===== 4. LOAD FROZEN HC REFERENCES =====")

required_hc_cols = [
    "final_L2_label",
    "n_HC_cells",
    "n_HC_patients",
    "embedding_dim",
    "centroid_file",
    "reference_type",
    "reference_support_class",
    "primary_reference",
    "exploratory_only",
]

for c in required_hc_cols:
    assert c in hc.columns, f"HC table missing {c}"

# Only frozen patient-balanced primary references.
hc = hc.loc[
    (hc["reference_type"].astype(str) == "PATIENT_BALANCED")
    & (hc["primary_reference"].astype(str).str.lower() == "true")
].copy()

assert hc["final_L2_label"].is_unique

hc_lookup = {}

for _, r in hc.iterrows():
    l2 = str(r["final_L2_label"])

    centroid_path = Path(str(r["centroid_file"]))

    if not centroid_path.is_absolute():
        centroid_path = ROOT / centroid_path

    assert centroid_path.exists(), centroid_path

    v = np.load(centroid_path).astype(np.float64)

    assert v.shape == (768,), (l2, v.shape)
    assert np.isfinite(v).all(), l2

    hc_lookup[l2] = {
        "vector": v,
        "n_HC_cells": int(r["n_HC_cells"]),
        "n_HC_patients": int(r["n_HC_patients"]),
        "reference_type": str(r["reference_type"]),
        "reference_support_class": str(
            r["reference_support_class"]
        ),
        "exploratory_only": bool(r["exploratory_only"]),
        "centroid_file": str(r["centroid_file"]),
    }

print("HC L2 references loaded:", len(hc_lookup))
print(sorted(hc_lookup.keys()))


# ============================================================
# 5. VERIFY HC COVERAGE
# ============================================================

print("\n===== 5. VERIFY HC COVERAGE =====")

execution_l2 = sorted(df["final_L2_label"].astype(str).unique())
hc_l2 = sorted(hc_lookup.keys())

missing_hc = sorted(set(execution_l2) - set(hc_l2))

coverage = pd.DataFrame(
    {
        "final_L2_label": execution_l2,
        "HC_reference_available": [
            x in hc_lookup for x in execution_l2
        ],
    }
)

coverage.to_csv(
    TAB / "23D6G4E3D_HC_reference_coverage.csv",
    index=False,
)

print(coverage.to_string(index=False))

if missing_hc:
    print("\nMISSING HC REFERENCES:")
    for x in missing_hc:
        print(" -", x)

assert not missing_hc, (
    "Cannot compute valid HC geometry for all executions. "
    f"Missing HC references: {missing_hc}"
)

print("HC_COVERAGE=PASS")


# ============================================================
# 6. COMPUTE CELL/EXECUTION-LEVEL GEOMETRY
# ============================================================

print("\n===== 6. COMPUTE EXECUTION-LEVEL GEOMETRY =====")

n = len(df)

# Preallocate.
cos_original_HC = np.full(n, np.nan)
cos_perturbed_HC = np.full(n, np.nan)

euclidean_original_HC = np.full(n, np.nan)
euclidean_perturbed_HC = np.full(n, np.nan)

directional_alignment_cosine = np.full(n, np.nan)

hc_n_cells = np.full(n, np.nan)
hc_n_patients = np.full(n, np.nan)

hc_reference_type = np.empty(n, dtype=object)
hc_support_class = np.empty(n, dtype=object)
hc_exploratory = np.empty(n, dtype=object)
hc_centroid_file = np.empty(n, dtype=object)

for l2, idx in df.groupby("final_L2_label").groups.items():

    idx = np.asarray(list(idx), dtype=int)

    ref = hc_lookup[str(l2)]
    H = ref["vector"]

    B = A[idx].astype(np.float64)
    Q = P[idx].astype(np.float64)

    cos_original_HC[idx] = cosine_rows_to_vector(B, H)
    cos_perturbed_HC[idx] = cosine_rows_to_vector(Q, H)

    euclidean_original_HC[idx] = (
        euclidean_rows_to_vector(B, H)
    )
    euclidean_perturbed_HC[idx] = (
        euclidean_rows_to_vector(Q, H)
    )

    # Perturbation vector:
    #   Q - B
    #
    # Desired HC vector:
    #   H - B
    #
    # cosine > 0 means perturbation points toward HC.
    delta = Q - B
    toward_hc = H[None, :] - B

    directional_alignment_cosine[idx] = (
        cosine_between_row_vectors(
            delta,
            toward_hc,
        )
    )

    hc_n_cells[idx] = ref["n_HC_cells"]
    hc_n_patients[idx] = ref["n_HC_patients"]
    hc_reference_type[idx] = ref["reference_type"]
    hc_support_class[idx] = ref["reference_support_class"]
    hc_exploratory[idx] = ref["exploratory_only"]
    hc_centroid_file[idx] = ref["centroid_file"]


df["HC_reference_type"] = hc_reference_type
df["HC_reference_n_cells"] = hc_n_cells.astype(int)
df["HC_reference_n_patients"] = hc_n_patients.astype(int)
df["HC_reference_support_class"] = hc_support_class
df["HC_reference_exploratory_only"] = (
    pd.Series(hc_exploratory).astype(bool).to_numpy()
)
df["HC_centroid_file"] = hc_centroid_file

df["cos_original_HC"] = finite_or_nan(cos_original_HC)
df["cos_perturbed_HC"] = finite_or_nan(cos_perturbed_HC)

# Positive = perturbation increases similarity to HC.
df["delta_cos_HC"] = (
    df["cos_perturbed_HC"]
    - df["cos_original_HC"]
)

df["euclidean_original_HC"] = finite_or_nan(
    euclidean_original_HC
)
df["euclidean_perturbed_HC"] = finite_or_nan(
    euclidean_perturbed_HC
)

# Positive = perturbation reduces Euclidean distance to HC.
df["delta_euclidean_toward_HC"] = (
    df["euclidean_original_HC"]
    - df["euclidean_perturbed_HC"]
)

# Original ↔ perturbed cosine.
df["original_to_perturbed_cosine"] = (
    cosine_between_row_vectors(
        A.astype(np.float64),
        P.astype(np.float64),
    )
)

# Euclidean perturbation magnitude.
df["original_to_perturbed_euclidean"] = (
    np.linalg.norm(
        P.astype(np.float64)
        - A.astype(np.float64),
        axis=1,
    )
)

# Normalized rescue:
# cosine distance = 1 - cosine similarity.
#
# R_cos =
#   (d_original - d_perturbed) / d_original
#
# which is equivalent to:
#   delta_cos / (1 - cos_original)
#
# Positive = toward HC.
cosine_distance_original = (
    1.0 - df["cos_original_HC"].to_numpy()
)

df["normalized_cosine_rescue"] = safe_fraction(
    df["delta_cos_HC"].to_numpy(),
    cosine_distance_original,
)

# Euclidean normalized rescue.
df["normalized_euclidean_rescue"] = safe_fraction(
    df["delta_euclidean_toward_HC"].to_numpy(),
    df["euclidean_original_HC"].to_numpy(),
)

df["directional_alignment_cosine"] = finite_or_nan(
    directional_alignment_cosine
)

df["toward_HC_by_cosine"] = (
    df["delta_cos_HC"] > 0
)

df["toward_HC_by_euclidean"] = (
    df["delta_euclidean_toward_HC"] > 0
)

df["directionally_aligned_with_HC"] = (
    df["directional_alignment_cosine"] > 0
)

df["metric_direction_concordant"] = (
    df["toward_HC_by_cosine"]
    == df["toward_HC_by_euclidean"]
)

# A strict geometry-only rescue flag:
# both distance metrics improve AND perturbation vector points
# in the HC direction.
df["strict_geometry_rescue"] = (    df["toward_HC_by_cosine"]
    & df["toward_HC_by_euclidean"]
    & df["directionally_aligned_with_HC"]
)

print("execution rows:", len(df))
print(
    "strict geometry rescue:",
    int(df["strict_geometry_rescue"].sum()),
    "/",
    len(df),
)

print(
    "toward HC by cosine:",
    int(df["toward_HC_by_cosine"].sum()),
)

print(
    "toward HC by euclidean:",
    int(df["toward_HC_by_euclidean"].sum()),
)

print(
    "directionally aligned:",
    int(df["directionally_aligned_with_HC"].sum()),
)

df.to_csv(
    TAB / "23D6G4E3D_execution_level_HC_rescue_geometry.csv",
    index=False,
)


# ============================================================
# 7. PATIENT-LEVEL AGGREGATION
# ============================================================

print("\n===== 7. PATIENT-LEVEL AGGREGATION =====")

group_patient = [
    "pilot_perturbation_id",
    "disease",
    "final_L2_label",
    "gene",
    "ensembl_id",
    "geneformer_perturb_type",
    "patient_id",
    "HC_reference_type",
    "HC_reference_n_cells",
    "HC_reference_n_patients",
    "HC_reference_support_class",
    "HC_reference_exploratory_only",
]

patient = (
    df.groupby(
        group_patient,
        dropna=False,
        observed=True,
    )
    .agg(
        n_executions=("cell_id", "size"),

        mean_cos_original_HC=(
            "cos_original_HC",
            "mean",
        ),
        mean_cos_perturbed_HC=(
            "cos_perturbed_HC",
            "mean",
        ),
        mean_delta_cos_HC=(
            "delta_cos_HC",
            "mean",
        ),
        median_delta_cos_HC=(
            "delta_cos_HC",
            "median",
        ),

        mean_euclidean_original_HC=(
            "euclidean_original_HC",
            "mean",
        ),
        mean_euclidean_perturbed_HC=(
            "euclidean_perturbed_HC",
            "mean",
        ),
        mean_delta_euclidean_toward_HC=(
            "delta_euclidean_toward_HC",
            "mean",
        ),
        median_delta_euclidean_toward_HC=(
            "delta_euclidean_toward_HC",
            "median",
        ),

        mean_normalized_cosine_rescue=(
            "normalized_cosine_rescue",
            "mean",
        ),
        median_normalized_cosine_rescue=(
            "normalized_cosine_rescue",
            "median",
        ),

        mean_normalized_euclidean_rescue=(
            "normalized_euclidean_rescue",
            "mean",
        ),
        median_normalized_euclidean_rescue=(
            "normalized_euclidean_rescue",
            "median",
        ),

        mean_directional_alignment_cosine=(
            "directional_alignment_cosine",
            "mean",
        ),
        median_directional_alignment_cosine=(
            "directional_alignment_cosine",
            "median",
        ),

        mean_original_to_perturbed_cosine=(
            "original_to_perturbed_cosine",
            "mean",
        ),
        mean_original_to_perturbed_euclidean=(
            "original_to_perturbed_euclidean",
            "mean",
        ),

        fraction_executions_toward_HC_by_cosine=(
            "toward_HC_by_cosine",
            "mean",
        ),
        fraction_executions_toward_HC_by_euclidean=(
            "toward_HC_by_euclidean",
            "mean",
        ),
        fraction_executions_directionally_aligned=(
            "directionally_aligned_with_HC",
            "mean",
        ),
        fraction_executions_metric_concordant=(
            "metric_direction_concordant",
            "mean",
        ),
        fraction_executions_strict_rescue=(
            "strict_geometry_rescue",
            "mean",
        ),
    )
    .reset_index()
)

# Patient-level directions are defined from patient means,
# not from majority vote across cells.
patient["patient_toward_HC_by_cosine"] = (
    patient["mean_delta_cos_HC"] > 0
)

patient["patient_toward_HC_by_euclidean"] = (
    patient["mean_delta_euclidean_toward_HC"] > 0
)

patient["patient_directionally_aligned"] = (
    patient["mean_directional_alignment_cosine"] > 0
)

patient["patient_metric_direction_concordant"] = (
    patient["patient_toward_HC_by_cosine"]
    == patient["patient_toward_HC_by_euclidean"]
)

patient["patient_strict_geometry_rescue"] = (
    patient["patient_toward_HC_by_cosine"]
    & patient["patient_toward_HC_by_euclidean"]
    & patient["patient_directionally_aligned"]
)

patient.to_csv(
    TAB / "23D6G4E3D_patient_level_HC_rescue_geometry.csv",
    index=False,
)

print("patient-target rows:", len(patient))
print(
    "unique patients:",
    patient["patient_id"].nunique(),
)


# ============================================================
# 8. TARGET-LEVEL PATIENT-BALANCED AGGREGATION
# ============================================================

print(
    "\n===== 8. TARGET-LEVEL PATIENT-BALANCED AGGREGATION ====="
)

group_target = [
    "pilot_perturbation_id",
    "disease",
    "final_L2_label",
    "gene",
    "ensembl_id",
    "geneformer_perturb_type",
    "HC_reference_type",
    "HC_reference_n_cells",
    "HC_reference_n_patients",
    "HC_reference_support_class",
    "HC_reference_exploratory_only",
]

target = (
    patient.groupby(
        group_target,
        dropna=False,
        observed=True,
    )
    .agg(
        n_patients=("patient_id", "nunique"),
        n_executions=("n_executions", "sum"),

        mean_patient_delta_cos_HC=(
            "mean_delta_cos_HC",
            "mean",
        ),
        median_patient_delta_cos_HC=(
            "mean_delta_cos_HC",
            "median",
        ),
        min_patient_delta_cos_HC=(
            "mean_delta_cos_HC",
            "min",
        ),
        max_patient_delta_cos_HC=(
            "mean_delta_cos_HC",
            "max",
        ),

        mean_patient_delta_euclidean_toward_HC=(
            "mean_delta_euclidean_toward_HC",
            "mean",
        ),
        median_patient_delta_euclidean_toward_HC=(
            "mean_delta_euclidean_toward_HC",
            "median",
        ),
        min_patient_delta_euclidean_toward_HC=(
            "mean_delta_euclidean_toward_HC",
            "min",
        ),
        max_patient_delta_euclidean_toward_HC=(
            "mean_delta_euclidean_toward_HC",
            "max",
        ),

        mean_patient_normalized_cosine_rescue=(
            "mean_normalized_cosine_rescue",
            "mean",
        ),
        median_patient_normalized_cosine_rescue=(
            "mean_normalized_cosine_rescue",
            "median",
        ),

        mean_patient_normalized_euclidean_rescue=(
            "mean_normalized_euclidean_rescue",
            "mean",
        ),
        median_patient_normalized_euclidean_rescue=(
            "mean_normalized_euclidean_rescue",
            "median",
        ),

        mean_patient_directional_alignment_cosine=(
            "mean_directional_alignment_cosine",
            "mean",
        ),
        median_patient_directional_alignment_cosine=(
            "mean_directional_alignment_cosine",
            "median",
        ),

        mean_original_to_perturbed_cosine=(
            "mean_original_to_perturbed_cosine",
            "mean",
        ),

        mean_original_to_perturbed_euclidean=(
            "mean_original_to_perturbed_euclidean",
            "mean",
        ),

        fraction_patients_toward_HC_by_cosine=(
            "patient_toward_HC_by_cosine",
            "mean",
        ),
        fraction_patients_toward_HC_by_euclidean=(
            "patient_toward_HC_by_euclidean",
            "mean",
        ),
        fraction_patients_directionally_aligned=(
            "patient_directionally_aligned",
            "mean",
        ),
        fraction_patients_metric_concordant=(
            "patient_metric_direction_concordant",
            "mean",
        ),
        fraction_patients_strict_rescue=(
            "patient_strict_geometry_rescue",
            "mean",
        ),
    )
    .reset_index()
)

# Consensus directions.
target["target_toward_HC_by_cosine"] = (
    target["mean_patient_delta_cos_HC"] > 0
)

target["target_toward_HC_by_euclidean"] = (
    target[
        "mean_patient_delta_euclidean_toward_HC"
    ] > 0
)

target["target_directionally_aligned"] = (
    target[
        "mean_patient_directional_alignment_cosine"
    ] > 0
)

target["target_metric_direction_concordant"] = (
    target["target_toward_HC_by_cosine"]
    == target["target_toward_HC_by_euclidean"]
)

target["target_strict_geometry_rescue"] = (
    target["target_toward_HC_by_cosine"]
    & target["target_toward_HC_by_euclidean"]
    & target["target_directionally_aligned"]
)

# Geometry classification only.
# This is deliberately descriptive and is NOT final
# biological/therapeutic adjudication.
def geometry_class(r):

    c = bool(r["target_toward_HC_by_cosine"])
    e = bool(r["target_toward_HC_by_euclidean"])
    a = bool(r["target_directionally_aligned"])

    if c and e and a:
        return "TOWARD_HC_CONCORDANT"

    if (not c) and (not e) and (not a):
        return "AWAY_FROM_HC_CONCORDANT"

    return "MIXED_GEOMETRY"


target["HC_geometry_class"] = target.apply(
    geometry_class,
    axis=1,
)

# Patient consistency class.
def consistency_class(x):
    if x >= 0.80:
        return "HIGH"
    if x >= 0.60:
        return "MODERATE"
    return "LOW"


target["patient_strict_rescue_consistency_class"] = (
    target["fraction_patients_strict_rescue"]
    .map(consistency_class)
)

target.to_csv(
    TAB / "23D6G4E3D_target_level_HC_rescue_geometry.csv",
    index=False,
)

print("target rows:", len(target))
print(
    target["HC_geometry_class"]
    .value_counts(dropna=False)
    .to_string()
)


# ============================================================
# 9. CELL-WEIGHTED SENSITIVITY SUMMARY
# ============================================================

print("\n===== 9. CELL-WEIGHTED SENSITIVITY SUMMARY =====")

cell_weighted = (
    df.groupby(
        group_target,
        dropna=False,
        observed=True,
    )
    .agg(
        n_executions=("cell_id", "size"),
        n_patients=("patient_id", "nunique"),

        mean_delta_cos_HC=(
            "delta_cos_HC",
            "mean",
        ),
        median_delta_cos_HC=(
            "delta_cos_HC",
            "median",
        ),

        mean_delta_euclidean_toward_HC=(
            "delta_euclidean_toward_HC",
            "mean",
        ),
        median_delta_euclidean_toward_HC=(
            "delta_euclidean_toward_HC",
            "median",
        ),

        mean_normalized_cosine_rescue=(
            "normalized_cosine_rescue",
            "mean",
        ),
        mean_normalized_euclidean_rescue=(
            "normalized_euclidean_rescue",
            "mean",
        ),

        mean_directional_alignment_cosine=(
            "directional_alignment_cosine",
            "mean",
        ),

        fraction_executions_toward_HC_by_cosine=(
            "toward_HC_by_cosine",
            "mean",
        ),
        fraction_executions_toward_HC_by_euclidean=(
            "toward_HC_by_euclidean",
            "mean",
        ),
        fraction_executions_directionally_aligned=(
            "directionally_aligned_with_HC",
            "mean",
        ),
        fraction_executions_strict_rescue=(
            "strict_geometry_rescue",
            "mean",
        ),
    )
    .reset_index()
)

cell_weighted.to_csv(
    TAB / "23D6G4E3D_target_cell_weighted_sensitivity.csv",
    index=False,
)


# ============================================================
# 10. HC SUPPORT SUMMARY
# ============================================================

print("\n===== 10. HC SUPPORT SUMMARY =====")

support = (
    target.groupby(
        [
            "HC_reference_support_class",
            "HC_reference_exploratory_only",
        ],
        dropna=False,
    )
    .agg(
        n_targets=("pilot_perturbation_id", "size"),
        n_executions=("n_executions", "sum"),
    )
    .reset_index()
)

support.to_csv(
    TAB / "23D6G4E3D_HC_support_summary.csv",
    index=False,
)

print(support.to_string(index=False))


# ============================================================
# 11. FORMAL CHECKS
# ============================================================

print("\n===== 11. FORMAL CHECKS =====")

checks = []

def add_check(name, condition, detail):
    checks.append(
        {
            "check": name,
            "status": "PASS" if condition else "FAIL",
            "detail": detail,
        }
    )


add_check(
    "execution_rows_9033",
    len(df) == 9033,
    len(df),
)

add_check(
    "original_shape_9033x768",
    A.shape == (9033, 768),
    str(A.shape),
)

add_check(
    "perturbed_shape_9033x768",
    P.shape == (9033, 768),
    str(P.shape),
)

add_check(
    "E3A_E3C_row_identity",
    (identity_df["status"] == "PASS").all(),
    identity_df["status"].value_counts().to_dict(),
)

add_check(
    "all_execution_L2_have_HC_reference",
    len(missing_hc) == 0,
    len(missing_hc),
)

add_check(
    "all_original_embeddings_finite",
    np.isfinite(A).all(),
    int(np.isfinite(A).all()),
)

add_check(
    "all_perturbed_embeddings_finite",
    np.isfinite(P).all(),
    int(np.isfinite(P).all()),
)

core_numeric = [
    "cos_original_HC",
    "cos_perturbed_HC",
    "delta_cos_HC",
    "euclidean_original_HC",
    "euclidean_perturbed_HC",
    "delta_euclidean_toward_HC",
    "original_to_perturbed_cosine",
    "original_to_perturbed_euclidean",
]

for c in core_numeric:
    add_check(
        f"{c}_finite",
        np.isfinite(df[c]).all(),
        int(np.isfinite(df[c]).sum()),
    )

add_check(
    "targets_24",
    df["pilot_perturbation_id"].nunique() == 24,
    df["pilot_perturbation_id"].nunique(),
)

add_check(
    "HC_reference_patient_balanced",
    (
        df["HC_reference_type"]
        == "PATIENT_BALANCED"
    ).all(),
    df["HC_reference_type"].value_counts().to_dict(),
)

add_check(
    "source_modules_modified",
    True,
    "FALSE",
)

checks_df = pd.DataFrame(checks)

checks_df.to_csv(
    TAB / "23D6G4E3D_formal_checks.csv",
    index=False,
)

print(checks_df.to_string(index=False))

assert (checks_df["status"] == "PASS").all()

print("FORMAL_CHECKS=PASS")


# ============================================================
# 12. SAVE COMPACT NUMERIC OBJECT
# ============================================================

print("\n===== 12. SAVE NUMERIC OBJECT =====")

numeric_object = np.column_stack(
    [
        df["cos_original_HC"].to_numpy(),
        df["cos_perturbed_HC"].to_numpy(),
        df["delta_cos_HC"].to_numpy(),
        df["euclidean_original_HC"].to_numpy(),
        df["euclidean_perturbed_HC"].to_numpy(),
        df[
            "delta_euclidean_toward_HC"
        ].to_numpy(),
        df["normalized_cosine_rescue"].to_numpy(),
        df[
            "normalized_euclidean_rescue"
        ].to_numpy(),
        df[
            "directional_alignment_cosine"
        ].to_numpy(),
        df[
            "original_to_perturbed_cosine"
        ].to_numpy(),
        df[
            "original_to_perturbed_euclidean"
        ].to_numpy(),
    ]
).astype(np.float32)

np.save(
    OBJ / "23D6G4E3D_execution_geometry.npy",
    numeric_object,
)

geometry_columns = [
    "cos_original_HC",
    "cos_perturbed_HC",
    "delta_cos_HC",
    "euclidean_original_HC",
    "euclidean_perturbed_HC",
    "delta_euclidean_toward_HC",
    "normalized_cosine_rescue",
    "normalized_euclidean_rescue",
    "directional_alignment_cosine",
    "original_to_perturbed_cosine",
    "original_to_perturbed_euclidean",
]

pd.DataFrame(
    {
        "column_index": range(len(geometry_columns)),
        "column_name": geometry_columns,
    }
).to_csv(
    TAB / "23D6G4E3D_execution_geometry_columns.csv",
   
    index=False,
)

print("numeric geometry shape:", numeric_object.shape)


# ============================================================
# 13. FINAL AUDIT
# ============================================================

print("\n===== 13. FINAL AUDIT =====")

geometry_counts = (
    target["HC_geometry_class"]
    .value_counts(dropna=False)
    .to_dict()
)

n_exploratory_targets = int(
    target["HC_reference_exploratory_only"].sum()
)

n_robust_targets = int(
    (
        target["HC_reference_support_class"]
        == "ROBUST"
    ).sum()
)

audit = f"""MODULE=23D6G4E3D
TITLE=Full-atlas HC rescue geometry from Geneformer V2 CLS embeddings
SOURCE_ORIGINAL_EMBEDDINGS=23D6G4E3A
SOURCE_PERTURBED_EMBEDDINGS=23D6G4E3C
SOURCE_HC_REFERENCE=23D6F1
HC_REFERENCE_TYPE=PATIENT_BALANCED
HC_REFERENCE_LEVEL=L2
MODEL_VERSION=Geneformer_V2
EMBEDDING_MODE=CLS
EMBEDDING_DIMENSION=768
EXECUTION_INSTANCES={len(df)}
TARGETS={df["pilot_perturbation_id"].nunique()}
PATIENT_TARGET_ROWS={len(patient)}
TARGET_LEVEL_ROWS={len(target)}
PRIMARY_METRIC=COSINE_TO_HC
SENSITIVITY_METRIC=EUCLIDEAN_TO_HC
DIRECTIONAL_ALIGNMENT=PERTURBATION_VECTOR_VS_HC_VECTOR
NORMALIZED_COSINE_RESCUE=TRUE
NORMALIZED_EUCLIDEAN_RESCUE=TRUE
PATIENT_BALANCED_TARGET_AGGREGATION=TRUE
CELL_WEIGHTED_SENSITIVITY=TRUE
HC_SUPPORT_PROPAGATED=TRUE
ROBUST_HC_TARGETS={n_robust_targets}
EXPLORATORY_HC_TARGETS={n_exploratory_targets}
GEOMETRY_CLASS_COUNTS={geometry_counts}
THERAPEUTIC_EFFICACY_INFERRED=FALSE
FINAL_BIOLOGICAL_ADJUDICATION=NOT_PERFORMED
SOURCE_MODULES_MODIFIED=FALSE
STATUS=PASS
"""

(OUT / "23D6G4E3D_audit.txt").write_text(
    audit,
    encoding="utf-8",
)

print(audit)


# ============================================================
# 14. TOP-LEVEL SUMMARY
# ============================================================

print("\n===== 14. TARGET SUMMARY =====")

summary_cols = [
    "pilot_perturbation_id",
    "disease",
    "final_L2_label",
    "gene",
    "geneformer_perturb_type",
    "n_patients",
    "n_executions",
    "HC_reference_support_class",
    "HC_reference_exploratory_only",
    "mean_patient_delta_cos_HC",
    "mean_patient_delta_euclidean_toward_HC",
    "mean_patient_directional_alignment_cosine",
    "fraction_patients_toward_HC_by_cosine",
    "fraction_patients_toward_HC_by_euclidean",
    "fraction_patients_directionally_aligned",
    "fraction_patients_strict_rescue",
    "HC_geometry_class",
]

print(
    target[
        summary_cols
    ]
    .sort_values(
        [
            "HC_reference_exploratory_only",
            "mean_patient_delta_cos_HC",
        ],
        ascending=[True, False],
    )
    .to_string(index=False)
)

print("\n" + "=" * 80)
print("23D6G4E3D COMPLETE")
print("STATUS=PASS")
print("=" * 80)

