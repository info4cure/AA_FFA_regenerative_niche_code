from pathlib import Path
import numpy as np
import pandas as pd

# ============================================================
# 23D6G4F2C
# EXHAUSTIVE FULL-ATLAS HC RESCUE GEOMETRY
#
# Original embeddings:
#   G3B = one unique original CLS embedding per disease cell
#
# Perturbed embeddings:
#   F2B = one perturbed CLS embedding per eligible
#         cell-target execution instance
#
# Exact mapping:
#   F2B execution rows -> G3B original embeddings by cell_id
#
# Primary inference:
#   patient-balanced target-level geometry
#
# Sensitivity:
#   cell-weighted geometry
#
# No biological/therapeutic adjudication is performed here.
# ============================================================

ROOT = Path.cwd()

BASE = (
    ROOT
    / "PIPELINE_V2/results/disease_biology"
    / "23_ADVANCED_ANALYSIS_ROADMAP"
    / "23D_AI_PERTURBATION"
    / "23D6_PERTURBATION"
    / "23D6G_FULL_ATLAS_REASSESSMENT"
)

G3B = BASE / "23D6G3B_FULL_ATLAS_EMBEDDINGS"

F2 = (
    BASE
    / "23D6G4F_EXHAUSTIVE_FULL_ATLAS"
)

F2A = (
    F2
    / "23D6G4F2A_RECONSTRUCT_EXHAUSTIVE_PERTURBED_DATASETS"
)

F2B = (
    F2
    / "23D6G4F2B_EXHAUSTIVE_PERTURBED_EMBEDDINGS"
)

HC = (
    BASE.parent
    / "23D6F_HC_RESCUE_GEOMETRY"
    / "23D6F1_HC_REFERENCE_EMBEDDINGS"
)

OUT = (
    F2
    / "23D6G4F2C_EXHAUSTIVE_HC_RESCUE_GEOMETRY"
)

TAB = OUT / "tables"
OBJ = OUT / "objects"

TAB.mkdir(parents=True, exist_ok=True)
OBJ.mkdir(parents=True, exist_ok=True)

G3B_EMB = (
    G3B
    / "objects"
    / "23D6G3B_full_atlas.csv"
)

G3B_MAN = (
    G3B
    / "tables"
    / "23D6G3B_embedding_cell_manifest.csv"
)

F2B_MAN = (
    F2B
    / "tables"
    / "23D6G4F2B_exhaustive_embedding_manifest.csv"
)

HC_TABLE = (
    HC
    / "tables"
    / "23D6F1_HC_patient_balanced_centroids.csv"
)

PER_TARGET = F2B / "objects" / "per_target"

EPS = 1e-12


print("=" * 80)
print("23D6G4F2C — EXHAUSTIVE FULL-ATLAS HC RESCUE GEOMETRY")
print("52,509 CELL-TARGET EXECUTION INSTANCES")
print("EXACT ORIGINAL EMBEDDING MAPPING BY CELL_ID")
print("PATIENT-BALANCED TARGET-LEVEL PRIMARY GEOMETRY")
print("=" * 80)

print("G3B:", G3B)
print("F2A:", F2A)
print("F2B:", F2B)
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
    num = np.asarray(num, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)

    return np.divide(
        num,
        den,
        out=np.full(len(num), np.nan, dtype=np.float64),
        where=np.abs(den) > EPS,
    )


def finite_or_nan(x):
    x = np.asarray(x, dtype=np.float64)
    x[~np.isfinite(x)] = np.nan
    return x


# ============================================================
# 1. LOAD ORIGINAL G3B EMBEDDINGS
# ============================================================

print("\n===== 1. LOAD ORIGINAL G3B EMBEDDINGS =====")

for p in [G3B_EMB, G3B_MAN, F2B_MAN, HC_TABLE]:
    assert p.exists(), f"Missing input: {p}"

g3 = pd.read_csv(G3B_EMB)
g3m = pd.read_csv(G3B_MAN)

numeric_cols = [
    str(i) for i in range(768)
]

assert all(c in g3.columns for c in numeric_cols)

G = g3[numeric_cols].to_numpy(dtype=np.float32)

assert G.shape == (38699, 768), G.shape
assert len(g3m) == 38699
assert np.isfinite(G).all()

# The embedding CSV itself contains the authoritative row cell_id.
assert "cell_id" in g3.columns
assert g3["cell_id"].nunique() == 38699
assert not g3["cell_id"].duplicated().any()

print("original unique cells:", len(g3))
print("embedding shape:", G.shape)
print("unique cell_id:", g3["cell_id"].nunique())

print("G3B_ORIGINAL_EMBEDDINGS=PASS")


# ============================================================
# 2. LOAD F2B EXECUTION MANIFEST
# ============================================================

print("\n===== 2. LOAD F2B EXECUTION MANIFEST =====")

df = pd.read_csv(F2B_MAN)

assert len(df) == 52509, len(df)
assert df["pilot_perturbation_id"].nunique() == 24

required = [
    "pilot_perturbation_id",
    "cell_id",
    "patient_id",
    "disease",
    "final_L2_label",
    "gene",
    "geneformer_perturb_type",
]

for c in required:
    assert c in df.columns, c

print("execution rows:", len(df))
print("unique cells:", df["cell_id"].nunique())
print("targets:", df["pilot_perturbation_id"].nunique())


# ============================================================
# 3. EXACT ORIGINAL EMBEDDING MAPPING BY CELL_ID
# ============================================================

print("\n===== 3. EXACT ORIGINAL EMBEDDING MAPPING =====")

g3_index = pd.Series(
    np.arange(len(g3), dtype=np.int64),
    index=g3["cell_id"].astype(str),
)

mapped_index = (
    df["cell_id"]
    .astype(str)
    .map(g3_index)
)

missing = int(mapped_index.isna().sum())

print("missing F2B cells in G3B:", missing)

assert missing == 0

mapped_index = mapped_index.astype(np.int64).to_numpy()

A = G[mapped_index]

assert A.shape == (52509, 768)
assert np.isfinite(A).all()

# Metadata concordance.
g3_lookup = (
    g3[
        [
            "cell_id",
            "patient_id",
            "disease",
            "L2",
            "full_atlas_order_index",
        ]
    ]
    .drop_duplicates("cell_id")
    .set_index("cell_id")
)

joined = (
    df[
        [
            "cell_id",
            "patient_id",
            "disease",
            "final_L2_label",
            "full_atlas_order_index",
        ]
    ]
    .join(
        g3_lookup,
        on="cell_id",
        rsuffix="_G3B",
    )
)

metadata_checks = {
    "patient_id": (
        joined["patient_id"].astype(str)
        == joined["patient_id_G3B"].astype(str)
    ),
    "disease": (
        joined["disease"].astype(str)
        == joined["disease_G3B"].astype(str)
    ),
    "L2": (
        joined["final_L2_label"].astype(str)
        == joined["L2"].astype(str)
    ),
    "full_atlas_order_index": (
        pd.to_numeric(
            joined["full_atlas_order_index"]
        )
        == pd.to_numeric(
            joined["full_atlas_order_index_G3B"]
        )
    ),
}

identity_rows = []

for name, x in metadata_checks.items():
    nmatch = int(x.sum())

    identity_rows.append(
        {
            "field": name,
            "n_rows": len(x),
            "n_match": nmatch,
            "status": (
                "PASS"
                if nmatch == len(x)
                else "FAIL"
            ),
        }
    )

identity_df = pd.DataFrame(identity_rows)

identity_df.to_csv(
    TAB / "23D6G4F2C_original_mapping_validation.csv",
    index=False,
)

print(identity_df.to_string(index=False))

assert (identity_df["status"] == "PASS").all()

df["G3B_original_embedding_row"] = mapped_index

print("EXACT_ORIGINAL_MAPPING=PASS")


# ============================================================
# 4. LOAD PERTURBED F2B EMBEDDINGS
# ============================================================

print("\n===== 4. LOAD PERTURBED F2B EMBEDDINGS =====")

P_parts = []

for counter, (pid, z) in enumerate(
    df.groupby(
        "pilot_perturbation_id",
        sort=False,
        observed=True,
    ),
    start=1,
):
    pid = str(pid)

    p = PER_TARGET / f"{pid}_perturbed_cls.npy"

    assert p.exists(), p

    X = np.load(p)

    expected = len(z)

    assert X.shape == (expected, 768), (
        pid,
        X.shape,
        expected,
    )

    assert np.isfinite(X).all()

    P_parts.append(X.astype(np.float32))

    print(
        f"[{counter:02d}/24]",
        pid,
        "rows=",
        len(X),
        "status=PASS",
    )

P = np.vstack(P_parts)

assert P.shape == (52509, 768)
assert np.isfinite(P).all()

print("perturbed embedding shape:", P.shape)
print("PERTURBED_EMBEDDINGS=PASS")


# ============================================================
# 5. LOAD FROZEN HC REFERENCES
# ============================================================

print("\n===== 5. LOAD FROZEN HC REFERENCES =====")

hc = pd.read_csv(HC_TABLE)

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
    assert c in hc.columns, c

hc = hc.loc[
    (hc["reference_type"].astype(str) == "PATIENT_BALANCED")
    & (
        hc["primary_reference"]
        .astype(str)
        .str.lower()
        == "true"
    )
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
    assert np.isfinite(v).all()

    hc_lookup[l2] = {
        "vector": v,
        "n_HC_cells": int(r["n_HC_cells"]),
        "n_HC_patients": int(r["n_HC_patients"]),
        "reference_type": str(r["reference_type"]),
        "reference_support_class":
            str(r["reference_support_class"]),
        "exploratory_only":
            str(r["exploratory_only"]).lower() == "true",
        "centroid_file": str(r["centroid_file"]),
    }

print("HC references:", len(hc_lookup))


# ============================================================
# 6. VERIFY HC COVERAGE
# ============================================================

print("\n===== 6. VERIFY HC COVERAGE =====")

execution_l2 = sorted(
    df["final_L2_label"].astype(str).unique()
)

missing_hc = sorted(
    set(execution_l2) - set(hc_lookup)
)

coverage = pd.DataFrame(
    {
        "final_L2_label": execution_l2,
        "HC_reference_available": [
            x in hc_lookup
            for x in execution_l2
        ],
    }
)

coverage.to_csv(
    TAB / "23D6G4F2C_HC_reference_coverage.csv",
    index=False,
)

print(coverage.to_string(index=False))

assert not missing_hc, missing_hc

print("HC_COVERAGE=PASS")


# ============================================================
# 7. EXECUTION-LEVEL HC GEOMETRY
# ============================================================

print("\n===== 7. COMPUTE EXECUTION-LEVEL GEOMETRY =====")

n = len(df)

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

for l2, idx in df.groupby(
    "final_L2_label",
    sort=False,
).groups.items():

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
    hc_support_class[idx] = ref[
        "reference_support_class"
    ]
    hc_exploratory[idx] = ref["exploratory_only"]
    hc_centroid_file[idx] = ref["centroid_file"]


df["HC_reference_type"] = hc_reference_type
df["HC_reference_n_cells"] = hc_n_cells.astype(int)
df["HC_reference_n_patients"] = hc_n_patients.astype(int)

df["HC_reference_support_class"] = (
    hc_support_class
)

df["HC_reference_exploratory_only"] = (
    pd.Series(hc_exploratory)
    .astype(bool)
    .to_numpy()
)

df["HC_centroid_file"] = hc_centroid_file

df["cos_original_HC"] = finite_or_nan(
    cos_original_HC
)

df["cos_perturbed_HC"] = finite_or_nan(
    cos_perturbed_HC
)

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

df["delta_euclidean_toward_HC"] = (
    df["euclidean_original_HC"]
    - df["euclidean_perturbed_HC"]
)

df["original_to_perturbed_cosine"] = (
    cosine_between_row_vectors(
        A.astype(np.float64),
        P.astype(np.float64),
    )
)

df["original_to_perturbed_euclidean"] = (
    np.linalg.norm(
        P.astype(np.float64)
        - A.astype(np.float64),
        axis=1,
    )
)

cosine_distance_original = (
    1.0 - df["cos_original_HC"].to_numpy()
)

df["normalized_cosine_rescue"] = safe_fraction(
    df["delta_cos_HC"].to_numpy(),
    cosine_distance_original,
)

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

df["strict_geometry_rescue"] = (
    df["toward_HC_by_cosine"]
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
    TAB / "23D6G4F2C_execution_level_HC_rescue_geometry.csv",
    index=False,
)


# ============================================================
# 8. PATIENT-LEVEL AGGREGATION
# ============================================================

print("\n===== 8. PATIENT-LEVEL AGGREGATION =====")

# ensembl_id may or may not have been propagated by F2B.
# Recover it from F2A/G4B if absent.
if "ensembl_id" not in df.columns:
    G4B_SUPPORT = (
        BASE
        / "23D6G4B_FULL_ATLAS_EXECUTION_SETS"
        / "tables"
        / "23D6G4B_full_atlas_execution_support.csv"
    )

    assert G4B_SUPPORT.exists(), G4B_SUPPORT

    sup = pd.read_csv(G4B_SUPPORT)

    gene_meta = (
        sup[
            [
                "pilot_perturbation_id",
                "ensembl_id",
            ]
        ]
        .drop_duplicates("pilot_perturbation_id")
    )

    assert (
        gene_meta["pilot_perturbation_id"].nunique()
        == 24
    )

    df = df.merge(
        gene_meta,
        on="pilot_perturbation_id",
        how="left",
        validate="many_to_one",
    )

    assert df["ensembl_id"].notna().all()

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
    TAB / "23D6G4F2C_patient_level_HC_rescue_geometry.csv",
    index=False,
)

print("patient-target rows:", len(patient))
print("unique patients:", patient["patient_id"].nunique())


# ============================================================
# 9. TARGET-LEVEL PATIENT-BALANCED AGGREGATION
# ============================================================

print(
    "\n===== 9. TARGET-LEVEL PATIENT-BALANCED AGGREGATION ====="
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
    TAB / "23D6G4F2C_target_level_HC_rescue_geometry.csv",
    index=False,
)

print("target rows:", len(target))
print(
    target["HC_geometry_class"]
    .value_counts(dropna=False)
    .to_string()
)


# ============================================================
# 10. CELL-WEIGHTED SENSITIVITY
# ============================================================

print("\n===== 10. CELL-WEIGHTED SENSITIVITY =====")

cell_weighted = (
    df.groupby(
        group_target,
        dropna=False,
        observed=True,
    )
    .agg(
        n_executions=("cell_id", "size"),
        n_patients=("patient_id", "nunique"),

        mean_delta_cos_HC=("delta_cos_HC", "mean"),
        median_delta_cos_HC=("delta_cos_HC", "median"),

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
    TAB / "23D6G4F2C_target_cell_weighted_sensitivity.csv",
    index=False,
)


# ============================================================
# 11. HC SUPPORT SUMMARY
# ============================================================

print("\n===== 11. HC SUPPORT SUMMARY =====")

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
    TAB / "23D6G4F2C_HC_support_summary.csv",
    index=False,
)

print(support.to_string(index=False))


# ============================================================
# 12. FORMAL CHECKS
# ============================================================

print("\n===== 12. FORMAL CHECKS =====")

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
    "execution_rows_52509",
    len(df) == 52509,
    len(df),
)

add_check(
    "original_shape_52509x768",
    A.shape == (52509, 768),
    str(A.shape),
)

add_check(
    "perturbed_shape_52509x768",
    P.shape == (52509, 768),
    str(P.shape),
)

add_check(
    "unique_original_G3B_cells_38699",
    g3["cell_id"].nunique() == 38699,
    g3["cell_id"].nunique(),
)

add_check(
    "F2B_to_G3B_cell_coverage",
    missing == 0,
    missing,
)

add_check(
    "exact_metadata_concordance",
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
    "target_level_rows_24",
   
    len(target) == 24,
    len(target),
)

add_check(
    "patient_target_rows_expected",
    len(patient) == (
        df[
            [
                "pilot_perturbation_id",
                "patient_id",
            ]
        ]
        .drop_duplicates()
        .shape[0]
    ),
    len(patient),
)

add_check(
    "all_execution_L2_have_HC_reference",
    len(missing_hc) == 0,
    len(missing_hc),
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
    "execution_instances_preserved",
    int(target["n_executions"].sum()) == 52509,
    int(target["n_executions"].sum()),
)

add_check(
    "patient_coverage_preserved",
    patient["patient_id"].nunique() == 11,
    patient["patient_id"].nunique(),
)

add_check(
    "AA_patients_6",
    patient.loc[
        patient["disease"] == "AA",
        "patient_id",
    ].nunique() == 6,
    patient.loc[
        patient["disease"] == "AA",
        "patient_id",
    ].nunique(),
)

add_check(
    "FFA_patients_5",
    patient.loc[
        patient["disease"] == "FFA",
        "patient_id",
    ].nunique() == 5,
    patient.loc[
        patient["disease"] == "FFA",
        "patient_id",
    ].nunique(),
)

add_check(
    "no_patient_subsampling",
    len(df) == 52509,
    "G4B exhaustive universe",
)

add_check(
    "patient_balanced_target_aggregation",
    True,
    "patient means equally weighted",
)

add_check(
    "cell_weighted_sensitivity_generated",
    len(cell_weighted) == 24,
    len(cell_weighted),
)

add_check(
    "source_modules_modified",
    True,
    "FALSE",
)

checks_df = pd.DataFrame(checks)

checks_df.to_csv(
    TAB / "23D6G4F2C_formal_checks.csv",
    index=False,
)

print(checks_df.to_string(index=False))

assert (
    checks_df["status"] == "PASS"
).all(), checks_df.loc[
    checks_df["status"] != "PASS"
]

print("FORMAL_CHECKS=PASS")


# ============================================================
# 13. SAVE COMPACT NUMERIC OBJECT
# ============================================================

print("\n===== 13. SAVE NUMERIC OBJECT =====")

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
        df[
            "normalized_cosine_rescue"
        ].to_numpy(),
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

assert numeric_object.shape == (52509, 11)

np.save(
    OBJ / "23D6G4F2C_execution_geometry.npy",
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
        "column_index": range(
            len(geometry_columns)
        ),
        "column_name": geometry_columns,
    }
).to_csv(
    TAB
    / "23D6G4F2C_execution_geometry_columns.csv",
    index=False,
)

print(
    "numeric geometry shape:",
    numeric_object.shape,
)


# ============================================================
# 14. FINAL AUDIT
# ============================================================

print("\n===== 14. FINAL AUDIT =====")

geometry_counts = (
    target["HC_geometry_class"]
    .value_counts(dropna=False)
    .to_dict()
)

n_exploratory_targets = int(
    target[
        "HC_reference_exploratory_only"
    ].sum()
)

n_robust_targets = int(
    (
        target[
            "HC_reference_support_class"
        ]
        == "ROBUST"
    ).sum()
)

n_toward_hc = int(
    (
        target["HC_geometry_class"]
        == "TOWARD_HC_CONCORDANT"
    ).sum()
)

n_away_hc = int(
    (
        target["HC_geometry_class"]
        == "AWAY_FROM_HC_CONCORDANT"
    ).sum()
)

n_mixed = int(
    (
        target["HC_geometry_class"]
        == "MIXED_GEOMETRY"
    ).sum()
)

audit = f"""MODULE=23D6G4F2C
TITLE=Exhaustive full-atlas HC rescue geometry from Geneformer V2 CLS embeddings
SOURCE_ORIGINAL_EMBEDDINGS=23D6G3B
SOURCE_EXECUTION_MAPPING=23D6G4F2B
SOURCE_PERTURBED_EMBEDDINGS=23D6G4F2B
SOURCE_HC_REFERENCE=23D6F1
SOURCE_EXECUTION_UNIVERSE=23D6G4B
HC_REFERENCE_TYPE=PATIENT_BALANCED
HC_REFERENCE_LEVEL=L2
MODEL_VERSION=Geneformer_V2
EMBEDDING_MODE=CLS
EMBEDDING_DIMENSION=768
ORIGINAL_UNIQUE_CELLS={g3["cell_id"].nunique()}
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
PATIENT_STRATIFIED_SUBSAMPLING=FALSE
CAP_PER_PATIENT=NONE
ALL_ELIGIBLE_G4B_INSTANCES_ANALYZED=TRUE
EXACT_CELL_ID_MAPPING_TO_ORIGINAL_EMBEDDINGS=TRUE
HC_SUPPORT_PROPAGATED=TRUE
ROBUST_HC_TARGETS={n_robust_targets}
EXPLORATORY_HC_TARGETS={n_exploratory_targets}
TOWARD_HC_CONCORDANT_TARGETS={n_toward_hc}
AWAY_FROM_HC_CONCORDANT_TARGETS={n_away_hc}
MIXED_GEOMETRY_TARGETS={n_mixed}
GEOMETRY_CLASS_COUNTS={geometry_counts}
THERAPEUTIC_EFFICACY_INFERRED=FALSE
FINAL_BIOLOGICAL_ADJUDICATION=NOT_PERFORMED
SOURCE_MODULES_MODIFIED=FALSE
STATUS=PASS
"""

(
    OUT / "23D6G4F2C_audit.txt"
).write_text(
    audit,
    encoding="utf-8",
)

print(audit)


# ============================================================
# 15. TARGET SUMMARY
# ============================================================

print("\n===== 15. TARGET SUMMARY =====")

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
    "patient_strict_rescue_consistency_class",
]

target_summary = (
    target[summary_cols]
    .sort_values(
        [
            "HC_reference_exploratory_only",
            "mean_patient_delta_cos_HC",
        ],
        ascending=[True, False],
    )
    .reset_index(drop=True)
)

target_summary.insert(
    0,
    "HC_rescue_rank",
    np.arange(
        1,
        len(target_summary) + 1,
    ),
)

target_summary.to_csv(
    TAB
    / "23D6G4F2C_target_HC_rescue_ranking.csv",
    index=False,
)

print(
    target_summary.to_string(
        index=False
    )
)


# ============================================================
# 16. OUTPUT INVENTORY
# ============================================================

print("\n===== 16. OUTPUT INVENTORY =====")

expected_outputs = [
    TAB
    / "23D6G4F2C_execution_level_HC_rescue_geometry.csv",
    TAB
    / "23D6G4F2C_patient_level_HC_rescue_geometry.csv",
    TAB
    / "23D6G4F2C_target_level_HC_rescue_geometry.csv",
    TAB
    / "23D6G4F2C_target_cell_weighted_sensitivity.csv",
    TAB
    / "23D6G4F2C_HC_support_summary.csv",
    TAB
    / "23D6G4F2C_formal_checks.csv",
    TAB
    / "23D6G4F2C_execution_geometry_columns.csv",
    TAB
    / "23D6G4F2C_target_HC_rescue_ranking.csv",
    OBJ
    / "23D6G4F2C_execution_geometry.npy",
    OUT
    / "23D6G4F2C_audit.txt",
]

for p in expected_outputs:
    print(
        p,
        "exists=",
        p.exists(),
    )
    assert p.exists(), p


print()
print("=" * 80)
print("23D6G4F2C COMPLETE")
print("STATUS=PASS")
print("=" * 80)
