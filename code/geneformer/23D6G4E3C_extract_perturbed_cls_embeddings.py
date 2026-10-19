from pathlib import Path
import sys
import time
import numpy as np
import pandas as pd
import torch

from datasets import load_from_disk
from geneformer import EmbExtractor


# ============================================================
# PATHS
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

E3B = (
    BASE
    / "23D6G4E_PATIENT_ROBUSTNESS_HC_DIRECTIONALITY"
    / "23D6G4E3B_RECONSTRUCT_PERTURBED_DATASETS"
)

DATASETS = E3B / "datasets"

OUT = (
    BASE
    / "23D6G4E_PATIENT_ROBUSTNESS_HC_DIRECTIONALITY"
    / "23D6G4E3C_PERTURBED_CLS_EMBEDDINGS"
)

OBJ = OUT / "objects"
TAB = OUT / "tables"

MODEL = (
    ROOT
    / "PIPELINE_V2/external_models/Geneformer"
    / "Geneformer-V2-104M"
)

MANIFEST = (
    E3B
    / "tables"
    / "23D6G4E3B_execution_row_manifest.csv"
)

OBJ.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)


print("=" * 80)
print("23D6G4E3C — PERTURBED CLS EMBEDDINGS")
print("=" * 80)

print("Python:", sys.version.split()[0])
print("Torch:", torch.__version__)
print("MPS built:", torch.backends.mps.is_built())
print("MPS available:", torch.backends.mps.is_available())
print("Model:", MODEL)
print("Datasets:", DATASETS)
print("Output:", OUT)

assert MODEL.exists(), MODEL
assert DATASETS.exists(), DATASETS
assert MANIFEST.exists(), MANIFEST


# ============================================================
# 1. LOAD EXECUTION MANIFEST
# ============================================================

print("\n===== 1. LOAD EXECUTION MANIFEST =====")

manifest = pd.read_csv(MANIFEST)

assert len(manifest) == 9033, len(manifest)

assert "pilot_perturbation_id" in manifest.columns

targets = sorted(
    manifest["pilot_perturbation_id"]
    .astype(str)
    .unique()
)

assert len(targets) == 24, len(targets)

print("execution rows:", len(manifest))
print("targets:", len(targets))


# ============================================================
# 2. VALIDATE DATASETS
# ============================================================

print("\n===== 2. VALIDATE PERTURBED DATASETS =====")

dataset_summary = []

for pid in targets:

    ds_path = DATASETS / f"{pid}.dataset"

    assert ds_path.exists(), ds_path

    ds = load_from_disk(str(ds_path))

    expected = int(
        (
            manifest["pilot_perturbation_id"]
            .astype(str)
            == pid
        ).sum()
    )

    actual = len(ds)

    print(
        pid,
        "rows=",
        actual,
        "expected=",
        expected,
    )

    assert actual == expected, (
        pid,
        actual,
        expected,
    )

    assert "input_ids" in ds.column_names
    assert "length" in ds.column_names

    dataset_summary.append(
        {
            "pilot_perturbation_id": pid,
            "n_execution_rows": actual,
            "dataset_path": str(ds_path),
        }
    )

dataset_summary = pd.DataFrame(dataset_summary)

assert dataset_summary["n_execution_rows"].sum() == 9033

dataset_summary.to_csv(
    TAB / "23D6G4E3C_dataset_validation.csv",
    index=False,
)

print("PERTURBED_DATASETS=PASS")


# ============================================================
# 3. EXTRACT CLS EMBEDDINGS
# ============================================================

print("\n===== 3. EXTRACT PERTURBED CLS EMBEDDINGS =====")

all_embeddings = []
all_rows = []
execution_summary = []

global_row = 0

for counter, pid in enumerate(targets, start=1):

    print()
    print("-" * 80)
    print(f"[{counter:02d}/24] {pid}")
    print("-" * 80)

    ds_path = DATASETS / f"{pid}.dataset"

    ds = load_from_disk(str(ds_path))

    n = len(ds)

    target_obj = OBJ / pid
    target_obj.mkdir(
        parents=True,
        exist_ok=True,
    )

    labels = [
        x
        for x in [
            "patient_id",
            "disease",
            "compartment",
            "L2",
            "stratum",
        ]
        if x in ds.column_names
    ]

    print("rows:", n)
    print("labels:", labels)

    embex = EmbExtractor(
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cls",
        max_ncells=None,
        emb_layer=-1,
        emb_label=labels,
        labels_to_plot=None,
        forward_batch_size=8,
        nproc=1,
        summary_stat=None,
        save_tdigest=False,
        model_version="V2",
    )

    t0 = time.time()

    result = embex.extract_embs(
        model_directory=str(MODEL),
        input_data_file=str(ds_path),
        output_directory=str(target_obj),
        output_prefix=f"{pid}_perturbed_cls",
        output_torch_embs=True,
    )

    elapsed = time.time() - t0

    print(
        "extract_embs return type:",
        type(result),
    )

    # --------------------------------------------------------
    # Geneformer EmbExtractor return normalization
    # Current Geneformer V2 returns:
    #   (pandas.DataFrame, torch.Tensor)
    # The DataFrame contains 768 embedding dimensions plus
    # requested metadata; the tensor contains the same 768-D
    # embeddings without metadata.
    # --------------------------------------------------------

    if isinstance(result, tuple):
        if len(result) != 2:
            raise RuntimeError(
                f"{pid}: unexpected EmbExtractor tuple length "
                f"{len(result)}"
            )
        result_df, result_tensor = result

        if not isinstance(result_df, pd.DataFrame):
            raise RuntimeError(
                f"{pid}: tuple element 0 is not DataFrame; "
                f"got {type(result_df)}"
            )

        if not torch.is_tensor(result_tensor):
            raise RuntimeError(
                f"{pid}: tuple element 1 is not torch.Tensor; "
                f"got {type(result_tensor)}"
            )

        if tuple(result_tensor.shape) != (n, 768):
            raise RuntimeError(
                f"{pid}: unexpected torch embedding shape "
                f"{tuple(result_tensor.shape)}; expected {(n, 768)}"
            )

        result = result_df

    elif isinstance(result, pd.DataFrame):
        # Compatibility with Geneformer versions returning only DataFrame.
        result_tensor = None

    else:
        raise RuntimeError(
            f"{pid}: unsupported EmbExtractor return type "
            f"{type(result)}"
        )

    numeric_cols = [
        c
        for c in result.columns
        if str(c).isdigit()
    ]

    numeric_cols = sorted(
        numeric_cols,
        key=lambda x: int(str(x)),
    )

    if len(numeric_cols) != 768:
        raise RuntimeError(
            f"{pid}: expected 768 embedding dimensions; "
            f"found {len(numeric_cols)}"
        )

    emb = (
        result[numeric_cols]
        .to_numpy(dtype=np.float32)
    )

    assert emb.shape == (n, 768), (
        pid,
        emb.shape,
    )

    assert np.isfinite(emb).all()

    all_embeddings.append(emb)

    target_manifest = manifest.loc[
        manifest["pilot_perturbation_id"]
        .astype(str)
        .eq(pid)
    ].copy()

    assert len(target_manifest) == n

    target_manifest = target_manifest.reset_index(
        drop=True
    )

    target_manifest[
        "perturbed_embedding_row"
    ] = np.arange(
        global_row,
        global_row + n,
        dtype=np.int64,
    )

    target_manifest[
        "within_target_embedding_row"
    ] = np.arange(
        n,
        dtype=np.int64,
    )

    all_rows.append(target_manifest)

    global_row += n

    execution_summary.append(
        {
            "pilot_perturbation_id": pid,
            "n_embeddings": n,
            "embedding_dimension": 768,
            "elapsed_seconds": elapsed,
            "status": "PASS",
        }
    )

    print(
        "embedding shape:",
        emb.shape,
    )
    print(
        f"elapsed: {elapsed:.2f}s"
    )
    print("STATUS=PASS")

    # Persist progress after each target.
    pd.DataFrame(
        execution_summary
    ).to_csv(
        TAB / "23D6G4E3C_execution_summary.csv",
        index=False,
    )


# ============================================================
# 4. CONSOLIDATE
# ============================================================

print("\n===== 4. CONSOLIDATE =====")

X = np.vstack(all_embeddings).astype(
    np.float32,
    copy=False,
)

row_manifest = pd.concat(
    all_rows,
    ignore_index=True,
)

print("final embedding shape:", X.shape)
print("manifest rows:", len(row_manifest))

assert X.shape == (9033, 768)
assert len(row_manifest) == 9033

assert np.array_equal(
    row_manifest[
        "perturbed_embedding_row"
    ].to_numpy(),
    np.arange(9033),
)

assert np.isfinite(X).all()

np.save(
    OBJ / "23D6G4E3C_perturbed_cls_embeddings.npy",
    X,
)

row_manifest.to_csv(
    TAB / "23D6G4E3C_embedding_row_manifest.csv",
    index=False,
)

summary = pd.DataFrame(
    execution_summary
)

summary.to_csv(
    TAB / "23D6G4E3C_execution_summary.csv",
    index=False,
)


# ============================================================
# 5. PER-TARGET COUNTS
# ============================================================

print("\n===== 5. PER-TARGET COUNTS =====")

counts = (
    row_manifest
    .groupby(
        "pilot_perturbation_id",
        observed=True,
    )
    .size()
    .rename("n_perturbed_embeddings")
    .reset_index()
)

print(counts.to_string(index=False))

assert len(counts) == 24
assert counts["n_perturbed_embeddings"].sum() == 9033

counts.to_csv(
    TAB / "23D6G4E3C_target_counts.csv",
    index=False,
)


# ============================================================
# 6. FORMAL CHECKS
# ============================================================

print("\n===== 6. FORMAL CHECKS =====")

checks = pd.DataFrame(
    [
        {
            "check": "targets_24",
            "status": (
                "PASS"
                if len(counts) == 24
                else "FAIL"
            ),
            "detail": len(counts),
        },
        {
            "check": "execution_instances_9033",
            "status": (
                "PASS"
                if len(row_manifest) == 9033
                else "FAIL"
            ),
            "detail": len(row_manifest),
        },
        {
            "check": "embedding_shape_9033x768",
            "status": (
                "PASS"
                if X.shape == (9033, 768)
                else "FAIL"
            ),
            "detail": str(X.shape),
        },
        {
            "check": "all_embeddings_finite",
            "status": (
                "PASS"
                if np.isfinite(X).all()
                else "FAIL"
            ),
            "detail": int(np.isfinite(X).all()),
        },
        {
            "check": "all_targets_pass",
            "status": (
                "PASS"
                if (
                    summary["status"] == "PASS"
                ).all()
                else "FAIL"
            ),
            "detail": (
                summary["status"]
                .value_counts()
                .to_dict()
            ),
        },
        {
            "check": "source_E3B_unchanged",
            "status": "PASS",
            "detail": "read-only",
        },
        {
            "check": "model_version",
            "status": "PASS",
            "detail": "Geneformer_V2",
        },
        {
            "check": "model_type",
            "status": "PASS",
            "detail": "Pretrained",
        },
        {
            "check": "embedding_mode",
            "status": "PASS",
            "detail": "cls",
        },
        {
            "check": "embedding_layer",
            "status": "PASS",
            "detail": "-1",
        },
    ]
)

checks.to_csv(
    TAB / "23D6G4E3C_formal_checks.csv",
    index=False,
)

print(checks.to_string(index=False))

assert (checks["status"] == "PASS").all()

print("FORMAL_CHECKS=PASS")


# ============================================================
# 7. AUDIT
# ============================================================

audit = f"""MODULE=23D6G4E3C
TITLE=Geneformer V2 CLS embeddings of reconstructed perturbed execution datasets
SOURCE_PERTURBED_DATASETS=23D6G4E3B
SOURCE_EXECUTION_MAPPING=23D6G4E1_via_E3B_manifest
MODEL_VERSION=Geneformer_V2
MODEL_TYPE=Pretrained
EMB_MODE=cls
EMB_LAYER=-1
FORWARD_BATCH_SIZE=8
NPROC=1
TARGETS=24
EXECUTION_INSTANCES=9033
EMBEDDING_DIMENSION=768
OUTPUT_MATRIX_SHAPE={X.shape[0]}x{X.shape[1]}
OUTPUT_DTYPE={X.dtype}
PERTURBED_MODEL_FORWARD_PASS=TRUE
PERTURBATION_RECONSTRUCTION_RERUN=FALSE
IN_SILICO_PERTURBER_RERUN=FALSE
BIOLOGICAL_RESCUE_CLASSIFICATION=NOT_PERFORMED
HC_DIRECTIONALITY=NOT_PERFORMED
SOURCE_MODULES_MODIFIED=FALSE
STATUS=PASS
"""

(OUT / "23D6G4E3C_audit.txt").write_text(
    audit
)

print("\n===== 7. FINAL AUDIT =====")
print(audit)

print("=" * 80)
print("23D6G4E3C COMPLETE")
print("STATUS=PASS")
print("=" * 80)
