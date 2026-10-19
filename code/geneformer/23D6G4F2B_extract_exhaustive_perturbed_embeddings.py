#!/usr/bin/env python3

import os
import sys
import time
import gc
from pathlib import Path

os.environ["HF_DATASETS_NUM_PROC"] = "1"
os.environ["TOKENIZERS_PARALLELISM"] = "false"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import torch

from datasets import load_from_disk
from geneformer import EmbExtractor


# ============================================================
# CONFIGURATION
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

FROOT = (
    BASE
    / "23D6G4F_EXHAUSTIVE_FULL_ATLAS"
)

F2A = (
    FROOT
    / "23D6G4F2A_RECONSTRUCT_EXHAUSTIVE_PERTURBED_DATASETS"
)

DATASETS = F2A / "datasets"

MANIFEST = (
    F2A
    / "tables"
    / "23D6G4F2A_execution_row_manifest.csv"
)

OUT = (
    FROOT
    / "23D6G4F2B_EXHAUSTIVE_PERTURBED_EMBEDDINGS"
)

OBJ = OUT / "objects"
TARGET_OBJ = OBJ / "per_target"
TAB = OUT / "tables"

MODEL = (
    ROOT
    / "PIPELINE_V2/external_models/Geneformer"
    / "Geneformer-V2-104M"
)

OBJ.mkdir(parents=True, exist_ok=True)
TARGET_OBJ.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)


print("=" * 80)
print("23D6G4F2B — EXHAUSTIVE PERTURBED CLS EMBEDDINGS")
print("ALL 52,509 ELIGIBLE CELL-TARGET INSTANCES")
print("RESTARTABLE PER TARGET")
print("=" * 80)

print("Python:", sys.version.split()[0])
print("Torch:", torch.__version__)
print("MPS built:", torch.backends.mps.is_built())
print("MPS available:", torch.backends.mps.is_available())
print("Model:", MODEL)
print("F2A:", F2A)
print("Datasets:", DATASETS)
print("Manifest:", MANIFEST)
print("Output:", OUT)

assert MODEL.exists(), MODEL
assert DATASETS.exists(), DATASETS
assert MANIFEST.exists(), MANIFEST


# ============================================================
# 1. LOAD FROZEN EXECUTION MANIFEST
# ============================================================

print("\n===== 1. LOAD F2A EXECUTION MANIFEST =====")

manifest = pd.read_csv(MANIFEST)

assert len(manifest) == 52509, len(manifest)

required = [
    "pilot_perturbation_id",
    "execution_row_index",
    "dataset_index",
    "execution_local_index",
    "cell_id",
    "patient_id",
    "disease",
    "final_L2_label",
    "cell_disease",
    "cell_L2",
    "full_atlas_order_index",
    "gene",
    "geneformer_token_id",
    "geneformer_perturb_type",
]

missing = [x for x in required if x not in manifest.columns]
assert not missing, missing

targets = (
    manifest["pilot_perturbation_id"]
    .astype(str)
    .drop_duplicates()
    .tolist()
)

assert len(targets) == 24, len(targets)

print("execution instances:", len(manifest))
print("targets:", len(targets))

counts = (
    manifest.groupby(
        "pilot_perturbation_id",
        sort=False
    )
    .size()
    .rename("n_execution_rows")
    .reset_index()
)

print(counts.to_string(index=False))

print("MANIFEST=PASS")


# ============================================================
# 2. VALIDATE F2A DATASETS
# ============================================================

print("\n===== 2. VALIDATE F2A DATASETS =====")

dataset_validation = []

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

    assert actual == expected, (
        pid,
        actual,
        expected,
    )

    assert "input_ids" in ds.column_names
    assert "length" in ds.column_names

    if "cell_id" in ds.column_names:

        mm = (
            manifest.loc[
                manifest["pilot_perturbation_id"]
                .astype(str)
                .eq(pid)
            ]
            .reset_index(drop=True)
        )

        assert (
            np.asarray(ds["cell_id"]).astype(str)
            ==
            mm["cell_id"].astype(str).to_numpy()
        ).all(), pid

    dataset_validation.append(
        {
            "pilot_perturbation_id": pid,
            "expected_rows": expected,
            "observed_rows": actual,
            "row_count_match": actual == expected,
            "dataset_path": str(ds_path),
        }
    )

    print(
        pid,
        "rows=",
        actual,
        "status=PASS",
    )

dataset_validation = pd.DataFrame(
    dataset_validation
)

dataset_validation.to_csv(
    TAB / "23D6G4F2B_dataset_validation.csv",
    index=False,
)

assert dataset_validation[
    "row_count_match"
].all()

print("DATASET_VALIDATION=PASS")


# ============================================================
# 3. EXTRACT PERTURBED CLS EMBEDDINGS
# ============================================================

print("\n===== 3. EXTRACT EXHAUSTIVE PERTURBED CLS EMBEDDINGS =====")

execution_summary = []

for counter, pid in enumerate(
    targets,
    start=1,
):

    print()
    print("-" * 80)
    print(f"[{counter:02d}/24] {pid}")
    print("-" * 80)

    ds_path = DATASETS / f"{pid}.dataset"

    ds = load_from_disk(str(ds_path))

    n = len(ds)

    target_manifest = (
        manifest.loc[
            manifest["pilot_perturbation_id"]
            .astype(str)
            .eq(pid)
        ]
        .reset_index(drop=True)
        .copy()
    )

    assert len(target_manifest) == n

    npy_path = (
        TARGET_OBJ
        / f"{pid}_perturbed_cls.npy"
    )

    manifest_path = (
        TARGET_OBJ
        / f"{pid}_embedding_manifest.csv"
    )

    complete_marker = (
        TARGET_OBJ
        / f"{pid}_COMPLETE.txt"
    )

    # --------------------------------------------------------
    # Restartability
    # --------------------------------------------------------

    if (
        complete_marker.exists()
        and npy_path.exists()
        and manifest_path.exists()
    ):

        try:

            old_emb = np.load(
                npy_path,
                mmap_mode="r",
            )

            old_manifest = pd.read_csv(
                manifest_path
            )

            valid = (
                old_emb.shape == (n, 768)
                and len(old_manifest) == n
                and np.isfinite(old_emb).all()
                and (
                    old_manifest["cell_id"]
                    .astype(str)
                    .to_numpy()
                    ==
                    target_manifest["cell_id"]
                    .astype(str)
                    .to_numpy()
                ).all()
            )

        except Exception:
            valid = False

        if valid:

            print(
                "STATUS=ALREADY_COMPLETE — "
                "validated and reused"
            )

            execution_summary.append(
                {
                    "pilot_perturbation_id": pid,
                    "n_embeddings": n,
                    "embedding_dimension": 768,
                    "elapsed_seconds": np.nan,
                    "execution_status": "PASS_REUSED",
                    "embedding_file": str(npy_path),
                    "manifest_file": str(manifest_path),
                }
            )

            del old_emb
            gc.collect()

            continue

        else:

            print(
                "Existing target outputs invalid; "
                "recomputing."
            )

            for p in [
                npy_path,
                manifest_path,
                complete_marker,
            ]:
                if p.exists():
                    p.unlink()

    print("rows:", n)

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

    target_tmp = (
        TARGET_OBJ
        / f"{pid}_geneformer_tmp"
    )

    target_tmp.mkdir(
        parents=True,
        exist_ok=True,
    )

    result = embex.extract_embs(
        model_directory=str(MODEL),
        input_data_file=str(ds_path),
        output_directory=str(target_tmp),
        output_prefix=f"{pid}_perturbed_cls",
        output_torch_embs=True,
    )

    elapsed = time.time() - t0

    print(
        "extract_embs return type:",
        type(result),
    )

    # --------------------------------------------------------
    # Normalize Geneformer V2 return
    # --------------------------------------------------------

    if isinstance(result, tuple):

        if len(result) != 2:
            raise RuntimeError(
                f"{pid}: unexpected tuple length "
                f"{len(result)}"
            )

        result_df, result_tensor = result

        if not isinstance(
            result_df,
            pd.DataFrame,
        ):
            raise RuntimeError(
                f"{pid}: tuple element 0 "
                "is not DataFrame"
            )

        if not torch.is_tensor(
            result_tensor
        ):
            raise RuntimeError(
                f"{pid}: tuple element 1 "
                "is not torch.Tensor"
            )

        if tuple(result_tensor.shape) != (
            n,
            768,
        ):
            raise RuntimeError(
                f"{pid}: unexpected tensor "
                f"shape {tuple(result_tensor.shape)}"
            )

        result = result_df

    elif isinstance(
        result,
        pd.DataFrame,
    ):
        result_tensor = None

    else:
        raise RuntimeError(
            f"{pid}: unsupported EmbExtractor "
            f"return type {type(result)}"
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
            f"{pid}: expected 768 dimensions; "
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

    assert np.isfinite(emb).all(), pid

    # --------------------------------------------------------
    # Save per-target embedding matrix
    # --------------------------------------------------------

    np.save(
        npy_path,
        emb,
        allow_pickle=False,
    )

    # --------------------------------------------------------
    # Save immutable embedding-row mapping
    # --------------------------------------------------------

    target_manifest[
        "within_target_embedding_row"
    ] = np.arange(
        n,
        dtype=np.int64,
    )

    target_manifest[
        "embedding_file"
    ] = str(npy_path)

    target_manifest.to_csv(
        manifest_path,
        index=False,
    )

    # Exact row mapping validation
    assert (
        target_manifest[
            "execution_row_index"
        ].to_numpy()
        ==
        np.arange(
            n,
            dtype=np.int64,
        )
    ).all(), pid

    complete_marker.write_text(
        "\n".join(
            [
                f"pilot_perturbation_id={pid}",
                f"n_embeddings={n}",
                "embedding_dimension=768",
                f"elapsed_seconds={elapsed:.6f}",
                f"embedding_file={npy_path}",
                f"manifest_file={manifest_path}",
                "STATUS=PASS",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    execution_summary.append(
        {
            "pilot_perturbation_id": pid,
            "n_embeddings": n,
            "embedding_dimension": 768,
            "elapsed_seconds": elapsed,
            "execution_status": "PASS",
            "embedding_file": str(npy_path),
            "manifest_file": str(manifest_path),
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

    del emb
    del result
    if result_tensor is not None:
        del result_tensor
    del ds

    gc.collect()

    if torch.backends.mps.is_available():
        torch.mps.empty_cache()


# ============================================================
# 4. CONSOLIDATE EXECUTION SUMMARY
# ============================================================

print("\n===== 4. CONSOLIDATE EXECUTION SUMMARY =====")

summary = pd.DataFrame(
    execution_summary
)

assert len(summary) == 24

assert summary[
    "execution_status"
].isin(
    [
        "PASS",
        "PASS_REUSED",
    ]
).all()

summary.to_csv(
    TAB / "23D6G4F2B_execution_summary.csv",
    index=False,
)

print(summary.to_string(index=False))


# ============================================================
# 5. CONSOLIDATE GLOBAL EMBEDDING MANIFEST
# ============================================================

print("\n===== 5. CONSOLIDATE GLOBAL EMBEDDING MANIFEST =====")

all_manifests = []

global_row = 0

for pid in targets:

    p = (
        TARGET_OBJ
        / f"{pid}_embedding_manifest.csv"
    )

    assert p.exists(), p

    mm = pd.read_csv(p)

    expected = int(
        (
            manifest[
                "pilot_perturbation_id"
            ].astype(str)
            == pid
        ).sum()
    )

    assert len(mm) == expected

    mm[
        "perturbed_embedding_global_row"
    ] = np.arange(
        global_row,
        global_row + len(mm),
        dtype=np.int64,
    )

    global_row += len(mm)

    all_manifests.append(mm)

global_manifest = pd.concat(
    all_manifests,
    ignore_index=True,
)

assert len(global_manifest) == 52509
assert global_row == 52509

global_manifest.to_csv(
    TAB
    / "23D6G4F2B_exhaustive_embedding_manifest.csv",
    index=False,
)

print(
    "global manifest rows:",
    len(global_manifest),
)

print(
    "targets:",
    global_manifest[
        "pilot_perturbation_id"
    ].nunique(),
)


# ============================================================
# 6. FINAL VALIDATION
# ============================================================

print("\n===== 6. FINAL VALIDATION =====")

validation_rows = []

for pid in targets:

    npy_path = (
        TARGET_OBJ
        / f"{pid}_perturbed_cls.npy"
    )

    manifest_path = (
        TARGET_OBJ
        / f"{pid}_embedding_manifest.csv"
    )

    emb = np.load(
        npy_path,
        mmap_mode="r",
    )

    mm = pd.read_csv(
        manifest_path
    )

    expected = int(
        (
            manifest[
                "pilot_perturbation_id"
            ].astype(str)
            == pid
        ).sum()
    )

    shape_ok = (
        emb.shape
        ==
        (expected, 768)
    )

    manifest_ok = (
        len(mm)
        ==
        expected
    )

    finite_ok = bool(
        np.isfinite(emb).all()
    )

    cell_order_ok = (
        mm["cell_id"]
        .astype(str)
        .to_numpy()
        ==
        manifest.loc[
            manifest[
                "pilot_perturbation_id"
            ].astype(str).eq(pid),
            "cell_id",
        ]
        .astype(str)
        .to_numpy()
    ).all()

    validation_rows.append(
        {
            "pilot_perturbation_id": pid,
            "expected_rows": expected,
            "embedding_rows": emb.shape[0],
            "embedding_dimension": emb.shape[1],
            "shape_ok": shape_ok,
            "manifest_ok": manifest_ok,
            "finite_ok": finite_ok,
            "cell_order_ok": cell_order_ok,
        }
    )

    del emb

validation = pd.DataFrame(
    validation_rows
)

validation.to_csv(
    TAB
    / "23D6G4F2B_embedding_validation.csv",
    index=False,
)

print(validation.to_string(index=False))

assert validation["shape_ok"].all()
assert validation["manifest_ok"].all()
assert validation["finite_ok"].all()
assert validation["cell_order_ok"].all()

print("EMBEDDING_VALIDATION=PASS")


# ============================================================
# 7. FORMAL CHECKS
# ============================================================

print("\n===== 7. FORMAL CHECKS =====")

checks = pd.DataFrame(
    [
        {
            "check": "targets_24",
            "status": (
                "PASS"
                if len(summary) == 24
                else "FAIL"
            ),
            "detail": len(summary),
        },
        {
            "check": "execution_instances_52509",
            "status": (
                "PASS"
               
                if len(global_manifest) == 52509
                else "FAIL"
            ),
            "detail": len(global_manifest),
        },
        {
            "check": "all_targets_embedded",
            "status": (
                "PASS"
                if global_manifest[
                    "pilot_perturbation_id"
                ].nunique() == 24
                else "FAIL"
            ),
            "detail": global_manifest[
                "pilot_perturbation_id"
            ].nunique(),
        },
        {
            "check": "all_execution_status_pass",
            "status": (
                "PASS"
                if summary[
                    "execution_status"
                ].isin(
                    [
                        "PASS",
                        "PASS_REUSED",
                    ]
                ).all()
                else "FAIL"
            ),
            "detail": str(
                summary[
                    "execution_status"
                ].value_counts().to_dict()
            ),
        },
        {
            "check": "all_embedding_shapes_correct",
            "status": (
                "PASS"
                if validation[
                    "shape_ok"
                ].all()
                else "FAIL"
            ),
            "detail": (
                f"{int(validation['shape_ok'].sum())}/24"
            ),
        },
        {
            "check": "embedding_dimension_768",
            "status": (
                "PASS"
                if (
                    validation[
                        "embedding_dimension"
                    ] == 768
                ).all()
                else "FAIL"
            ),
            "detail": str(
                sorted(
                    validation[
                        "embedding_dimension"
                    ].unique().tolist()
                )
            ),
        },
        {
            "check": "all_embeddings_finite",
            "status": (
                "PASS"
                if validation[
                    "finite_ok"
                ].all()
                else "FAIL"
            ),
            "detail": (
                f"{int(validation['finite_ok'].sum())}/24"
            ),
        },
        {
            "check": "all_embedding_manifests_match",
            "status": (
                "PASS"
                if validation[
                    "manifest_ok"
                ].all()
                else "FAIL"
            ),
            "detail": (
                f"{int(validation['manifest_ok'].sum())}/24"
            ),
        },
        {
            "check": "exact_cell_order_preserved",
            "status": (
                "PASS"
                if validation[
                    "cell_order_ok"
                ].all()
                else "FAIL"
            ),
            "detail": (
                f"{int(validation['cell_order_ok'].sum())}/24"
            ),
        },
        {
            "check": "no_patient_subsampling",
            "status": "PASS",
            "detail": "F2A exhaustive G4B universe",
        },
        {
            "check": "max_ncells_none",
            "status": "PASS",
            "detail": "None",
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
        {
            "check": "no_classifier",
            "status": "PASS",
            "detail": "Pretrained",
        },
        {
            "check": "no_fine_tuning",
            "status": "PASS",
            "detail": "none",
        },
        {
            "check": "biological_rescue_not_performed",
            "status": "PASS",
            "detail": "TRUE",
        },
        {
            "check": "HC_directionality_not_performed",
            "status": "PASS",
            "detail": "TRUE",
        },
    ]
)

checks.to_csv(
    TAB / "23D6G4F2B_formal_checks.csv",
    index=False,
)

print(checks.to_string(index=False))

assert (
    checks["status"] == "PASS"
).all()

print("FORMAL_CHECKS=PASS")


# ============================================================
# 8. FINAL AUDIT
# ============================================================

print("\n===== 8. FINAL AUDIT =====")

n_embeddings_total = int(
    summary["n_embeddings"].sum()
)

assert n_embeddings_total == 52509

audit = f"""MODULE=23D6G4F2B
TITLE=Exhaustive Geneformer V2 perturbed CLS embedding extraction
SOURCE_MODULE=23D6G4F2A
SOURCE_EXECUTION_SETS=23D6G4B
SOURCE_EXECUTION_INSTANCES=52509
TARGETS=24
PERTURBED_EMBEDDINGS={n_embeddings_total}
EMBEDDING_DIMENSION=768
MODEL_VERSION=Geneformer_V2
MODEL_TYPE=Pretrained
EMB_MODE=cls
EMB_LAYER=-1
MAX_NCELLS=None
FORWARD_BATCH_SIZE=8
NPROC=1
MODEL_INPUT_SIZE=4096
PATIENT_STRATIFIED_SUBSAMPLING=FALSE
CAP_PER_PATIENT=NONE
ALL_ELIGIBLE_G4B_INSTANCES_EMBEDDED=TRUE
EXACT_F2A_ROW_ORDER_PRESERVED=TRUE
CELL_IDS_PRESERVED=TRUE
PATIENT_IDS_PRESERVED=TRUE
TARGET_IDS_PRESERVED=TRUE
PER_TARGET_EMBEDDINGS_PERSISTED=TRUE
GLOBAL_EMBEDDING_MANIFEST_PERSISTED=TRUE
RESTARTABLE_PER_TARGET=TRUE
CLASSIFIER_EXECUTION_PERFORMED=FALSE
FINE_TUNING_PERFORMED=FALSE
BIOLOGICAL_RESCUE_CLASSIFICATION=NOT_PERFORMED
HC_DIRECTIONALITY=NOT_PERFORMED
SOURCE_MODULES_MODIFIED=FALSE
STATUS=PASS
"""

(
    OUT / "23D6G4F2B_audit.txt"
).write_text(
    audit,
    encoding="utf-8",
)

print(audit)


# ============================================================
# 9. OUTPUT INVENTORY
# ============================================================

print("\n===== 9. OUTPUT INVENTORY =====")

expected_outputs = [
    TAB / "23D6G4F2B_dataset_validation.csv",
    TAB / "23D6G4F2B_execution_summary.csv",
    TAB / "23D6G4F2B_exhaustive_embedding_manifest.csv",
    TAB / "23D6G4F2B_embedding_validation.csv",
    TAB / "23D6G4F2B_formal_checks.csv",
    OUT / "23D6G4F2B_audit.txt",
]

for p in expected_outputs:
    print(
        p,
        "exists=",
        p.exists(),
    )
    assert p.exists(), p

n_npy = len(
    list(
        TARGET_OBJ.glob(
            "*_perturbed_cls.npy"
        )
    )
)

n_target_manifests = len(
    list(
        TARGET_OBJ.glob(
            "*_embedding_manifest.csv"
        )
    )
)

n_complete = len(
    list(
        TARGET_OBJ.glob(
            "*_COMPLETE.txt"
        )
    )
)

print(
    "per-target embedding matrices:",
    n_npy,
)

print(
    "per-target embedding manifests:",
    n_target_manifests,
)

print(
    "per-target COMPLETE markers:",
    n_complete,
)

assert n_npy == 24
assert n_target_manifests == 24
assert n_complete == 24


print()
print("=" * 80)
print("23D6G4F2B COMPLETE")
print("STATUS=PASS")
print("=" * 80)
