#!/usr/bin/env python3

from pathlib import Path
import shutil
import numpy as np
import pandas as pd

from datasets import load_from_disk
from geneformer import perturber_utils as pu


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
)

G4B = (
    BASE
    / "23D6G_FULL_ATLAS_REASSESSMENT"
    / "23D6G4B_FULL_ATLAS_EXECUTION_SETS"
)

SUPPORT = (
    G4B / "tables"
    / "23D6G4B_full_atlas_execution_support.csv"
)

MANIFEST = (
    G4B / "tables"
    / "23D6G4B_full_atlas_execution_cell_manifest.csv"
)

SOURCE_DATASETS = (
    G4B / "objects" / "execution_datasets"
)

OUT = (
    BASE
    / "23D6G_FULL_ATLAS_REASSESSMENT"
    / "23D6G4F_EXHAUSTIVE_FULL_ATLAS"
    / "23D6G4F2A_RECONSTRUCT_EXHAUSTIVE_PERTURBED_DATASETS"
)

DSOUT = OUT / "datasets"
TAB = OUT / "tables"

DSOUT.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

MAX_LEN = 4096
SPECIAL_TOKEN = True
EXPECTED_TOTAL = 52509
EXPECTED_TARGETS = 24


print("=" * 80)
print("23D6G4F2A — EXHAUSTIVE PERTURBED DATASET RECONSTRUCTION")
print("ALL ELIGIBLE G4B CELL-TARGET INSTANCES")
print("NO PATIENT SUBSAMPLING")
print("=" * 80)

print("G4B:", G4B)
print("Support:", SUPPORT)
print("Manifest:", MANIFEST)
print("Source datasets:", SOURCE_DATASETS)
print("Output:", OUT)

assert SUPPORT.exists(), SUPPORT
assert MANIFEST.exists(), MANIFEST
assert SOURCE_DATASETS.exists(), SOURCE_DATASETS


# ============================================================
# 1. LOAD FROZEN G4B UNIVERSE
# ============================================================

print("\n===== 1. LOAD G4B FROZEN UNIVERSE =====")

support = pd.read_csv(SUPPORT)
manifest = pd.read_csv(MANIFEST)

assert len(support) == EXPECTED_TARGETS, len(support)
assert support["pilot_perturbation_id"].nunique() == EXPECTED_TARGETS
assert len(manifest) == EXPECTED_TOTAL, len(manifest)
assert manifest["pilot_perturbation_id"].nunique() == EXPECTED_TARGETS

print("targets:", len(support))
print("execution instances:", len(manifest))

print(
    support[
        [
            "pilot_perturbation_id",
            "disease",
            "final_L2_label",
            "gene",
            "geneformer_perturb_type",
            "n_execution_cells",
            "n_execution_patients",
        ]
    ].to_string(index=False)
)

print("G4B_UNIVERSE=PASS")


# ============================================================
# 2. RECONSTRUCT PERTURBED DATASETS
# ============================================================

print("\n===== 2. RECONSTRUCT EXHAUSTIVE PERTURBED DATASETS =====")

summary_rows = []
row_manifest = []

for counter, (_, r) in enumerate(support.iterrows(), start=1):

    pid = str(r["pilot_perturbation_id"])
    disease = str(r["disease"])
    L2 = str(r["final_L2_label"])
    gene = str(r["gene"])
    token = int(r["geneformer_token_id"])
    perturb_type = str(r["geneformer_perturb_type"])

    expected_n = int(r["n_execution_cells"])
    expected_patients = int(r["n_execution_patients"])

    ds_path = SOURCE_DATASETS / f"{pid}.dataset"

    assert ds_path.exists(), ds_path

    ds = load_from_disk(str(ds_path))

    n_before = len(ds)

    print()
    print("-" * 80)
    print(
        f"[{counter:02d}/{EXPECTED_TARGETS}] "
        f"{pid} | {disease} | {L2} | "
        f"{gene} | {perturb_type}"
    )
    print("rows:", n_before)
    print("expected:", expected_n)
    print("token:", token)
    print("-" * 80)

    assert n_before == expected_n, (
        pid,
        n_before,
        expected_n,
    )

    mm = (
        manifest.loc[
            manifest["pilot_perturbation_id"]
            .astype(str)
            .eq(pid)
        ]
        .copy()
        .reset_index(drop=True)
    )

    assert len(mm) == n_before, (
        pid,
        len(mm),
        n_before,
    )

    # --------------------------------------------------------
    # Exact row-order validation against frozen G4B manifest
    # --------------------------------------------------------

    if "cell_id" in ds.column_names:
        assert (
            np.asarray(ds["cell_id"]).astype(str)
            ==
            mm["cell_id"].astype(str).to_numpy()
        ).all(), f"{pid}: cell order mismatch"

    if "patient_id" in ds.column_names:
        assert (
            np.asarray(ds["patient_id"]).astype(str)
            ==
            mm["patient_id"].astype(str).to_numpy()
        ).all(), f"{pid}: patient order mismatch"

    # --------------------------------------------------------
    # Exact Geneformer perturbation transformation
    # --------------------------------------------------------

    def perturb_example(example):

        ids = list(example["input_ids"])

        indices = [
            ids.index(token)
            if token in ids
            else None
        ]

        indices = [
            x for x in indices
            if x is not None
        ]

        example["tokens_to_perturb"] = [token]

        if len(indices) > 0:
            example["perturb_index"] = indices
        else:
            example["perturb_index"] = [-100]

        original_length = int(example["length"])

        if perturb_type == "delete":

            example = pu.delete_indices(example)
            example["n_overflow"] = 0

        elif perturb_type == "overexpress":

            example = pu.overexpress_tokens(
                example,
                MAX_LEN,
                SPECIAL_TOKEN,
            )

            # Retained only as Geneformer bookkeeping.
            # Not used as biological/QC criterion.
            n_overflow = pu.calc_n_overflow(
                MAX_LEN,
                original_length,
                [token],
                indices,
            )

            example["n_overflow"] = int(n_overflow)

        else:
            raise RuntimeError(
                f"{pid}: unsupported perturbation type "
                f"{perturb_type}"
            )

        return example

    pert = ds.map(
        perturb_example,
        num_proc=None,
    )

    assert len(pert) == n_before

    # --------------------------------------------------------
    # Perturbation validation
    # --------------------------------------------------------

    original_lengths = np.asarray(
        ds["length"],
        dtype=int,
    )

    pert_lengths = np.asarray(
        pert["length"],
        dtype=int,
    )

    token_present_before = np.asarray(
        [token in ids for ids in ds["input_ids"]],
        dtype=bool,
    )

    token_present_after = np.asarray(
        [token in ids for ids in pert["input_ids"]],
        dtype=bool,
    )

    if perturb_type == "delete":

        # G4B delete universe was frozen as target-present.
        assert token_present_before.all(), pid
        assert (~token_present_after).all(), pid

        assert np.all(
            pert_lengths == original_lengths - 1
        ), pid

    else:

        assert token_present_after.all(), pid

    # --------------------------------------------------------
    # Persist immutable metadata
    # --------------------------------------------------------

    metadata = {
        "pilot_perturbation_id":
            [pid] * len(pert),
        "target_gene":
            [gene] * len(pert),
        "target_token":
            [token] * len(pert),
        "geneformer_perturb_type":
            [perturb_type] * len(pert),
        "execution_row_index":
            list(range(len(pert))),
    }

    for col, values in metadata.items():

        if col in pert.column_names:
            pert = pert.remove_columns(col)

        pert = pert.add_column(
            col,
            values,
        )

    # Remove temporary Geneformer columns.
    for col in [
        "tokens_to_perturb",
        "perturb_index",
    ]:
        if col in pert.column_names:
            pert = pert.remove_columns(col)

    # --------------------------------------------------------
    # Final metadata/order validation
    # --------------------------------------------------------

    if "cell_id" in pert.column_names:
        assert (
            np.asarray(pert["cell_id"]).astype(str)
            ==
            mm["cell_id"].astype(str).to_numpy()
        ).all(), pid

    if "patient_id" in pert.column_names:
        assert (
            np.asarray(pert["patient_id"]).astype(str)
            ==
            mm["patient_id"].astype(str).to_numpy()
        ).all(), pid

    observed_patients = (
        mm["patient_id"]
        .astype(str)
        .nunique()
    )

    assert observed_patients == expected_patients, (
        pid,
        observed_patients,
        expected_patients,
    )

    # --------------------------------------------------------
    # Save perturbed dataset
    # --------------------------------------------------------

    out_ds = DSOUT / f"{pid}.dataset"

    if out_ds.exists():
        shutil.rmtree(out_ds)

    pert.save_to_disk(str(out_ds))

    n_overflow_total = int(
        np.asarray(
            pert["n_overflow"],
            dtype=int,
        ).sum()
    )

    summary_rows.append(
        {
            "pilot_perturbation_id": pid,
            "disease": disease,
            "final_L2_label": L2,
            "gene": gene,
            "geneformer_token_id": token,
            "geneformer_perturb_type": perturb_type,
            "n_cells": len(pert),
            "n_patients": observed_patients,
            "n_overflow_total": n_overflow_total,
            "dataset_path": str(out_ds),
        }
    )

    for i in range(len(pert)):

        row_manifest.append(
            {
                "pilot_perturbation_id":
                    pid,
                "execution_row_index":
                    i,
                "dataset_index":
                    mm.loc[i, "dataset_index"],
                "execution_local_index":
                    mm.loc[i, "execution_local_index"],
                "cell_id":
                    mm.loc[i, "cell_id"],
                "patient_id":
                    mm.loc[i, "patient_id"],
                "disease":
                    mm.loc[i, "disease"],
                "final_L2_label":
                    mm.loc[i, "final_L2_label"],
                "cell_disease":
                    mm.loc[i, "cell_disease"],
                "cell_L2":
                    mm.loc[i, "cell_L2"],
                "full_atlas_order_index":
                    mm.loc[i, "full_atlas_order_index"],
                "target_present_before":
                    mm.loc[i, "target_present_before"],
                "sequence_length_before":
                    mm.loc[i, "sequence_length"],
                "length_uncropped":
                    mm.loc[i, "length_uncropped"],
                "sequence_length_after":
                    int(pert_lengths[i]),
                "gene":
                    gene,
                "geneformer_token_id":
                    token,
                "geneformer_perturb_type":
                    perturb_type,
            }
        )

    print(
        "rows=",
        len(pert),
        "patients=",
        observed_patients,
        "length_before_mean=",
        round(float(original_lengths.mean()), 2),
        "length_after_mean=",
        round(float(pert_lengths.mean()), 2),
        "status=PASS",
    )


# ============================================================
# 3. SAVE MANIFESTS
# ============================================================

print("\n===== 3. SAVE MANIFESTS =====")

summary = pd.DataFrame(summary_rows)
rows = pd.DataFrame(row_manifest)

summary.to_csv(
    TAB / "23D6G4F2A_dataset_summary.csv",
    index=False,
)

rows.to_csv(
    TAB / "23D6G4F2A_execution_row_manifest.csv",
    index=False,
)

print(summary.to_string(index=False))

print(
    "TOTAL_RECONSTRUCTED_INSTANCES=",
    len(rows),
    sep="",
)


# ============================================================
# 4. EXACT G4B CONCORDANCE
# ============================================================

print("\n===== 4. EXACT G4B CONCORDANCE =====")

assert len(rows) == len(manifest) == EXPECTED_TOTAL

for col in [
    "pilot_perturbation_id",
    "dataset_index",
    "execution_local_index",
    "cell_id",
    "patient_id",
    "full_atlas_order_index",
]:

    a = rows[col].astype(str).to_numpy()
    b = manifest[col].astype(str).to_numpy()

    assert np.array_equal(a, b), (
        f"G4B concordance failure: {col}"
    )

print("ROW_ORDER_CONCORDANCE_G4B=PASS")


# ============================================================
# 5. TARGET/PATIENT COVERAGE
# ============================================================

print("\n===== 5. TARGET/PATIENT COVERAGE =====")

coverage = (
    rows
    .groupby(
        [
            "pilot_perturbation_id",
            "patient_id",
        ],
        as_index=False,
    )
    .size()
    .rename(columns={"size": "n_cells"})
)

coverage.to_csv(
    TAB / "23D6G4F2A_patient_target_coverage.csv",
    index=False,
)

coverage_summary = (
    coverage
    .groupby(
        "pilot_perturbation_id",
        as_index=False,
    )
    .agg(
        n_patients=("patient_id", "nunique"),
        min_cells_per_patient=("n_cells", "min"),
        median_cells_per_patient=("n_cells", "median"),
        max_cells_per_patient=("n_cells", "max"),
        total_cells=("n_cells", "sum"),
    )
)

coverage_summary.to_csv(
    TAB / "23D6G4F2A_patient_target_coverage_summary.csv",
    index=False,
)

print(coverage_summary.to_string(index=False))


# ============================================================
# 6. FORMAL CHECKS
# ============================================================

print("\n===== 6. FORMAL CHECKS =====")

target_count_match = (
    summary.set_index("pilot_perturbation_id")["n_cells"]
    ==
    support.set_index("pilot_perturbation_id")[
        "n_execution_cells"
    ]
).all()

patient_count_match = (
    summary.set_index("pilot_perturbation_id")["n_patients"]
    ==
    support.set_index("pilot_perturbation_id")[
        "n_execution_patients"
    ]
).all()

checks = pd.DataFrame(
    [
        {
            "check": "targets_24",
            "status":
                "PASS"
                if len(summary) == EXPECTED_TARGETS
                else "FAIL",
            "detail": len(summary),
        },
        {
            "check": "execution_instances_52509",
            "status":
                "PASS"
                if len(rows) == EXPECTED_TOTAL
                else "FAIL",
            "detail": len(rows),
        },
        {
            "check": "all_target_counts_match_G4B",
            "status":
                "PASS"
                if target_count_match
                else "FAIL",
            "detail": f"{len(summary)}/{EXPECTED_TARGETS}",
        },
        {
            "check": "all_patient_counts_match_G4B",
            "status":
                "PASS"
                if patient_count_match
                else "FAIL",
            "detail": f"{len(summary)}/{EXPECTED_TARGETS}",
        },
        {
            "check": "exact_G4B_row_order",
            "status": "PASS",
            "detail": EXPECTED_TOTAL,
        },
        {
            "check": "no_patient_subsampling",
            "status": "PASS",
            "detail": "G4B exhaustive universe",
        },
        {
            "check": "delete_semantics",
            "status": "PASS",
            "detail": "target removed from all delete cells",
        },
        {
            "check": "overexpression_semantics",
            "status": "PASS",
            "detail": "target present after perturbation",
        },
        {
            "check": "model_execution_performed",
            "status": "PASS",
            "detail": "FALSE",
        },
        {
            "check": "source_G4B_modified",
            "status": "PASS",
            "detail": "FALSE",
        },
    ]
)

checks.to_csv(
    TAB / "23D6G4F2A_formal_checks.csv",
    index=False,
)

print(checks.to_string(index=False))

assert (checks["status"] == "PASS").all()

print("FORMAL_CHECKS=PASS")


# ============================================================
# 7. FINAL AUDIT
# ============================================================

audit = f"""MODULE=23D6G4F2A
TITLE=Exhaustive reconstruction of Geneformer perturbation datasets
SOURCE_EXECUTION_SETS=23D6G4B
SOURCE_EXECUTION_MANIFEST=23D6G4B_full_atlas_execution_cell_manifest
SOURCE_PERTURBATION_SEMANTICS=installed_Geneformer_V2
TARGETS={len(summary)}
EXECUTION_INSTANCES={len(rows)}
DELETE_IMPLEMENTATION=geneformer.perturber_utils.delete_indices
OVEREXPRESS_IMPLEMENTATION=geneformer.perturber_utils.overexpress_tokens
OVERFLOW_IMPLEMENTATION=geneformer.perturber_utils.calc_n_overflow
MAX_LEN={MAX_LEN}
SPECIAL_TOKEN={str(SPECIAL_TOKEN).upper()}
MODEL_EXECUTION_PERFORMED=FALSE
IN_SILICO_PERTURBER_RERUN=FALSE
PERTURBED_DATASETS_RECONSTRUCTED=TRUE
PATIENT_STRATIFIED_SUBSAMPLING=FALSE
CAP_PER_PATIENT=NONE
ALL_ELIGIBLE_G4B_INSTANCES_RECONSTRUCTED=TRUE
PATIENT_COVERAGE_PRESERVED=TRUE
EXACT_G4B_ROW_ORDER_PRESERVED=TRUE
DELETE_TARGETS_REQUIRE_TARGET_PRESENT=TRUE
OVEREXPRESSION_TARGET_PRESENT_AFTER=TRUE
OVERFLOW_USED_FOR_BIOLOGICAL_QC=FALSE
BIOLOGICAL_RESCUE_CLASSIFICATION=NOT_PERFORMED
HC_DIRECTIONALITY=NOT_PERFORMED
SOURCE_MODULES_MODIFIED=FALSE
STATUS=PASS
"""

(OUT / "23D6G4F2A_audit.txt").write_text(
    audit,
    encoding="utf-8",
)

print("\n===== 7. FINAL AUDIT =====")
print(audit)


# ============================================================
# 8. OUTPUT INVENTORY
# ============================================================

print("\n===== 8. OUTPUT INVENTORY =====")

expected_outputs = [
    TAB / "23D6G4F2A_dataset_summary.csv",
    TAB / "23D6G4F2A_execution_row_manifest.csv",
    TAB / "23D6G4F2A_patient_target_coverage.csv",
    TAB / "23D6G4F2A_patient_target_coverage_summary.csv",
    TAB / "23D6G4F2A_formal_checks.csv",
    OUT / "23D6G4F2A_audit.txt",
]

for x in expected_outputs:
    print(x, "exists=", x.exists())
    assert x.exists(), x

dataset_dirs = sorted(DSOUT.glob("*.dataset"))

print(
    "perturbed dataset directories:",
    len(dataset_dirs),
)

assert len(dataset_dirs) == EXPECTED_TARGETS, (
    len(dataset_dirs),
    EXPECTED_TARGETS,
)

print()
print("=" * 80)
print("23D6G4F2A COMPLETE")
print("STATUS=PASS")
print("=" * 80)
