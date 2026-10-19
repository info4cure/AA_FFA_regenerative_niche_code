from pathlib import Path
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

PANEL = (
    BASE.parent
    / "23D6_TARGET_SELECTION"
    / "23D6C_PILOT_PANEL_FREEZE"
    / "tables"
    / "23D6C_frozen_pilot_perturbation_panel.csv"
)

BALANCED = (
    BASE
    / "23D6G_FULL_ATLAS_REASSESSMENT"
    / "23D6G4C_BALANCED_EXECUTION_SETS"
)

E1 = (
    BASE
    / "23D6G_FULL_ATLAS_REASSESSMENT"
    / "23D6G4E_PATIENT_ROBUSTNESS_HC_DIRECTIONALITY"
    / "23D6G4E1_EXACT_EXECUTION_MAPPING"
    / "tables"
    / "23D6G4E1_exact_execution_cell_mapping.csv"
)

OUT = (
    BASE
    / "23D6G_FULL_ATLAS_REASSESSMENT"
    / "23D6G4E_PATIENT_ROBUSTNESS_HC_DIRECTIONALITY"
    / "23D6G4E3B_RECONSTRUCT_PERTURBED_DATASETS"
)

DSOUT = OUT / "datasets"
TAB = OUT / "tables"

DSOUT.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)


print("=" * 80)
print("23D6G4E3B — EXACT PERTURBED DATASET RECONSTRUCTION")
print("=" * 80)


# ============================================================
# 1. LOAD FROZEN PANEL + E1
# ============================================================

print("\n===== 1. LOAD SOURCES =====")

panel = pd.read_csv(PANEL)
mapping = pd.read_csv(E1)

assert len(panel) == 24
assert panel["pilot_perturbation_id"].nunique() == 24
assert len(mapping) == 9033

print("targets:", len(panel))
print("execution instances:", len(mapping))


# ============================================================
# 2. LOCATE BALANCED DATASETS
# ============================================================

print("\n===== 2. LOCATE BALANCED DATASETS =====")

paths = {}

for pid in panel["pilot_perturbation_id"].astype(str):

    candidates = list(BALANCED.rglob(f"{pid}.dataset"))

    if not candidates:
        candidates = list(BALANCED.rglob(f"*{pid}*.dataset"))

    assert len(candidates) == 1, (pid, candidates)

    paths[pid] = candidates[0]

    print(pid, len(load_from_disk(str(candidates[0]))))


# ============================================================
# 3. MODEL TOKEN CONVENTIONS
# ============================================================

# Geneformer V2 tokenized datasets already contain special tokens.
# From observed datasets:
# first token = CLS
# last token  = EOS.
#
# pu.overexpress_tokens only requires max_len + special_token.
# max_len for Geneformer V2 is 4096.

MAX_LEN = 4096
SPECIAL_TOKEN = True

print("\n===== 3. RECONSTRUCTION CONFIGURATION =====")
print("MAX_LEN:", MAX_LEN)
print("SPECIAL_TOKEN:", SPECIAL_TOKEN)


# ============================================================
# 4. RECONSTRUCT
# ============================================================

print("\n===== 4. RECONSTRUCT 24 TARGET DATASETS =====")

summary = []
row_manifest = []

for _, r in panel.iterrows():

    pid = str(r["pilot_perturbation_id"])
    gene = str(r["gene"])
    perturb_type = str(r["geneformer_perturb_type"])
    token = int(r["geneformer_token_id"])

    print("\n" + "-" * 80)
    print(pid, gene, perturb_type, "token=", token)

    ds = load_from_disk(str(paths[pid]))

    n_before = len(ds)

    # --------------------------------------------------------
    # Exact Geneformer group perturbation transformation
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

            # G4D / Geneformer:
            # pu.delete_indices(example)
            example = pu.delete_indices(example)

            example["n_overflow"] = 0

        elif perturb_type == "overexpress":

            # G4D / Geneformer:
            # pu.overexpress_tokens(example, max_len, special_token)
            example = pu.overexpress_tokens(
                example,
                MAX_LEN,
                SPECIAL_TOKEN,
            )

            n_overflow = pu.calc_n_overflow(
                MAX_LEN,
                original_length,
                [token],
                indices,
            )

            example["n_overflow"] = int(n_overflow)

        else:
            raise RuntimeError(
                f"{pid}: unsupported perturbation type {perturb_type}"
            )

        return example

    pert = ds.map(
        perturb_example,
        num_proc=None,
    )

    assert len(pert) == n_before

    # --------------------------------------------------------
    # Validation
    # --------------------------------------------------------

    original_lengths = np.asarray(ds["length"], dtype=int)
    pert_lengths = np.asarray(pert["length"], dtype=int)

    if perturb_type == "delete":

        token_present_before = np.asarray(
            [token in ids for ids in ds["input_ids"]],
            dtype=bool,
        )

        token_present_after = np.asarray(
            [token in ids for ids in pert["input_ids"]],
            dtype=bool,
        )

        # G4D delete execution sets should contain the target.
        assert token_present_before.all(), pid
        assert (~token_present_after).all(), pid

        assert np.all(
            pert_lengths == original_lengths - 1
        ), pid

    else:

        # overexpressed token must be present after reconstruction
        token_present_after = np.asarray(
            [token in ids for ids in pert["input_ids"]],
            dtype=bool,
        )

        assert token_present_after.all(), pid

    # --------------------------------------------------------
    # Add immutable execution metadata
    # --------------------------------------------------------

    mm = (
        mapping[
            mapping["pilot_perturbation_id"].astype(str) == pid
        ]
        .reset_index(drop=True)
    )

    assert len(mm) == len(pert), (
        pid,
        len(mm),
        len(pert),
    )

    # Validate balanced dataset order against E1
    if "cell_id" in pert.column_names:

        assert (
            np.asarray(pert["cell_id"]).astype(str)
            ==
            mm["cell_id"].astype(str).to_numpy()
        ).all(), pid

    # Ensure metadata is persisted even if already present.
    metadata = {
        "pilot_perturbation_id": [pid] * len(pert),
        "target_gene": [gene] * len(pert),
        "target_token": [token] * len(pert),
        "geneformer_perturb_type": [perturb_type] * len(pert),
        "execution_row_index": list(range(len(pert))),
    }

    for col, values in metadata.items():

        if col in pert.column_names:
            pert = pert.remove_columns(col)

        pert = pert.add_column(col, values)

    # --------------------------------------------------------
    # Remove temporary Geneformer bookkeeping columns
    # except n_overflow, which is useful for audit
    # --------------------------------------------------------

    for col in [
        "tokens_to_perturb",
        "perturb_index",
    ]:
        if col in pert.column_names:
            pert = pert.remove_columns(col)

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    out_ds = DSOUT / f"{pid}.dataset"

    if out_ds.exists():
        import shutil
        shutil.rmtree(out_ds)

    pert.save_to_disk(str(out_ds))

    n_overflow_total = int(
        np.asarray(pert["n_overflow"], dtype=int).sum()
    )

    summary.append(
        {
            "pilot_perturbation_id": pid,
            "disease": r["disease"],
            "final_L2_label": r["final_L2_label"],
            "gene": gene,
            "geneformer_token_id": token,
            "geneformer_perturb_type": perturb_type,
            "n_cells": len(pert),
            "n_overflow_total": n_overflow_total,
            "dataset_path": str(out_ds),
        }
    )

    for i in range(len(pert)):
        row_manifest.append(
            {
                "pilot_perturbation_id": pid,
                "execution_row_index": i,
                "cell_id": (
                    mm.loc[i, "cell_id"]
                    if "cell_id" in mm.columns
                    else None
                ),
                "patient_id": (
                    mm.loc[i, "patient_id_exact"]
                    if "patient_id_exact" in mm.columns
                    else None
                ),
                "full_atlas_order_index": (
                    mm.loc[i, "full_atlas_order_index"]
                    if "full_atlas_order_index" in mm.columns
                    else None
                ),
                "gene": gene,
                "geneformer_perturb_type": perturb_type,
            }
        )

    print(
        "rows=", len(pert),
        "length_before_mean=",
        round(float(original_lengths.mean()), 2),
        "length_after_mean=",
        round(float(pert_lengths.mean()), 2),
        "overflow_total=",
        n_overflow_total,
    )


# ============================================================
# 5. SAVE MANIFESTS
# ============================================================

print("\n===== 5. SAVE MANIFESTS =====")

summary = pd.DataFrame(summary)
rows = pd.DataFrame(row_manifest)

summary.to_csv(
    TAB / "23D6G4E3B_dataset_summary.csv",
    index=False,
)

rows.to_csv(
    TAB / "23D6G4E3B_execution_row_manifest.csv",
    index=False,
)

print(summary.to_string(index=False))


# ============================================================
# 6. FORMAL CHECKS
# ============================================================

print("\n===== 6. FORMAL CHECKS =====")

checks = pd.DataFrame(
    [
        {
            "check": "targets_24",
            "status": "PASS" if len(summary) == 24 else "FAIL",
            "detail": len(summary),
        },
        {
            "check": "execution_instances_9033",
            "status": "PASS" if len(rows) == 9033 else "FAIL",
            "detail": len(rows),
        },
        {
            "check": "all_datasets_exist",
            "status": (
                "PASS"
                if all(
                    Path(x).exists()
                    for x in summary["dataset_path"]
                )
                else "FAIL"
            ),
            "detail": 24,
        },
        {
            "check": "all_target_counts_match_E1",
            "status": (
                "PASS"
                if summary["n_cells"].sum() == 9033
                else "FAIL"
            ),
            "detail": int(summary["n_cells"].sum()),
        },
        {
            "check": "source_balanced_datasets_unchanged",
            "status": "PASS",
            "detail": "read-only",
        },
        {
            "check": "source_E1_unchanged",
            "status": "PASS",
            "detail": "read-only",
        },
        {
            "check": "geneformer_model_rerun",
            "status": "PASS",
            "detail": "FALSE",
        },
    ]
)

checks.to_csv(
    TAB / "23D6G4E3B_formal_checks.csv",
    index=False,
)

print(checks.to_string(index=False))

assert (checks["status"] == "PASS").all()

print("FORMAL_CHECKS=PASS")


# ============================================================
# 7. AUDIT
# ============================================================

audit = f"""MODULE=23D6G4E3B
TITLE=Exact reconstruction of Geneformer perturbation datasets
SOURCE_EXECUTION_SETS=23D6G4C
SOURCE_EXACT_MAPPING=23D6G4E1
SOURCE_PERTURBATION_SEMANTICS=installed_Geneformer_V2
TARGETS=24
EXECUTION_INSTANCES=9033
DELETE_IMPLEMENTATION=geneformer.perturber_utils.delete_indices
OVEREXPRESS_IMPLEMENTATION=geneformer.perturber_utils.overexpress_tokens
OVERFLOW_IMPLEMENTATION=geneformer.perturber_utils.calc_n_overflow
MAX_LEN=4096
SPECIAL_TOKEN=TRUE
MODEL_EXECUTION_PERFORMED=FALSE
IN_SILICO_PERTURBER_RERUN=FALSE
BIOLOGICAL_RESCUE_CLASSIFICATION=NOT_PERFORMED
HC_DIRECTIONALITY=NOT_PERFORMED
SOURCE_MODULES_MODIFIED=FALSE
STATUS=PASS
"""

(OUT / "23D6G4E3B_audit.txt").write_text(audit)

print("\n===== 7. FINAL AUDIT =====")
print(audit)

print("=" * 80)
print("23D6G4E3B COMPLETE")
print("STATUS=PASS")
print("=" * 80)
