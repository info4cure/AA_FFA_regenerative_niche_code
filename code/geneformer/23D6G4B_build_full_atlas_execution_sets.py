from pathlib import Path
from collections import Counter
import pandas as pd
import numpy as np
from datasets import load_from_disk

print("=" * 80)
print("23D6G4B — BUILD FULL-ATLAS PERTURBATION EXECUTION SETS")
print("24-TARGET FROZEN PANEL — NO MODEL EXECUTION")
print("=" * 80)

ROOT = Path.cwd()

PANEL = ROOT / (
    "PIPELINE_V2/results/disease_biology/"
    "23_ADVANCED_ANALYSIS_ROADMAP/"
    "23D_AI_PERTURBATION/"
    "23D6_TARGET_SELECTION/"
    "23D6C_PILOT_PANEL_FREEZE/"
    "tables/23D6C_frozen_pilot_perturbation_panel.csv"
)

DATASET = ROOT / (
    "PIPELINE_V2/results/disease_biology/"
    "23_ADVANCED_ANALYSIS_ROADMAP/"
    "23D_AI_PERTURBATION/"
    "23D6_PERTURBATION/"
    "23D6G_FULL_ATLAS_REASSESSMENT/"
    "23D6G2C_FULL_ATLAS_TOKENIZATION/"
    "dataset/23D6G2C_geneformer_full_atlas.dataset"
)

OUT = ROOT / (
    "PIPELINE_V2/results/disease_biology/"
    "23_ADVANCED_ANALYSIS_ROADMAP/"
    "23D_AI_PERTURBATION/"
    "23D6_PERTURBATION/"
    "23D6G_FULL_ATLAS_REASSESSMENT/"
    "23D6G4B_FULL_ATLAS_EXECUTION_SETS"
)

OBJ = OUT / "objects" / "execution_datasets"
TAB = OUT / "tables"

OBJ.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------
# 1. LOAD
# ----------------------------------------------------------------------

panel = pd.read_csv(PANEL)
ds = load_from_disk(str(DATASET))

assert len(panel) == 24
assert panel["pilot_perturbation_id"].nunique() == 24
assert len(ds) == 38699

required = [
    "cell_id",
    "patient_id",
    "disease",
    "L2",
    "input_ids",
    "full_atlas_order_index",
]

missing = [x for x in required if x not in ds.column_names]
assert not missing, missing

print("\n===== SOURCE =====")
print("panel rows:", len(panel))
print("dataset rows:", len(ds))
print("patients:", len(set(ds["patient_id"])))

# ----------------------------------------------------------------------
# 2. MASTER MANIFEST
# ----------------------------------------------------------------------

manifest = pd.DataFrame({
    "dataset_index": np.arange(len(ds), dtype=int),
    "cell_id": ds["cell_id"],
    "patient_id": ds["patient_id"],
    "disease": ds["disease"],
    "L2": ds["L2"],
    "compartment": ds["compartment"],
    "full_atlas_order_index": ds["full_atlas_order_index"],
})

assert manifest["cell_id"].is_unique
assert manifest["full_atlas_order_index"].is_unique

# ----------------------------------------------------------------------
# 3. BUILD EXECUTION SETS
# ----------------------------------------------------------------------

summary = []
cell_manifest_rows = []

for _, r in panel.iterrows():

    pid = str(r["pilot_perturbation_id"])
    disease = str(r["disease"])
    l2 = str(r["final_L2_label"])
    gene = str(r["gene"])
    ensembl = str(r["ensembl_id"])
    token = int(r["geneformer_token_id"])
    ptype = str(r["geneformer_perturb_type"])

    print("\n" + "-" * 80)
    print(pid, "|", disease, "|", l2, "|", gene, "|", ptype)

    state_idx = manifest.index[
        (manifest["disease"] == disease) &
        (manifest["L2"] == l2)
    ].to_numpy(dtype=int)

    target_present = np.asarray(
        [token in ds[int(i)]["input_ids"] for i in state_idx],
        dtype=bool,
    )

    # Match Geneformer semantics used in the pilot:
    # delete requires the target to be present;
    # overexpression can introduce a target that is absent.
    if ptype == "delete":
        exec_idx = state_idx[target_present]
    elif ptype == "overexpress":
        exec_idx = state_idx
    else:
        raise ValueError(f"Unknown perturbation type: {ptype}")

    exec_idx = np.asarray(exec_idx, dtype=int)

    n_state = len(state_idx)
    n_present = int(target_present.sum())
    n_exec = len(exec_idx)

    if n_exec:
        pats = manifest.loc[exec_idx, "patient_id"].astype(str)
        n_patients = pats.nunique()
        patient_counts = Counter(pats)

        sub = ds.select(exec_idx.tolist())

        # Persist exact full-atlas cells for this frozen perturbation.
        outfile = OBJ / f"{pid}.dataset"

        if outfile.exists():
            import shutil
            shutil.rmtree(outfile)

        sub.save_to_disk(str(outfile))

        for local_i, global_i in enumerate(exec_idx):
            row = manifest.loc[int(global_i)]

            tokens = ds[int(global_i)]["input_ids"]
            present = token in tokens

            cell_manifest_rows.append({
                "pilot_perturbation_id": pid,
                "disease": disease,
                "final_L2_label": l2,
                "gene": gene,
                "ensembl_id": ensembl,
                "geneformer_token_id": token,
                "geneformer_perturb_type": ptype,
                "dataset_index": int(global_i),
                "execution_local_index": int(local_i),
                "cell_id": row["cell_id"],
                "patient_id": row["patient_id"],
                "cell_disease": row["disease"],
                "cell_L2": row["L2"],
                "full_atlas_order_index":
                    int(row["full_atlas_order_index"]),
                "target_present_before": bool(present),
                "sequence_length":
                    int(ds[int(global_i)]["length"]),
                "length_uncropped":
                    int(ds[int(global_i)]["length_uncropped"]),
            })

    else:
        n_patients = 0
        patient_counts = {}

    # Simple support classification for execution reporting only.
    if n_patients >= 4 and n_exec >= 20:
        support = "ROBUST"
    elif n_exec > 0:
        support = "LIMITED"
    else:
        support = "NONE"

    summary.append({
        "pilot_perturbation_id": pid,
        "disease": disease,
        "final_L2_label": l2,
        "gene": gene,
        "ensembl_id": ensembl,
        "geneformer_token_id": token,
        "geneformer_perturb_type": ptype,
        "rescue_direction": r["rescue_direction"],
        "pre_geneformer_evidence_class":
            r["pre_geneformer_evidence_class"],
        "n_state_cells": n_state,
        "n_target_present": n_present,
        "fraction_target_present":
            (n_present / n_state if n_state else np.nan),
        "n_execution_cells": n_exec,
        "n_execution_patients": n_patients,
        "execution_support_class_full_atlas": support,
        "patient_cell_counts":
            ";".join(
                f"{k}:{v}"
                for k, v in sorted(patient_counts.items())
            ),
        "execution_dataset":
            str(OBJ / f"{pid}.dataset")
            if n_exec else "",
    })

    print("state cells:", n_state)
    print("target present:", n_present)
    print(
        "target fraction:",
        round(n_present / n_state, 4)
        if n_state else "NA"
    )
    print("execution cells:", n_exec)
    print("execution patients:", n_patients)
    print("support:", support)

# ----------------------------------------------------------------------
# 4. OUTPUT TABLES
# ----------------------------------------------------------------------

summary = pd.DataFrame(summary)

cell_manifest = pd.DataFrame(cell_manifest_rows)

summary.to_csv(
    TAB / "23D6G4B_full_atlas_execution_support.csv",
    index=False,
)

cell_manifest.to_csv(
    TAB / "23D6G4B_full_atlas_execution_cell_manifest.csv",
    index=False,
)

print("\n===== FULL PANEL =====")

cols = [
    "pilot_perturbation_id",
    "disease",
    "final_L2_label",
    "gene",
    "geneformer_perturb_type",
    "n_state_cells",
    "n_target_present",
    "fraction_target_present",
    "n_execution_cells",
    "n_execution_patients",
    "execution_support_class_full_atlas",
]

print(summary[cols].to_string(index=False))

# ----------------------------------------------------------------------
# 5. KEY QUESTION: IRF7
# ----------------------------------------------------------------------

print("\n===== 23D6C_006 / IRF7 REASSESSMENT =====")

irf7 = summary[
    summary["pilot_perturbation_id"] == "23D6C_006"
]

assert len(irf7) == 1

print(irf7[cols].to_string(index=False))

# ----------------------------------------------------------------------
# 6. INVARIANTS
# ----------------------------------------------------------------------

print("\n===== FORMAL CHECKS =====")

checks = []

def check(name, condition, detail):
    status = "PASS" if bool(condition) else "FAIL"
    checks.append({
        "check": name,
        "status": status,
        "detail": detail,
    })
    print(name, status, detail)

check(
    "frozen_panel_n24",
    len(summary) == 24,
    len(summary),
)

check(
    "unique_perturbation_ids",
    summary["pilot_perturbation_id"].nunique() == 24,
    summary["pilot_perturbation_id"].nunique(),
)

check(
    "all_state_universes_nonzero",
    (summary["n_state_cells"] > 0).all(),
    int((summary["n_state_cells"] > 0).sum()),
)

delete_rows = summary[
    summary["geneformer_perturb_type"] == "delete"
]

check(
    "delete_execution_equals_target_present",
    (
        delete_rows["n_execution_cells"] ==
        delete_rows["n_target_present"]
    ).all(),
    f"{len(delete_rows)}/{len(delete_rows)}",
)

over_rows = summary[
    summary["geneformer_perturb_type"] == "overexpress"
]

check(
    "overexpress_execution_equals_state_universe",
    (
        over_rows["n_execution_cells"] ==
        over_rows["n_state_cells"]
    ).all(),
    f"{len(over_rows)}/{len(over_rows)}",
)

if len(cell_manifest):
    check(
        "execution_cell_ids_unique_within_perturbation",
        not cell_manifest.duplicated(
            ["pilot_perturbation_id", "cell_id"]
        ).any(),
        len(cell_manifest),
    )

    check(
        "all_execution_cells_state_specific",
        (
            (cell_manifest["disease"] ==
             cell_manifest["cell_disease"]) &
            (cell_manifest["final_L2_label"] ==
             cell_manifest["cell_L2"])
        ).all(),
        len(cell_manifest),
    )

checks = pd.DataFrame(checks)

checks.to_csv(
    TAB / "23D6G4B_formal_checks.csv",
    index=False,
)

assert (checks["status"] == "PASS").all()

# ----------------------------------------------------------------------
# 7. AUDIT
# ----------------------------------------------------------------------

n_exec_targets = int(
    (summary["n_execution_cells"] > 0).sum()
)

n_nonexec_targets = int(
    (summary["n_execution_cells"] == 0).sum()
)

total_instances = int(
    summary["n_execution_cells"].sum()
)

irf7_n = int(
    irf7.iloc[0]["n_execution_cells"]
)

audit = f"""MODULE=23D6G4B
TITLE=Full-atlas frozen perturbation execution-set construction
SOURCE_PANEL=23D6C
SOURCE_DATASET=23D6G2C
FROZEN_TARGETS=24
EXECUTABLE_TARGETS={n_exec_targets}
NONEXECUTABLE_TARGETS={n_nonexec_targets}
TOTAL_CELL_PERTURBATION_INSTANCES={total_instances}
IRF7_DELETE_EXECUTION_CELLS={irf7_n}
MODEL_EXECUTION_PERFORMED=FALSE
PERTURBATION_PERFORMED=FALSE
STATUS=PASS
"""

(OUT / "23D6G4B_audit.txt").write_text(audit)

print("\n" + audit)

print("=" * 80)
print("23D6G4B COMPLETE")
print("=" * 80)
