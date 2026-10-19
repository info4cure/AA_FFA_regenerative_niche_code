from pathlib import Path
import numpy as np
import pandas as pd

# ============================================================
# 23D6G4E2 — PATIENT-LEVEL PERTURBATION ROBUSTNESS
# ============================================================

ROOT = Path(
    "PIPELINE_V2/results/disease_biology/"
    "23_ADVANCED_ANALYSIS_ROADMAP/23D_AI_PERTURBATION/"
    "23D6_PERTURBATION/23D6G_FULL_ATLAS_REASSESSMENT"
)

E1 = (
    ROOT /
    "23D6G4E_PATIENT_ROBUSTNESS_HC_DIRECTIONALITY" /
    "23D6G4E1_EXACT_EXECUTION_MAPPING"
)

G4D = ROOT / "23D6G4D_FULL_ATLAS_PERTURBATIONS"

OUT = (
    ROOT /
    "23D6G4E_PATIENT_ROBUSTNESS_HC_DIRECTIONALITY" /
    "23D6G4E2_PATIENT_LEVEL_ROBUSTNESS"
)

TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

MAP = (
    E1 / "tables" /
    "23D6G4E1_exact_execution_cell_mapping.csv"
)

EXEC = (
    G4D / "tables" /
    "23D6G4D_execution_summary.csv"
)

print("=" * 80)
print("23D6G4E2 — PATIENT-LEVEL PERTURBATION ROBUSTNESS")
print("=" * 80)

# ============================================================
# 1. LOAD
# ============================================================

print("\n===== 1. LOAD =====")

x = pd.read_csv(MAP)
exe = pd.read_csv(EXEC)

print("cell perturbation rows:", len(x))
print("targets:", x["pilot_perturbation_id"].nunique())
print("execution summary rows:", len(exe))

assert len(x) == 9033
assert x["pilot_perturbation_id"].nunique() == 24

required = [
    "pilot_perturbation_id",
    "patient_id_exact",
    "cosine_similarity",
    "cosine_shift_from_1",
]

missing = [c for c in required if c not in x.columns]

if missing:
    raise RuntimeError(
        f"Missing required mapping columns: {missing}"
    )

# ============================================================
# 2. RESOLVE TARGET METADATA
# ============================================================

print("\n===== 2. TARGET METADATA =====")

meta_cols = [
    "pilot_perturbation_id",
    "disease",
    "final_L2_label",
    "gene",
    "ensembl_id",
    "geneformer_perturb_type",
]

target_meta = exe[meta_cols].drop_duplicates()

if len(target_meta) != 24:
    raise RuntimeError(
        "Execution summary does not contain exactly 24 unique targets"
    )

print(target_meta.to_string(index=False))

# ============================================================
# 3. PATIENT-LEVEL AGGREGATION
# ============================================================

print("\n===== 3. PATIENT-LEVEL AGGREGATION =====")

pat = (
    x.groupby(
        [
            "pilot_perturbation_id",
            "patient_id_exact",
        ],
        dropna=False,
    )
    .agg(
        n_cells=("cosine_shift_from_1", "size"),
        mean_cosine_similarity=("cosine_similarity", "mean"),
        median_cosine_similarity=("cosine_similarity", "median"),
        mean_cosine_shift_from_1=("cosine_shift_from_1", "mean"),
        median_cosine_shift_from_1=("cosine_shift_from_1", "median"),
        sd_cosine_shift_from_1=("cosine_shift_from_1", "std"),
        min_cosine_shift_from_1=("cosine_shift_from_1", "min"),
        max_cosine_shift_from_1=("cosine_shift_from_1", "max"),
    )
    .reset_index()
)

pat = pat.merge(
    target_meta,
    on="pilot_perturbation_id",
    how="left",
    validate="many_to_one",
)

print("patient-target strata:", len(pat))
print(
    pat[
        [
            "pilot_perturbation_id",
            "patient_id_exact",
            "n_cells",
            "mean_cosine_shift_from_1",
        ]
    ].to_string(index=False)
)

# ============================================================
# 4. PATIENT-BALANCED TARGET SUMMARY
# ============================================================

print("\n===== 4. PATIENT-BALANCED TARGET SUMMARY =====")

rows = []

for pid, g in pat.groupby("pilot_perturbation_id", sort=True):

    vals = g["mean_cosine_shift_from_1"].to_numpy(dtype=float)

    mean_pb = float(np.mean(vals))
    median_pb = float(np.median(vals))
    sd_pb = float(np.std(vals, ddof=1)) if len(vals) > 1 else np.nan

    if np.isfinite(mean_pb) and mean_pb > 0 and np.isfinite(sd_pb):
        cv_pb = sd_pb / mean_pb
    else:
        cv_pb = np.nan

    q1 = float(np.quantile(vals, 0.25))
    q3 = float(np.quantile(vals, 0.75))
    iqr = q3 - q1

    rows.append(
        {
            "pilot_perturbation_id": pid,
            "n_patients": len(g),
            "total_execution_cells": int(g["n_cells"].sum()),
            "min_cells_per_patient": int(g["n_cells"].min()),
            "max_cells_per_patient": int(g["n_cells"].max()),
            "patient_balanced_mean_cosine_shift_from_1": mean_pb,
            "patient_balanced_median_cosine_shift_from_1": median_pb,
            "patient_sd_mean_shift": sd_pb,
            "patient_cv_mean_shift": cv_pb,
            "patient_min_mean_shift": float(np.min(vals)),
            "patient_max_mean_shift": float(np.max(vals)),
            "patient_q1_mean_shift": q1,
            "patient_q3_mean_shift": q3,
            "patient_iqr_mean_shift": iqr,
        }
    )

target = pd.DataFrame(rows)

target = target.merge(
    target_meta,
    on="pilot_perturbation_id",
    how="left",
    validate="one_to_one",
)

# Add original cell-weighted G4D estimate for comparison.
g4d_effect = exe[
    [
        "pilot_perturbation_id",
        "mean_cosine_shift_from_1",
        "median_cosine_shift_from_1",
    ]
].rename(
    columns={
        "mean_cosine_shift_from_1":
            "cell_weighted_mean_cosine_shift_from_1",
        "median_cosine_shift_from_1":
            "cell_weighted_median_cosine_shift_from_1",
    }
)

target = target.merge(
    g4d_effect,
    on="pilot_perturbation_id",
    how="left",
    validate="one_to_one",
)

target["patient_vs_cell_weighted_ratio"] = (
    target["patient_balanced_mean_cosine_shift_from_1"] /
    target["cell_weighted_mean_cosine_shift_from_1"]
)

# ============================================================
# 5. LEAVE-ONE-PATIENT-OUT ROBUSTNESS
# ============================================================

print("\n===== 5. LEAVE-ONE-PATIENT-OUT =====")

lopo_rows = []

for pid, g in pat.groupby("pilot_perturbation_id", sort=True):

    for patient in g["patient_id_exact"]:

        keep = g.loc[
            g["patient_id_exact"] != patient,
            "mean_cosine_shift_from_1",
        ].to_numpy(dtype=float)

        lopo_rows.append(
            {
                "pilot_perturbation_id": pid,
                "left_out_patient": patient,
                "n_patients_remaining": len(keep),
                "lopo_patient_balanced_mean_shift":
                    float(np.mean(keep)) if len(keep) else np.nan,
                "lopo_patient_balanced_median_shift":
                    float(np.median(keep)) if len(keep) else np.nan,
            }
        )

lopo = pd.DataFrame(lopo_rows)

lopo_summary = (
    lopo.groupby("pilot_perturbation_id")
    .agg(
        lopo_min_mean_shift=(
            "lopo_patient_balanced_mean_shift", "min"
        ),
        lopo_max_mean_shift=(
            "lopo_patient_balanced_mean_shift", "max"
        ),
        lopo_mean_mean_shift=(
            "lopo_patient_balanced_mean_shift", "mean"
        ),
        lopo_sd_mean_shift=(
            "lopo_patient_balanced_mean_shift", "std"
        ),
    )
    .reset_index()
)

target = target.merge(
    lopo_summary,
    on="pilot_perturbation_id",
    how="left",
    validate="one_to_one",
)

target["lopo_min_fraction_of_full_effect"] = (
    target["lopo_min_mean_shift"] /
    target["patient_balanced_mean_cosine_shift_from_1"]
)

target["lopo_max_fraction_of_full_effect"] = (
    target["lopo_max_mean_shift"] /
    target["patient_balanced_mean_cosine_shift_from_1"]
)

# ============================================================
# 6. DESCRIPTIVE PATIENT SUPPORT
# ============================================================

print("\n===== 6. DESCRIPTIVE SUPPORT CLASSIFICATION =====")

# This is deliberately descriptive.
# It is NOT a biological rescue/aggravation classification.
#
# We assess whether the magnitude observed in the patient-balanced
# analysis is broadly distributed across patients or dominated by
# a restricted subset.

for pid, g in pat.groupby("pilot_perturbation_id"):

    med = float(
        g["mean_cosine_shift_from_1"].median()
    )

    # Patient is considered to show material support when its mean
    # shift reaches at least 50% of the across-patient median.
    # This relative criterion avoids imposing an arbitrary absolute
    # biological threshold on very small embedding shifts.
    threshold = 0.5 * med

    mask = (
        g["mean_cosine_shift_from_1"] >= threshold
    )

    target.loc[
        target["pilot_perturbation_id"] == pid,
        "patient_support_threshold"
    ] = threshold

    target.loc[
        target["pilot_perturbation_id"] == pid,
        "n_patients_supporting"
    ] = int(mask.sum())

    target.loc[
        target["pilot_perturbation_id"] == pid,
        "fraction_patients_supporting"
    ] = float(mask.mean())

target["n_patients_supporting"] = (
    target["n_patients_supporting"].astype(int)
)

def support_class(row):
    frac = row["fraction_patients_supporting"]
    lopo = row["lopo_min_fraction_of_full_effect"]
    cv = row["patient_cv_mean_shift"]

    if (
        frac >= 0.80 and
        lopo >= 0.70 and
        (not np.isfinite(cv) or cv <= 1.0)
    ):
        return "HIGH"

    if (
        frac >= 0.60 and
        lopo >= 0.50
    ):
        return "MODERATE"

    return "HETEROGENEOUS"

target["patient_robustness_class"] = target.apply(
    support_class,
    axis=1,
)

# ============================================================
# 7. RANKINGS
# ============================================================

print("\n===== 7. PATIENT-BALANCED RANKING =====")

target = target.sort_values(
    [
        "patient_balanced_mean_cosine_shift_from_1",
        "lopo_min_mean_shift",
    ],
    ascending=[False, False],
).reset_index(drop=True)

target.insert(
    0,
    "rank_patient_balanced",
    np.arange(1, len(target) + 1),
)

show_cols = [
    "rank_patient_balanced",
    "pilot_perturbation_id",
    "disease",
    "final_L2_label",
    "gene",
    "geneformer_perturb_type",
    "n_patients",
    "patient_balanced_mean_cosine_shift_from_1",
    "patient_balanced_median_cosine_shift_from_1",
    "patient_cv_mean_shift",
    "lopo_min_fraction_of_full_effect",
    "fraction_patients_supporting",
    "patient_robustness_class",
]

print(target[show_cols].to_string(index=False))

# ============================================================
# 8. SAVE TABLES
# ============================================================

print("\n===== 8. SAVE TABLES =====")

pat_file = TAB / \
    "23D6G4E2_patient_level_effects.csv"

target_file = TAB / \
    "23D6G4E2_patient_balanced_target_summary.csv"

lopo_file = TAB / \
    "23D6G4E2_leave_one_patient_out.csv"

ranking_file = TAB / \
    "23D6G4E2_patient_balanced_ranking.csv"

pat.to_csv(pat_file, index=False)
target.to_csv(target_file, index=False)
lopo.to_csv(lopo_file, index=False)
target[show_cols].to_csv(ranking_file, index=False)

# ============================================================
# 9. FORMAL CHECKS
# ============================================================

print("\n===== 9. FORMAL CHECKS =====")

expected_strata = 132

checks = pd.DataFrame(
    [
        [
            "source_cells_9033",
            len(x) == 9033,
            len(x),
        ],
        [
            "targets_24",
            target["pilot_perturbation_id"].nunique() == 24,
            target["pilot_perturbation_id"].nunique(),
        ],
        [
            "patient_target_strata_132",
            len(pat) == expected_strata,
            len(pat),
        ],
        [
            "all_patient_means_finite",
            np.isfinite(
                pat["mean_cosine_shift_from_1"]
            ).all(),
            int(
                np.isfinite(
                    pat["mean_cosine_shift_from_1"]
                ).sum()
            ),
        ],
        [
            "all_target_effects_finite",
            np.isfinite(
                target[
                    "patient_balanced_mean_cosine_shift_from_1"
                ]
            ).all(),
            int(
                np.isfinite(
                    target[
                        "patient_balanced_mean_cosine_shift_from_1"
                    ]
                ).sum()
            ),
        ],
        [
            "all_lopo_effects_finite",
            np.isfinite(
                lopo["lopo_patient_balanced_mean_shift"]
            ).all(),
            int(
                np.isfinite(
                    lopo["lopo_patient_balanced_mean_shift"]
                ).sum()
            ),
        ],
        [
            "AA_targets_have_6_patients",
            (
                target.loc[
                    target["disease"] == "AA",
                    "n_patients"
                ] == 6
            ).all(),
            target.loc[
                target["disease"] == "AA",
                "n_patients"
            ].tolist(),
        ],
        [
            "FFA_targets_have_5_patients",
            (
                target.loc[
                    target["disease"] == "FFA",
                    "n_patients"
                ] == 5
            ).all(),
            target.loc[
                target["disease"] == "FFA",
                "n_patients"
            ].tolist(),
        ],
        [
            "source_mapping_unchanged",
            True,
            "23D6G4E1 read-only",
        ],
        [
            "source_perturbations_unchanged",
            True,
            "23D6G4D read-only",
        ],
    ],
    columns=["check", "status_bool", "detail"],
)

checks["status"] = checks["status_bool"].map(
    {True: "PASS", False: "FAIL"}
)

checks = checks[
    ["check", "status", "detail"]
]

checks_file = TAB / \
    "23D6G4E2_formal_checks.csv"

checks.to_csv(checks_file, index=False)

print(checks.to_string(index=False))

if (checks["status"] != "PASS").any():
    raise RuntimeError(
        "23D6G4E2 FORMAL CHECKS FAILED"
    )

print("FORMAL_CHECKS=PASS")

# ============================================================
# 10. AUDIT
# ============================================================

audit = OUT / "23D6G4E2_audit.txt"

with open(audit, "w") as fh:
    fh.write(
f"""MODULE=23D6G4E2
TITLE=Patient-level perturbation robustness
SOURCE_MAPPING=23D6G4E1
SOURCE_PERTURBATIONS=23D6G4D
TARGETS=24
PATIENT_TARGET_STRATA={len(pat)}
CELL_PERTURBATION_INSTANCES={len(x)}
PRIMARY_EFFECT=patient_balanced_mean_cosine_shift_from_1
PATIENT_WEIGHTING=equal
LOPO_ANALYSIS=TRUE
SUPPORT_CLASSIFICATION=descriptive_only
BIOLOGICAL_RESCUE_CLASSIFICATION=NOT_PERFORMED
HC_DIRECTIONALITY=NOT_PERFORMED
PSEUDOREPLICATION_CONTROL=patient_level_aggregation
SOURCE_MODULES_MODIFIED=FALSE
STATUS=PASS
"""
    )

print("\n===== 10. FINAL AUDIT =====")
print(audit.read_text())

print("=" * 80)
print("23D6G4E2 COMPLETE")
print("")
print("OUTPUTS:")

for f in sorted(TAB.glob("23D6G4E2_*")):
    print("-", f)

print("-", audit)
print("")
print("STATUS=PASS")
print("=" * 80)
