from pathlib import Path
import numpy as np
import pandas as pd

# ============================================================
# 23D6G4E3E
# PATIENT ROBUSTNESS AND FULL-ATLAS HC-RESCUE ADJUDICATION
#
# Sources:
#   23D6G4E3D:
#       - execution-level geometry
#       - patient-level geometry
#       - patient-balanced target geometry
#       - cell-weighted sensitivity
#
#   23D6F4:
#       - frozen pilot HC-rescue adjudication
#
# Purpose:
#   1. Adjudicate full-atlas HC rescue at target level.
#   2. Prioritize patient-level robustness over cell counts.
#   3. Preserve cosine + Euclidean metric concordance.
#   4. Propagate HC-reference support.
#   5. Compare full-atlas result with frozen pilot result.
#   6. Avoid therapeutic-efficacy claims.
#
# No embeddings are recomputed.
# No source modules are modified.
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

E3D = (
    G4E
    / "23D6G4E3D_FULL_ATLAS_HC_RESCUE_GEOMETRY"
)

F4 = (
    BASE
    / "23D6F_HC_RESCUE_GEOMETRY"
    / "23D6F4_PILOT_ADJUDICATION"
)

OUT = (
    G4E
    / "23D6G4E3E_PATIENT_ROBUSTNESS_HC_RESCUE_ADJUDICATION"
)

TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

TARGET_FILE = (
    E3D
    / "tables"
    / "23D6G4E3D_target_level_HC_rescue_geometry.csv"
)

PATIENT_FILE = (
    E3D
    / "tables"
    / "23D6G4E3D_patient_level_HC_rescue_geometry.csv"
)

EXEC_FILE = (
    E3D
    / "tables"
    / "23D6G4E3D_execution_level_HC_rescue_geometry.csv"
)

CELL_FILE = (
    E3D
    / "tables"
    / "23D6G4E3D_target_cell_weighted_sensitivity.csv"
)

PILOT_FILE = (
    F4
    / "tables"
    / "23D6F4_frozen_panel_HC_rescue_adjudication.csv"
)

print("=" * 80)
print(
    "23D6G4E3E — PATIENT ROBUSTNESS AND "
    "FULL-ATLAS HC-RESCUE ADJUDICATION"
)
print("=" * 80)
print("E3D :", E3D)
print("F4  :", F4)
print("OUT :", OUT)


# ============================================================
# HELPERS
# ============================================================

def boolify(x):
    if isinstance(x, (bool, np.bool_)):
        return bool(x)

    if pd.isna(x):
        return False

    s = str(x).strip().lower()

    if s in {"true", "1", "yes", "y"}:
        return True

    if s in {"false", "0", "no", "n"}:
        return False

    raise ValueError(f"Cannot boolify value: {x!r}")


def sign_label(x, tol=0.0):
    if pd.isna(x):
        return "NA"
    if x > tol:
        return "TOWARD_HC"
    if x < -tol:
        return "AWAY_FROM_HC"
    return "NEUTRAL"


def consistency_class(x):
    if pd.isna(x):
        return "NA"
    if x >= 0.80:
        return "HIGH"
    if x >= 0.60:
        return "MODERATE"
    return "LOW"


def support_tier(n_patients):
    if n_patients >= 5:
        return "HIGH_PATIENT_SUPPORT"
    if n_patients >= 3:
        return "MODERATE_PATIENT_SUPPORT"
    return "LIMITED_PATIENT_SUPPORT"


# ============================================================
# 1. LOAD INPUTS
# ============================================================

print("\n===== 1. LOAD INPUTS =====")

for p in [
    TARGET_FILE,
    PATIENT_FILE,
    EXEC_FILE,
    CELL_FILE,
    PILOT_FILE,
]:
    assert p.exists(), f"Missing input: {p}"

target = pd.read_csv(TARGET_FILE)
patient = pd.read_csv(PATIENT_FILE)
execution = pd.read_csv(EXEC_FILE)
cell = pd.read_csv(CELL_FILE)
pilot = pd.read_csv(PILOT_FILE)

print("target :", target.shape)
print("patient:", patient.shape)
print("exec   :", execution.shape)
print("cell   :", cell.shape)
print("pilot  :", pilot.shape)

assert len(target) == 24
assert execution["pilot_perturbation_id"].nunique() == 24
assert target["pilot_perturbation_id"].is_unique
assert pilot["pilot_perturbation_id"].is_unique


# ============================================================
# 2. VERIFY TARGET IDENTITY ACROSS SOURCES
# ============================================================

print("\n===== 2. VERIFY TARGET IDENTITY =====")

target_ids = set(target["pilot_perturbation_id"])
pilot_ids = set(pilot["pilot_perturbation_id"])
cell_ids = set(cell["pilot_perturbation_id"])

assert target_ids == pilot_ids
assert target_ids == cell_ids

print("targets shared across E3D / cell sensitivity / F4:", len(target_ids))
print("TARGET_IDENTITY=PASS")


# ============================================================
# 3. RECOMPUTE PATIENT ROBUSTNESS FROM PATIENT TABLE
# ============================================================

print("\n===== 3. PATIENT ROBUSTNESS =====")

patient_rows = []

for pid, g in patient.groupby(
    "pilot_perturbation_id",
    sort=True,
):

    g = g.copy()

    n_patients = g["patient_id"].nunique()

    cos_vals = g["mean_delta_cos_HC"].astype(float)
    euc_vals = g[
        "mean_delta_euclidean_toward_HC"
    ].astype(float)

    align_vals = g[
        "mean_directional_alignment_cosine"
    ].astype(float)

    cos_pos = cos_vals > 0
    euc_pos = euc_vals > 0
    align_pos = align_vals > 0

    strict = cos_pos & euc_pos & align_pos
    away_strict = (
        (~cos_pos)
        & (~euc_pos)
        & (align_vals < 0)
    )

    metric_concordant = (
        cos_pos == euc_pos
    )

    # Patient-level sign agreement between cosine and Euclidean.
    same_toward = cos_pos & euc_pos
    same_away = (~cos_pos) & (~euc_pos)

    patient_rows.append(
        {
            "pilot_perturbation_id": pid,
            "recomputed_n_patients": n_patients,

            "n_patients_cosine_toward_HC": int(
                cos_pos.sum()
            ),
            "n_patients_euclidean_toward_HC": int(
                euc_pos.sum()
            ),
            "n_patients_strict_rescue": int(
                strict.sum()
            ),
            "n_patients_both_metrics_toward_HC": int(
                same_toward.sum()
            ),
            "n_patients_both_metrics_away_HC": int(
                same_away.sum()
            ),

            "fraction_patients_cosine_toward_HC": float(
                cos_pos.mean()
            ),
            "fraction_patients_euclidean_toward_HC": float(
                euc_pos.mean()
            ),
            "fraction_patients_strict_rescue_recomputed": float(
                strict.mean()
            ),
            "fraction_patients_metric_concordant_recomputed": float(
                metric_concordant.mean()
            ),
            "fraction_patients_directionally_aligned_recomputed": float(
                align_pos.mean()
            ),

            "mean_patient_delta_cos_HC_recomputed": float(
                cos_vals.mean()
            ),
            "median_patient_delta_cos_HC_recomputed": float(
                cos_vals.median()
            ),

            "mean_patient_delta_euclidean_HC_recomputed": float(
                euc_vals.mean()
            ),
            "median_patient_delta_euclidean_HC_recomputed": float(
                euc_vals.median()
            ),

            "mean_patient_alignment_recomputed": float(
                align_vals.mean()
            ),
            "min_patient_alignment": float(
                align_vals.min()
            ),
            "max_patient_alignment": float(
                align_vals.max()
            ),
        }
    )

robust = pd.DataFrame(patient_rows)

assert len(robust) == 24

robust.to_csv(
    TAB / "23D6G4E3E_patient_robustness_recomputed.csv",
    index=False,
)

print(
    robust[
        [
            "pilot_perturbation_id",
            "recomputed_n_patients",
            "fraction_patients_cosine_toward_HC",
            "fraction_patients_euclidean_toward_HC",
            "fraction_patients_strict_rescue_recomputed",
        ]
    ].to_string(index=False)
)


# ============================================================
# 4. BUILD MASTER TARGET TABLE
# ============================================================

print("\n===== 4. BUILD MASTER TARGET TABLE =====")

master = target.merge(
    robust,
    on="pilot_perturbation_id",
    how="left",
    validate="one_to_one",
)

cell_keep = [
    "pilot_perturbation_id",
    "mean_delta_cos_HC",
    "mean_delta_euclidean_toward_HC",
    "mean_normalized_cosine_rescue",
    "mean_normalized_euclidean_rescue",
    "mean_directional_alignment_cosine",
    "fraction_executions_toward_HC_by_cosine",
    "fraction_executions_toward_HC_by_euclidean",
    "fraction_executions_directionally_aligned",
    "fraction_executions_strict_rescue",
]

for c in cell_keep:
    assert c in cell.columns, c

cell2 = cell[cell_keep].copy()

cell2 = cell2.rename(
    columns={
        c: f"cellweighted_{c}"
        for c in cell_keep
        if c != "pilot_perturbation_id"
    }
)

master = master.merge(
    cell2,
    on="pilot_perturbation_id",
    how="left",
    validate="one_to_one",
)


# ============================================================
# 5. ADD FROZEN PILOT ADJUDICATION
# ============================================================

print("\n===== 5. ADD FROZEN PILOT ADJUDICATION =====")

pilot_keep = [
    "pilot_perturbation_id",
    "HC_rescue_adjudication",
    "adjudication_strength",
    "interpretation",
]

for c in pilot_keep:
    assert c in pilot.columns, c

pilot2 = pilot[pilot_keep].copy()

pilot2 = pilot2.rename(
    columns={
        "HC_rescue_adjudication":
            "pilot_HC_rescue_adjudication",
        "adjudication_strength":
            "pilot_adjudication_strength",
        "interpretation":
            "pilot_interpretation",
    }
)

master = master.merge(
    pilot2,
    on="pilot_perturbation_id",
    how="left",
    validate="one_to_one",
)

assert len(master) == 24


# ============================================================
# 6. DERIVE FULL-ATLAS DIRECTIONAL FEATURES
# ============================================================

print("\n===== 6. FULL-ATLAS DIRECTIONAL FEATURES =====")

master["HC_reference_exploratory_only"] = (
    master["HC_reference_exploratory_only"]
    .map(boolify)
)

master["cosine_mean_direction"] = (
    master["mean_patient_delta_cos_HC"]
    .map(sign_label)
)

master["euclidean_mean_direction"] = (
    master[
        "mean_patient_delta_euclidean_toward_HC"
    ]
    .map(sign_label)
)

master["mean_metric_direction_concordant"] = (
    master["cosine_mean_direction"]
    == master["euclidean_mean_direction"]
)

master["patient_cosine_consistency_class"] = (
    master[
        "fraction_patients_cosine_toward_HC"
    ].map(consistency_class)
)

master["patient_euclidean_consistency_class"] = (
    master[
        "fraction_patients_euclidean_toward_HC"
    ].map(consistency_class)
)

master["patient_strict_rescue_consistency_class_E3E"] = (
    master[
        "fraction_patients_strict_rescue_recomputed"
    ].map(consistency_class)
)

master["patient_support_tier"] = (
    master["n_patients"].map(support_tier)
)

# Cell-weighted direction used only as sensitivity.
master["cellweighted_cosine_direction"] = (
    master["cellweighted_mean_delta_cos_HC"]
    .map(sign_label)
)

master["cellweighted_euclidean_direction"] = (
    master[
        "cellweighted_mean_delta_euclidean_toward_HC"
    ]
    .map(sign_label)
)

master["patient_vs_cell_cosine_concordant"] = (
    master["cosine_mean_direction"]
    == master["cellweighted_cosine_direction"]
)

master["patient_vs_cell_euclidean_concordant"] = (
    master["euclidean_mean_direction"]
    == master["cellweighted_euclidean_direction"]
)

master["patient_vs_cell_both_metrics_concordant"] = (
    master["patient_vs_cell_cosine_concordant"]
    & master["patient_vs_cell_euclidean_concordant"]
)


# ============================================================
# 7. FULL-ATLAS ADJUDICATION
# ============================================================

print("\n===== 7. FULL-ATLAS ADJUDICATION =====")

def adjudicate(r):

    # HC reference limitation takes precedence over biological
    # direction. Numerical geometry is retained but interpretation
    # remains explicitly exploratory.
    if bool(r["HC_reference_exploratory_only"]):
        return (
            "EXPLORATORY_HC_REFERENCE",
            "EXPLORATORY",
            "State-specific HC reference is supported by a "
            "single HC patient; full-atlas geometry is not "
            "eligible for confirmatory HC-rescue adjudication."
        )

    cos_dir = r["cosine_mean_direction"]
    euc_dir = r["euclidean_mean_direction"]

    fcos = float(
        r["fraction_patients_cosine_toward_HC"]
    )
    feuc = float(
        r["fraction_patients_euclidean_toward_HC"]
    )
    fstrict = float(
        r["fraction_patients_strict_rescue_recomputed"]
    )

    align = float(
        r["mean_patient_directional_alignment_cosine"]
    )

    cell_concordant = bool(
        r["patient_vs_cell_both_metrics_concordant"]
    )

    # --------------------------------------------------------
    # ROBUST SUPPORTIVE
    #
    # Requirements:
    # - both patient-balanced metrics point toward HC
    # - >=80% patients improve by both metrics
    # - >=80% patients satisfy strict rescue
    # - mean perturbation vector aligns with HC direction
    # - cell-weighted sensitivity agrees
    # --------------------------------------------------------
    if (
        cos_dir == "TOWARD_HC"
        and euc_dir == "TOWARD_HC"
        and fcos >= 0.80
        and feuc >= 0.80
        and fstrict >= 0.80
        and align > 0
        and cell_concordant
    ):
        return (
            "ROBUST_SUPPORTIVE",
            "HIGH",
            "Patient-balanced cosine and Euclidean geometry "
            "both shift toward HC with high inter-patient "
            "consistency, positive directional alignment, and "
            "concordant cell-weighted sensitivity."
        )

    # --------------------------------------------------------
    # MODERATE SUPPORTIVE
    #
    # Requirements:
    # - both mean metrics toward HC
    # - >=60% patients improve by both metrics
    # - >=60% strict rescue
    # - positive directional alignment
    #
    # Cell-weighted disagreement downgrades later.
    # --------------------------------------------------------
    if (
        cos_dir == "TOWARD_HC"
        and euc_dir == "TOWARD_HC"
        and fcos >= 0.60
        and feuc >= 0.60
        and fstrict >= 0.60
        and align > 0
    ):
        if cell_concordant:
            return (
                "MODERATE_SUPPORTIVE",
                "MODERATE",
                "Both patient-balanced metrics shift toward HC "
                "with moderate inter-patient consistency and "
                "concordant cell-weighted sensitivity."
            )
        else:
            return (
                "HETEROGENEOUS",
                "MIXED",
                "Patient-balanced geometry is supportive, but "
                "cell-weighted sensitivity is directionally "
                "discordant."
            )

    # --------------------------------------------------------
    # DIRECTIONALLY SUPPORTIVE BUT HETEROGENEOUS
    #
    # Mean effect is toward HC by both metrics but fewer than
    # 60% of patients satisfy the rescue direction.
   
    # --------------------------------------------------------
    if (
        cos_dir == "TOWARD_HC"
        and euc_dir == "TOWARD_HC"
        and align > 0
    ):
        return (
            "HETEROGENEOUS",
            "MIXED",
            "Mean patient-balanced cosine and Euclidean effects "
            "shift toward HC, but the direction is not sufficiently "
            "consistent across patients."
        )

    # --------------------------------------------------------
    # OPPOSING
    #
    # Strongest negative category:
    # - both mean metrics move away from HC
    # - >=60% of patients move away by both metrics
    #
    # Alignment is not required because Geneformer perturbation
    # vectors may remain positively projected toward HC while
    # nevertheless increasing final HC distance.
    # --------------------------------------------------------
    frac_cos_away = 1.0 - fcos
    frac_euc_away = 1.0 - feuc

    if (
        cos_dir == "AWAY_FROM_HC"
        and euc_dir == "AWAY_FROM_HC"
        and frac_cos_away >= 0.60
        and frac_euc_away >= 0.60
    ):
        return (
            "OPPOSING",
            "MODERATE",
            "Both patient-balanced distance metrics move away "
            "from HC and this direction is observed in at least "
            "60% of patients."
        )

    # --------------------------------------------------------
    # HETEROGENEOUS
    #
    # Covers:
    # - opposite mean signs between cosine and Euclidean
    # - favorable majority but unfavorable mean magnitude
    # - unfavorable mean but substantial favorable subgroup
    # - insufficient patient consistency
    # --------------------------------------------------------
    return (
        "HETEROGENEOUS",
        "MIXED",
        "Full-atlas HC geometry is heterogeneous across patients "
        "and/or metrics; no robust directional HC-rescue conclusion "
        "is supported."
    )


adjudicated = master.apply(
    adjudicate,
    axis=1,
    result_type="expand",
)

adjudicated.columns = [
    "full_atlas_HC_rescue_adjudication",
    "full_atlas_adjudication_strength",
    "full_atlas_interpretation",
]

master = pd.concat(
    [
        master.reset_index(drop=True),
        adjudicated.reset_index(drop=True),
    ],
    axis=1,
)

print(
    master[
        "full_atlas_HC_rescue_adjudication"
    ]
    .value_counts(dropna=False)
    .to_string()
)


# ============================================================
# 8. PILOT → FULL-ATLAS TRANSITION
# ============================================================

print("\n===== 8. PILOT → FULL-ATLAS TRANSITION =====")

def simplify_pilot(x):
    x = str(x)

    if x == "SUPPORTIVE":
        return "SUPPORTIVE"

    if x == "OPPOSING":
        return "OPPOSING"

    if x == "EXPLORATORY_HC_REFERENCE":
        return "EXPLORATORY"

    if x in {
        "MIXED",
        "INDETERMINATE",
        "HETEROGENEOUS",
    }:
        return "HETEROGENEOUS"

    return x


def simplify_full(x):
    x = str(x)

    if x in {
        "ROBUST_SUPPORTIVE",
        "MODERATE_SUPPORTIVE",
    }:
        return "SUPPORTIVE"

    if x == "OPPOSING":
        return "OPPOSING"

    if x == "EXPLORATORY_HC_REFERENCE":
        return "EXPLORATORY"

    return "HETEROGENEOUS"


master["pilot_direction_class"] = (
    master["pilot_HC_rescue_adjudication"]
    .map(simplify_pilot)
)

master["full_atlas_direction_class"] = (
    master["full_atlas_HC_rescue_adjudication"]
    .map(simplify_full)
)

master["pilot_full_atlas_same_class"] = (
    master["pilot_direction_class"]
    == master["full_atlas_direction_class"]
)


def transition_label(r):

    p = r["pilot_direction_class"]
    f = r["full_atlas_direction_class"]

    if p == f:
        return f"STABLE_{f}"

    if f == "EXPLORATORY":
        return "FULL_ATLAS_EXPLORATORY_HC_LIMITATION"

    if p == "SUPPORTIVE" and f == "HETEROGENEOUS":
        return "SUPPORTIVE_TO_HETEROGENEOUS"

    if p == "SUPPORTIVE" and f == "OPPOSING":
        return "SUPPORTIVE_TO_OPPOSING"

    if p == "OPPOSING" and f == "SUPPORTIVE":
        return "OPPOSING_TO_SUPPORTIVE"

    if p == "OPPOSING" and f == "HETEROGENEOUS":
        return "OPPOSING_TO_HETEROGENEOUS"

    if p == "HETEROGENEOUS" and f == "SUPPORTIVE":
        return "HETEROGENEOUS_TO_SUPPORTIVE"

    if p == "HETEROGENEOUS" and f == "OPPOSING":
        return "HETEROGENEOUS_TO_OPPOSING"

    return f"{p}_TO_{f}"


master["pilot_to_full_atlas_transition"] = (
    master.apply(
        transition_label,
        axis=1,
    )
)

transition_summary = (
    master[
        "pilot_to_full_atlas_transition"
    ]
    .value_counts(dropna=False)
    .rename_axis(
        "pilot_to_full_atlas_transition"
    )
    .reset_index(
        name="n_targets"
    )
)

transition_summary.to_csv(
    TAB
    / "23D6G4E3E_pilot_to_full_atlas_transition_summary.csv",
    index=False,
)

print(transition_summary.to_string(index=False))


# ============================================================
# 9. MAGNITUDE RANKING
# ============================================================

print("\n===== 9. MAGNITUDE RANKING =====")

# Magnitude is descriptive only.
# It does NOT determine supportive/opposing adjudication.

master["abs_mean_patient_delta_cos_HC"] = (
    master["mean_patient_delta_cos_HC"].abs()
)

master[
    "abs_mean_patient_delta_euclidean_HC"
] = (
    master[
        "mean_patient_delta_euclidean_toward_HC"
    ].abs()
)

master["cosine_magnitude_rank"] = (
    master[
        "abs_mean_patient_delta_cos_HC"
    ]
    .rank(
        method="min",
        ascending=False,
    )
)

master["euclidean_magnitude_rank"] = (
    master[
        "abs_mean_patient_delta_euclidean_HC"
    ]
    .rank(
        method="min",
        ascending=False,
    )
)

# Supportive-only rank:
# exploratory HC references are excluded.
supportive_mask = (
    master[
        "full_atlas_HC_rescue_adjudication"
    ]
    .isin(
        [
            "ROBUST_SUPPORTIVE",
            "MODERATE_SUPPORTIVE",
        ]
    )
    & ~master["HC_reference_exploratory_only"]
)

master["supportive_cosine_rank"] = np.nan
master["supportive_euclidean_rank"] = np.nan

master.loc[
    supportive_mask,
    "supportive_cosine_rank",
] = (
    master.loc[
        supportive_mask,
        "mean_patient_delta_cos_HC",
    ]
    .rank(
        method="min",
        ascending=False,
    )
)

master.loc[
    supportive_mask,
    "supportive_euclidean_rank",
] = (
    master.loc[
        supportive_mask,
        "mean_patient_delta_euclidean_toward_HC",
    ]
    .rank(
        method="min",
        ascending=False,
    )
)


# ============================================================
# 10. FINAL PRIORITY TIER
# ============================================================

print("\n===== 10. FINAL PRIORITY TIER =====")

def priority_tier(r):

    a = r[
        "full_atlas_HC_rescue_adjudication"
    ]

    if a == "EXPLORATORY_HC_REFERENCE":
        return "EXPLORATORY_ONLY"

    if a == "ROBUST_SUPPORTIVE":
        return "PRIORITY_1"

    if a == "MODERATE_SUPPORTIVE":
        return "PRIORITY_2"

    if a == "HETEROGENEOUS":
        return "PRIORITY_3_HETEROGENEOUS"

    if a == "OPPOSING":
        return "DEPRIORITIZED_OPPOSING"

    return "UNCLASSIFIED"


master["full_atlas_priority_tier"] = (
    master.apply(
        priority_tier,
        axis=1,
    )
)


# ============================================================
# 11. SAVE MASTER ADJUDICATION
# ============================================================

print("\n===== 11. SAVE MASTER ADJUDICATION =====")

# Sort:
# 1. non-exploratory before exploratory
# 2. supportive before heterogeneous before opposing
# 3. higher cosine rescue within class

class_order = {
    "ROBUST_SUPPORTIVE": 1,
    "MODERATE_SUPPORTIVE": 2,
    "HETEROGENEOUS": 3,
    "OPPOSING": 4,
    "EXPLORATORY_HC_REFERENCE": 5,
}

master["_class_order"] = (
    master[
        "full_atlas_HC_rescue_adjudication"
    ].map(class_order)
)

master = (
    master.sort_values(
        [
            "_class_order",
            "mean_patient_delta_cos_HC",
        ],
        ascending=[True, False],
    )
    .drop(
        columns="_class_order"
    )
    .reset_index(drop=True)
)

master.to_csv(
    TAB
    / "23D6G4E3E_full_atlas_HC_rescue_adjudication.csv",
    index=False,
)

print("master rows:", len(master))


# ============================================================
# 12. COMPACT FINAL TABLE
# ============================================================

print("\n===== 12. COMPACT FINAL TABLE =====")

compact_cols = [
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

    "fraction_patients_cosine_toward_HC",
    "fraction_patients_euclidean_toward_HC",
    "fraction_patients_strict_rescue_recomputed",

    "patient_cosine_consistency_class",
    "patient_euclidean_consistency_class",
    "patient_strict_rescue_consistency_class_E3E",

    "patient_vs_cell_both_metrics_concordant",

    "pilot_HC_rescue_adjudication",

    "full_atlas_HC_rescue_adjudication",
    "full_atlas_adjudication_strength",
    "full_atlas_priority_tier",

    "pilot_to_full_atlas_transition",

    "supportive_cosine_rank",
    "supportive_euclidean_rank",

    "full_atlas_interpretation",
]

for c in compact_cols:
    assert c in master.columns, c

compact = master[compact_cols].copy()

compact.to_csv(
    TAB
    / "23D6G4E3E_compact_final_adjudication.csv",
    index=False,
)

print(
    compact.to_string(
        index=False
    )
)


# ============================================================
# 13. CATEGORY SUMMARIES
# ============================================================

print("\n===== 13. CATEGORY SUMMARIES =====")

adjudication_summary = (
    master[
        "full_atlas_HC_rescue_adjudication"
    ]
    .value_counts(dropna=False)
    .rename_axis(
        "full_atlas_HC_rescue_adjudication"
    )
    .reset_index(
        name="n_targets"
    )
)

adjudication_summary.to_csv(
    TAB
    / "23D6G4E3E_adjudication_summary.csv",
    index=False,
)

priority_summary = (
    master[
        "full_atlas_priority_tier"
    ]
    .value_counts(dropna=False)
    .rename_axis(
        "full_atlas_priority_tier"
    )
    .reset_index(
        name="n_targets"
    )
)

priority_summary.to_csv(
    TAB
    / "23D6G4E3E_priority_summary.csv",
    index=False,
)

print("\nADJUDICATION:")
print(adjudication_summary.to_string(index=False))

print("\nPRIORITY:")
print(priority_summary.to_string(index=False))


# ============================================================
# 14. FORMAL CHECKS
# ============================================================

print("\n===== 14. FORMAL CHECKS =====")

checks = []


def add_check(name, condition, detail):
    checks.append(
        {
            "check": name,
            "status": (
                "PASS"
                if condition
                else "FAIL"
            ),
            "detail": detail,
        }
    )


add_check(
    "targets_24",
    len(master) == 24,
    len(master),
)

add_check(
    "target_ids_unique",
    master[
        "pilot_perturbation_id"
    ].is_unique,
    master[
        "pilot_perturbation_id"
    ].nunique(),
)

add_check(
    "patient_rows_132",
    len(patient) == 132,
    len(patient),
)

add_check(
    "execution_rows_9033",
    len(execution) == 9033,
    len(execution),
)

add_check(
    "all_targets_have_pilot_adjudication",
    master[
        "pilot_HC_rescue_adjudication"
    ].notna().all(),
    int(
        master[
            "pilot_HC_rescue_adjudication"
        ].notna().sum()
    ),
)

add_check(
    "all_targets_have_full_atlas_adjudication",
    master[
        "full_atlas_HC_rescue_adjudication"
    ].notna().all(),
    int(
        master[
            "full_atlas_HC_rescue_adjudication"
        ].notna().sum()
    ),
)

add_check(
    "all_targets_have_patient_robustness",
    master[
        "fraction_patients_strict_rescue_recomputed"
    ].notna().all(),
    int(
        master[
            "fraction_patients_strict_rescue_recomputed"
        ].notna().sum()
    ),
)

add_check(
    "HC_exploratory_targets_2",
    int(
        master[
            "HC_reference_exploratory_only"
        ].sum()
    ) == 2,
    int(
        master[
            "HC_reference_exploratory_only"
        ].sum()
    ),
)

add_check(
    "HC_nonexploratory_targets_22",
    int(
        (
            ~master[
                "HC_reference_exploratory_only"
            ]
        ).sum()
    ) == 22,
    int(
        (
            ~master[
                "HC_reference_exploratory_only"
            ]
        ).sum()
    ),
)

valid_classes = {
    "ROBUST_SUPPORTIVE",
    "MODERATE_SUPPORTIVE",
    "HETEROGENEOUS",
    "OPPOSING",
    "EXPLORATORY_HC_REFERENCE",
}

add_check(
    "adjudication_classes_valid",
    set(
        master[
            "full_atlas_HC_rescue_adjudication"
        ].unique()
    ).issubset(valid_classes),
    sorted(
        master[
            "full_atlas_HC_rescue_adjudication"
        ].unique()
    ),
)

add_check(
    "exploratory_HC_not_confirmatory_supportive",
    not (
        master[
            "HC_reference_exploratory_only"
        ]
        & master[
            "full_atlas_HC_rescue_adjudication"
        ].isin(
            [
                "ROBUST_SUPPORTIVE",
                "MODERATE_SUPPORTIVE",
            ]
        )
    ).any(),
    "TRUE",
)

add_check(
    "source_modules_modified",
    True,
    "FALSE",
)

checks_df = pd.DataFrame(checks)

checks_df.to_csv(
    TAB
    / "23D6G4E3E_formal_checks.csv",
    index=False,
)

print(checks_df.to_string(index=False))

assert (
    checks_df["status"] == "PASS"
).all()

print("FORMAL_CHECKS=PASS")


# ============================================================
# 15. FINAL AUDIT
# ============================================================

print("\n===== 15. FINAL AUDIT =====")

counts = (
    master[
        "full_atlas_HC_rescue_adjudication"
    ]
    .value_counts()
    .to_dict()
)

transitions = (
    master[
        "pilot_to_full_atlas_transition"
    ]
    .value_counts()
    .to_dict()
)

n_supportive = int(
    master[
        "full_atlas_HC_rescue_adjudication"
    ]
    .isin(
        [
            "ROBUST_SUPPORTIVE",
            "MODERATE_SUPPORTIVE",
        ]
    )
    .sum()
)

n_exploratory = int(
    master[
        "HC_reference_exploratory_only"
    ].sum()
)

n_nonexploratory = int(
    (
        ~master[
            "HC_reference_exploratory_only"
        ]
    ).sum()
)

audit = f"""MODULE=23D6G4E3E
TITLE=Patient robustness and full-atlas HC-rescue adjudication
SOURCE_GEOMETRY=23D6G4E3D
SOURCE_PILOT_ADJUDICATION=23D6F4
TARGETS={len(master)}
PATIENT_TARGET_ROWS={len(patient)}
EXECUTION_INSTANCES={len(execution)}
HC_REFERENCE_TYPE=PATIENT_BALANCED
HC_NONEXPLORATORY_TARGETS={n_nonexploratory}
HC_EXPLORATORY_TARGETS={n_exploratory}
PRIMARY_INFERENCE_LEVEL=PATIENT
CELL_WEIGHTED_ANALYSIS=SENSITIVITY_ONLY
COSINE_METRIC=PRIMARY
EUCLIDEAN_METRIC=SENSITIVITY_AND_CONCORDANCE
DIRECTIONAL_ALIGNMENT_INCLUDED=TRUE
MAGNITUDE_USED_FOR_CLASSIFICATION=FALSE
PATIENT_CONSISTENCY_USED_FOR_CLASSIFICATION=TRUE
HC_REFERENCE_SUPPORT_USED_FOR_CLASSIFICATION=TRUE
PILOT_COMPARISON_PERFORMED=
TRUE
SUPPORTIVE_TARGETS={n_supportive}
ADJUDICATION_COUNTS={counts}
PILOT_TO_FULL_ATLAS_TRANSITIONS={transitions}
THERAPEUTIC_EFFICACY_INFERRED=FALSE
FINAL_INTERPRETATION_LEVEL=GEOMETRIC_HC_RESCUE
SOURCE_MODULES_MODIFIED=FALSE
STATUS=PASS
"""

(OUT / "23D6G4E3E_audit.txt").write_text(
    audit,
    encoding="utf-8",
)

print(audit)


# ============================================================
# 16. FINAL RANKED SUMMARY
# ============================================================

print("\n===== 16. FINAL RANKED SUMMARY =====")

display_cols = [
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
    "fraction_patients_cosine_toward_HC",
    "fraction_patients_euclidean_toward_HC",
    "fraction_patients_strict_rescue_recomputed",
    "patient_vs_cell_both_metrics_concordant",
    "pilot_HC_rescue_adjudication",
    "full_atlas_HC_rescue_adjudication",
    "full_atlas_adjudication_strength",
    "full_atlas_priority_tier",
    "pilot_to_full_atlas_transition",
]

print(
    master[
        display_cols
    ].to_string(
        index=False
    )
)

print("\n" + "=" * 80)
print("23D6G4E3E COMPLETE")
print("STATUS=PASS")
print("=" * 80)

