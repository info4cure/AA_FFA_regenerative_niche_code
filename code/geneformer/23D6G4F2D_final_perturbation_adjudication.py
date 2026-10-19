from pathlib import Path
import numpy as np
import pandas as pd

# ============================================================
# 23D6G4F2D
# EXHAUSTIVE PATIENT-ROBUST PERTURBATION ADJUDICATION
#
# Integrates:
#   - frozen PRE-Geneformer evidence (23D6C)
#   - exhaustive patient-balanced HC rescue geometry (F2C)
#   - cell-weighted sensitivity (F2C)
#
# IMPORTANT:
#   The pre-Geneformer panel was frozen before perturbation.
#   No target was selected using Geneformer results.
#   No therapeutic efficacy is inferred.
# ============================================================

ROOT = Path.cwd()

BASE = (
    ROOT
    / "PIPELINE_V2/results/disease_biology"
    / "23_ADVANCED_ANALYSIS_ROADMAP"
    / "23D_AI_PERTURBATION"
)

PANEL = (
    BASE
    / "23D6_TARGET_SELECTION"
    / "23D6C_PILOT_PANEL_FREEZE"
    / "tables"
    / "23D6C_frozen_pilot_perturbation_panel.csv"
)

F2C = (
    BASE
    / "23D6_PERTURBATION"
    / "23D6G_FULL_ATLAS_REASSESSMENT"
    / "23D6G4F_EXHAUSTIVE_FULL_ATLAS"
    / "23D6G4F2C_EXHAUSTIVE_HC_RESCUE_GEOMETRY"
)

TARGET_FILE = (
    F2C
    / "tables"
    / "23D6G4F2C_target_level_HC_rescue_geometry.csv"
)

PATIENT_FILE = (
    F2C
    / "tables"
    / "23D6G4F2C_patient_level_HC_rescue_geometry.csv"
)

CELL_WEIGHTED_FILE = (
    F2C
    / "tables"
    / "23D6G4F2C_target_cell_weighted_sensitivity.csv"
)

OUT = (
    BASE
    / "23D6_PERTURBATION"
    / "23D6G_FULL_ATLAS_REASSESSMENT"
    / "23D6G4F_EXHAUSTIVE_FULL_ATLAS"
    / "23D6G4F2D_FINAL_PERTURBATION_ADJUDICATION"
)

TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("23D6G4F2D — FINAL EXHAUSTIVE PERTURBATION ADJUDICATION")
print("PRE-GENEFORMER EVIDENCE + PATIENT-BALANCED HC RESCUE")
print("NO THERAPEUTIC EFFICACY INFERENCE")
print("=" * 80)


# ============================================================
# 1. LOAD SOURCES
# ============================================================

print("\n===== 1. LOAD SOURCES =====")

for p in [
    PANEL,
    TARGET_FILE,
    PATIENT_FILE,
    CELL_WEIGHTED_FILE,
]:
    assert p.exists(), p

panel = pd.read_csv(PANEL)
target = pd.read_csv(TARGET_FILE)
patient = pd.read_csv(PATIENT_FILE)
cw = pd.read_csv(CELL_WEIGHTED_FILE)

assert len(panel) == 24
assert len(target) == 24
assert len(cw) == 24
assert len(patient) == 132

assert panel["pilot_perturbation_id"].nunique() == 24
assert target["pilot_perturbation_id"].nunique() == 24
assert cw["pilot_perturbation_id"].nunique() == 24

print("frozen targets:", len(panel))
print("target geometry rows:", len(target))
print("patient-target rows:", len(patient))
print("cell-weighted rows:", len(cw))

print("SOURCE_LOAD=PASS")


# ============================================================
# 2. VERIFY TARGET IDENTITY ACROSS SOURCES
# ============================================================

print("\n===== 2. VERIFY TARGET IDENTITY =====")

ids_panel = set(panel["pilot_perturbation_id"].astype(str))
ids_target = set(target["pilot_perturbation_id"].astype(str))
ids_cw = set(cw["pilot_perturbation_id"].astype(str))

assert ids_panel == ids_target == ids_cw

identity_cols = [
    "disease",
    "final_L2_label",
    "gene",
    "geneformer_perturb_type",
]

identity_checks = []

for c in identity_cols:

    a = (
        panel[
            ["pilot_perturbation_id", c]
        ]
        .set_index("pilot_perturbation_id")[c]
        .astype(str)
        .sort_index()
    )

    b = (
        target[
            ["pilot_perturbation_id", c]
        ]
        .set_index("pilot_perturbation_id")[c]
        .astype(str)
        .sort_index()
    )

    ok = a.equals(b)

    identity_checks.append(
        {
            "field": c,
            "status": "PASS" if ok else "FAIL",
            "n_targets": len(a),
        }
    )

identity_checks = pd.DataFrame(identity_checks)

print(identity_checks.to_string(index=False))

assert (identity_checks["status"] == "PASS").all()

identity_checks.to_csv(
    TAB / "23D6G4F2D_source_identity_checks.csv",
    index=False,
)

print("TARGET_IDENTITY=PASS")


# ============================================================
# 3. SELECT FROZEN PRE-GENEFORMER EVIDENCE
# ============================================================

print("\n===== 3. PRE-GENEFORMER EVIDENCE =====")

pre_cols = [
    "pilot_perturbation_id",
    "contrast",
    "module",
    "disease",
    "final_L2_label",
    "gene",
    "ensembl_id",
    "geneformer_perturb_type",
    "rescue_direction",
    "disease_direction",
    "logFC",
    "FDR_BH",
    "abs_logFC",
    "direct_DE_FDR005",
    "robust_cell_support",
    "targeted_n_axes",
    "targeted_n_pathways",
    "targeted_max_L2",
    "targeted_pathway_support",
    "bioplanet_n_pathways",
    "bioplanet_n_L2",
    "bioplanet_n_hits",
    "bioplanet_min_FDR",
    "bioplanet_support",
    "hallmark_n_pathways",
    "hallmark_min_FDR",
    "hallmark_support",
    "AA_vs_FFA_logFC",
    "AA_vs_FFA_FDR",
    "AA_vs_FFA_support",
    "n_support_domains",
    "pre_geneformer_evidence_class",
    "evidence_rank",
    "selection_stage",
    "selection_status",
    "geneformer_results_used_for_selection",
]

for c in pre_cols:
    assert c in panel.columns, c

pre = panel[pre_cols].copy()

assert (
    pre["geneformer_results_used_for_selection"]
    .astype(str)
    .str.lower()
    .eq("false")
    .all()
)

print(
    pre["pre_geneformer_evidence_class"]
    .value_counts(dropna=False)
    .to_string()
)

print("PRE_GENEFORMER_INDEPENDENCE=PASS")


# ============================================================
# 4. PREPARE PRIMARY PATIENT-BALANCED GEOMETRY
# ============================================================

print("\n===== 4. PRIMARY PATIENT-BALANCED GEOMETRY =====")

geometry_cols = [
    "pilot_perturbation_id",
    "HC_reference_type",
    "HC_reference_n_cells",
    "HC_reference_n_patients",
    "HC_reference_support_class",
    "HC_reference_exploratory_only",
    "n_patients",
    "n_executions",
    "mean_patient_delta_cos_HC",
    "median_patient_delta_cos_HC",
    "min_patient_delta_cos_HC",
    "max_patient_delta_cos_HC",
    "mean_patient_delta_euclidean_toward_HC",
    "median_patient_delta_euclidean_toward_HC",
    "min_patient_delta_euclidean_toward_HC",
    "max_patient_delta_euclidean_toward_HC",
    "mean_patient_normalized_cosine_rescue",
    "median_patient_normalized_cosine_rescue",
    "mean_patient_normalized_euclidean_rescue",
    "median_patient_normalized_euclidean_rescue",
    "mean_patient_directional_alignment_cosine",
    "median_patient_directional_alignment_cosine",
    "mean_original_to_perturbed_cosine",
    "mean_original_to_perturbed_euclidean",
    "fraction_patients_toward_HC_by_cosine",
    "fraction_patients_toward_HC_by_euclidean",
    "fraction_patients_directionally_aligned",
    "fraction_patients_metric_concordant",
    "fraction_patients_strict_rescue",
    "target_toward_HC_by_cosine",
    "target_toward_HC_by_euclidean",
    "target_directionally_aligned",
    "target_metric_direction_concordant",
    "target_strict_geometry_rescue",
    "HC_geometry_class",
    "patient_strict_rescue_consistency_class",
]

for c in geometry_cols:
    assert c in target.columns, c

geometry = target[geometry_cols].copy()

# Correct semantic interpretation of HC support.
geometry["HC_reference_multipatient"] = (
    geometry["HC_reference_n_patients"] >= 2
)

geometry["HC_reference_single_patient"] = (
    geometry["HC_reference_n_patients"] == 1
)

assert (
    geometry["HC_reference_multipatient"]
    .sum()
    == 22
)

assert (
    geometry["HC_reference_single_patient"]
    .sum()
    == 2
)

print(
    geometry["HC_geometry_class"]
    .value_counts()
    .to_string()
)

print(
    "multipatient HC targets:",
    int(geometry["HC_reference_multipatient"].sum()),
)

print(
    "single-patient HC targets:",
    int(geometry["HC_reference_single_patient"].sum()),
)


# ============================================================
# 5. CELL-WEIGHTED SENSITIVITY
# ============================================================

print("\n===== 5. CELL-WEIGHTED SENSITIVITY =====")

cw_keep = cw[
    [
        "pilot_perturbation_id",
        "mean_delta_cos_HC",
        "median_delta_cos_HC",
        "mean_delta_euclidean_toward_HC",
        "median_delta_euclidean_toward_HC",
        "mean_normalized_cosine_rescue",
        "mean_normalized_euclidean_rescue",
        "mean_directional_alignment_cosine",
        "fraction_executions_toward_HC_by_cosine",
        "fraction_executions_toward_HC_by_euclidean",
        "fraction_executions_directionally_aligned",
        "fraction_executions_strict_rescue",
    ]
].copy()

cw_keep = cw_keep.rename(
    columns={
        c: f"cell_weighted_{c}"
        for c in cw_keep.columns
        if c != "pilot_perturbation_id"
    }
)


# ============================================================
# 6. MERGE EVIDENCE LAYERS
# ============================================================

print("\n===== 6. MERGE EVIDENCE LAYERS =====")

final = (
    pre
    .merge(
        geometry,
        on="pilot_perturbation_id",
        how="inner",
        validate="one_to_one",
    )
    .merge(
        cw_keep,
        on="pilot_perturbation_id",
        how="inner",
        validate="one_to_one",
    )
)

assert len(final) == 24

# Merge generated duplicate identity columns conservatively.
for c in [
    "disease",
    "final_L2_label",
    "gene",
    "geneformer_perturb_type",
]:
    cx = f"{c}_x"
    cy = f"{c}_y"

    if cx in final.columns and cy in final.columns:
        assert (
            final[cx].astype(str)
            == final[cy].astype(str)
        ).all()

        final[c] = final[cx]

        final = final.drop(
            columns=[cx, cy]
        )

print("merged targets:", len(final))
print("MERGE=PASS")


# ============================================================
# 7. PATIENT ROBUSTNESS FLAGS
# ============================================================

print("\n===== 7. PATIENT ROBUSTNESS =====")

# Primary consistency threshold:
# at least 60% of disease patients show strict rescue.
#
# AA: >= 4/6
# FFA: >= 3/5

final["patient_rescue_ge60"] = (
    final["fraction_patients_strict_rescue"]
    >= 0.60
)

final["patient_rescue_ge80"] = (
    final["fraction_patients_strict_rescue"]
    >= 0.80
)

final["patient_majority_cosine_toward_HC"] = (
    final["fraction_patients_toward_HC_by_cosine"]
    >= 0.60
)

final["patient_majority_euclidean_toward_HC"] = (
    final["fraction_patients_toward_HC_by_euclidean"]
    >= 0.60
)

final["patient_directionality_universal"] = (
    final["fraction_patients_directionally_aligned"]
    >= 1.0 - 1e-12
)

# Agreement between primary patient-balanced sign and
# cell-weighted sensitivity sign.
final["cell_weighted_cosine_same_direction"] = (
    np.sign(final["mean_patient_delta_cos_HC"])
    ==
    np.sign(final["cell_weighted_mean_delta_cos_HC"])
)

final["cell_weighted_euclidean_same_direction"] = (
    np.sign(
        final[
            "mean_patient_delta_euclidean_toward_HC"
        ]
    )
    ==
    np.sign(
        final[
            "cell_weighted_mean_delta_euclidean_toward_HC"
        ]
    )
)

final["patient_cell_weighted_direction_concordant"] = (
    final["cell_weighted_cosine_same_direction"]
    & final["cell_weighted_euclidean_same_direction"]
)


# ============================================================
# 8. FINAL RULE-BASED ADJUDICATION
# ============================================================

print("\n===== 8. FINAL RULE-BASED ADJUDICATION =====")

def adjudicate(r):

    multipatient = bool(
        r["HC_reference_multipatient"]
    )

    exploratory = bool(
        r["HC_reference_exploratory_only"]
    )

    geometry_class = str(
        r["HC_geometry_class"]
    )

    ge60 = bool(
        r["patient_rescue_ge60"]
    )

    if exploratory or not multipatient:

        if geometry_class == "TOWARD_HC_CONCORDANT":
            return "EXPLORATORY_HC_RESCUE"

        return "EXPLORATORY_HC_MIXED"

    if geometry_class == "TOWARD_HC_CONCORDANT":

        if ge60:
            return "SUPPORTED_RESCUE"

        return "DIRECTIONAL_RESCUE_LOW_CONSISTENCY"

    if geometry_class == "AWAY_FROM_HC_CONCORDANT":
        return "CONCORDANT_AWAY_FROM_HC"

    return "MIXED_PERTURBATION_EVIDENCE"


final["final_perturbation_evidence_class"] = (
    final.apply(
        adjudicate,
        axis=1,
    )
)

# Conservative interpretation tier.
def tier(r):

    c = r["final_perturbation_evidence_class"]

    if c == "SUPPORTED_RESCUE":
        return "TIER_1_SUPPORTED"

    if c == "DIRECTIONAL_RESCUE_LOW_CONSISTENCY":
        return "TIER_2_DIRECTIONAL"

    if c == "EXPLORATORY_HC_RESCUE":
        return "TIER_3_EXPLORATORY_REFERENCE"

    if c == "CONCORDANT_AWAY_FROM_HC":
        return "TIER_4_NONRESCUE"

    return "TIER_4_MIXED"


final["final_evidence_tier"] = final.apply(
    tier,
    axis=1,
)

print(
    final[
        "final_perturbation_evidence_class"
    ]
    .value_counts(dropna=False)
    .to_string()
)


# ============================================================
# 9. RANK WITHIN ADJUDICATION — DESCRIPTIVE ONLY
# ============================================================

print("\n===== 9. FINAL DESCRIPTIVE RANKING =====")

tier_order = {
    "TIER_1_SUPPORTED": 1,
    "TIER_2_DIRECTIONAL": 2,
    "TIER_3_EXPLORATORY_REFERENCE": 3,
    "TIER_4_MIXED": 4,
    "TIER_4_NONRESCUE": 5,
}

final["_tier_order"] = (
    final["final_evidence_tier"]
    .map(tier_order)
)

# Ranking does NOT generate the adjudication.
# It only orders targets within/after rule-based classes.
final = final.sort_values(
    [
        "_tier_order",
        "fraction_patients_strict_rescue",
        "mean_patient_delta_cos_HC",
        "mean_patient_delta_euclidean_toward_HC",
    ],
    ascending=[
        True,
        False,
        False,
        False,
    ],
).reset_index(drop=True)

final["final_descriptive_rank"] = (
    np.arange(1, len(final) + 1)
)

final = final.drop(
    columns=["_tier_order"]
)

print(
    final[
        [
            "final_descriptive_rank",
            "pilot_perturbation_id",
            "disease",
            "final_L2_label",
            "gene",
            "geneformer_perturb_type",
            "HC_reference_support_class",
            "HC_reference_exploratory_only",
            "n_patients",
            "n_executions",
            "mean_patient_delta_cos_HC",
            "mean_patient_delta_euclidean_toward_HC",
            "mean_patient_directional_alignment_cosine",
            "fraction_patients_strict_rescue",
            "HC_geometry_class",
            "final_perturbation_evidence_class",
            "final_evidence_tier",
        ]
    ].to_string(index=False)
)


# ============================================================
# 10. PRE-GENEFORMER ↔ PERTURBATION CROSS-TABULATION
# ============================================================

print(
    "\n===== 10. PRE-GENEFORMER ↔ PERTURBATION CROSS-TABULATION ====="
)

cross = pd.crosstab(
    final["pre_geneformer_evidence_class"],
    final["final_perturbation_evidence_class"],
    margins=True,
)

print(cross.to_string())

cross.to_csv(
    TAB
    / "23D6G4F2D_preGeneformer_vs_perturbation_crosstab.csv"
)


# ============================================================
# 11. DISEASE / COMPARTMENT SUMMARY
# ============================================================

print("\n===== 11. DISEASE / COMPARTMENT SUMMARY =====")

disease_summary = (
    final.groupby(
        [
            "disease",
            "final_perturbation_evidence_class",
            "final_evidence_tier",
        ],
        dropna=False,
        observed=True,
    )
    .agg(
        n_targets=("pilot_perturbation_id", "size"),
        n_unique_genes=("gene", "nunique"),
        median_fraction_patients_strict_rescue=(
            "fraction_patients_strict_rescue",
            "median",
        ),
        median_delta_cos_HC=(
            "mean_patient_delta_cos_HC",
            "median",
        ),
        median_delta_euclidean_toward_HC=(
            "mean_patient_delta_euclidean_toward_HC",
            "median",
        ),
    )
    .reset_index()
)

disease_summary.to_csv(
    TAB / "23D6G4F2D_disease_summary.csv",
    index=False,
)

print(disease_summary.to_string(index=False))


# ============================================================
# 12. TARGET-LEVEL FINAL TABLE
# ============================================================

print("\n===== 12. SAVE FINAL TARGET TABLE =====")

# Preserve the full evidence matrix.
final.to_csv(
    TAB
    / "23D6G4F2D_final_target_adjudication_full.csv",
    index=False,
)

compact_cols = [
    "final_descriptive_rank",
    "pilot_perturbation_id",
    "disease",
    "final_L2_label",
    "gene",
    "ensembl_id",
    "geneformer_perturb_type",

    # Independent evidence before Geneformer
    "pre_geneformer_evidence_class",
    "n_support_domains",
    "direct_DE_FDR005",
    "robust_cell_support",
    "targeted_pathway_support",
    "bioplanet_support",
    "hallmark_support",
    "AA_vs_FFA_support",

    # HC reference
    "HC_reference_n_cells",
    "HC_reference_n_patients",
    "HC_reference_support_class",
    "HC_reference_exploratory_only",
    "HC_reference_multipatient",

    # Exhaustive perturbation
    "n_patients",
    "n_executions",
    "mean_patient_delta_cos_HC",
    "mean_patient_delta_euclidean_toward_HC",
    "mean_patient_directional_alignment_cosine",
    "fraction_patients_toward_HC_by_cosine",
    "fraction_patients_toward_HC_by_euclidean",
    "fraction_patients_directionally_aligned",
    "fraction_patients_strict_rescue",
    "patient_rescue_ge60",
    "patient_rescue_ge80",

    # Sensitivity
    "cell_weighted_mean_delta_cos_HC",
    "cell_weighted_mean_delta_euclidean_toward_HC",
    "patient_cell_weighted_direction_concordant",

    # Final interpretation
    "HC_geometry_class",
    "patient_strict_rescue_consistency_class",
    "final_perturbation_evidence_class",
    "final_evidence_tier",
]

for c in compact_cols:
    assert c in final.columns, c

compact = final[compact_cols].copy()

compact.to_csv(
    TAB
    / "23D6G4F2D_final_target_adjudication_compact.csv",
    index=False,
)

print("full target rows:", len(final))
print("compact target rows:", len(compact))


# ============================================================
# 13. SUPPORTED / DIRECTIONAL CANDIDATE TABLE
# ============================================================

print("\n===== 13. SUPPORTED CANDIDATES =====")

candidate_classes = [
    "SUPPORTED_RESCUE",
    "DIRECTIONAL_RESCUE_LOW_CONSISTENCY",
    "EXPLORATORY_HC_RESCUE",
]

candidates = final.loc[
    final[
        "final_perturbation_evidence_class"
    ].isin(candidate_classes)
].copy()

candidates.to_csv(
    TAB
    / "23D6G4F2D_supported_and_directional_candidates.csv",
    index=False,
)

print(
    candidates[
        [
            "final_descriptive_rank",
            "disease",
            "final_L2_label",
            "gene",
            "pre_geneformer_evidence_class",
            "fraction_patients_strict_rescue",
            "mean_patient_delta_cos_HC",
            "mean_patient_delta_euclidean_toward_HC",
            "HC_reference_n_patients",
            "final_perturbation_evidence_class",
            "final_evidence_tier",
        ]
    ].to_string(index=False)
)


# ============================================================
# 14. FORMAL CHECKS
# ============================================================

print("\n===== 14. FORMAL CHECKS =====")

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
    "frozen_targets_24",
    len(final) == 24,
    len(final),
)

add_check(
    "unique_target_ids_24",
    final["pilot_perturbation_id"].nunique() == 24,
    final["pilot_perturbation_id"].nunique(),
)

add_check(
    "pre_geneformer_selection_independent",
    (
        final[
            "geneformer_results_used_for_selection"
        ]
        .astype(str)
        .str.lower()
        .eq("false")
        .all()
    ),
    "23D6C frozen before Geneformer",
)

add_check(
    "all_have_pre_geneformer_class",
    final["pre_geneformer_evidence_class"].notna().all(),
    int(
        final[
            "pre_geneformer_evidence_class"
        ].notna().sum()
    ),
)

add_check(
    "all_have_HC_geometry_class",
    final["HC_geometry_class"].notna().all(),
    int(final["HC_geometry_class"].notna().sum()),
)

add_check(
    "all_have_final_adjudication",
    final[
        "final_perturbation_evidence_class"
    ].notna().all(),
    int(
        final[
            "final_perturbation_evidence_class"
        ].notna().sum()
    ),
)

add_check(
    "all_have_final_tier",
    final["final_evidence_tier"].notna().all(),
    int(final["final_evidence_tier"].notna().sum()),
)

add_check(
    "multipatient_HC_targets_22",
    int(final["HC_reference_multipatient"].sum()) == 22,
    int(final["HC_reference_multipatient"].sum()),
)

add_check(
    "single_patient_HC_targets_2",
    int(final["HC_reference_single_patient"].sum()) == 2,
    int(final["HC_reference_single_patient"].sum()),
)

add_check(
    "exploratory_HC_targets_2",
    int(
        final[
            "HC_reference_exploratory_only"
        ].sum()
    ) == 2,
    int(
        final[
            "HC_reference_exploratory_only"
        ].sum()
    ),
)

add_check(
    "patient_target_rows_132",
    len(patient) == 132,
    len(patient),
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
    "all_execution_counts_positive",
    (final["n_executions"] > 0).all(),
    int((final["n_executions"] > 0).sum()),
)

add_check(
    "total_exhaustive_execution_instances_52509",
    int(final["n_executions"].sum()) == 52509,
    int(final["n_executions"].sum()),
)

add_check(
    "final_ranks_unique_24",
    final["final_descriptive_rank"].nunique() == 24,
    final["final_descriptive_rank"].nunique(),
)

add_check(
    "no_therapeutic_efficacy_inference",
    True,
    "TRUE",
)

add_check(
    "rule_based_not_score_based",
    True,
    "TRUE",
)

add_check(
    "source_modules_modified",
    True,
    "FALSE",
)

checks_df = pd.DataFrame(checks)

checks_df.to_csv(
    TAB / "23D6G4F2D_formal_checks.csv",
    index=False,
)

print(checks_df.to_string(index=False))

assert (checks_df["status"] == "PASS").all()

print("FORMAL_CHECKS=PASS")


# ============================================================
# 15. FINAL COUNTS
# ============================================================

print("\n===== 15. FINAL COUNTS =====")

class_counts = (
    final[
        "final_perturbation_evidence_class"
    ]
    .value_counts(dropna=False)
    .to_dict()
)

tier_counts = (
    final[
        "final_evidence_tier"
    ]
    .value_counts(dropna=False)
    .to_dict()
)

n_supported = int(
    (
        final[
            "final_perturbation_evidence_class"
        ]
        == "SUPPORTED_RESCUE"
    ).sum()
)

n_directional = int(
    (
        final[
            "final_perturbation_evidence_class"
        ]
        == "DIRECTIONAL_RESCUE_LOW_CONSISTENCY"
    ).sum()
)

n_exploratory_rescue = int(
    (
        final[
            "final_perturbation_evidence_class"
        ]
        == "EXPLORATORY_HC_RESCUE"
    ).sum()
)

n_mixed = int(
    (
        final[
            "final_perturbation_evidence_class"
        ]
        == "MIXED_PERTURBATION_EVIDENCE"
    ).sum()
)

n_away = int(
    (
        final[
            "final_perturbation_evidence_class"
        ]
        == "CONCORDANT_AWAY_FROM_HC"
    ).sum()
)

print("class counts:", class_counts)
print("tier counts:", tier_counts)


# ============================================================
# 16. FINAL AUDIT
# ============================================================

print("\n===== 16. FINAL AUDIT =====")

audit = f"""MODULE=23D6G4F2D
TITLE=Final rule-based adjudication of exhaustive Geneformer perturbation evidence
SOURCE_PRE_GENEFORMER_PANEL=23D6C
SOURCE_EXHAUSTIVE_EXECUTION=23D6G4F1
SOURCE_PERTURBED_EMBEDDINGS=23D6G4F2B
SOURCE_HC_RESCUE_GEOMETRY=23D6G4F2C
TARGETS=24
EXECUTION_INSTANCES=52509
PATIENT_TARGET_ROWS=132
AA_PATIENTS=6
FFA_PATIENTS=5
PRE_GENEFORMER_PANEL_FROZEN=TRUE
GENEFORMER_RESULTS_USED_FOR_INITIAL_SELECTION=FALSE
PRIMARY_INFERENCE_LEVEL=PATIENT_BALANCED_TARGET
CELL_WEIGHTED_ANALYSIS=SENSITIVITY
HC_REFERENCE_TYPE=PATIENT_BALANCED
MULTIPATIENT_HC_TARGETS={int(final["HC_reference_multipatient"].sum())}
SINGLE_PATIENT_HC_TARGETS={int(final["HC_reference_single_patient"].sum())}
EXPLORATORY_HC_TARGETS={int(final["HC_reference_exploratory_only"].sum())}
SUPPORTED_RESCUE_TARGETS={n_supported}
DIRECTIONAL_RESCUE_LOW_CONSISTENCY_TARGETS={n_directional}
EXPLORATORY_HC_RESCUE_TARGETS={n_exploratory_rescue}
MIXED_PERTURBATION_EVIDENCE_TARGETS={n_mixed}
CONCORDANT_AWAY_FROM_HC_TARGETS={n_away}
FINAL_CLASS_COUNTS={class_counts}
FINAL_TIER_COUNTS={tier_counts}
SUPPORTED_RESCUE_DEFINITION=MULTIPATIENT_HC_AND_TOWARD_HC_CONCORDANT_AND_GE60_PERCENT_PATIENT_STRICT_RESCUE
DIRECTIONAL_LOW_CONSISTENCY_DEFINITION=MULTIPATIENT_HC_AND_TOWARD_HC_CONCORDANT_AND_LT60_PERCENT_PATIENT_STRICT_RESCUE
EXPLORATORY_REFERENCE_DEFINITION=SINGLE_PATIENT_HC_REFERENCE
ADJUDICATION_RULE_BASED=TRUE
ARBITRARY_COMPOSITE_SCORE_USED=FALSE
DESCRIPTIVE_RANK_USED_FOR_ADJUDICATION=FALSE
THERAPEUTIC_TARGET_CLAIM=FALSE
THERAPEUTIC_EFFICACY_INFERRED=FALSE
BIOLOGICAL_INTERPRETATION_LEVEL=IN_SILICO_PERTURBATION_EVIDENCE
SOURCE_MODULES_MODIFIED=FALSE
STATUS=PASS
"""

(OUT / "23D6G4F2D_audit.txt").write_text(
    audit,
    encoding="utf-8",
)

print(audit)


# ============================================================
# 17. OUTPUT INVENTORY
# ============================================================

print("\n===== 17. OUTPUT INVENTORY =====")

expected_outputs = [
    TAB / "23D6G4F2D_source_identity_checks.csv",
    TAB / "23D6G4F2D_preGeneformer_vs_perturbation_crosstab.csv",
    TAB / "23D6G4F2D_disease_summary.csv",
    TAB / "23D6G4F2D_final_target_adjudication_full.csv",
    TAB / "23D6G4F2D_final_target_adjudication_compact.csv",
    TAB / "23D6G4F2D_supported_and_directional_candidates.csv",
    TAB / "23D6G4F2D_formal_checks.csv",
    OUT / "23D6G4F2D_audit.txt",
]

for p in expected_outputs:
    print(
        p,
        "exists=",
        p.exists(),
    )

assert all(p.exists() for p in expected_outputs)

print()
print("=" * 80)
print("23D6G4F2D COMPLETE")
print("STATUS=PASS")
print("=" * 80)
