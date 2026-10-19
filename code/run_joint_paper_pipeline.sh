#!/bin/bash
set -euo pipefail

ATLAS="04_ANALYSIS/12_scrnaseq/04_integration/AA_FFA_HC_joint_Harmony_v0.1.rds"
ANNOT="04_ANALYSIS/12_scrnaseq/06_broad_annotation/AA_FFA_HC_joint_annotated_v0.1.rds"
EPI="04_ANALYSIS/12_scrnaseq/07_fine_annotation/AA_FFA_HC_epithelial_atlas_v0.1.rds"

echo "============================================================"
echo "AA–FFA–HC PAPER PIPELINE"
echo "============================================================"

# ------------------------------------------------------------
# 03
# ------------------------------------------------------------

if [ ! -s "$ATLAS" ]; then
    echo "03 atlas not available."
    echo "Run/finish 03 first."
    exit 1
fi

echo "03 PASS — joint atlas exists"

# ------------------------------------------------------------
# 04
# ------------------------------------------------------------

echo
echo "RUNNING 04 — broad annotation"

Rscript \
08_REPRODUCIBILITY/scripts/R/scrnaseq/04_annotate_joint_atlas.R \
2>&1 | tee \
08_REPRODUCIBILITY/logs/scrnaseq/04_annotate_joint_atlas.log

if [ ! -s "$ANNOT" ]; then
    echo "ERROR: 04 did not generate annotated atlas"
    exit 1
fi

echo "04 PASS"

# ------------------------------------------------------------
# 05
# ------------------------------------------------------------

echo
echo "RUNNING 05 — patient-level composition"

Rscript \
08_REPRODUCIBILITY/scripts/R/scrnaseq/05_celltype_composition.R \
2>&1 | tee \
08_REPRODUCIBILITY/logs/scrnaseq/05_celltype_composition.log

COMP="05_RESULTS/scrnaseq/composition/broad_celltype_proportions_by_patient.tsv"

if [ ! -s "$COMP" ]; then
    echo "ERROR: 05 did not generate composition table"
    exit 1
fi

echo "05 PASS"

# ------------------------------------------------------------
# 06
# ------------------------------------------------------------

echo
echo "RUNNING 06 — follicular/epithelial atlas"

Rscript \
08_REPRODUCIBILITY/scripts/R/scrnaseq/06_epithelial_follicular_atlas.R \
2>&1 | tee \
08_REPRODUCIBILITY/logs/scrnaseq/06_epithelial_follicular_atlas.log

if [ ! -s "$EPI" ]; then
    echo "ERROR: 06 did not generate epithelial atlas"
    exit 1
fi

echo "06 PASS"

echo
echo "============================================================"
echo "PIPELINE COMPLETE THROUGH EPITHELIAL ATLAS"
echo "============================================================"
