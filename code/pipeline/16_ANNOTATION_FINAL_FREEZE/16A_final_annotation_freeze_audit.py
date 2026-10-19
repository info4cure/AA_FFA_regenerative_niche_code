#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib
from datetime import datetime

ROOT = Path(".")
MAN = Path("08_REPRODUCIBILITY/manuscript")

METHODS = MAN / "SCRNASEQ_Methods_frozen_modules.md"
RESULTS = MAN / "SCRNASEQ_Results_frozen_modules.md"
STATUS = MAN / "SCRNASEQ_MANUSCRIPT_STATUS.tsv"
AUDIT = MAN / "SCRNASEQ_MANUSCRIPT_AUDIT.txt"

OUT = Path(
    "PIPELINE_V2/results/annotation/"
    "16_ANNOTATION_FINAL_FREEZE"
)
OUT.mkdir(parents=True, exist_ok=True)

for p in (METHODS, RESULTS, STATUS, AUDIT):
    if not p.exists():
        raise SystemExit(f"ERROR: missing {p}")


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


methods = METHODS.read_text(encoding="utf-8")
results = RESULTS.read_text(encoding="utf-8")
status_text = STATUS.read_text(encoding="utf-8")

with STATUS.open(
    encoding="utf-8",
    newline=""
) as f:
    status_rows = list(
        csv.DictReader(f, delimiter="\t")
    )

expected_modules = [
    "06_GLOBAL_L1",
    "07_EPITHELIAL_L2",
    "08_IMMUNE_L2",
    "09_FIBROBLAST_L2",
    "10_VASCULAR_L2",
    "11_MURAL_L2",
    "12_SWEAT_GLAND_L2",
    "13_MELANOCYTE_L2",
    "14_SEBOCYTE_L2",
]

status_map = {
    r["module"]: r["status"]
    for r in status_rows
}

checks = {}

# ------------------------------------------------------------
# STATUS completeness
# ------------------------------------------------------------

checks["status_nine_modules"] = (
    len(status_rows) == 9
)

checks["status_exact_module_set"] = (
    set(status_map) == set(expected_modules)
)

for module in expected_modules:
    checks[f"{module}_FROZEN"] = (
        status_map.get(module) == "FROZEN"
    )

checks["15_not_independent_module"] = (
    not any(
        r["module"].startswith("15_")
        for r in status_rows
    )
)

# ------------------------------------------------------------
# Key manuscript section uniqueness
# ------------------------------------------------------------

method_headings = [
    "## Epithelial compartment — FROZEN",
    "## Immune-compartment refinement and L2 annotation",
    "## Fibroblast-compartment refinement and hierarchical L2 annotation",
    "## Vascular-compartment refinement and L2 annotation",
    "## Mural and smooth-muscle compartment refinement and L2 annotation",
    "## Sweat-gland L2 annotation",
    "## Melanocyte compartment refinement and L2 annotation",
    "## Sebocyte compartment refinement and L2 annotation",
]

result_headings = [
    "## Epithelial compartment — FROZEN",
    "## Immune-compartment refinement resolves nine L2 identities",
    "## Fibroblast refinement resolves consensus families and specialized states",
    "## Vascular refinement resolves major endothelial lineages and specialized states",
    "## Mural refinement resolves pericyte, contractile, APM-like and follicular-interface states",
    "## Sweat-gland refinement resolves secretory, ductal and basal/myoepithelial states",
    "## Melanocyte refinement resolves four recurrent transcriptional states",
    "## Sebocyte refinement resolves progenitor-like and differentiated lipid states",
]

for i, heading in enumerate(method_headings, start=7):
    checks[f"methods_module_{i:02d}_once"] = (
        methods.count(heading) == 1
    )

for i, heading in enumerate(result_headings, start=7):
    checks[f"results_module_{i:02d}_once"] = (
        results.count(heading) == 1
    )

# ------------------------------------------------------------
# Critical biological invariants
# ------------------------------------------------------------

checks["global_82893_documented"] = (
    "82,893" in status_text
)

checks["immune_4542_documented"] = (
    "4,542 refined immune cells" in status_text
)

checks["fibroblast_4685_documented"] = (
    "4,685 refined fibroblast cells" in status_text
)

checks["vascular_3023_documented"] = (
    "3,023 refined vascular cells" in status_text
)

checks["mural_1603_documented"] = (
    "1,603 mural/smooth-muscle cells" in status_text
)

checks["sweat_1843_documented"] = (
    "1,843 sweat-gland cells" in status_text
)

checks["melanocyte_1517_documented"] = (
    "1,517 refined melanocytes" in status_text
)

checks["sebocyte_1388_documented"] = (
    "1,388 sebocyte-compartment cells" in status_text
)

checks["vascular_lymphatic_integrated"] = (
    "182 lymphatic" in results
    and "185 cells" in results
    and "lymphatic endothelium" in results
)

checks["melanocyte_39_return_documented"] = (
    "39-cell boundary population" in results
    and "1,517" in results
    and "after return of 39 basal epithelial boundary cells" in status_text
)

checks["disease_not_used_annotation"] = (
    methods.count("Disease status was not used") >= 5
    and results.count("Disease status was not used") >= 5
)

# ------------------------------------------------------------
# Final result
# ------------------------------------------------------------

all_pass = all(checks.values())

print("=" * 68)
print("PIPELINE_V2 — FINAL ANNOTATION FREEZE AUDIT")
print("=" * 68)

print("\n=== FORMAL CHECKS ===")

for name, passed in checks.items():
    print(f"{name:<55} {passed}")

print("\n=== DOCUMENT HASHES ===")

hashes = {
    "METHODS": sha256(METHODS),
    "RESULTS": sha256(RESULTS),
    "STATUS": sha256(STATUS),
    "AUDIT": sha256(AUDIT),
}

for name, value in hashes.items():
    print(f"{name:<10} {value}")

# ------------------------------------------------------------
# Write machine-readable manifest
# ------------------------------------------------------------

manifest = OUT / "ANNOTATION_FINAL_FREEZE_manifest.tsv"

with manifest.open(
    "w",
    encoding="utf-8",
    newline=""
) as f:

    w = csv.writer(f, delimiter="\t")

    w.writerow([
        "module",
        "status",
        "manuscript_state",
    ])

    for r in status_rows:
        w.writerow([
            r["module"],
            r["status"],
            r["manuscript_state"],
        ])

checks_file = OUT / "ANNOTATION_FINAL_FREEZE_checks.tsv"

with checks_file.open(
    "w",
    encoding="utf-8",
    newline=""
) as f:

    w = csv.writer(f, delimiter="\t")
    w.writerow(["check", "pass"])

    for name, passed in checks.items():
        w.writerow([name, str(passed).upper()])

hash_file = OUT / "ANNOTATION_FINAL_FREEZE_hashes.tsv"

with hash_file.open(
    "w",
    encoding="utf-8",
    newline=""
) as f:

    w = csv.writer(f, delimiter="\t")
    w.writerow(["artifact", "sha256"])

    for name, value in hashes.items():
        w.writerow([name, value])

freeze_txt = OUT / "ANNOTATION_FINAL_FREEZE.txt"

freeze_txt.write_text(
    "\n".join([
        "PIPELINE_V2 ANNOTATION FINAL FREEZE",
        f"Timestamp: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "Global atlas: 82,893 cells",
        "Patients: 16",
        "Batches: 2",
        "",
        "Frozen modules: 06-14",
        "15_LYMPHATIC_EC: validation audit only; not an independent L2 module",
        "",
        "Annotation identities are immutable for downstream disease analyses.",
        "Disease labels must not be used to redefine frozen identities.",
        "",
        f"All formal checks PASS: {all_pass}",
        "",
        *[
            f"{k}\t{v}"
            for k, v in hashes.items()
        ],
        "",
    ]),
    encoding="utf-8",
)

print("\n=== OUTPUTS ===")
print(manifest)
print(checks_file)
print(hash_file)
print(freeze_txt)

print("\n" + "=" * 68)

if not all_pass:
    print("FINAL ANNOTATION FREEZE: FAILED")
    print("No annotation state should be considered closed.")
    raise SystemExit(1)

print("GLOBAL ATLAS: 82,893 CELLS")
print("MODULES 06-14: FROZEN")
print("15_LYMPHATIC_EC: VALIDATION ONLY")
print("ANNOTATION IDENTITIES: IMMUTABLE")
print("DISEASE LABELS USED FOR ANNOTATION: NO")
print("ANNOTATION PHASE: CLOSED")
print("READY FOR DISEASE BIOLOGY")
print("PASS")
print("=" * 68)
