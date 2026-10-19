#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib
import shutil
from datetime import datetime

# ============================================================
# 21G3B4 — CROSS-MODAL INTEGRATION MANUSCRIPT FREEZE
# ============================================================

ROOT = Path(".")

METHODS_FROZEN = Path(
    "08_REPRODUCIBILITY/manuscript/SCRNASEQ_Methods_frozen_modules.md"
)
RESULTS_FROZEN = Path(
    "08_REPRODUCIBILITY/manuscript/SCRNASEQ_Results_frozen_modules.md"
)
STATUS_FILE = Path(
    "08_REPRODUCIBILITY/manuscript/SCRNASEQ_MANUSCRIPT_STATUS.tsv"
)
AUDIT_FILE = Path(
    "08_REPRODUCIBILITY/manuscript/SCRNASEQ_MANUSCRIPT_AUDIT.txt"
)

MAIN_METHODS = Path(
    "07_MANUSCRIPT/main/sections/02_methods.md"
)
MAIN_RESULTS = Path(
    "07_MANUSCRIPT/main/sections/03_results_skeleton.md"
)

BASE = Path(
    "PIPELINE_V2/results/disease_biology/"
    "21_CELL_CELL_COMMUNICATION/"
    "21G_MULTICOMPARTMENT_COMMUNICATION"
)

CROSS = BASE / "21G4_CROSSMODAL_INTEGRATION"

B2_MATRIX = CROSS / "tables/21G3B2_crossmodal_concordance_matrix.csv"
B2_THRESHOLD = CROSS / "tables/21G3B2_communication_threshold_signals_crossmodal.csv"
B2_FDR = CROSS / "tables/21G3B2_communication_FDR_signals_crossmodal.csv"

B3_AUDIT = CROSS / "21G3B3_audit.txt"

# allow audit to be located elsewhere under 21G tree if needed
if not B3_AUDIT.exists():
    hits = list(BASE.rglob("21G3B3_audit.txt"))
    if len(hits) == 1:
        B3_AUDIT = hits[0]

REQUIRED = [
    METHODS_FROZEN,
    RESULTS_FROZEN,
    STATUS_FILE,
    AUDIT_FILE,
    MAIN_METHODS,
    MAIN_RESULTS,
    B2_MATRIX,
    B2_THRESHOLD,
    B2_FDR,
    B3_AUDIT,
]

for p in REQUIRED:
    if not p.exists():
        raise FileNotFoundError(f"Required file missing: {p}")

# ============================================================
# Helpers
# ============================================================

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def backup(path):
    dst = path.with_suffix(path.suffix + ".pre_21G3B4.bak")
    shutil.copy2(path, dst)
    return dst

def append_once(path, marker, block):
    txt = path.read_text(encoding="utf-8")
    if marker in txt:
        print(f"Already integrated: {path}")
        return False

    if not txt.endswith("\n"):
        txt += "\n"

    txt += "\n" + block.strip() + "\n"
    path.write_text(txt, encoding="utf-8")
    return True

# ============================================================
# Recover frozen counts
# ============================================================

matrix = read_csv(B2_MATRIX)
threshold = read_csv(B2_THRESHOLD)
fdr = read_csv(B2_FDR)

assert len(matrix) == 317, len(matrix)
assert len(threshold) == 16, len(threshold)
assert len(fdr) == 8, len(fdr)

receivers = {r["receiver"] for r in matrix}
contrasts = {r["contrast"] for r in matrix}

assert len(receivers) == 28
assert contrasts == {"AA_vs_HC", "FFA_vs_HC", "AA_vs_FFA"}

def evidence_count(rows, contrast, cls):
    return sum(
        r["contrast"] == contrast
        and r["communication_evidence_class"] == cls
        for r in rows
    )

aa_hc_fdr005 = evidence_count(fdr, "AA_vs_HC", "FDR_LT_0.05")
ffa_hc_fdr010 = evidence_count(fdr, "FFA_vs_HC", "FDR_LT_0.10")

assert aa_hc_fdr005 == 4
assert ffa_hc_fdr010 == 4

nominal = [
    r for r in threshold
    if r["communication_evidence_class"] == "NOMINAL_P_LT_0.05"
]
assert len(nominal) == 8

# Decomposition summary
decomp_present = [
    r for r in threshold
    if r.get("decomposition_pattern", "")
]

density_supported = [
    r for r in threshold
    if r.get("decomposition_pattern") == "DENSITY_SUPPORTED"
]

# ============================================================
# Backup manuscript state
# ============================================================

for p in [
    METHODS_FROZEN,
    RESULTS_FROZEN,
    STATUS_FILE,
    AUDIT_FILE,
    MAIN_METHODS,
    MAIN_RESULTS,
]:
    backup(p)

# ============================================================
# Frozen Methods
# ============================================================

METHODS_MARKER = "## 21G3B4 — Cross-modal communication integration"

methods_block = r"""
## 21G3B4 — Cross-modal communication integration

Cross-modal integration was performed only after freezing the
multicompartment communication analysis. The 317 prespecified primary
incoming communication records spanning 28 receiver identities and three
disease contrasts were retained as the anchor universe. Differential
abundance (DA), pseudobulk differential expression (DE), pathway-level
evidence, and communication decomposition were attached without refitting
any upstream analysis.

DA evidence was joined by exact receiver identity and disease contrast.
DE was incorporated as evidence for the presence or absence of genes
differentially expressed at FDR < 0.05 within the corresponding frozen
L2 identity; because the DE summary did not provide a single state-level
direction, no aggregate DE direction was inferred.

Communication-to-pathway comparison was restricted to a bridge specified
before concordance analysis. Communication axes were mapped to frozen
20F2B pathway axes only when a direct or biologically related mapping had
been prespecified. DIRECT and RELATED mappings were retained separately.
Communication axes without a sufficiently specific pathway counterpart
were explicitly classified as having no prespecified bridge, and absence
of a matching 20F2B state-axis record was not interpreted as negative
pathway evidence.

Communication evidence retained the multiplicity-adjusted classes frozen
in 21G2. Signals were summarized as FDR < 0.05, additional FDR < 0.10,
nominal P < 0.05, or below-threshold evidence. Cross-modal triangulation
was descriptive: DA and DE were treated as orthogonal evidence domains,
without construction of a cross-modal numerical score and without
combination or recalculation of P values or q values.

For communication signals with a frozen decomposition result, changes in
incoming communication were additionally described according to the
previously defined sender-active-intensity and pooled-density components.
These decomposition analyses were interpreted descriptively and were not
used to infer causal sender-receiver mechanisms.

No DE, DA, pathway, communication, or annotation model was refitted during
cross-modal integration.
"""

append_once(
    METHODS_FROZEN,
    METHODS_MARKER,
    methods_block
)

# ============================================================
# Frozen Results
# ============================================================

RESULTS_MARKER = "## 21G3B4 — Cross-modal integration of communication signals"

results_block = r"""
## 21G3B4 — Cross-modal integration of communication signals

Cross-modal integration retained all 317 frozen communication records
across 28 receiver identities and three disease contrasts. Sixteen
communication records reached a prespecified reporting threshold. Eight
were FDR-supported: four AA-versus-HC signals reached FDR < 0.05 and four
additional FFA-versus-HC signals reached FDR < 0.10. The remaining eight
signals were nominal (P < 0.05) and were retained as contextual evidence.

In AA versus HC, FDR-supported communication differences were concentrated
in follicular epithelial states. Higher IGF2R-associated
clearance/trafficking communication was observed in Bulge/Isthmus,
ORS-basal, and ORS receiver states, whereas the Bulge/Isthmus
growth-survival communication signal was lower in AA. These signals were
accompanied by differential transcriptional structure in the corresponding
receiver states, and several also had orthogonal differential-abundance
evidence. The IGF2R-associated results were retained as contextual
clearance/trafficking signals and were not interpreted as equivalent to
increased canonical IGF trophic signaling.

In FFA versus HC, four additional FDR < 0.10 communication signals involved
reduced growth-survival communication in basal or infundibular epithelial
states. Additional nominal growth-survival reductions extended across
follicular epithelial receivers. These FFA communication differences were
therefore treated as suggestive rather than FDR < 0.05 findings.

The only AA-versus-FFA threshold signal was nominal and involved higher
IGF-trophic communication toward ORS-basal in AA; it did not survive
multiplicity correction.

Communication decomposition was available for all 16 threshold signals.
The strongest FDR-supported signals were predominantly accompanied by
density-supported decomposition, whereas sender-active-intensity evidence
did not show comparable multiplicity-adjusted support. This pattern is
consistent with a substantial contribution from the abundance or
representation of communicating cellular components, but does not establish
a causal cellular mechanism.

DA and DE were used as orthogonal triangulation rather than as proof of
communication direction. No cross-modal score, combined significance test,
or causal signaling inference was generated.
"""

append_once(
    RESULTS_FROZEN,
    RESULTS_MARKER,
    results_block
)

# ============================================================
# Main manuscript — condensed Methods
# ============================================================

MAIN_METHODS_MARKER = "<!-- 21G3B4_CROSSMODAL_METHODS -->"

main_methods_block = r"""
<!-- 21G3B4_CROSSMODAL_METHODS -->

### Cross-modal integration of cell-cell communication

Frozen receiver-level communication results were integrated descriptively
with differential abundance, pseudobulk differential expression and
pathway-level analyses. Communication served as the anchor modality.
Orthogonal evidence was linked by exact frozen L2 identity and disease
contrast. Communication-to-pathway comparisons were restricted to
prespecified direct or related biological-axis mappings; absence of a
prespecified bridge or matching pathway record was not interpreted as
negative evidence. Differential expression was summarized by the presence
of FDR-significant genes without inferring an aggregate state-level
direction. No cross-modal score or combined significance statistic was
calculated, and no upstream model was refitted.
"""

append_once(
    MAIN_METHODS,
    MAIN_METHODS_MARKER,
    main_methods_block
)

# ============================================================
# Main manuscript — condensed Results
# ============================================================

MAIN_RESULTS_MARKER = "<!-- 21G3B4_CROSSMODAL_RESULTS -->"

main_results_block = r"""
<!-- 21G3B4_CROSSMODAL_RESULTS -->

### Cross-modal integration identifies epithelial communication remodeling

Cross-modal integration identified 16 communication signals meeting the
prespecified reporting threshold, including eight FDR-supported signals.
AA versus HC contained four signals at FDR < 0.05, concentrated in
follicular epithelial receivers. These comprised increased
IGF2R-associated clearance/trafficking communication toward
Bulge/Isthmus, ORS-basal and ORS states and reduced growth-survival
communication toward Bulge/Isthmus. Several of these receiver states also
showed differential abundance and/or transcriptional remodeling.

FFA versus HC showed four additional communication signals at FDR < 0.10,
characterized by reduced growth-survival communication across basal and
infundibular epithelial receivers, with additional nominal reductions in
follicular epithelial states. The AA-versus-FFA comparison yielded only
one nominal communication signal and no FDR-supported difference.

Decomposition of the threshold signals indicated that the strongest
communication differences were predominantly supported by the density
component rather than sender-active intensity. Thus, communication
remodeling was interpreted together with cellular abundance and
transcriptional changes, without inferring causal intercellular signaling.
"""

append_once(
    MAIN_RESULTS,
    MAIN_RESULTS_MARKER,
    main_results_block
)

# ============================================================
# Status
# ============================================================

status_txt = STATUS_FILE.read_text(encoding="utf-8")

if "21G3_CROSSMODAL_INTEGRATION" not in status_txt:
    if not status_txt.endswith("\n"):
        status_txt += "\n"

    status_txt += (
        "21G3_CROSSMODAL_INTEGRATION\tFROZEN\t"
        "317 communication-anchor rows; 28 receivers; 3 contrasts; "
        "16 threshold signals; 8 FDR signals; "
        "4 AA-HC q<0.05; 4 additional FFA-HC q<0.10; "
        "DA/DE/pathway/decomposition integrated descriptively; "
        "no cross-modal score; no combined p/q values\n"
    )

    STATUS_FILE.write_text(status_txt, encoding="utf-8")

# ============================================================
# Audit
# ============================================================

timestamp = datetime.now().isoformat(timespec="seconds")

audit_block = f"""

========================================================================
21G3B4 CROSS-MODAL INTEGRATION — MANUSCRIPT FREEZE
========================================================================
Date: {timestamp}

STATUS
------
21G3_CROSSMODAL_INTEGRATION = FROZEN

FROZEN COUNTS
-------------
Communication anchor rows = {len(matrix)}
Receivers = {len(receivers)}
Disease contrasts = {len(contrasts)}
Threshold communication signals = {len(threshold)}
FDR communication signals = {len(fdr)}
AA-vs-HC FDR<0.05 signals = {aa_hc_fdr005}
FFA-vs-HC additional FDR<0.10 signals = {ffa_hc_fdr010}
Nominal-only communication signals = {len(nominal)}
Threshold signals with decomposition = {len(decomp_present)}
Density-supported decompositions = {len(density_supported)}

INTERPRETIVE BOUNDARIES
-----------------------
- communication retained as anchor modality
- DA treated as orthogonal abundance evidence
- DE treated as transcriptional-remodeling evidence
- no aggregate DE direction inferred
- DIRECT and RELATED pathway bridges distinguished
- no-prespecified-bridge is not negative pathway evidence
- missing 20F2B state-axis record is not negative pathway evidence
- IGF2R_CONTEXTUAL denotes clearance/trafficking context and is not
  interpreted as canonical IGF trophic signaling
- decomposition does not establish causal sender-receiver mechanisms
- no cross-modal numerical score
- no combined p-value
- no combined q-value
- no causal signaling inference

SOURCE SHA-256
--------------
{sha256(B2_MATRIX)}  {B2_MATRIX}
{sha256(B2_THRESHOLD)}  {B2_THRESHOLD}
{sha256(B2_FDR)}  {B2_FDR}
{sha256(B3_AUDIT)}  {B3_AUDIT}

MANUSCRIPT FILES UPDATED
------------------------
{METHODS_FROZEN}
{RESULTS_FROZEN}
{STATUS_FILE}
{MAIN_METHODS}
{MAIN_RESULTS}

ANALYTICAL INTEGRITY
--------------------
- no communication refit
- no DE refit
- no DA refit
- no pathway refit
- no p-value recalculation
- no q-value recalculation
- no BH recalculation
- no annotation modification
- manuscript integration only

POST-INTEGRATION SHA-256
-----------------------
{sha256(METHODS_FROZEN)}  {METHODS_FROZEN}
{sha256(RESULTS_FROZEN)}  {RESULTS_FROZEN}
{sha256(STATUS_FILE)}  {STATUS_FILE}
{sha256(MAIN_METHODS)}  {MAIN_METHODS}
{sha256(MAIN_RESULTS)}  {MAIN_RESULTS}

21G3B4 COMPLETE
PASS
========================================================================
"""

with open(AUDIT_FILE, "a", encoding="utf-8") as f:
    f.write(audit_block)

# ============================================================
# Final report
# ============================================================

print("=" * 72)
print("SCRNASEQ MANUSCRIPT UPDATE — 21G3B4 COMPLETE")
print("=" * 72)
print("21G3_CROSSMODAL_INTEGRATION -> FROZEN")
print("Frozen Methods integrated")
print("Frozen Results integrated")
print("Main Methods integrated")
print("Main Results integrated")
print(f"{len(matrix)} communication-anchor rows")
print(f"{len(receivers)} receivers")
print(f"{len(threshold)} threshold signals")
print(f"{len(fdr)} FDR signals")
print(f"{aa_hc_fdr005} AA-vs-HC q<0.05")
print(f"{ffa_hc_fdr010} additional FFA-vs-HC q<0.10")
print("No statistical recalculation")
print("No analytical refit")
print("PASS")
