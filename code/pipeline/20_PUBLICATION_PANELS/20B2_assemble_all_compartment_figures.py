#!/usr/bin/env python3

from pathlib import Path
import re
import copy
import xml.etree.ElementTree as ET
import hashlib
import csv

# ============================================================
# 20B2 — ASSEMBLE ALL COMPARTMENT FIGURES
#
# A = L2 UMAP
# B = cell-marker heatmap
# C = marker violins
# D = patient composition
# E = disease composition
#
# SVG master output only.
# No analysis, clustering, annotation or UMAP recalculation.
# ============================================================

ROOT = Path(
    "PIPELINE_V2/results/figure_specification/"
    "19D4_final_publication_L2_suite/figures"
)

OUT = Path(
    "PIPELINE_V2/results/figure_specification/"
    "19D5_assembled_compartment_figures"
)

SCRIPT_AUDIT = OUT / "manifest"

OUT.mkdir(parents=True, exist_ok=True)
SCRIPT_AUDIT.mkdir(parents=True, exist_ok=True)

GROUPS = [
    "Epithelial_lineage",
    "Fibroblast",
    "Immune",
    "Endothelial",
    "Mural",
    "Melanocyte",
]

# ============================================================
# MASTER GEOMETRY
# ============================================================

W = 1800
H = 1120

# Row 1:
# A UMAP | B heatmap | C violins
#
# Row 2:
# D patient composition | E disease composition

BOXES = {
    "A": (35,   55,  500, 500),
    "B": (555,  55,  600, 500),
    "C": (1175, 55,  590, 500),

    "D": (35,  600, 1080, 475),
    "E": (1140,600, 625, 475),
}

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"

ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)

# ============================================================
# HELPERS
# ============================================================

def resolve_one(directory, pattern):
    hits = sorted(directory.glob(pattern))

    if len(hits) != 1:
        raise RuntimeError(
            f"Expected exactly one file for:\n"
            f"  directory = {directory}\n"
            f"  pattern   = {pattern}\n"
            f"Found {len(hits)}:\n" +
            "\n".join(str(x) for x in hits)
        )

    return hits[0]


def panel_files(group):

    g = ROOT / group

    if not g.exists():
        raise RuntimeError(f"Missing group directory: {g}")

    return {
        "A": resolve_one(
            g / "A_L2_UMAP/SVG",
            "*.svg"
        ),

        "B": resolve_one(
            g / "B_cell_marker_heatmap/SVG",
            "*.svg"
        ),

        "C": resolve_one(
            g / "C_marker_violins/SVG",
            "*.svg"
        ),

        "D": resolve_one(
            g / "D_patient_composition/SVG",
            "*.svg"
        ),

        "E": resolve_one(
            g / "D2_disease_composition/SVG",
            "*.svg"
        ),
    }


def parse_viewbox(root):

    vb = root.attrib.get("viewBox")

    if vb:
        vals = [float(x) for x in vb.replace(",", " ").split()]
        if len(vals) != 4:
            raise RuntimeError(f"Invalid viewBox: {vb}")
        return vals

    def numeric(x):
        return float(
            re.sub(
                r"[A-Za-z%]+$",
                "",
                x
            )
        )

    return [
        0.0,
        0.0,
        numeric(root.attrib["width"]),
        numeric(root.attrib["height"]),
    ]


def import_svg(master, path, box, letter):

    x, y, bw, bh = box

    tree = ET.parse(path)
    src = tree.getroot()

    vx, vy, vw, vh = parse_viewbox(src)

    if vw <= 0 or vh <= 0:
        raise RuntimeError(f"Invalid SVG dimensions: {path}")

    # Preserve original aspect ratio.
    scale = min(
        bw / vw,
        bh / vh
    )

    rendered_w = vw * scale
    rendered_h = vh * scale

    tx = x + (bw - rendered_w) / 2
    ty = y + (bh - rendered_h) / 2

    wrapper = ET.SubElement(
        master,
        f"{{{SVG_NS}}}g",
        {
            "id": f"panel_{letter}",
            "transform":
                f"translate({tx:.4f},{ty:.4f})"
        }
    )

    inner = ET.SubElement(
        wrapper,
        f"{{{SVG_NS}}}g",
        {
            "transform":
                f"scale({scale:.10f}) "
                f"translate({-vx:.4f},{-vy:.4f})"
        }
    )

    # Copy vector elements directly.
    for child in list(src):
        inner.append(
            copy.deepcopy(child)
        )

    # Panel letter.
    lab = ET.SubElement(
        master,
        f"{{{SVG_NS}}}text",
        {
            "x": str(x),
            "y": str(y + 2),
            "font-family":
                "Arial, Helvetica, sans-serif",
            "font-size": "32",
            "font-weight": "700",
            "fill": "#000000",
            "dominant-baseline": "hanging",
        }
    )

    lab.text = letter


def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)

    return h.hexdigest()


# ============================================================
# PREFLIGHT
# ============================================================

print("=" * 68)
print("20B2 — ALL COMPARTMENT SVG ASSEMBLY")
print("=" * 68)

all_inputs = {}

for group in GROUPS:

    print(f"\n===== {group} =====")

    files = panel_files(group)
    all_inputs[group] = files

    for letter in ["A", "B", "C", "D", "E"]:
        print(
            letter,
            files[letter]
        )

assert len(all_inputs) == 6

# ============================================================
# ASSEMBLY
# ============================================================

inventory = []

for group in GROUPS:

    print(f"\nASSEMBLING: {group}")

    files = all_inputs[group]

    group_out = OUT / group / "SVG"
    group_out.mkdir(
        parents=True,
        exist_ok=True
    )

    master = ET.Element(
        f"{{{SVG_NS}}}svg",
        {
            "width": f"{W}pt",
            "height": f"{H}pt",
            "viewBox": f"0 0 {W} {H}",
            "version": "1.1",
        }
    )

    # White background.
    ET.SubElement(
        master,
        f"{{{SVG_NS}}}rect",
        {
            "x": "0",
            "y": "0",
            "width": str(W),
            "height": str(H),
            "fill": "#FFFFFF",
        }
    )

    for letter in ["A", "B", "C", "D", "E"]:

        import_svg(
            master,
            files[letter],
            BOXES[letter],
            letter
        )

    outfile = (
        group_out /
        f"19D5_{group}_assembled.svg"
    )

    ET.ElementTree(master).write(
        outfile,
        encoding="utf-8",
        xml_declaration=True
    )

    if not outfile.exists():
        raise RuntimeError(
            f"Output missing: {outfile}"
        )

    size = outfile.stat().st_size

    if size == 0:
        raise RuntimeError(
            f"Zero-byte output: {outfile}"
        )

    digest = sha256(outfile)

    inventory.append({
        "group": group,
        "n_panels": 5,
        "width_pt": W,
        "height_pt": H,
        "output": str(outfile),
        "bytes": size,
        "sha256": digest,
    })

    print(
        "  ->",
        outfile,
        f"({size / 1024 / 1024:.2f} MB)"
    )

# ============================================================
# MANIFEST
# ============================================================

manifest = (
    SCRIPT_AUDIT /
    "19D5_assembled_SVG_manifest.tsv"
)

with manifest.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "group",
            "n_panels",
            "width_pt",
            "height_pt",
            "output",
            "bytes",
            "sha256",
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(inventory)

# ============================================================
# CHECKS
# ============================================================

outputs = sorted(
    OUT.glob(
        "*/SVG/19D5_*_assembled.svg"
    )
)

assert len(outputs) == 6
assert all(x.stat().st_size > 0 for x in outputs)

print("\n" + "=" * 68)
print("OUTPUT INVENTORY")
print("=" * 68)

for x in outputs:
    print(x)

print()
print("Groups assembled:", len(outputs))
print("Panels/group: 5")
print("Total embedded panels:", len(outputs) * 5)
print("Master format: SVG")
print("Rasterization during assembly: NO")
print("UMAP recalculated: NO")
print("Clustering recalculated: NO")
print("Annotation modified: NO")
print("STATUS=PASS")
