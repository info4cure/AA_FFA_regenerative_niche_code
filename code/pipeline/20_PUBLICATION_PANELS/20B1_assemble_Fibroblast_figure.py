#!/usr/bin/env python3

from pathlib import Path
import re
import copy
import xml.etree.ElementTree as ET

# ============================================================
# PATHS
# ============================================================

ROOT = Path(
    "PIPELINE_V2/results/figure_specification/"
    "19D4_final_publication_L2_suite/figures/Fibroblast"
)

OUT = Path(
    "PIPELINE_V2/results/figure_specification/"
    "19D5_assembled_compartment_figures/Fibroblast"
)

SVG_OUT = OUT / "SVG"
SVG_OUT.mkdir(parents=True, exist_ok=True)

PANELS = {
    "A": ROOT / "A_L2_UMAP/SVG/19D4_Fibroblast_A_L2_UMAP.svg",
    "B": ROOT / "B_cell_marker_heatmap/SVG/19D4_Fibroblast_B_cell_marker_heatmap.svg",
    "C": ROOT / "C_marker_violins/SVG/19D4_Fibroblast_C_marker_violins.svg",
    "D": ROOT / "D_patient_composition/SVG/19D4_Fibroblast_D_patient_composition.svg",
    "E": ROOT / "D2_disease_composition/SVG/19D4_Fibroblast_D2_disease_composition.svg",
}

for p in PANELS.values():
    assert p.exists(), f"Missing panel: {p}"

# ============================================================
# MASTER CANVAS
# ============================================================

# Wide journal-style canvas.
W = 1800
H = 1120

MARGIN = 35
GAP = 22

# Layout:
#
# A/B occupy upper-left/middle.
# C is tall and occupies upper-right through row 2.
# D occupies row 2 under A+B.
# E/F occupy bottom row.

BOXES = {
    # ROW 1 — molecular identity
    "A": (35,   55,  500, 500),
    "B": (555,  55,  600, 500),
    "C": (1175, 55,  590, 500),

    # ROW 2 — composition
    "D": (35,  600, 1080, 475),
    "E": (1140,600, 625, 475),
}

# ============================================================
# SVG HELPERS
# ============================================================

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"

ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)

def parse_viewbox(root):
    vb = root.attrib.get("viewBox")
    if vb:
        vals = [float(x) for x in vb.split()]
        assert len(vals) == 4
        return vals

    def val(x):
        return float(re.sub(r"[A-Za-z]+$", "", x))

    w = val(root.attrib["width"])
    h = val(root.attrib["height"])
    return [0.0, 0.0, w, h]


def import_svg(master, path, box, panel_letter):
    x, y, bw, bh = box

    tree = ET.parse(path)
    src = tree.getroot()

    vx, vy, vw, vh = parse_viewbox(src)

    # Preserve aspect ratio.
    scale = min(bw / vw, bh / vh)

    rw = vw * scale
    rh = vh * scale

    # Center inside allocated box.
    tx = x + (bw - rw) / 2
    ty = y + (bh - rh) / 2

    wrapper = ET.SubElement(
        master,
        f"{{{SVG_NS}}}g",
        {
            "id": f"panel_{panel_letter}",
            "transform": f"translate({tx:.3f},{ty:.3f})"
        }
    )

    inner = ET.SubElement(
        wrapper,
        f"{{{SVG_NS}}}g",
        {
            "transform":
                f"scale({scale:.8f}) "
                f"translate({-vx:.3f},{-vy:.3f})"
        }
    )

    # Copy children rather than nesting an external SVG.
    for child in list(src):
        inner.append(copy.deepcopy(child))

    # Panel letter.
    label = ET.SubElement(
        master,
        f"{{{SVG_NS}}}text",
        {
            "x": str(x),
            "y": str(y + 4),
            "font-family": "Arial, Helvetica, sans-serif",
            "font-size": "32",
            "font-weight": "700",
            "fill": "#000000",
            "dominant-baseline": "hanging",
        }
    )
    label.text = panel_letter


# ============================================================
# BUILD MASTER SVG
# ============================================================

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
        "fill": "#FFFFFF"
    }
)

for letter in ["A", "B", "C", "D", "E"]:
    import_svg(
        master,
        PANELS[letter],
        BOXES[letter],
        letter
    )

outfile = SVG_OUT / "19D5_Fibroblast_assembled.svg"

ET.ElementTree(master).write(
    outfile,
    encoding="utf-8",
    xml_declaration=True
)

print("============================================================")
print("19D5 — FIBROBLAST ASSEMBLY")
print("============================================================")
print("canvas:", W, "x", H)
print("panels:", len(PANELS))
print("output:", outfile)
print("STATUS=PASS")
