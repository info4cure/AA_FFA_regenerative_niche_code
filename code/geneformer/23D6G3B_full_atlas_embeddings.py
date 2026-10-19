
from pathlib import Path
from geneformer import EmbExtractor
from datasets import load_from_disk
import torch
import pandas as pd
import numpy as np

ROOT = Path.cwd()

BASE = (
    ROOT /
    "PIPELINE_V2/results/disease_biology/"
    "23_ADVANCED_ANALYSIS_ROADMAP/"
    "23D_AI_PERTURBATION/"
    "23D6_PERTURBATION/"
    "23D6G_FULL_ATLAS_REASSESSMENT"
)

DATASET = (
    BASE /
    "23D6G2C_FULL_ATLAS_TOKENIZATION/"
    "dataset/"
    "23D6G2C_geneformer_full_atlas.dataset"
)

MODEL = (
    ROOT /
    "PIPELINE_V2/external_models/Geneformer/"
    "Geneformer-V2-104M"
)

OUT = BASE / "23D6G3B_FULL_ATLAS_EMBEDDINGS"
OBJ = OUT / "objects"
TAB = OUT / "tables"

OBJ.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("23D6G3B — FULL-ATLAS GENEFORMER V2 EMBEDDINGS")
print("=" * 80)

# ------------------------------------------------------------------
# Minimal input checks
# ------------------------------------------------------------------

assert DATASET.exists(), DATASET
assert (MODEL / "config.json").exists()
assert (MODEL / "model.safetensors").exists()

ds = load_from_disk(str(DATASET))

assert len(ds) == 38699
assert len(set(ds["cell_id"])) == 38699
assert max(ds["length"]) <= 4096

print("cells:", len(ds))
print("model:", MODEL)
print("MPS available:", torch.backends.mps.is_available())

# ------------------------------------------------------------------
# Geneformer embeddings
# ------------------------------------------------------------------

labels = [
    "cell_id",
    "full_atlas_order_index",
    "patient_id",
    "disease",
    "compartment",
    "L2",
    "stratum",
]

embex = EmbExtractor(
    model_type="Pretrained",
    num_classes=0,
    emb_mode="cls",
    max_ncells=None,
    emb_layer=-1,
    emb_label=labels,
    labels_to_plot=None,
    forward_batch_size=8,
    nproc=1,
    summary_stat=None,
    save_tdigest=False,
    model_version="V2",
)

print()
print("===== EXTRACT EMBEDDINGS =====")

result = embex.extract_embs(
    model_directory=str(MODEL),
    input_data_file=str(DATASET),
    output_directory=str(OBJ),
    output_prefix="23D6G3B_full_atlas",
    output_torch_embs=True,
)

print()
print("extract_embs return type:", type(result))

# ------------------------------------------------------------------
# Inspect outputs
# ------------------------------------------------------------------

print()
print("===== GENERATED FILES =====")

for p in sorted(OBJ.rglob("*")):
    if p.is_file():
        print(p)

# ------------------------------------------------------------------
# Save execution manifest
# ------------------------------------------------------------------

manifest = pd.DataFrame({
    "cell_id": ds["cell_id"],
    "full_atlas_order_index": ds["full_atlas_order_index"],
    "patient_id": ds["patient_id"],
    "disease": ds["disease"],
    "compartment": ds["compartment"],
    "L2": ds["L2"],
    "stratum": ds["stratum"],
    "length": ds["length"],
})

manifest = manifest.sort_values(
    "full_atlas_order_index"
).reset_index(drop=True)

manifest.to_csv(
    TAB / "23D6G3B_embedding_cell_manifest.csv",
    index=False,
)

assert len(manifest) == 38699
assert manifest["cell_id"].is_unique
assert manifest["full_atlas_order_index"].tolist() == list(
    range(1, 38700)
)

# ------------------------------------------------------------------
# Audit
# ------------------------------------------------------------------

audit = f"""MODULE=23D6G3B
TITLE=Full-atlas Geneformer V2 CLS embedding extraction
SOURCE_MODULE=23D6G2C
SOURCE_CELLS={len(ds)}
MODEL_VERSION=Geneformer_V2
MODEL_TYPE=Pretrained
EMB_MODE=cls
EMB_LAYER=-1
MAX_NCELLS=None
FORWARD_BATCH_SIZE=8
NPROC=1
MODEL_INPUT_SIZE=4096
CELL_IDS_PRESERVED=TRUE
FULL_ATLAS_ORDER_PRESERVED=TRUE
PERTURBATION_PERFORMED=FALSE
CLASSIFIER_EXECUTION_PERFORMED=FALSE
FINE_TUNING_PERFORMED=FALSE
STATUS=PASS
"""

(OUT / "23D6G3B_audit.txt").write_text(audit)

print()
print(audit)

print("=" * 80)
print("23D6G3B COMPLETE")
print("=" * 80)

