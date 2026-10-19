#!/usr/bin/env python3

import os
import sys
import time
import pickle
import shutil
import warnings
from pathlib import Path

# ----------------------------------------------------------------------
# SERIAL GUARDS — important on macOS / MPS
# ----------------------------------------------------------------------

os.environ["HF_DATASETS_NUM_PROC"] = "1"
os.environ["TOKENIZERS_PARALLELISM"] = "false"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
import torch

from datasets import load_from_disk
from geneformer import InSilicoPerturber


# ======================================================================
# CONFIGURATION
# ======================================================================

ROOT = Path.cwd()

BASE = (
    ROOT
    / "PIPELINE_V2/results/disease_biology"
    / "23_ADVANCED_ANALYSIS_ROADMAP"
    / "23D_AI_PERTURBATION"
    / "23D6_PERTURBATION"
)

PANEL = (
    BASE.parent
    / "23D6_TARGET_SELECTION"
    / "23D6C_PILOT_PANEL_FREEZE"
    / "tables"
    / "23D6C_frozen_pilot_perturbation_panel.csv"
)

BALANCED = (
    BASE
    / "23D6G_FULL_ATLAS_REASSESSMENT"
    / "23D6G4C_BALANCED_EXECUTION_SETS"
)

OUT = (
    BASE
    / "23D6G_FULL_ATLAS_REASSESSMENT"
    / "23D6G4D_FULL_ATLAS_PERTURBATIONS"
)

MODEL = (
    ROOT
    / "PIPELINE_V2/external_models/Geneformer"
    / "Geneformer-V2-104M"
)

OBJ = OUT / "objects"
RUNS = OBJ / "individual_runs"
TAB = OUT / "tables"

RUNS.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)



def main():
    print("=" * 80)
    print("23D6G4D — FULL-ATLAS BALANCED GENEFORMER PERTURBATIONS")
    print("24-TARGET FROZEN PANEL")
    print("PATIENT-STRATIFIED BALANCED EXECUTION SETS")
    print("RESTARTABLE EXECUTION")
    print("=" * 80)

    print()
    print("Python:", sys.version.split()[0])
    print("Torch:", torch.__version__)
    print("MPS built:", torch.backends.mps.is_built())
    print("MPS available:", torch.backends.mps.is_available())
    print("Model:", MODEL)
    print("Panel:", PANEL)
    print("Balanced source:", BALANCED)
    print("Output:", OUT)

    assert MODEL.exists(), MODEL
    assert PANEL.exists(), PANEL
    assert BALANCED.exists(), BALANCED

    if not torch.backends.mps.is_available():
        print(
            "WARNING: MPS unavailable. Geneformer may execute on CPU."
        )


    # ======================================================================
    # 1. LOAD FROZEN PANEL
    # ======================================================================

    print("\n===== 1. LOAD FROZEN PANEL =====")

    panel = pd.read_csv(PANEL)

    assert len(panel) == 24
    assert panel["pilot_perturbation_id"].nunique() == 24

    required = [
        "pilot_perturbation_id",
        "disease",
        "final_L2_label",
        "gene",
        "ensembl_id",
        "geneformer_token_id",
        "geneformer_perturb_type",
        "rescue_direction",
        "pre_geneformer_evidence_class",
    ]

    missing = [x for x in required if x not in panel.columns]
    assert not missing, missing

    print("targets:", len(panel))
    print(
        panel[
            [
                "pilot_perturbation_id",
                "disease",
                "final_L2_label",
                "gene",
                "geneformer_perturb_type",
            ]
        ].to_string(index=False)
    )

    print("FROZEN_PANEL=PASS")


    # ======================================================================
    # 2. LOCATE BALANCED DATASETS
    # ======================================================================

    print("\n===== 2. LOCATE BALANCED EXECUTION DATASETS =====")

    dataset_paths = {}

    for pid in panel["pilot_perturbation_id"]:

        candidates = list(
            BALANCED.rglob(f"{pid}.dataset")
        )

        if len(candidates) == 0:
            candidates = list(
                BALANCED.rglob(f"*{pid}*.dataset")
            )

        if len(candidates) != 1:
            raise RuntimeError(
                f"{pid}: expected exactly one balanced dataset; "
                f"found {len(candidates)}: {candidates}"
            )

        dataset_paths[pid] = candidates[0]

        ds_tmp = load_from_disk(str(candidates[0]))

        print(
            pid,
            "rows=",
            len(ds_tmp),
            "path=",
            candidates[0],
        )

    assert len(dataset_paths) == 24

    print("BALANCED_DATASETS=PASS")


    # ======================================================================
    # 3. HELPERS
    # ======================================================================

    def find_result_pickles(path):
        """Return candidate Geneformer perturbation result files."""
        path = Path(path)

        out = []

        for pattern in [
            "*.pickle",
            "*.pkl",
        ]:
            out.extend(path.rglob(pattern))

        return sorted(set(out))


    def flatten_numeric(obj):
        """
        Recursively recover scalar cosine-like numerical values from
        Geneformer result objects.

        This intentionally does not assume one exact pickle nesting layout.
        """
        vals = []

        if obj is None:
            return vals

        if torch.is_tensor(obj):
            arr = obj.detach().cpu().numpy()
            vals.extend(
                np.asarray(arr, dtype=float).ravel().tolist()
            )
            return vals

        if isinstance(obj, np.ndarray):
            try:
                vals.extend(
                    np.asarray(obj, dtype=float).ravel().tolist()
                )
            except Exception:
                pass
            return vals

        if isinstance(
            obj,
            (
                float,
                int,
                np.floating,
                np.integer,
            ),
        ):
            vals.append(float(obj))
            return vals

        if isinstance(obj, dict):
            for _, v in obj.items():
                vals.extend(flatten_numeric(v))
            return vals

        if isinstance(obj, (list, tuple)):
            for v in obj:
                vals.extend(flatten_numeric(v))
            return vals

        return vals


    def recover_cosines_from_pickle(p):
        """
        Recover finite cosine values from a Geneformer output pickle.

        Geneformer ISP result structures may vary slightly by perturbation
        mode/version; preserve the raw object and extract numerical values
        conservatively.
        """

        with open(p, "rb") as f:
            obj = pickle.load(f)

        rows = []

        if isinstance(obj, dict):

            for key, value in obj.items():

                vals = flatten_numeric(value)

                vals = [
                    float(v)
                    for v in vals
                    if np.isfinite(v)
                ]

                for j, v in enumerate(vals):
                    rows.append(
                        {
                            "source_file": str(p),
                            "key": repr(key),
                            "value_index": j,
                            "cosine_similarity": v,
                        }
                    )

        else:

            vals = flatten_numeric(obj)

            vals = [
                float(v)
                for v in vals
                if np.isfinite(v)
            ]

            for j, v in enumerate(vals):
                rows.append(
                    {
                        "source_file": str(p),
                        "key": "ROOT",
                        "value_index": j,
                        "cosine_similarity": v,
                    }
                )

        return rows


    # ======================================================================
    # 4. EXECUTION
    # ======================================================================

    print("\n===== 3. RUN FULL-ATLAS BALANCED PERTURBATIONS =====")

    execution_summary = []
    raw_cosine_rows = []

    n_total = len(panel)

    for counter, (_, row) in enumerate(panel.iterrows(), start=1):

        pid = str(row["pilot_perturbation_id"])
        disease = str(row["disease"])
        L2 = str(row["final_L2_label"])
        gene = str(row["gene"])
        ensembl = str(row["ensembl_id"])
        perturb_type = str(row["geneformer_perturb_type"])

        ds_path = dataset_paths[pid]

        ds = load_from_disk(str(ds_path))
        n_cells = len(ds)

        run_dir = RUNS / pid
        run_dir.mkdir(parents=True, exist_ok=True)

        complete_marker = run_dir / "COMPLETE.txt"

        print()
        print("-" * 80)
        print(
            f"[{counter:02d}/{n_total}] "
            f"{pid} | {disease} | {L2} | "
            f"{gene} | {perturb_type}"
        )
        print("cells:", n_cells)
        print("dataset:", ds_path)
        print("-" * 80)

        # ------------------------------------------------------------------
        # Restartability
        # ------------------------------------------------------------------

        if complete_marker.exists():

            print("STATUS=ALREADY_COMPLETE — skipping model execution")

            existing_csv = (
                TAB / f"{pid}_raw_cosine_values.csv"
            )

            if existing_csv.exists():
                old = pd.read_csv(existing_csv)
                raw_cosine_rows.extend(
                    old.to_dict("records")
                )

            execution_summary.append(
                {
                    "pilot_perturbation_id": pid,
                    "disease": disease,
                    "final_L2_label": L2,
                    "gene": gene,
                    "ensembl_id": ensembl,
                    "geneformer_perturb_type": perturb_type,
                    "n_execution_cells": n_cells,
                    "execution_status": "PASS_REUSED",
                    "elapsed_seconds": np.nan,
                }
            )

            continue

        # clean incomplete prior run
        for old_file in run_dir.iterdir():
            if old_file.name != "COMPLETE.txt":
                if old_file.is_dir():
                    shutil.rmtree(old_file)
                else:
                    old_file.unlink()

        # ------------------------------------------------------------------
        # Validate state specificity
        # ------------------------------------------------------------------

        assert all(
            str(x) == disease
            for x in ds["disease"]
        ), pid

        assert all(
            str(x) == L2
            for x in ds["L2"]
        ), pid

        # ------------------------------------------------------------------
        # Construct ISP
        # ------------------------------------------------------------------

        kwargs = dict(
            perturb_type=perturb_type,
            genes_to_perturb=[ensembl],
            combos=0,
            anchor_gene=None,
            model_type="Pretrained",
            num_classes=0,
            emb_mode="cls",
            cell_emb_style="mean_pool",
            filter_data=None,
            cell_states_to_model=None,
            state_embs_dict=None,
            max_ncells=None,
            cell_inds_to_perturb="all",
            emb_layer=-1,
            forward_batch_size=8,
            nproc=1,
            model_version="V2",
            clear_mem_ncells=100,
        )

        if perturb_type == "overexpress":
            kwargs["perturb_rank_shift"] = None

        isp = InSilicoPerturber(**kwargs)

        print("genes_to_perturb:", [ensembl])
        print("forward_batch_size:", 8)
        print("nproc:", 1)
        print("emb_mode: cls")
        print("emb_layer: -1")

        # ------------------------------------------------------------------
        # Execute
        # ------------------------------------------------------------------

        t0 = time.time()

        output_prefix = f"{pid}_{gene}_{perturb_type}"

        try:

            isp.perturb_data(
                model_directory=str(MODEL),
                input_data_file=str(ds_path),
                output_directory=str(run_dir),
                output_prefix=output_prefix,
            )

            elapsed = time.time() - t0

            print(
                f"MODEL_EXECUTION_FINISHED elapsed={elapsed:.2f}s"
            )

            # --------------------------------------------------------------
            # Recover result files
            # --------------------------------------------------------------

            pickles = find_result_pickles(run_dir)

            print("result pickle files:", len(pickles))

            if len(pickles) == 0:
                raise RuntimeError(
                    f"{pid}: Geneformer completed but no result "
                    "pickle was found."
                )

            this_raw = []

            for p in pickles:

                print(" result:", p)

                extracted = recover_cosines_from_pickle(p)

                for z in extracted:
                    z.update(
                        {
                            "pilot_perturbation_id": pid,
                            "disease": disease,
                            "final_L2_label": L2,
                            "gene": gene,
                            "ensembl_id": ensembl,
                            "geneformer_perturb_type": perturb_type,
                        }
                    )

                this_raw.extend(extracted)

            if len(this_raw) == 0:
                raise RuntimeError(
                    f"{pid}: no finite numerical perturbation "
                    "values recovered from Geneformer result."
                )

            this_raw_df = pd.DataFrame(this_raw)

            # Geneformer perturbation outputs are cosine similarities.
            # Allow only minimal floating-point tolerance.
            if (
                this_raw_df["cosine_similarity"].min() < -1.0001
                or
                this_raw_df["cosine_similarity"].max() > 1.0001
            ):
                raise RuntimeError(
                    f"{pid}: recovered numerical values outside "
                    "expected cosine-similarity range."
                )

            this_raw_df["cosine_shift_from_1"] = (
                1.0 - this_raw_df["cosine_similarity"]
            )

            raw_file = TAB / f"{pid}_raw_cosine_values.csv"

            this_raw_df.to_csv(
                raw_file,
                index=False,
            )

            raw_cosine_rows.extend(
                this_raw_df.to_dict("records")
            )

            # --------------------------------------------------------------
            # Per-target audit
            # --------------------------------------------------------------

            n_values = len(this_raw_df)

            target_audit = f"""MODULE=23D6G4D
    PERTURBATION_ID={pid}
    DISEASE={disease}
    L2={L2}
    GENE={gene}
    ENSEMBL_ID={ensembl}
    PERTURB_TYPE={perturb_type}
    EXECUTION_CELLS={n_cells}
    RECOVERED_COSINE_VALUES={n_values}
    MEAN_COSINE_SIMILARITY={this_raw_df['cosine_similarity'].mean()}
    MEDIAN_COSINE_SIMILARITY={this_raw_df['cosine_similarity'].median()}
    MEAN_COSINE_SHIFT_FROM_1={this_raw_df['cosine_shift_from_1'].mean()}
    MEDIAN_COSINE_SHIFT_FROM_1={this_raw_df['cosine_shift_from_1'].median()}
    MIN_COSINE_SIMILARITY={this_raw_df['cosine_similarity'].min()}
    MAX_COSINE_SIMILARITY={this_raw_df['cosine_similarity'].max()}
    ELAPSED_SECONDS={elapsed}
    MODEL_VERSION=Geneformer_V2
    MODEL_TYPE=Pretrained
    EMB_MODE=cls
    EMB_LAYER=-1
    FORWARD_BATCH_SIZE=8
    NPROC=1
    STATUS=PASS
    """

            (run_dir / "audit.txt").write_text(target_audit)

            complete_marker.write_text(
                f"{pid}\tPASS\n"
            )

            execution_summary.append(
                {
                    "pilot_perturbation_id": pid,
                    "disease": disease,
                    "final_L2_label": L2,
                    "gene": gene,
                    "ensembl_id": ensembl,
                    "geneformer_perturb_type": perturb_type,
                    "n_execution_cells": n_cells,
                    "n_recovered_cosine_values": n_values,
                    "mean_cosine_similarity": float(
                        this_raw_df["cosine_similarity"].mean()
                    ),
                    "sd_cosine_similarity": float(
                        this_raw_df["cosine_similarity"].std()
                    ),
                    "median_cosine_similarity": float(
                        this_raw_df["cosine_similarity"].median()
                    ),
                    "min_cosine_similarity": float(
                        this_raw_df["cosine_similarity"].min()
                    ),
                    "max_cosine_similarity": float(
                        this_raw_df["cosine_similarity"].max()
                    ),
                    "mean_cosine_shift_from_1": float(
                        this_raw_df["cosine_shift_from_1"].mean()
                    ),
                    "median_cosine_shift_from_1": float(
                        this_raw_df["cosine_shift_from_1"].median()
                    ),
                    "elapsed_seconds": elapsed,
                    "execution_status": "PASS",
                }
            )

            # Persist after every successful perturbation.
            pd.DataFrame(
                execution_summary
            ).to_csv(
                TAB / "23D6G4D_execution_summary.csv",
                index=False,
            )

            print("STATUS=PASS")

        except Exception as e:

            elapsed = time.time() - t0

            print()
            print("STATUS=FAIL")
            print("ERROR:", repr(e))

            execution_summary.append(
                {
                    "pilot_perturbation_id": pid,
                    "disease": disease,
                    "final_L2_label": L2,
                    "gene": gene,
                    "ensembl_id": ensembl,
                    "geneformer_perturb_type": perturb_type,
                    "n_execution_cells": n_cells,
                    "elapsed_seconds": elapsed,
                    "execution_status": "FAIL",
                    "error": repr(e),
                }
            )

            pd.DataFrame(
                execution_summary
            ).to_csv(
                TAB / "23D6G4D_execution_summary.csv",
                index=False,
            )

            raise


    # ======================================================================
    # 5. CONSOLIDATE RESULTS
    # ======================================================================

    print("\n===== 4. CONSOLIDATE RESULTS =====")

    summary = pd.DataFrame(execution_summary)

    summary.to_csv(
        TAB / "23D6G4D_execution_summary.csv",
        index=False,
    )

    print("execution summary rows:", len(summary))

    print(
        summary[
            [
                "pilot_perturbation_id",
                "disease",
                "final_L2_label",
                "gene",
                "geneformer_perturb_type",
                "n_execution_cells",
                "execution_status",
            ]
        ].to_string(index=False)
    )

    assert len(summary) == 24
    assert summary["pilot_perturbation_id"].nunique() == 24
    assert summary["execution_status"].isin(
        ["PASS", "PASS_REUSED"]
    ).all()

    print("EXECUTION_SUMMARY=PASS")


    # ======================================================================
    # 6. CONSOLIDATE RAW COSINE VALUES
    # ======================================================================

    print("\n===== 5. CONSOLIDATE COSINE RESULTS =====")

    if len(raw_cosine_rows) == 0:
        raise RuntimeError(
            "No perturbation cosine results available."
        )

    raw = pd.DataFrame(raw_cosine_rows)

    required_raw = [
        "pilot_perturbation_id",
        "disease",
        "final_L2_label",
        "gene",
        "ensembl_id",
        "geneformer_perturb_type",
        "cosine_similarity",
    ]

    missing_raw = [
        c for c in required_raw
        if c not in raw.columns
    ]

    assert not missing_raw, missing_raw

    raw["cosine_similarity"] = pd.to_numeric(
        raw["cosine_similarity"],
        errors="raise",
    )

    raw["cosine_shift_from_1"] = (
        1.0 - raw["cosine_similarity"]
    )

    assert np.isfinite(
        raw["cosine_similarity"].to_numpy()
    ).all()

    assert (
        raw["cosine_similarity"].to_numpy() >= -1.0001
    ).all()

    assert (
        raw["cosine_similarity"].to_numpy() <= 1.0001
    ).all()

    raw.to_csv(
        TAB / "23D6G4D_all_raw_cosine_values.csv",
        index=False,
    )

    print("raw cosine rows:", len(raw))
    print(
        "perturbations represented:",
        raw["pilot_perturbation_id"].nunique(),
    )

    assert raw["pilot_perturbation_id"].nunique() == 24

    print("RAW_COSINE_CONSOLIDATION=PASS")


    # ======================================================================
    # 7. PERTURBATION-LEVEL SUMMARY
    # ======================================================================

    print("\n===== 6. PERTURBATION-LEVEL SUMMARY =====")

    group_cols = [
        "pilot_perturbation_id",
        "disease",
        "final_L2_label",
        "gene",
        "ensembl_id",
        "geneformer_perturb_type",
    ]

    perturb_summary = (
        raw
        .groupby(
            group_cols,
            dropna=False,
            observed=True,
        )
        .agg(
            n_cosine_values=(
                "cosine_similarity",
                "size",
            ),
            mean_cosine_similarity=(
                "cosine_similarity",
                "mean",
            ),
            sd_cosine_similarity=(
                "cosine_similarity",
                "std",
            ),
            median_cosine_similarity=(
                "cosine_similarity",
                "median",
            ),
            min_cosine_similarity=(
                "cosine_similarity",
                "min",
            ),
            max_cosine_similarity=(
                "cosine_similarity",
                "max",
            ),
            mean_cosine_shift_from_1=(
                "cosine_shift_from_1",
                "mean",
            ),
            median_cosine_shift_from_1=(
                "cosine_shift_from_1",
                "median",
            ),
            max_cosine_shift_from_1=(
                "cosine_shift_from_1",
                "max",
            ),
        )
        .reset_index()
    )

    perturb_summary = perturb_summary.merge(
        summary[
            [
                "pilot_perturbation_id",
                "n_execution_cells",
            ]
        ],
        on="pilot_perturbation_id",
        how="left",
        validate="one_to_one",
    )

    perturb_summary = perturb_summary.sort_values(
        [
            "disease",
            "geneformer_perturb_type",
            "mean_cosine_shift_from_1",
        ],
        ascending=[
            True,
            True,
            False,
        ],
    )

    perturb_summary.to_csv(
        TAB / "23D6G4D_perturbation_level_summary.csv",
        index=False,
    )

    print(
        perturb_summary[
            [
                "pilot_perturbation_id",
                "disease",
                "final_L2_label",
                "gene",
                "geneformer_perturb_type",
                "n_execution_cells",
                "n_cosine_values",
                "mean_cosine_similarity",
                "mean_cosine_shift_from_1",
            ]
        ].to_string(index=False)
    )

    assert len(perturb_summary) == 24

    print("PERTURBATION_LEVEL_SUMMARY=PASS")


    # ======================================================================
    # 8. RANK PERTURBATIONS BY EMBEDDING DISPLACEMENT
    # ======================================================================

    print("\n===== 7. PERTURBATION RANKING =====")

    ranking = perturb_summary.copy()

    ranking["rank_global"] = (
        ranking["mean_cosine_shift_from_1"]
        .rank(
            method="min",
            ascending=False,
        )
        .astype(int)
    )

    ranking["rank_within_disease"] = (
        ranking
        .groupby("disease")[
            "mean_cosine_shift_from_1"
        ]
        .rank(
            method="min",
            ascending=False,
        )
        .astype(int)
    )

    ranking["rank_within_disease_type"] = (
        ranking
        .groupby(
            [
                "disease",
                "geneformer_perturb_type",
            ]
        )[
            "mean_cosine_shift_from_1"
        ]
        .rank(
            method="min",
            ascending=False,
        )
        .astype(int)
    )

    ranking = ranking.sort_values(
        "rank_global"
    )

    ranking.to_csv(
        TAB / "23D6G4D_perturbation_ranking.csv",
        index=False,
    )

    print(
        ranking[
            [
                "rank_global",
                "pilot_perturbation_id",
                "disease",
                "final_L2_label",
                "gene",
                "geneformer_perturb_type",
                "mean_cosine_shift_from_1",
                "median_cosine_shift_from_1",
            ]
        ].to_string(index=False)
    )

    print("PERTURBATION_RANKING=PASS")


    # ======================================================================
    # 9. FORMAL CHECKS
    # ======================================================================

    print("\n===== 8. FORMAL CHECKS =====")

    checks = pd.DataFrame(
        [
            {
                "check": "frozen_panel_n24",
                "status": (
                    "PASS"
                    if len(panel) == 24
                    else "FAIL"
                ),
                "detail": len(panel),
            },
            {
                "check": "all_24_targets_executed",
                "status": (
                    "PASS"
                    if len(summary) == 24
                    else "FAIL"
                ),
                "detail": len(summary),
            },
            {
                "check": "all_execution_status_pass",
                "status": (
                    "PASS"
                    if summary["execution_status"]
                    .isin(["PASS", "PASS_REUSED"])
                    .all()
                    else "FAIL"
                ),
                "detail": str(
                    summary["execution_status"]
                    .value_counts()
                    .to_dict()
                ),
            },
            {
                "check": "all_24_targets_have_cosine_results",
                "status": (
                    "PASS"
                    if raw["pilot_perturbation_id"]
                    .nunique() == 24
                    else "FAIL"
                ),
                "detail": raw[
                    "pilot_perturbation_id"
                ].nunique(),
            },
            {
                "check": "all_cosines_finite",
                "status": (
                    "PASS"
                    if np.isfinite(
                        raw["cosine_similarity"]
                        .to_numpy()
                    ).all()
                    else "FAIL"
                ),
                "detail": len(raw),
            },
            {
                "check": "cosines_within_valid_range",
                "status": (
                    "PASS"
                    if (
                        (
                            raw["cosine_similarity"]
                            >= -1.0001
                        ).all()
                        and
                        (
                            raw["cosine_similarity"]
                            <= 1.0001
                        ).all()
                    )
                    else "FAIL"
                ),
                "detail": (
                    f"{raw['cosine_similarity'].min()} "
                    f"to "
                    f"{raw['cosine_similarity'].max()}"
                ),
            },
            {
                "check": "irf7_full_atlas_rescued",
                "status": (
                    "PASS"
                    if (
                        "23D6C_006"
                        in set(
                            raw[
                                "pilot_perturbation_id"
                            ]
                        )
                    )
                    else "FAIL"
                ),
                "detail": "23D6C_006",
            },
            {
                "check": "no_classifier",
                "status": "PASS",
                "detail": "Pretrained",
            },
            {
                "check": "no_fine_tuning",
                "status": "PASS",
                "detail": "none",
            },
        ]
    )

    checks.to_csv(
        TAB / "23D6G4D_formal_checks.csv",
        index=False,
    )

    print(checks.to_string(index=False))

    assert (
        checks["status"] == "PASS"
    ).all()

    print("FORMAL_CHECKS=PASS")


    # ======================================================================
    # 10. FINAL AUDIT
    # ======================================================================

    print("\n===== 10. FINAL AUDIT =====")

    n_cells_total = int(
        summary["n_execution_cells"].sum()
    )

    n_raw = int(len(raw))

    top = ranking.iloc[0]

    audit = f"""MODULE=23D6G4D
    TITLE=Balanced full-atlas Geneformer V2 in-silico perturbation execution
    SOURCE_PANEL=23D6C
    SOURCE_EXECUTION_SETS=23D6G4C
    MODEL_VERSION=Geneformer_V2
    MODEL_TYPE=Pretrained
    EMB_MODE=cls
    EMB_LAYER=-1
    FORWARD_BATCH_SIZE=8
    NPROC=1
    FROZEN_TARGETS=24
    EXECUTED_TARGETS=24
    CELL_PERTURBATION_INPUT_INSTANCES={n_cells_total}
    RECOVERED_COSINE_VALUES={n_raw}
    IRF7_DELETE_FULL_ATLAS_EXECUTED=TRUE
    PATIENT_STRATIFIED_BALANCED_EXECUTION=TRUE
    MAX_CELLS_PER_PATIENT_PER_TARGET=100
    CLASSIFIER_EXECUTION_PERFORMED=FALSE
    FINE_TUNING_PERFORMED=FALSE
    TOP_EMBEDDING_SHIFT_PERTURBATION={top['pilot_perturbation_id']}
    TOP_EMBEDDING_SHIFT_GENE={top['gene']}
    TOP_EMBEDDING_SHIFT_DISEASE={top['disease']}
    TOP_EMBEDDING_SHIFT_L2={top['final_L2_label']}
    TOP_MEAN_COSINE_SHIFT_FROM_1={top['mean_cosine_shift_from_1']}
    STATUS=PASS
    """

    (OUT / "23D6G4D_audit.txt").write_text(audit)

    print(audit)

    print("=" * 80)
    print("23D6G4D COMPLETE")
    print()
    print("OUTPUTS:")
    print("- tables/23D6G4D_execution_summary.csv")
    print("- tables/23D6G4D_all_raw_cosine_values.csv")
    print("- tables/23D6G4D_perturbation_level_summary.csv")
    print("- tables/23D6G4D_perturbation_ranking.csv")
    print("- tables/23D6G4D_formal_checks.csv")
    print("- objects/individual_runs/")
    print("- 23D6G4D_audit.txt")
    print()
    print("STATUS=PASS")
    print("=" * 80)


if __name__ == "__main__":
    from multiprocess import freeze_support
    freeze_support()
    main()
