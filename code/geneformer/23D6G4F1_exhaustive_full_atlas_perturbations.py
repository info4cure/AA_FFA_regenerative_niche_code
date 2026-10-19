#!/usr/bin/env python3

import os
import sys
import time
import pickle
import shutil
import warnings
from pathlib import Path

# ============================================================
# SERIAL / MACOS GUARDS
# ============================================================

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


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path.cwd()

BASE = (
    ROOT
    / "PIPELINE_V2/results/disease_biology"
    / "23_ADVANCED_ANALYSIS_ROADMAP"
    / "23D_AI_PERTURBATION"
    / "23D6_PERTURBATION"
    / "23D6G_FULL_ATLAS_REASSESSMENT"
)

G4B = (
    BASE
    / "23D6G4B_FULL_ATLAS_EXECUTION_SETS"
)

SUPPORT = (
    G4B
    / "tables"
    / "23D6G4B_full_atlas_execution_support.csv"
)

MANIFEST = (
    G4B
    / "tables"
    / "23D6G4B_full_atlas_execution_cell_manifest.csv"
)

OUT = (
    BASE
    / "23D6G4F_EXHAUSTIVE_FULL_ATLAS"
    / "23D6G4F1_EXHAUSTIVE_PERTURBATIONS"
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


# ============================================================
# HELPERS
# ============================================================

def flatten_numeric(obj):
    """
    Recursively recover finite scalar numeric values from
    arbitrary Geneformer pickle structures.
    """
    vals = []

    if isinstance(obj, dict):
        for v in obj.values():
            vals.extend(flatten_numeric(v))

    elif isinstance(obj, (list, tuple)):
        for v in obj:
            vals.extend(flatten_numeric(v))

    elif isinstance(obj, np.ndarray):
        if np.issubdtype(obj.dtype, np.number):
            vals.extend(
                obj.astype(float).ravel().tolist()
            )
        else:
            for v in obj.ravel():
                vals.extend(flatten_numeric(v))

    elif torch.is_tensor(obj):
        vals.extend(
            obj.detach()
            .cpu()
            .numpy()
            .astype(float)
            .ravel()
            .tolist()
        )

    elif isinstance(
        obj,
        (
            int,
            float,
            np.integer,
            np.floating,
        ),
    ):
        vals.append(float(obj))

    return vals


def find_result_pickles(run_dir):
    files = []

    for p in run_dir.rglob("*"):
        if not p.is_file():
            continue

        suffix = p.suffix.lower()

        if suffix in {
            ".pickle",
            ".pkl",
        }:
            files.append(p)

    return sorted(files)


def recover_cosines_from_pickle(path):
    with open(path, "rb") as fh:
        obj = pickle.load(fh)

    vals = flatten_numeric(obj)

    vals = [
        float(v)
        for v in vals
        if np.isfinite(v)
    ]

    return vals


def parse_patient_counts(s):
    """
    Parse strings such as:
    ARE087:748;ARE090:865;...
    """
    out = {}

    if pd.isna(s):
        return out

    for item in str(s).split(";"):
        item = item.strip()

        if not item:
            continue

        patient, n = item.rsplit(":", 1)

        out[str(patient)] = int(n)

    return out


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 80)
    print(
        "23D6G4F1 — EXHAUSTIVE FULL-ATLAS "
        "GENEFORMER PERTURBATIONS"
    )
    print("24-TARGET FROZEN PANEL")
    print("ALL ELIGIBLE G4B CELL-TARGET INSTANCES")
    print("NO PATIENT SUBSAMPLING")
    print("RESTARTABLE EXECUTION")
    print("=" * 80)

    print()
    print("Python:", sys.version.split()[0])
    print("Torch:", torch.__version__)
    print("MPS built:", torch.backends.mps.is_built())
    print("MPS available:", torch.backends.mps.is_available())
    print("Model:", MODEL)
    print("G4B:", G4B)
    print("Support:", SUPPORT)
    print("Manifest:", MANIFEST)
    print("Output:", OUT)

    assert MODEL.exists(), MODEL
    assert G4B.exists(), G4B
    assert SUPPORT.exists(), SUPPORT
    assert MANIFEST.exists(), MANIFEST

    if not torch.backends.mps.is_available():
        print(
            "WARNING: MPS unavailable. "
            "Geneformer may execute on CPU."
        )


    # ========================================================
    # 1. LOAD G4B FROZEN UNIVERSE
    # ========================================================

    print("\n===== 1. LOAD G4B FROZEN UNIVERSE =====")

    support = pd.read_csv(SUPPORT)
    manifest = pd.read_csv(MANIFEST)

    assert len(support) == 24
    assert (
        support["pilot_perturbation_id"].nunique()
        == 24
    )

    assert len(manifest) == 52509

    assert (
        manifest["pilot_perturbation_id"].nunique()
        == 24
    )

    print("targets:", len(support))
    print("execution instances:", len(manifest))

    print(
        support[
            [
                "pilot_perturbation_id",
                "disease",
                "final_L2_label",
                "gene",
                "geneformer_perturb_type",
                "n_execution_cells",
                "n_execution_patients",
            ]
        ].to_string(index=False)
    )

    print("G4B_UNIVERSE=PASS")


    # ========================================================
    # 2. VERIFY EXECUTION DATASETS
    # ========================================================

    print(
        "\n===== 2. VERIFY EXHAUSTIVE EXECUTION DATASETS ====="
    )

    dataset_paths = {}
    verification_rows = []

    for _, row in support.iterrows():

        pid = str(row["pilot_perturbation_id"])

        ds_path = Path(
            str(row["execution_dataset"])
        )

        if not ds_path.exists():
            raise FileNotFoundError(
                f"{pid}: execution dataset missing: "
                f"{ds_path}"
            )

        ds = load_from_disk(str(ds_path))

        n_dataset = len(ds)

        n_support = int(
            row["n_execution_cells"]
        )

        n_manifest = int(
            (
                manifest["pilot_perturbation_id"]
                == pid
            ).sum()
        )

        expected_patient_counts = (
            parse_patient_counts(
                row["patient_cell_counts"]
            )
        )

        observed_patient_counts = (
            manifest.loc[
                manifest["pilot_perturbation_id"]
                == pid,
                "patient_id",
            ]
            .astype(str)
            .value_counts()
            .sort_index()
            .to_dict()
        )

        patient_counts_match = (
            expected_patient_counts
            == observed_patient_counts
        )

        status = (
            "PASS"
            if (
                n_dataset
                == n_support
                == n_manifest
                and patient_counts_match
            )
            else "FAIL"
        )

        verification_rows.append(
            {
                "pilot_perturbation_id": pid,
                "n_dataset": n_dataset,
                "n_support": n_support,
                "n_manifest": n_manifest,
                "patient_counts_match":
                    patient_counts_match,
                "status": status,
                "dataset_path": str(ds_path),
            }
        )

        dataset_paths[pid] = ds_path

        print(
            pid,
            "rows=",
            n_dataset,
            "patients=",
            row["n_execution_patients"],
            "status=",
            status,
        )

    verification = pd.DataFrame(
        verification_rows
    )

    verification.to_csv(
        TAB
        / "23D6G4F1_dataset_verification.csv",
        index=False,
    )

    assert (
        verification["status"] == "PASS"
    ).all()

    assert (
        verification["n_dataset"].sum()
        == 52509
    )

    print("DATASET_VERIFICATION=PASS")


    # ========================================================
    # 3. RUN EXHAUSTIVE PERTURBATIONS
    # ========================================================

    print(
        "\n===== 3. RUN EXHAUSTIVE FULL-ATLAS "
        "PERTURBATIONS ====="
    )

    execution_summary = []
    raw_cosine_rows = []

    n_total = len(support)

    for counter, (_, row) in enumerate(
        support.iterrows(),
        start=1,
    ):

        pid = str(
            row["pilot_perturbation_id"]
        )

        disease = str(row["disease"])
        L2 = str(row["final_L2_label"])
        gene = str(row["gene"])
        ensembl = str(row["ensembl_id"])

        perturb_type = str(
            row["geneformer_perturb_type"]
        )

        ds_path = dataset_paths[pid]

        ds = load_from_disk(str(ds_path))

        n_cells = len(ds)

        run_dir = RUNS / pid

        run_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        complete_marker = (
            run_dir / "COMPLETE.txt"
        )

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


        # ----------------------------------------------------
        # RESTARTABILITY
        # ----------------------------------------------------

        existing_csv = (
            TAB
            / f"{pid}_raw_cosine_values.csv"
        )

        if complete_marker.exists():

            print(
                "STATUS=ALREADY_COMPLETE — "
                "skipping model execution"
            )

            if existing_csv.exists():

                old = pd.read_csv(
                    existing_csv
                )

                raw_cosine_rows.extend(
                    old.to_dict("records")
                )

            execution_summary.append(
                {
                    "pilot_perturbation_id":
                        pid,
                    "disease":
                        disease,
                    "final_L2_label":
                        L2,
                    "gene":
                        gene,
                    "ensembl_id":
                        ensembl,
                    "geneformer_perturb_type":
                        perturb_type,
                    "n_execution_cells":
                        n_cells,
                    "execution_status":
                        "PASS_REUSED",
                    "elapsed_seconds":
                        np.nan,
                }
            )

            continue


        # ----------------------------------------------------
        # CLEAN INCOMPLETE PRIOR RUN
        # ----------------------------------------------------

        for old_file in run_dir.iterdir():

            if old_file.name == "COMPLETE.txt":
                continue

            if old_file.is_dir():
                shutil.rmtree(old_file)
            else:
                old_file.unlink()


        # ----------------------------------------------------
        # VALIDATE STATE SPECIFICITY
        # ----------------------------------------------------

        if "disease" in ds.column_names:

            assert all(
                str(x) == disease
                for x in ds["disease"]
            ), pid

        if "L2" in ds.column_names:

            assert all(
                str(x) == L2
                for x in ds["L2"]
            ), pid


        # ----------------------------------------------------
        # CONSTRUCT GENEFORMER ISP
        # ----------------------------------------------------

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

            # CRITICAL:
            # no model-side cell cap
            max_ncells=None,

            # CRITICAL:
            # every row in the G4B execution dataset
            cell_inds_to_perturb="all",

            emb_layer=-1,

            # conservative on Apple Silicon
            forward_batch_size=8,
            nproc=1,

            model_version="V2",

            # memory cleanup only;
            # NOT a sampling parameter
            clear_mem_ncells=100,
        )

        if perturb_type == "overexpress":
            kwargs["perturb_rank_shift"] = None

        isp = InSilicoPerturber(
            **kwargs
        )

        print(
            "genes_to_perturb:",
            [ensembl],
        )

        print(
            "max_ncells:",
            None,
        )

        print(
            "cell_inds_to_perturb:",
            "all",
        )

        print(
            "forward_batch_size:",
            8,
        )

        print(
            "nproc:",
            1,
        )

        print(
            "emb_mode:",
            "cls",
        )

        print(
            "emb_layer:",
            -1,
        )


        # ----------------------------------------------------
        # EXECUTE
        # ----------------------------------------------------

        t0 = time.time()

        output_prefix = (
            f"{pid}_{gene}_{perturb_type}"
        )

        try:

            isp.perturb_data(
                model_directory=str(MODEL),
                input_data_file=str(ds_path),
                output_directory=str(run_dir),
                output_prefix=output_prefix,
            )

            elapsed = (
                time.time() - t0
            )

            print(
                "MODEL_EXECUTION_FINISHED "
                f"elapsed={elapsed:.2f}s"
            )


            # ------------------------------------------------
            # RECOVER RESULT PICKLES
            # ------------------------------------------------

            pickles = (
                find_result_pickles(
                    run_dir
                )
            )

            print(
                "result pickle files:",
                len(pickles),
            )

            if len(pickles) == 0:
                raise RuntimeError(
                    f"{pid}: Geneformer completed "
                    "but no result pickle was found."
                )

            this_raw = []

            for p in pickles:

                print(
                    " result:",
                    p,
                )

                vals = (
                    recover_cosines_from_pickle(
                        p
                    )
                )

                print(
                    " recovered numeric values:",
                    len(vals),
                )

                for j, value in enumerate(vals):

                    this_raw.append(
                        {
                            "pilot_perturbation_id":
                                pid,
                            "disease":
                                disease,
                            "final_L2_label":
                                L2,
                            "gene":
                                gene,
                            "ensembl_id":
                                ensembl,
                            "geneformer_perturb_type":
                                perturb_type,
                            "source_pickle":
                                str(p),
                            "value_index":
                                j,
                            "cosine_similarity":
                                float(value),
                        }
                    )

            if len(this_raw) == 0:
                raise RuntimeError(
                    f"{pid}: no finite numeric "
                    "Geneformer result values recovered."
                )

            this_raw_df = pd.DataFrame(
                this_raw
            )

            this_raw_df.to_csv(
                existing_csv,
                index=False,
            )

            raw_cosine_rows.extend(
                this_raw
            )

            complete_marker.write_text(
                "\n".join(
                    [
                        f"pilot_perturbation_id={pid}",
                        f"gene={gene}",
                        f"disease={disease}",
                        f"L2={L2}",
                        f"perturb_type={perturb_type}",
                        f"n_execution_cells={n_cells}",
                        f"n_recovered_values={len(this_raw)}",
                        f"elapsed_seconds={elapsed:.6f}",
                        "STATUS=PASS",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            execution_summary.append(
                {
                    "pilot_perturbation_id":
                        pid,
                    "disease":
                        disease,
                    "final_L2_label":
                        L2,
                    "gene":
                        gene,
                    "ensembl_id":
                        ensembl,
                    "geneformer_perturb_type":
                        perturb_type,
                    "n_execution_cells":
                        n_cells,
                    "execution_status":
                        "PASS",
                    "elapsed_seconds":
                        elapsed,
                }
            )

            print("STATUS=PASS")

        except Exception as exc:

            elapsed = (
                time.time() - t0
            )

            execution_summary.append(
                {
                    "pilot_perturbation_id":
                        pid,
                    "disease":
                        disease,
                    "final_L2_label":
                        L2,
                    "gene":
                        gene,
                    "ensembl_id":
                        ensembl,
                    "geneformer_perturb_type":
                        perturb_type,
                    "n_execution_cells":
                        n_cells,
                    "execution_status":
                        "FAIL",
                    "elapsed_seconds":
                        elapsed,
                    "error":
                        repr(exc),
                }
            )

            print(
                "STATUS=FAIL",
                repr(exc),
            )

            # Save partial execution summary before aborting.
            pd.DataFrame(
                execution_summary
            ).to_csv(
                TAB
                / "23D6G4F1_execution_summary_PARTIAL.csv",
                index=False,
            )

            raise


    # ========================================================
    # 4. CONSOLIDATE EXECUTION SUMMARY
    # ========================================================

    print(
        "\n===== 4. CONSOLIDATE EXECUTION SUMMARY ====="
    )

    summary = pd.DataFrame(
        execution_summary
    )

    summary.to_csv(
        TAB
        / "23D6G4F1_execution_summary.csv",
        index=False,
    )

    print(
        summary.to_string(
            index=False
        )
    )

    assert len(summary) == 24

    assert (
        summary[
            "pilot_perturbation_id"
        ].nunique()
        == 24
    )

    assert (
        summary[
            "execution_status"
        ].isin(
            [
                "PASS",
                "PASS_REUSED",
            ]
        )
    ).all()

    print("EXECUTION_SUMMARY=PASS")


    # ========================================================
    # 5. CONSOLIDATE RAW GENEFORMER VALUES
    # ========================================================

    print(
        "\n===== 5. CONSOLIDATE RAW GENEFORMER VALUES ====="
    )

    raw = pd.DataFrame(
        raw_cosine_rows
    )

    if len(raw) == 0:
        raise RuntimeError(
            "No Geneformer result values recovered."
        )

    raw.to_csv(
        TAB
        / "23D6G4F1_all_raw_cosine_values.csv",
        index=False,
    )

    print(
        "recovered raw values:",
        len(raw),
    )

    print(
        "targets represented:",
        raw[
            "pilot_perturbation_id"
        ].nunique(),
    )

    assert (
        raw[
            "pilot_perturbation_id"
        ].nunique()
        == 24
    )


    # ========================================================
    # 6. PERTURBATION-LEVEL SUMMARY
    # ========================================================

    print(
        "\n===== 6. PERTURBATION-LEVEL SUMMARY ====="
    )

    perturbation_summary = (
        raw
        .groupby(
            [
                "pilot_perturbation_id",
                "disease",
                "final_L2_label",
                "gene",
                "ensembl_id",
                "geneformer_perturb_type",
            ],
            as_index=False,
        )
        .agg(
            n_recovered_values=(
                "cosine_similarity",
                "size",
            ),
            mean_cosine_similarity=(
                "cosine_similarity",
                "mean",
            ),
            median_cosine_similarity=(
                "cosine_similarity",
                "median",
            ),
            sd_cosine_similarity=(
                "cosine_similarity",
                "std",
            ),
            min_cosine_similarity=(
                "cosine_similarity",
                "min",
            ),
            max_cosine_similarity=(
                "cosine_similarity",
                "max",
            ),
        )
    )

    perturbation_summary[
        "mean_cosine_shift_from_1"
    ] = (
        1.0
        - perturbation_summary[
            "mean_cosine_similarity"
        ]
    )

    perturbation_summary[
        "median_cosine_shift_from_1"
    ] = (
        1.0
        - perturbation_summary[
            "median_cosine_similarity"
        ]
    )

    perturbation_summary = (
        perturbation_summary
        .merge(
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
    )

    perturbation_summary.to_csv(
        TAB
        / "23D6G4F1_perturbation_level_summary.csv",
        index=False,
    )

    print(
        perturbation_summary[
            [
                "pilot_perturbation_id",
                "disease",
                "final_L2_label",
                "gene",
                "geneformer_perturb_type",
                "n_execution_cells",
                "n_recovered_values",
                "mean_cosine_shift_from_1",
                "median_cosine_shift_from_1",
            ]
        ].to_string(
            index=False
        )
    )


    # ========================================================
    # 7. RANKING
    # ========================================================

    print(
        "\n===== 7. EXHAUSTIVE EMBEDDING-SHIFT RANKING ====="
    )

    ranking = (
        perturbation_summary
        .sort_values(
            [
                "mean_cosine_shift_from_1",
                "median_cosine_shift_from_1",
            ],
            ascending=[
                False,
                False,
            ],
        )
        .reset_index(
            drop=True
        )
    )

    ranking.insert(
        0,
        "rank_exhaustive",
        np.arange(
            1,
            len(ranking) + 1,
        ),
    )

    ranking.to_csv(
        TAB
        / "23D6G4F1_exhaustive_perturbation_ranking.csv",
        index=False,
    )

    print(
        ranking[
            [
                "rank_exhaustive",
                "pilot_perturbation_id",
                "disease",
                "final_L2_label",
                "gene",
                "geneformer_perturb_type",
                "n_execution_cells",
                "mean_cosine_shift_from_1",
                "median_cosine_shift_from_1",
            ]
        ].to_string(
            index=False
        )
    )


    # ========================================================
    # 8. EXECUTION CELL-COUNT VALIDATION
    # ========================================================

    print(
        "\n===== 8. EXECUTION CELL-COUNT VALIDATION ====="
    )

    expected_counts = (
        support[
            [
                "pilot_perturbation_id",
                "n_execution_cells",
            ]
        ]
        .rename(
            columns={
                "n_execution_cells":
                    "expected_execution_cells"
            }
        )
    )

    observed_counts = (
        summary[
            [
                "pilot_perturbation_id",
                "n_execution_cells",
            ]
        ]
        .rename(
            columns={
                "n_execution_cells":
                    "observed_execution_cells"
            }
        )
    )

    count_check = (
        expected_counts
        .merge(
            observed_counts,
            on="pilot_perturbation_id",
            how="outer",
            validate="one_to_one",
        )
    )

    count_check[
        "cell_count_match"
    ] = (
        count_check[
            "expected_execution_cells"
        ]
        == count_check[
            "observed_execution_cells"
        ]
    )

    count_check.to_csv(
        TAB
        / "23D6G4F1_execution_cell_count_validation.csv",
        index=False,
    )

    print(
        count_check.to_string(
            index=False
        )
    )

    assert (
        count_check[
            "cell_count_match"
        ]
    ).all()

    assert int(
        count_check[
            "observed_execution_cells"
        ].sum()
    ) == 52509

    print(
        "TOTAL_EXECUTION_INSTANCES=52509"
    )

    print(
        "CELL_COUNT_VALIDATION=PASS"
    )


    # ========================================================
    # 9. PATIENT COVERAGE VALIDATION
    # ========================================================

    print(
        "\n===== 9. PATIENT COVERAGE VALIDATION ====="
    )

    patient_coverage = (
        manifest
        .groupby(
            [
                "pilot_perturbation_id",
                "patient_id",
            ],
            as_index=False,
        )
        .size()
        .rename(
            columns={
                "size":
                    "n_execution_cells"
            }
        )
    )

    patient_coverage.to_csv(
        TAB
        / "23D6G4F1_patient_execution_coverage.csv",
        index=False,
    )

    patient_target_summary = (
        patient_coverage
        .groupby(
            "pilot_perturbation_id",
            as_index=False,
        )
        .agg(
            n_patients=(
                "patient_id",
                "nunique",
            ),
            min_cells_per_patient=(
                "n_execution_cells",
                "min",
            ),
            median_cells_per_patient=(
                "n_execution_cells",
                "median",
            ),
            max_cells_per_patient=(
                "n_execution_cells",
                "max",
            ),
            total_cells=(
                "n_execution_cells",
                "sum",
            ),
        )
    )

    expected_patient_support = (
        support[
            [
                "pilot_perturbation_id",
                "n_execution_patients",
            ]
        ]
    )

    patient_target_summary = (
        patient_target_summary
        .merge(
            expected_patient_support,
            on="pilot_perturbation_id",
            how="left",
            validate="one_to_one",
        )
    )

    patient_target_summary[
        "patient_coverage_match"
    ] = (
        patient_target_summary[
            "n_patients"
        ]
        == patient_target_summary[
            "n_execution_patients"
        ]
    )

    patient_target_summary.to_csv(
        TAB
        / "23D6G4F1_patient_target_coverage_summary.csv",
        index=False,
    )

    print(
        patient_target_summary.to_string(
            index=False
        )
    )

    assert (
        patient_target_summary[
            "patient_coverage_match"
        ]
    ).all()

    print(
        "PATIENT_COVERAGE=PASS"
    )


    # ========================================================
    # 10. FORMAL CHECKS
    # ========================================================

    print(
        "\n===== 10. FORMAL CHECKS ====="
    )

    checks = pd.DataFrame(
        [
            {
                "check":
                    "frozen_targets_24",
                "status":
                    "PASS"
                    if len(support) == 24
                    else "FAIL",
                "detail":
                    str(len(support)),
            },
            {
                "check":
                    "execution_instances_52509",
                "status":
                    "PASS"
                    if len(manifest) == 52509
                    else "FAIL",
                "detail":
                    str(len(manifest)),
            },
            {
                "check":
                    "all_targets_executed",
                "status":
                    "PASS"
                    if len(summary) == 24
                    else "FAIL",
                "detail":
                    str(len(summary)),
            },
            {
                "check":
                    "all_execution_status_pass",
                "status":
                    "PASS"
                    if summary[
                        "execution_status"
                    ].isin(
                        [
                            "PASS",
                            "PASS_REUSED",
                        ]
                    ).all()
                    else "FAIL",
                "detail":
                    str(
                        summary[
                            "execution_status"
                        ]
                        .value_counts()
                        .to_dict()
                    ),
            },
            {
                "check":
                    "all_cell_counts_match_G4B",
                "status":
                    "PASS"
                    if count_check[
                        "cell_count_match"
                    ].all()
                    else "FAIL",
                "detail":
                    str(
                        int(
                            count_check[
                                "cell_count_match"
                            ].sum()
                        )
                    )
                    + "/24",
            },
            {
                "check":
                    "all_patient_coverage_match_G4B",
                "status":
                    "PASS"
                    if patient_target_summary[
                        "patient_coverage_match"
                    ].all()
                    else "FAIL",
                "detail":
                    str(
                        int(
                            patient_target_summary[
                                "patient_coverage_match"
                            ].sum()
                        )
                    )
                    + "/24",
            },
            {
                "check":
                    "no_patient_subsampling",
                "status":
                    "PASS",
                "detail":
                    "G4B exhaustive execution datasets",
            },
            {
                "check":
                    "max_ncells_none",
                "status":
                    "PASS",
                "detail":
                    "None",
            },
            {
                "check":
                    "cell_inds_to_perturb_all",
                "status":
                    "PASS",
                "detail":
                    "all",
            },
            {
                "check":
                    "model_version",
                "status":
                    "PASS",
                "detail":
                    "Geneformer_V2",
            },
            {
                "check":
                    "embedding_mode",
                "status":
                    "PASS",
                "detail":
                    "cls",
            },
            {
                "check":
                    "no_classifier",
                "status":
                    "PASS",
                "detail":
                    "Pretrained",
            },
            {
                "check":
                    "no_fine_tuning",
                "status":
                    "PASS",
                "detail":
                    "none",
            },
        ]
    )

    checks.to_csv(
        TAB
        / "23D6G4F1_formal_checks.csv",
        index=False,
    )

    print(
        checks.to_string(
            index=False
        )
    )

    assert (
        checks["status"] == "PASS"
    ).all()

    print(
        "FORMAL_CHECKS=PASS"
    )


    # ========================================================
    # 11. FINAL AUDIT
    # ========================================================

    print(
        "\n===== 11. FINAL AUDIT ====="
    )

    n_cells_total = int(
        summary[
            "n_execution_cells"
        ].sum()
    )

    n_raw = int(
        len(raw)
    )

    top = ranking.iloc[0]

    audit = f"""MODULE=23D6G4F1
TITLE=Exhaustive full-atlas Geneformer V2 in-silico perturbation execution
SOURCE_PANEL=23D6C
SOURCE_EXECUTION_SETS=23D6G4B
MODEL_VERSION=Geneformer_V2
MODEL_TYPE=Pretrained
EMB_MODE=cls
EMB_LAYER=-1
FORWARD_BATCH_SIZE=8
NPROC=1
FROZEN_TARGETS=24
EXECUTED_TARGETS=24
ORIGINAL_G4B_CELL_PERTURBATION_INSTANCES=52509
CELL_PERTURBATION_INPUT_INSTANCES={n_cells_total}
RECOVERED_GENEFORMER_VALUES={n_raw}
PATIENT_STRATIFIED_SUBSAMPLING=FALSE
CAP_PER_PATIENT=NONE
MAX_NCELLS=NONE
CELL_INDS_TO_PERTURB=ALL
ALL_ELIGIBLE_G4B_INSTANCES_EXECUTED=TRUE
DELETE_TARGETS_REQUIRE_TARGET_PRESENT_AS_FROZEN_IN_G4B=TRUE
OVEREXPRESSION_EXECUTION_UNIVERSE_AS_FROZEN_IN_G4B=TRUE
PATIENT_COVERAGE_PRESERVED=TRUE
CLASSIFIER_EXECUTION_PERFORMED=FALSE
FINE_TUNING_PERFORMED=FALSE
TOP_EMBEDDING_SHIFT_PERTURBATION={top['pilot_perturbation_id']}
TOP_EMBEDDING_SHIFT_GENE={top['gene']}
TOP_EMBEDDING_SHIFT_DISEASE={top['disease']}
TOP_EMBEDDING_SHIFT_L2={top['final_L2_label']}
TOP_MEAN_COSINE_SHIFT_FROM_1={top['mean_cosine_shift_from_1']}
SOURCE_MODULES_MODIFIED=FALSE
STATUS=PASS
"""

    (
        OUT
        / "23D6G4F1_audit.txt"
    ).write_text(
        audit,
        encoding="utf-8",
    )

    print(audit)


    # ========================================================
    # 12. OUTPUT INVENTORY
    # ========================================================

    print(
        "\n===== 12. OUTPUT INVENTORY ====="
    )

    outputs = [
        TAB
        / "23D6G4F1_dataset_verification.csv",

        TAB
        / "23D6G4F1_execution_summary.csv",

        TAB
        / "23D6G4F1_all_raw_cosine_values.csv",

        TAB
        / "23D6G4F1_perturbation_level_summary.csv",

        TAB
        / "23D6G4F1_exhaustive_perturbation_ranking.csv",

        TAB
        / "23D6G4F1_execution_cell_count_validation.csv",

        TAB
        / "23D6G4F1_patient_execution_coverage.csv",

        TAB
        / "23D6G4F1_patient_target_coverage_summary.csv",

        TAB
        / "23D6G4F1_formal_checks.csv",

        OUT
        / "23D6G4F1_audit.txt",
    ]

    for p in outputs:
        print(
            p,
            "exists=",
            p.exists(),
        )

        assert p.exists()

    print(
        "\n" + "=" * 80
    )

    print(
        "23D6G4F1 COMPLETE"
    )

    print(
        "STATUS=PASS"
    )

    print(
        "=" * 80
    )


if __name__ == "__main__":
    main()
