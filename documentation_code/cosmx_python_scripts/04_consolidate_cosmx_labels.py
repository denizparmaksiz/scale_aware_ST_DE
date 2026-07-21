#!/usr/bin/env python3
"""
Consolidate MapMyCells and scANVI label-transfer results for CosMx cells.

The script reproduces the hybrid labeling logic used in the original notebook:

1. Use class-level labels for most taxonomy classes.
2. Replace the broad Astro-Epen, OPC-Oligo, and Vascular classes with
   subclass-level labels.
3. Calculate a confidence-weighted consensus across MapMyCells runs.
4. Retain high-confidence MapMyCells assignments.
5. Rescue lower-confidence assignments when scANVI agrees exactly or when
   both methods agree at the broad Glut/GABA neurotransmitter-class level.
6. Mark all remaining cells as ambiguous.

Unlike the exploratory notebook, this script does not attempt to restore gene
names on the gene-reduced mapping objects. Mapping outputs are used only as
sources of cell labels. Final labels are joined to the original, full
preprocessed CosMx AnnData object, preserving its complete expression matrix
and gene metadata.

Expected inputs
---------------
Expression object:
- Full preprocessed CosMx H5AD.
- Raw counts and gene metadata are preserved as prepared upstream.
- ``adata.obs_names`` contains unique cell identifiers.

scANVI object:
- H5AD produced by ``run_scanvi_cosmx.py``.
- ``adata.obs`` contains ``scanvi_class`` and ``scanvi_subclass``.

MapMyCells results:
- CSV files produced by ``run_mapmycells_cosmx.py``.
- Each CSV contains cell IDs, class/subclass labels, and corresponding
  bootstrapping probabilities.

Example
-------
python consolidate_cosmx_labels.py \
    --expression-query data/processed/cosmx/cosmx_data_rawx_ensembl.h5ad \
    --scanvi-query results/scanvi/cosmx_scanvi_class_subclass.h5ad \
    --mapmycells-dir results/mapmycells \
    --output results/cosmx/cosmx_data_ct_final.h5ad
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path
import sys
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd


# For these broad Allen classes, subclass labels are more biologically useful
# than the class label and were used in the original consensus analysis.
SUBCLASS_REPLACEMENT_CLASSES: tuple[str, ...] = (
    "30 Astro-Epen",
    "31 OPC-Oligo",
    "33 Vascular",
)

# Mapping from the taxonomy label after removal of its numeric prefix to the
# broad analysis annotations used downstream. This dictionary was present in
# the original notebook but was not applied there; the script stores the
# mapped value in ``final_annotation`` while retaining the unmodified
# consensus taxonomy label in ``final_label``.
CLASS_TO_ANNOTATION: dict[str, str] = {
    "Oligo": "OLG",
    "Astro-TE": "Astro",
    "IT-ET": "IT-ET Glut",
    "NP-CT-L6b": "NP-CT-L6b Glut",
    "CTX-CGE": "GABA",
    "CTX-MGE": "GABA",
    "Peri": "PC",
    "Immune": "Immune",
    "Astro-NT": "Astro",
    "VLMC": "VLMC",
    "OPC": "OPC",
    "Astroependymal": "Astro",
    "Endo": "EC",
    "OB-CR": "OB-CR Glut",
    "LSX": "GABA",
    "CNU-LGE": "GABA",
    "TH": "TH Glut",
    "CNU-MGE": "GABA",
    "Ependymal": "Ependymal",
    "SMC": "VSMC",
    "MH-LH": "MH-LH Glut",
    "OB-IMN": "GABA",
    "DG-IMN": "DG-IMN Glut",
    "Tanycyte": "Ependymal",
    "ABC": "ABC",
    "CHOR": "Ependymal",
    "Astro-OLF": "Astro",
    "Hypendymal": "Ependymal",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Consolidate MapMyCells and scANVI labels for CosMx cells."
    )
    parser.add_argument(
        "--expression-query",
        type=Path,
        required=True,
        help=(
            "Full preprocessed CosMx H5AD to which final labels will be added. "
            "This object supplies the final expression matrix and gene metadata."
        ),
    )
    parser.add_argument(
        "--scanvi-query",
        type=Path,
        required=True,
        help="CosMx H5AD containing scanvi_class and scanvi_subclass.",
    )
    parser.add_argument(
        "--mapmycells-dir",
        type=Path,
        required=True,
        help="Directory containing MapMyCells CSV result files.",
    )
    parser.add_argument(
        "--mapmycells-pattern",
        default="cosmx_mmc_*.csv",
        help="Glob used to locate MapMyCells CSV files. Default: cosmx_mmc_*.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output H5AD containing final consensus labels.",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=None,
        help="Optional path for a cell-level consensus summary CSV.",
    )
    parser.add_argument(
        "--weighted-threshold",
        type=float,
        default=0.70,
        help="Minimum weighted agreement for high-confidence MMC. Default: 0.70.",
    )
    parser.add_argument(
        "--dominant-threshold",
        type=float,
        default=0.70,
        help="Minimum mean confidence of dominant MMC label. Default: 0.70.",
    )
    parser.add_argument(
        "--retain-ambiguous",
        action="store_true",
        help=(
            "Retain cells without an accepted final label. By default, ambiguous "
            "cells are removed to match the original notebook."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing output files.",
    )
    return parser.parse_args()


def require_file(path: Path, description: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{description} does not exist: {path}")
    return path


def package_version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def validate_threshold(value: float, name: str) -> None:
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must be between 0 and 1; received {value}.")


def choose_class_or_subclass(
    frame: pd.DataFrame,
    *,
    class_col: str,
    subclass_col: str,
    class_probability_col: str,
    subclass_probability_col: str,
) -> pd.DataFrame:
    """Select class or subclass label/probability using the original rule."""
    required = {
        class_col,
        subclass_col,
        class_probability_col,
        subclass_probability_col,
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(
            "MapMyCells result is missing required columns: " + ", ".join(missing)
        )

    use_subclass = frame[class_col].isin(SUBCLASS_REPLACEMENT_CLASSES)
    labels = frame[class_col].where(~use_subclass, frame[subclass_col])
    probabilities = frame[class_probability_col].where(
        ~use_subclass,
        frame[subclass_probability_col],
    )
    return pd.DataFrame(
        {
            "label": labels.astype("string"),
            "probability": pd.to_numeric(probabilities, errors="coerce"),
        },
        index=frame.index,
    )


def read_mapmycells_csv(path: Path) -> pd.DataFrame:
    """Read a verbose MapMyCells CSV and index it by cell ID."""
    frame = pd.read_csv(path, skiprows=3)
    if "cell_id" not in frame.columns:
        # Some versions may omit the three-line metadata preamble.
        frame = pd.read_csv(path)
    if "cell_id" not in frame.columns:
        raise ValueError(f"MapMyCells CSV has no cell_id column: {path}")
    if frame["cell_id"].duplicated().any():
        raise ValueError(f"MapMyCells CSV contains duplicate cell IDs: {path}")
    return frame.set_index("cell_id")


def config_name_from_path(path: Path) -> str:
    stem = path.stem
    prefix = "cosmx_mmc_"
    return stem[len(prefix):] if stem.startswith(prefix) else stem


def load_mapmycells_runs(
    *,
    result_dir: Path,
    pattern: str,
    target_cells: pd.Index,
) -> tuple[pd.DataFrame, list[Path]]:
    paths = sorted(result_dir.glob(pattern))
    if not paths:
        raise FileNotFoundError(
            f"No MapMyCells CSV files matched '{pattern}' in {result_dir}"
        )

    pieces: list[pd.DataFrame] = []
    seen_names: set[str] = set()

    for path in paths:
        name = config_name_from_path(path)
        if name in seen_names:
            raise ValueError(f"Duplicate MapMyCells configuration name: {name}")
        seen_names.add(name)

        frame = read_mapmycells_csv(path)
        selected = choose_class_or_subclass(
            frame,
            class_col="class_name",
            subclass_col="subclass_name",
            class_probability_col="class_bootstrapping_probability",
            subclass_probability_col="subclass_bootstrapping_probability",
        )
        selected = selected.reindex(target_cells)
        selected.columns = [f"{name}_final_class", f"{name}_bp"]
        pieces.append(selected)

    combined = pd.concat(pieces, axis=1)
    return combined, paths


def construct_scanvi_label(scanvi: ad.AnnData) -> pd.Series:
    required = {"scanvi_class", "scanvi_subclass"}
    missing = sorted(required - set(scanvi.obs.columns))
    if missing:
        raise ValueError(
            "scANVI query is missing required obs columns: " + ", ".join(missing)
        )

    use_subclass = scanvi.obs["scanvi_class"].isin(SUBCLASS_REPLACEMENT_CLASSES)
    label = scanvi.obs["scanvi_class"].where(
        ~use_subclass,
        scanvi.obs["scanvi_subclass"],
    )
    return label.astype("string").rename("scanvi_label")


def calculate_mmc_consensus(
    mmc: pd.DataFrame,
    *,
    weighted_threshold: float,
    dominant_threshold: float,
) -> pd.DataFrame:
    label_cols = [column for column in mmc.columns if column.endswith("_final_class")]
    probability_cols = [column for column in mmc.columns if column.endswith("_bp")]

    label_roots = [column.removesuffix("_final_class") for column in label_cols]
    probability_roots = [column.removesuffix("_bp") for column in probability_cols]
    if label_roots != probability_roots:
        raise ValueError(
            "MapMyCells label and probability columns could not be paired by run."
        )

    results: list[dict[str, Any]] = []

    for _, row in mmc.iterrows():
        labels: list[str] = []
        weights: list[float] = []

        for label_col, probability_col in zip(label_cols, probability_cols):
            label = row[label_col]
            weight = row[probability_col]
            if pd.notna(label) and pd.notna(weight):
                labels.append(str(label))
                weights.append(float(weight))

        if not labels:
            results.append(
                {
                    "dominant_label": pd.NA,
                    "raw_agreement": np.nan,
                    "weighted_agreement": np.nan,
                    "dominant_confidence": np.nan,
                    "mmc_high_conf": False,
                }
            )
            continue

        label_array = np.asarray(labels, dtype=object)
        weight_array = np.asarray(weights, dtype=float)
        if np.any(weight_array < 0):
            raise ValueError("MapMyCells probabilities must be nonnegative.")
        if weight_array.sum() == 0:
            weighted_agreement = np.nan
        else:
            weighted_support = {
                label: float(weight_array[label_array == label].sum())
                for label in np.unique(label_array)
            }
            dominant_label = max(weighted_support, key=weighted_support.get)
            dominant_mask = label_array == dominant_label
            weighted_agreement = float(
                weight_array[dominant_mask].sum() / weight_array.sum()
            )

        if weight_array.sum() == 0:
            # Preserve a deterministic dominant label even when all supplied
            # probabilities are zero.
            values, counts = np.unique(label_array, return_counts=True)
            dominant_label = str(values[np.argmax(counts)])
            dominant_mask = label_array == dominant_label

        raw_agreement = float(dominant_mask.mean())
        dominant_confidence = float(weight_array[dominant_mask].mean())
        high_confidence = bool(
            pd.notna(weighted_agreement)
            and weighted_agreement >= weighted_threshold
            and dominant_confidence >= dominant_threshold
        )

        results.append(
            {
                "dominant_label": dominant_label,
                "raw_agreement": raw_agreement,
                "weighted_agreement": weighted_agreement,
                "dominant_confidence": dominant_confidence,
                "mmc_high_conf": high_confidence,
            }
        )

    return pd.DataFrame(results, index=mmc.index)


def has_broad_neurotransmitter_agreement(
    dominant_label: object,
    scanvi_label: object,
) -> bool:
    dominant = str(dominant_label)
    scanvi = str(scanvi_label)
    return ("Glut" in dominant and "Glut" in scanvi) or (
        "GABA" in dominant and "GABA" in scanvi
    )


def assign_final_labels(consensus: pd.DataFrame) -> pd.DataFrame:
    consensus = consensus.copy()
    consensus["scanvi_agrees"] = (
        consensus["scanvi_label"].notna()
        & consensus["dominant_label"].notna()
        & consensus["scanvi_label"].eq(consensus["dominant_label"])
    )

    final_labels: list[object] = []
    final_statuses: list[str] = []

    for _, row in consensus.iterrows():
        dominant = row["dominant_label"]
        scanvi_label = row["scanvi_label"]

        if bool(row["mmc_high_conf"]):
            final_labels.append(dominant)
            final_statuses.append("high_conf_mmc")
        elif bool(row["scanvi_agrees"]):
            final_labels.append(dominant)
            final_statuses.append("rescued_by_scanvi")
        elif (
            pd.notna(dominant)
            and pd.notna(scanvi_label)
            and has_broad_neurotransmitter_agreement(dominant, scanvi_label)
        ):
            final_labels.append(dominant)
            final_statuses.append("rescued_by_broad_class")
        else:
            final_labels.append(pd.NA)
            final_statuses.append("ambiguous")

    consensus["final_label"] = pd.Series(
        final_labels,
        index=consensus.index,
        dtype="string",
    )
    consensus["final_status"] = pd.Categorical(final_statuses)

    # Match the notebook's extraction of the taxonomy term after the leading
    # numeric class code, while handling labels with or without that prefix.
    consensus["final_class"] = consensus["final_label"].str.replace(
        r"^\d+\s+",
        "",
        regex=True,
    )
    consensus["final_annotation"] = consensus["final_class"].map(
        CLASS_TO_ANNOTATION
    )
    return consensus


def main() -> None:
    args = parse_args()
    validate_threshold(args.weighted_threshold, "--weighted-threshold")
    validate_threshold(args.dominant_threshold, "--dominant-threshold")

    expression_path = require_file(args.expression_query, "Expression query H5AD")
    scanvi_path = require_file(args.scanvi_query, "scANVI query H5AD")
    mapmycells_dir = args.mapmycells_dir.expanduser().resolve()
    if not mapmycells_dir.is_dir():
        raise NotADirectoryError(
            f"MapMyCells result directory does not exist: {mapmycells_dir}"
        )

    output_path = args.output.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and not args.overwrite:
        raise FileExistsError(
            f"Output already exists: {output_path}. Use --overwrite to replace it."
        )

    summary_path = (
        args.summary_csv.expanduser().resolve()
        if args.summary_csv is not None
        else output_path.with_name(f"{output_path.stem}_label_summary.csv")
    )
    if summary_path.exists() and not args.overwrite:
        raise FileExistsError(
            f"Summary already exists: {summary_path}. Use --overwrite to replace it."
        )
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    expression = ad.read_h5ad(expression_path)
    scanvi = ad.read_h5ad(scanvi_path)

    if not expression.obs_names.is_unique:
        raise ValueError("Expression query contains duplicate cell identifiers.")
    if not scanvi.obs_names.is_unique:
        raise ValueError("scANVI query contains duplicate cell identifiers.")

    missing_scanvi = expression.obs_names.difference(scanvi.obs_names)
    extra_scanvi = scanvi.obs_names.difference(expression.obs_names)
    if len(missing_scanvi) or len(extra_scanvi):
        raise ValueError(
            "Expression and scANVI objects do not contain the same cells. "
            f"Missing from scANVI: {len(missing_scanvi):,}; "
            f"extra in scANVI: {len(extra_scanvi):,}."
        )

    scanvi_label = construct_scanvi_label(scanvi).reindex(expression.obs_names)
    mmc, mmc_paths = load_mapmycells_runs(
        result_dir=mapmycells_dir,
        pattern=args.mapmycells_pattern,
        target_cells=expression.obs_names,
    )

    missing_by_run = {
        column.removesuffix("_final_class"): int(mmc[column].isna().sum())
        for column in mmc.columns
        if column.endswith("_final_class")
    }

    consensus = calculate_mmc_consensus(
        mmc,
        weighted_threshold=args.weighted_threshold,
        dominant_threshold=args.dominant_threshold,
    )
    consensus["scanvi_label"] = scanvi_label
    consensus = assign_final_labels(consensus)

    overall_scanvi_agreement = float(consensus["scanvi_agrees"].mean())
    status_counts = {
        str(key): int(value)
        for key, value in consensus["final_status"].value_counts().items()
    }

    # Join mapping metrics and labels to the original full expression object.
    columns_to_replace = consensus.columns.intersection(expression.obs.columns)
    if len(columns_to_replace):
        expression.obs = expression.obs.drop(columns=list(columns_to_replace))
    expression.obs = expression.obs.join(consensus)

    n_input_cells = expression.n_obs
    if not args.retain_ambiguous:
        expression = expression[expression.obs["final_label"].notna()].copy()

    expression.uns["consensus_labeling"] = {
        "method": "confidence-weighted MapMyCells consensus with scANVI rescue",
        "subclass_replacement_classes": list(SUBCLASS_REPLACEMENT_CLASSES),
        "weighted_threshold": args.weighted_threshold,
        "dominant_threshold": args.dominant_threshold,
        "broad_rescue_categories": ["Glut", "GABA"],
        "ambiguous_cells_retained": bool(args.retain_ambiguous),
        "mapmycells_runs": [config_name_from_path(path) for path in mmc_paths],
        "scanvi_columns": ["scanvi_class", "scanvi_subclass"],
        "final_columns": [
            "dominant_label",
            "raw_agreement",
            "weighted_agreement",
            "dominant_confidence",
            "mmc_high_conf",
            "scanvi_label",
            "scanvi_agrees",
            "final_label",
            "final_status",
            "final_class",
            "final_annotation",
        ],
    }
    expression.write_h5ad(output_path)
    consensus.to_csv(summary_path, index_label="cell_id")

    provenance: dict[str, Any] = {
        "expression_query": str(expression_path),
        "scanvi_query": str(scanvi_path),
        "mapmycells_files": [str(path.resolve()) for path in mmc_paths],
        "output_h5ad": str(output_path),
        "summary_csv": str(summary_path),
        "parameters": {
            "weighted_threshold": args.weighted_threshold,
            "dominant_threshold": args.dominant_threshold,
            "retain_ambiguous": bool(args.retain_ambiguous),
            "subclass_replacement_classes": list(SUBCLASS_REPLACEMENT_CLASSES),
        },
        "dimensions": {
            "input_cells": int(n_input_cells),
            "output_cells": int(expression.n_obs),
            "output_genes": int(expression.n_vars),
        },
        "overall_scanvi_exact_agreement": overall_scanvi_agreement,
        "final_status_counts": status_counts,
        "missing_mapmycells_labels_by_run": missing_by_run,
        "python_version": sys.version,
        "package_versions": {
            "anndata": package_version("anndata"),
            "numpy": package_version("numpy"),
            "pandas": package_version("pandas"),
        },
    }
    provenance_path = output_path.with_name(
        f"{output_path.stem}_run_configuration.json"
    )
    provenance_path.write_text(
        json.dumps(provenance, indent=2),
        encoding="utf-8",
    )

    print("\nConsensus labeling completed successfully.")
    print(f"Overall exact scANVI/MMC agreement: {overall_scanvi_agreement:.3f}")
    print("Final status counts:")
    for status, count in status_counts.items():
        print(f"  {status}: {count:,}")
    print(f"Output: {output_path}")
    print(f"Cell-level summary: {summary_path}")
    print(f"Final dimensions: {expression.n_obs:,} cells x {expression.n_vars:,} genes")


if __name__ == "__main__":
    main()
