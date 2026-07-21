#!/usr/bin/env python3
"""
Map preprocessed CosMx cells to the Allen Whole Mouse Brain taxonomy with scANVI.

The script reconstructs a panel-matched WMB-10Xv2 scRNA-seq reference from the
Allen Brain Cell Atlas cache, trains or reloads class- and subclass-level scANVI
models, and writes the predicted labels to a CosMx AnnData object.

Expected query
--------------
- ``adata.X`` contains raw transcript counts.
- ``adata.var_names`` contains mouse Ensembl gene identifiers.
- ``adata.obs_names`` contains unique cell identifiers.

Allen reference
---------------
The reference uses the WMB-10Xv2 raw expression matrices and CCN20230722
class/subclass annotations distributed through ``abc_atlas_access``. Large
Allen files are downloaded once into ``--cache-dir`` and reused on later runs.

Example
-------
python run_scanvi_cosmx.py \
    --query data/processed/cosmx/cosmx_data_rawx_ensembl.h5ad \
    --cache-dir references/abc_cache \
    --output-dir results/scanvi
"""

from __future__ import annotations

import argparse
import gc
import importlib.metadata
import json
from pathlib import Path
import sys
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scvi

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache


# Stable CCN20230722 class accessions excluded from both the Allen reference and
# the anatomically restricted MapMyCells analysis.
CLASS_ACCESSIONS_TO_DROP: tuple[str, ...] = (
    "CS20230722_CLAS_11",
    "CS20230722_CLAS_12",
    "CS20230722_CLAS_13",
    "CS20230722_CLAS_14",
    "CS20230722_CLAS_15",
    "CS20230722_CLAS_16",
    "CS20230722_CLAS_19",
    "CS20230722_CLAS_20",
    "CS20230722_CLAS_21",
    "CS20230722_CLAS_22",
    "CS20230722_CLAS_23",
    "CS20230722_CLAS_24",
    "CS20230722_CLAS_25",
    "CS20230722_CLAS_26",
    "CS20230722_CLAS_27",
    "CS20230722_CLAS_28",
    "CS20230722_CLAS_29",
    "CS20230722_CLAS_32",
)

# WMB-10Xv2 regional matrices retained in the original CosMx scANVI analysis.
MATRIX_PREFIXES: tuple[str, ...] = (
    "WMB-10Xv2-CTXsp",
    "WMB-10Xv2-HPF",
    "WMB-10Xv2-Isocortex-1",
    "WMB-10Xv2-Isocortex-2",
    "WMB-10Xv2-Isocortex-3",
    "WMB-10Xv2-Isocortex-4",
    "WMB-10Xv2-TH",
)

LABEL_LEVELS: tuple[str, ...] = ("class", "subclass")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run class- and subclass-level scANVI mapping on CosMx data."
    )
    parser.add_argument(
        "--query",
        type=Path,
        required=True,
        help="CosMx H5AD with raw counts in X and mouse Ensembl IDs in var_names.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        required=True,
        help="Local ABC Atlas cache; downloaded Allen files are reused here.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for the reference, models, predictions, and provenance.",
    )
    parser.add_argument(
        "--output-name",
        default="cosmx_scanvi_class_subclass.h5ad",
        help="Filename for the mapped CosMx AnnData.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=512,
        help="Training batch size. Default: 512.",
    )
    parser.add_argument(
        "--max-epochs",
        type=int,
        default=5,
        help="Training epochs for each model. Default: 5, matching the analysis.",
    )
    parser.add_argument(
        "--accelerator",
        default="cpu",
        choices=("cpu", "gpu", "auto"),
        help="Training accelerator. Default: cpu.",
    )
    parser.add_argument(
        "--rebuild-reference",
        action="store_true",
        help="Rebuild the panel-matched reference even if it already exists.",
    )
    parser.add_argument(
        "--retrain-models",
        action="store_true",
        help="Retrain models even if saved model directories already exist.",
    )
    parser.add_argument(
        "--overwrite-output",
        action="store_true",
        help="Replace an existing mapped H5AD output.",
    )
    return parser.parse_args()


def package_version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def require_file(path: Path, description: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{description} does not exist: {path}")
    return path


def resolve_dropped_class_names(
    cache: AbcProjectCache,
    accessions: tuple[str, ...],
) -> list[str]:
    """Resolve stable class accessions to display names used in cell metadata."""
    terms = cache.get_metadata_dataframe(
        directory="WMB-taxonomy",
        file_name="cluster_annotation_term",
    )
    class_terms = terms.loc[
        terms["cluster_annotation_term_set_name"].eq("class"),
        ["label", "name"],
    ].drop_duplicates("label")
    accession_to_name = class_terms.set_index("label")["name"]

    missing = sorted(set(accessions) - set(accession_to_name.index))
    if missing:
        raise ValueError(
            "Class accessions were not found in WMB-taxonomy: " + ", ".join(missing)
        )

    return accession_to_name.loc[list(accessions)].astype(str).tolist()


def locate_cell_metadata_csv(
    cache: AbcProjectCache,
    cache_dir: Path,
) -> Path:
    """Download WMB-10X metadata if needed and return the cell-metadata CSV."""
    cache.get_directory_metadata("WMB-10X")
    matches = sorted(cache_dir.rglob("cell_metadata_with_cluster_annotation.csv"))
    if not matches:
        raise FileNotFoundError(
            "ABC cache did not contain cell_metadata_with_cluster_annotation.csv "
            "after downloading WMB-10X metadata."
        )
    if len(matches) > 1:
        print(f"Found {len(matches)} metadata versions; using {matches[-1]}")
    return matches[-1]


def create_filtered_cell_metadata(
    metadata_csv: Path,
    output_parquet: Path,
    dropped_class_names: list[str],
    matrix_prefixes: tuple[str, ...],
) -> pd.DataFrame:
    """Stream the 4-million-cell CSV and retain only required cells/columns."""
    usecols = ["cell_label", "class", "subclass", "feature_matrix_label"]
    writer: pq.ParquetWriter | None = None

    try:
        for chunk in pd.read_csv(
            metadata_csv,
            usecols=usecols,
            chunksize=250_000,
            dtype={column: "string" for column in usecols},
        ):
            keep = (
                ~chunk["class"].isin(dropped_class_names)
                & chunk["feature_matrix_label"].isin(matrix_prefixes)
            )
            chunk = chunk.loc[keep].copy()
            if chunk.empty:
                continue

            table = pa.Table.from_pandas(chunk, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(output_parquet, table.schema)
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()

    if writer is None or not output_parquet.is_file():
        raise RuntimeError("Class/matrix filtering retained no Allen reference cells.")

    return pd.read_parquet(output_parquet)


def get_expression_paths(
    cache: AbcProjectCache,
    matrix_prefixes: tuple[str, ...],
) -> list[Path]:
    paths: list[Path] = []
    for prefix in matrix_prefixes:
        file_name = f"{prefix}/raw"
        path = Path(cache.get_file_path("WMB-10Xv2", file_name))
        paths.append(path)
    return paths


def common_genes(query_path: Path, first_reference_path: Path) -> np.ndarray:
    query = ad.read_h5ad(query_path, backed="r")
    reference = ad.read_h5ad(first_reference_path, backed="r")
    try:
        genes = np.asarray(
            sorted(set(query.var_names.astype(str)) & set(reference.var_names.astype(str)))
        )
    finally:
        query.file.close()
        reference.file.close()

    if len(genes) == 0:
        raise ValueError("The CosMx query and Allen reference share no gene identifiers.")
    return genes


def build_reference(
    expression_paths: list[Path],
    genes: np.ndarray,
    metadata: pd.DataFrame,
    output_path: Path,
) -> ad.AnnData:
    """Build and immediately save the panel-matched Allen reference."""
    pieces: list[ad.AnnData] = []

    for path in expression_paths:
        print(f"Reading and subsetting {path.name}")
        backed = ad.read_h5ad(path, backed="r")
        try:
            piece = backed[:, genes].to_memory()
        finally:
            backed.file.close()
        pieces.append(piece)

    reference = ad.concat(pieces, axis=0, merge="same")
    del pieces
    gc.collect()

    labels = (
        metadata[["cell_label", "class", "subclass"]]
        .drop_duplicates("cell_label")
        .set_index("cell_label")
    )
    reference.obs["class"] = labels.reindex(reference.obs_names)["class"]
    reference.obs["subclass"] = labels.reindex(reference.obs_names)["subclass"]

    missing = reference.obs[list(LABEL_LEVELS)].isna().any(axis=1)
    if missing.any():
        print(f"Dropping {int(missing.sum()):,} reference cells without taxonomy labels")
        reference = reference[~missing].copy()

    for level in LABEL_LEVELS:
        reference.obs[level] = reference.obs[level].astype("category")

    reference.uns["scanvi_reference"] = {
        "atlas_directory": "WMB-10Xv2",
        "taxonomy": "CCN20230722",
        "matrix_prefixes": list(MATRIX_PREFIXES),
        "class_accessions_dropped": list(CLASS_ACCESSIONS_TO_DROP),
        "n_shared_genes": int(reference.n_vars),
    }
    reference.write_h5ad(output_path)
    return reference


def load_or_build_reference(
    *,
    cache: AbcProjectCache,
    cache_dir: Path,
    query_path: Path,
    output_dir: Path,
    rebuild: bool,
) -> tuple[ad.AnnData, Path, Path, list[Path], list[str]]:
    reference_path = output_dir / "wmb10xv2_cosmx_panel_reference.h5ad"
    metadata_parquet = output_dir / "wmb10xv2_cosmx_reference_metadata.parquet"

    dropped_names = resolve_dropped_class_names(cache, CLASS_ACCESSIONS_TO_DROP)
    expression_paths = get_expression_paths(cache, MATRIX_PREFIXES)

    if reference_path.exists() and not rebuild:
        print(f"Loading existing panel-matched reference: {reference_path}")
        return (
            ad.read_h5ad(reference_path),
            reference_path,
            metadata_parquet,
            expression_paths,
            dropped_names,
        )

    metadata_csv = locate_cell_metadata_csv(cache, cache_dir)
    print(f"Filtering Allen cell metadata: {metadata_csv}")
    metadata = create_filtered_cell_metadata(
        metadata_csv,
        metadata_parquet,
        dropped_names,
        MATRIX_PREFIXES,
    )

    genes = common_genes(query_path, expression_paths[0])
    print(f"Shared CosMx/WMB-10Xv2 genes: {len(genes):,}")
    reference = build_reference(expression_paths, genes, metadata, reference_path)
    return reference, reference_path, metadata_parquet, expression_paths, dropped_names


def load_query_for_mapping(query_path: Path, reference_genes: pd.Index) -> ad.AnnData:
    backed = ad.read_h5ad(query_path, backed="r")
    try:
        missing = reference_genes.difference(backed.var_names)
        if len(missing):
            raise ValueError(
                "The saved reference contains genes absent from the query: "
                + ", ".join(map(str, missing[:10]))
            )
        query = backed[:, reference_genes].to_memory()
    finally:
        backed.file.close()
    return query


def train_or_load_model(
    *,
    reference: ad.AnnData,
    level: str,
    model_dir: Path,
    retrain: bool,
    accelerator: str,
    batch_size: int,
    max_epochs: int,
) -> scvi.model.SCANVI:
    scvi.model.SCANVI.setup_anndata(
        reference,
        labels_key=level,
        unlabeled_category="unlabeled",
    )

    if model_dir.exists() and not retrain:
        print(f"Loading saved {level} model: {model_dir}")
        model = scvi.model.SCANVI.load(model_dir, adata=reference)
    else:
        print(f"Training {level} model")
        model = scvi.model.SCANVI(reference)
        model.train(
            accelerator=accelerator,
            devices=1,
            batch_size=batch_size,
            max_epochs=max_epochs,
            enable_progress_bar=True,
            datasplitter_kwargs={"num_workers": 0},
        )
        model.save(model_dir, overwrite=True)

    return model


def predict_level(
    *,
    model: scvi.model.SCANVI,
    query: ad.AnnData,
    level: str,
) -> np.ndarray:
    query.obs[level] = "unlabeled"
    scvi.model.SCANVI.setup_anndata(
        query,
        labels_key=level,
        unlabeled_category="unlabeled",
    )

    # The original analysis used CPU inference because the available Quadro
    # P5000 was incompatible with the installed PyTorch CUDA build.
    model.module.cpu()
    return np.asarray(model.predict(query))


def main() -> None:
    args = parse_args()

    query_path = require_file(args.query, "CosMx query H5AD")
    cache_dir = args.cache_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    mapped_path = output_dir / args.output_name
    if mapped_path.exists() and not args.overwrite_output:
        raise FileExistsError(
            f"Mapped output already exists: {mapped_path}. "
            "Use --overwrite-output to replace it."
        )

    cache = AbcProjectCache.from_s3_cache(cache_dir=cache_dir)

    reference, reference_path, metadata_path, expression_paths, dropped_names = (
        load_or_build_reference(
            cache=cache,
            cache_dir=cache_dir,
            query_path=query_path,
            output_dir=output_dir,
            rebuild=args.rebuild_reference,
        )
    )

    query = load_query_for_mapping(query_path, reference.var_names)

    model_paths: dict[str, str] = {}
    prediction_counts: dict[str, dict[str, int]] = {}

    for level in LABEL_LEVELS:
        model_dir = output_dir / f"scanvi_wmb10xv2_cosmx_{level}"
        model = train_or_load_model(
            reference=reference,
            level=level,
            model_dir=model_dir,
            retrain=args.retrain_models,
            accelerator=args.accelerator,
            batch_size=args.batch_size,
            max_epochs=args.max_epochs,
        )
        predictions = predict_level(model=model, query=query, level=level)
        output_key = f"scanvi_{level}"
        query.obs[output_key] = pd.Categorical(predictions)
        prediction_counts[level] = {
            str(key): int(value)
            for key, value in query.obs[output_key].value_counts().items()
        }
        model_paths[level] = str(model_dir)

        del predictions, model
        gc.collect()

    query.uns["scanvi_mapping"] = {
        "reference": "Allen WMB-10Xv2",
        "taxonomy": "CCN20230722",
        "reference_file": reference_path.name,
        "class_accessions_dropped": list(CLASS_ACCESSIONS_TO_DROP),
        "matrix_prefixes": list(MATRIX_PREFIXES),
        "shared_gene_count": int(reference.n_vars),
        "prediction_columns": ["scanvi_class", "scanvi_subclass"],
    }
    query.write_h5ad(mapped_path)

    provenance: dict[str, Any] = {
        "query_path": str(query_path),
        "mapped_output": str(mapped_path),
        "abc_cache_dir": str(cache_dir),
        "reference_path": str(reference_path),
        "filtered_metadata_path": str(metadata_path),
        "expression_matrix_paths": [str(path) for path in expression_paths],
        "matrix_prefixes": list(MATRIX_PREFIXES),
        "class_accessions_dropped": list(CLASS_ACCESSIONS_TO_DROP),
        "class_names_dropped": dropped_names,
        "model_paths": model_paths,
        "parameters": {
            "accelerator": args.accelerator,
            "batch_size": args.batch_size,
            "max_epochs": args.max_epochs,
            "num_workers": 0,
            "batch_covariate": None,
        },
        "dimensions": {
            "reference_cells": int(reference.n_obs),
            "query_cells": int(query.n_obs),
            "shared_genes": int(reference.n_vars),
        },
        "prediction_counts": prediction_counts,
        "python_version": sys.version,
        "package_versions": {
            "abc_atlas_access": package_version("abc_atlas_access"),
            "anndata": package_version("anndata"),
            "numpy": package_version("numpy"),
            "pandas": package_version("pandas"),
            "pyarrow": package_version("pyarrow"),
            "scvi-tools": package_version("scvi-tools"),
            "torch": package_version("torch"),
        },
    }
    (output_dir / "scanvi_run_configuration.json").write_text(
        json.dumps(provenance, indent=2),
        encoding="utf-8",
    )

    print("\nscANVI mapping completed successfully.")
    print(f"Mapped output: {mapped_path}")
    print(f"Dimensions: {query.n_obs:,} cells x {query.n_vars:,} shared genes")


if __name__ == "__main__":
    main()
