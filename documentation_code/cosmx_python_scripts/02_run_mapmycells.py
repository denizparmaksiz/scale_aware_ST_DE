#!/usr/bin/env python3
"""
Run multiple MapMyCells configurations on a preprocessed CosMx dataset.

Expected input
--------------
An H5AD file in which:
1. adata.X contains raw transcript counts.
2. adata.var_names contains mouse Ensembl gene IDs.
3. adata.obs_names contains unique cell identifiers.

Reference files
---------------
- mouse_markers_230821.json
- precomputed_stats_ABC_revision_230821.h5

These reference files correspond to the Allen Brain Cell Atlas Whole Mouse
Brain (WMB-10X) taxonomy (CCN20230722).

The script maps cells to an anatomically restricted subset of the Allen
Whole Mouse Brain 10x taxonomy using several MapMyCells configurations.

Outputs
-------
For each mapping configuration, the script writes:
- CSV file of cell-type assignments and confidence metrics.
- HDF5 file containing extended mapping results.
- Mapping results stored in adata.obsm under a unique key.

Example
-------
python run_mapmycells_cosmx.py \
    --query data/processed/cosmx_data_rawx_ensembl.h5ad \
    --markers references/mapmycells/mouse_markers_230821.json \
    --precomputed-stats references/mapmycells/precomputed_stats_ABC_revision_230821.h5 \
    --output-dir results/mapmycells \
    --n-processors 16
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any


# WMB taxonomy classes excluded because they were not anatomically plausible
# for the sagittal CosMx fields of view analyzed in this study.
NODES_TO_DROP: list[tuple[str, str]] = [
    ("class", "CS20230722_CLAS_11"),
    ("class", "CS20230722_CLAS_12"),
    ("class", "CS20230722_CLAS_13"),
    ("class", "CS20230722_CLAS_14"),
    ("class", "CS20230722_CLAS_15"),
    ("class", "CS20230722_CLAS_16"),
    ("class", "CS20230722_CLAS_19"),
    ("class", "CS20230722_CLAS_20"),
    ("class", "CS20230722_CLAS_21"),
    ("class", "CS20230722_CLAS_22"),
    ("class", "CS20230722_CLAS_23"),
    ("class", "CS20230722_CLAS_24"),
    ("class", "CS20230722_CLAS_25"),
    ("class", "CS20230722_CLAS_26"),
    ("class", "CS20230722_CLAS_27"),
    ("class", "CS20230722_CLAS_28"),
    ("class", "CS20230722_CLAS_29"),
    ("class", "CS20230722_CLAS_32"),
]


# The four MapMyCells specifications retained for comparison.
#
# flatten=False:
#   Traverse the hierarchical taxonomy.
#
# flatten=True, bootstrap_iteration=1, bootstrap_factor=1.0:
#   Flat correlation mapping without marker bootstrapping.
#
# flatten=True, bootstrap_iteration=100, bootstrap_factor=0.9:
#   Legacy bootstrapped flat mapping.
#
# collapse_markers=True:
#   Compile the supplied marker sets into one marker set.
RUN_CONFIGS: list[dict[str, Any]] = [
    {
        "name": "hier_b100_f0p5_coll_subset",
        "flatten": False,
        "collapse_markers": True,
        "bootstrap_iteration": 100,
        "bootstrap_factor": 0.5,
    },
    {
        "name": "flat_b100_f0p9_subset",
        "flatten": True,
        "collapse_markers": False,
        "bootstrap_iteration": 100,
        "bootstrap_factor": 0.9,
    },
    {
        "name": "flat_b100_f0p9_coll_subset",
        "flatten": True,
        "collapse_markers": True,
        "bootstrap_iteration": 100,
        "bootstrap_factor": 0.9,
    },
    {
        "name": "flat_b1_f1p0_coll_subset",
        "flatten": True,
        "collapse_markers": True,
        "bootstrap_iteration": 1,
        "bootstrap_factor": 1.0,
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run multiple MapMyCells specifications on CosMx data."
    )

    parser.add_argument(
        "--query",
        type=Path,
        required=True,
        help="Preprocessed H5AD with raw counts in X and Ensembl IDs in var_names.",
    )
    parser.add_argument(
        "--markers",
        type=Path,
        required=True,
        help="Allen WMB marker lookup JSON.",
    )
    parser.add_argument(
        "--precomputed-stats",
        type=Path,
        required=True,
        help="Allen WMB precomputed-statistics HDF5 file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory in which mapping outputs will be written.",
    )
    parser.add_argument(
        "--n-processors",
        type=int,
        default=4,
        help="Number of MapMyCells worker processes. Default: 4.",
    )
    parser.add_argument(
        "--overwrite-obsm",
        action="store_true",
        help="Allow an existing MapMyCells obsm key to be overwritten.",
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


def build_command(
    *,
    query_path: Path,
    marker_path: Path,
    stats_path: Path,
    output_dir: Path,
    run_config: dict[str, Any],
    n_processors: int,
    overwrite_obsm: bool,
) -> list[str]:
    name = str(run_config["name"])

    command = [
        sys.executable,
        "-m",
        "cell_type_mapper.cli.from_specified_markers",
        "--query_path",
        str(query_path),
        "--precomputed_stats.path",
        str(stats_path),
        "--query_markers.serialized_lookup",
        str(marker_path),
        "--type_assignment.normalization",
        "raw",
        "--drop_level",
        "CCN20230722_SUPT",
        "--nodes_to_drop",
        repr(NODES_TO_DROP),
        "--cloud_safe",
        "False",
        "--verbose_csv",
        "True",
        "--flatten",
        str(run_config["flatten"]),
        "--query_markers.collapse_markers",
        str(run_config["collapse_markers"]),
        "--type_assignment.bootstrap_iteration",
        str(run_config["bootstrap_iteration"]),
        "--type_assignment.bootstrap_factor",
        str(run_config["bootstrap_factor"]),
        "--type_assignment.n_processors",
        str(n_processors),
        "--hdf5_result_path",
        str(output_dir / f"cosmx_mmc_{name}.h5"),
        "--csv_result_path",
        str(output_dir / f"cosmx_mmc_{name}.csv"),
        "--log_path",
        str(output_dir / f"cosmx_mmc_{name}.log"),
        "--obsm_key",
        f"cosmx_mmc_{name}",
    ]

    if overwrite_obsm:
        command.extend(["--obsm_clobber", "True"])

    return command


def main() -> None:
    args = parse_args()

    query_path = require_file(args.query, "Query H5AD")
    marker_path = require_file(args.markers, "Marker lookup")
    stats_path = require_file(
        args.precomputed_stats,
        "Precomputed-statistics file",
    )

    if args.n_processors < 1:
        raise ValueError("--n-processors must be at least 1.")

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Prevent each worker process from spawning additional BLAS threads.
    run_environment = os.environ.copy()
    run_environment["OMP_NUM_THREADS"] = "1"
    run_environment["MKL_NUM_THREADS"] = "1"
    run_environment["NUMEXPR_NUM_THREADS"] = "1"

    provenance = {
        "query_path": str(query_path),
        "marker_path": str(marker_path),
        "precomputed_stats_path": str(stats_path),
        "nodes_to_drop": NODES_TO_DROP,
        "n_processors": args.n_processors,
        "python_version": sys.version,
        "cell_type_mapper_version": package_version("cell-type-mapper"),
        "anndata_version": package_version("anndata"),
        "numpy_version": package_version("numpy"),
        "run_configurations": RUN_CONFIGS,
    }

    provenance_path = output_dir / "mapmycells_run_configuration.json"
    provenance_path.write_text(
        json.dumps(provenance, indent=2),
        encoding="utf-8",
    )

    print(f"Query: {query_path}")
    print(f"Output directory: {output_dir}")
    print(f"MapMyCells version: {provenance['cell_type_mapper_version']}")
    print(f"Number of configurations: {len(RUN_CONFIGS)}")

    for run_number, config in enumerate(RUN_CONFIGS, start=1):
        name = str(config["name"])

        command = build_command(
            query_path=query_path,
            marker_path=marker_path,
            stats_path=stats_path,
            output_dir=output_dir,
            run_config=config,
            n_processors=args.n_processors,
            overwrite_obsm=args.overwrite_obsm,
        )

        command_path = output_dir / f"cosmx_mmc_{name}_command.txt"
        command_path.write_text(
            subprocess.list2cmdline(command),
            encoding="utf-8",
        )

        print(
            f"\n[{run_number}/{len(RUN_CONFIGS)}] "
            f"Running MapMyCells configuration: {name}"
        )

        try:
            subprocess.run(
                command,
                env=run_environment,
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"MapMyCells configuration '{name}' failed "
                f"with exit code {exc.returncode}."
            ) from exc

    print("\nAll MapMyCells configurations completed successfully.")


if __name__ == "__main__":
    main()