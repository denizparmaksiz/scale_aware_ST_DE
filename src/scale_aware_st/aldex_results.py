"""Load ALDEx results from either manuscript runs or pooled Additional Files.

The central contract is a dictionary mapping an exact stratum label to the raw
ALDEx result table for that stratum. Both input modes preserve column order,
gene order, numeric values, and stratum spelling.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Mapping

import pandas as pd


@dataclass(frozen=True)
class PooledALDExSpec:
    additional_file: str
    sheet_name: str
    stratum_column: str


POOLED_ALDEX_SPECS: Mapping[str, PooledALDExSpec] = {
    "primary_celltype_region": PooledALDExSpec(
        "Additional_File_4.xlsx",
        "Celltype_region_ALDEx_results",
        "Celltype x region",
    ),
    "primary_celltype": PooledALDExSpec(
        "Additional_File_4.xlsx", "Celltype_ALDEx_results", "Cell type"
    ),
    "primary_anterior_celltype_region": PooledALDExSpec(
        "Additional_File_4.xlsx",
        "Anterior_plane_ALDEx_results",
        "Celltype x region",
    ),
    "sun_binary_celltype_region": PooledALDExSpec(
        "Additional_File_5.xlsx",
        "Binary_age_ct_region_ALDEx",
        "Celltype x region",
    ),
    "sun_continuous_celltype": PooledALDExSpec(
        "Additional_File_5.xlsx", "Cont_age_ct_ALDEx", "Cell type"
    ),
    "sun_continuous_celltype_subregion": PooledALDExSpec(
        "Additional_File_5.xlsx",
        "Cont_age_ct_subregion_ALDEx",
        "Celltype x subregion",
    ),
    "cosmx_celltype": PooledALDExSpec(
        "Additional_File_6.xlsx", "Cosmx_ALDEx_results", "Celltype"
    ),
}


def _validate_result_table(table: pd.DataFrame, context: str) -> None:
    if "gene" not in table.columns:
        raise ValueError(f"{context} is missing the required 'gene' column.")
    if table["gene"].isna().any():
        raise ValueError(f"{context} contains missing gene identifiers.")
    duplicated = table["gene"].astype(str).duplicated()
    if duplicated.any():
        examples = table.loc[duplicated, "gene"].astype(str).head(5).tolist()
        raise ValueError(f"{context} contains duplicate genes: {examples}")


def load_pooled_aldex_results(
    workbook: str | Path,
    *,
    sheet_name: str,
    stratum_column: str,
) -> dict[str, pd.DataFrame]:
    """Split one pooled Additional File worksheet into per-stratum tables."""

    workbook = Path(workbook)
    if not workbook.is_file():
        raise FileNotFoundError(f"Pooled ALDEx workbook not found: {workbook}")
    pooled = pd.read_excel(workbook, sheet_name=sheet_name, engine="openpyxl")
    if stratum_column not in pooled.columns:
        raise ValueError(
            f"Sheet {sheet_name!r} in {workbook.name} is missing stratum column "
            f"{stratum_column!r}. Available columns: {pooled.columns.tolist()}"
        )
    if pooled[stratum_column].isna().any():
        raise ValueError(
            f"Sheet {sheet_name!r} contains rows with missing stratum labels."
        )

    results: dict[str, pd.DataFrame] = {}
    # sort=False and a boolean mask preserve first-seen stratum and row order.
    for label in pooled[stratum_column].astype(str).drop_duplicates():
        table = pooled.loc[
            pooled[stratum_column].astype(str).eq(label)
        ].drop(columns=[stratum_column])
        table = table.reset_index(drop=True)
        _validate_result_table(table, f"{workbook.name}:{sheet_name}:{label}")
        results[label] = table
    return results


def load_pooled_strata(
    workbook: str | Path,
    *,
    sheet_name: str,
    stratum_column: str,
    gene_column: str,
) -> dict[str, pd.DataFrame]:
    """Split any pooled manuscript table while preserving its published schema."""

    workbook = Path(workbook)
    if not workbook.is_file():
        raise FileNotFoundError(f"Pooled workbook not found: {workbook}")
    pooled = pd.read_excel(workbook, sheet_name=sheet_name, engine="openpyxl")
    missing = {stratum_column, gene_column}.difference(pooled.columns)
    if missing:
        raise ValueError(
            f"{workbook.name}:{sheet_name} is missing columns: {sorted(missing)}"
        )
    results = {}
    for label in pooled[stratum_column].astype(str).drop_duplicates():
        table = pooled.loc[pooled[stratum_column].astype(str).eq(label)].copy()
        table = table.drop(columns=[stratum_column]).reset_index(drop=True)
        if table[gene_column].isna().any() or table[gene_column].astype(str).duplicated().any():
            raise ValueError(
                f"{workbook.name}:{sheet_name}:{label} has missing or duplicate genes."
            )
        results[label] = table
    return results


def index_result_tables(
    results: Mapping[str, pd.DataFrame],
    *,
    gene_column: str = "gene",
) -> dict[str, pd.DataFrame]:
    """Return copies indexed by gene for downstream comparison code."""

    indexed = {}
    for label, table in results.items():
        if gene_column not in table.columns:
            raise ValueError(f"Result {label!r} is missing {gene_column!r}.")
        current = table.copy().set_index(gene_column)
        current.index.name = None
        indexed[label] = current
    return indexed


def load_pooled_aldex_spec(
    additional_data_dir: str | Path,
    analysis: str,
) -> dict[str, pd.DataFrame]:
    """Load a named manuscript analysis from its pooled Additional File."""

    try:
        spec = POOLED_ALDEX_SPECS[analysis]
    except KeyError as exc:
        raise KeyError(
            f"Unknown pooled ALDEx analysis {analysis!r}. Choose from "
            f"{sorted(POOLED_ALDEX_SPECS)}."
        ) from exc
    return load_pooled_aldex_results(
        Path(additional_data_dir) / spec.additional_file,
        sheet_name=spec.sheet_name,
        stratum_column=spec.stratum_column,
    )


def load_individual_aldex_results(
    paths: Iterable[str | Path],
    *,
    label_from_path: Callable[[Path], str],
) -> dict[str, pd.DataFrame]:
    """Load individual ALDEx Excel files into the same per-stratum contract."""

    results: dict[str, pd.DataFrame] = {}
    for value in paths:
        path = Path(value)
        if not path.is_file():
            raise FileNotFoundError(f"Individual ALDEx result not found: {path}")
        label = str(label_from_path(path))
        if not label:
            raise ValueError(f"Empty stratum label parsed from {path}")
        if label in results:
            raise ValueError(f"Duplicate stratum label {label!r} from {path}")
        table = pd.read_excel(path, engine="openpyxl")
        _validate_result_table(table, str(path))
        results[label] = table.reset_index(drop=True)
    if not results:
        raise ValueError("No individual ALDEx result files were supplied.")
    return results
