"""Published aging-gene set comparison used for manuscript Figure 1D."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class AgingGeneSourceFiles:
    """Published supplementary workbooks used by the Figure 1D analysis."""

    tabula_muris_senis: Path
    ximerakis: Path
    jin_allen: Path


CELL_TYPE_MAPPINGS: Mapping[str, Mapping[str, str]] = {
    "Oligodendrocytes": {
        "tabula_muris_senis": "Brain_Non-Myeloid.oligodendrocyte",
        "ximerakis": "OLG",
        "jin_allen": "MOL NN_4",
    },
    "Endothelial cells": {
        "tabula_muris_senis": "Brain_Non-Myeloid.endothelial cell",
        "ximerakis": "EC",
        "jin_allen": "Endo NN_1",
    },
    "Microglia": {
        "tabula_muris_senis": "Brain_Myeloid.microglial cell",
        "ximerakis": "MG",
        "jin_allen": "Microglia NN_1",
    },
}


def empirical_effect_cutoff(values: pd.Series) -> float:
    """Reproduce the original symmetric empirical-tail cutoff.

    The cutoff is the smaller absolute value of the 10th and 90th percentiles.
    Genes must have an absolute effect strictly greater than this value.
    """

    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(numeric).any():
        raise ValueError("Cannot derive an effect cutoff from all-missing values.")
    lower = abs(float(np.nanquantile(numeric, 0.1)))
    upper = abs(float(np.nanquantile(numeric, 0.9)))
    return min(lower, upper)


def _clean_gene_set(values) -> set[str]:
    return {
        str(value).strip()
        for value in values
        if pd.notna(value) and str(value).strip()
    }


def load_figure1d_aging_gene_sets(
    sources: AgingGeneSourceFiles,
    *,
    adjusted_p_threshold: float = 0.01,
) -> tuple[dict[str, dict[str, set[str]]], pd.DataFrame]:
    """Load the three published studies using the original selection rules.

    Returns
    -------
    gene_sets
        ``{matched_cell_type: {study_name: genes}}`` for plotting.
    selection_summary
        Auditable counts and empirical cutoffs for each study/cell type.
    """

    for path in (
        sources.tabula_muris_senis,
        sources.ximerakis,
        sources.jin_allen,
    ):
        if not Path(path).is_file():
            raise FileNotFoundError(f"Figure 1D source workbook not found: {path}")

    tms_padj = pd.read_excel(
        sources.tabula_muris_senis,
        sheet_name="DGE_result.tissue_cell.bh_p",
        index_col=0,
        engine="openpyxl",
    )
    tms_effect = pd.read_excel(
        sources.tabula_muris_senis,
        sheet_name="DGE_result.tissue_cell.age_coef",
        index_col=0,
        engine="openpyxl",
    )
    jin = pd.read_excel(
        sources.jin_allen,
        sheet_name="supertype table",
        engine="openpyxl",
    )
    required_jin = {"supertype_label", "gene"}
    if missing := required_jin.difference(jin.columns):
        raise ValueError(f"Jin et al. workbook is missing columns: {sorted(missing)}")

    gene_sets: dict[str, dict[str, set[str]]] = {}
    summary_rows: list[dict[str, object]] = []

    for matched_cell_type, labels in CELL_TYPE_MAPPINGS.items():
        tms_label = labels["tabula_muris_senis"]
        if tms_label not in tms_padj or tms_label not in tms_effect:
            raise KeyError(f"Missing Tabula Muris Senis column: {tms_label}")
        tms_cutoff = empirical_effect_cutoff(tms_effect[tms_label])
        tms_mask = (
            pd.to_numeric(tms_padj[tms_label], errors="coerce")
            <= adjusted_p_threshold
        ) & (
            pd.to_numeric(tms_effect[tms_label], errors="coerce").abs()
            > tms_cutoff
        )
        tms_genes = _clean_gene_set(tms_padj.index[tms_mask])

        ximerakis_label = labels["ximerakis"]
        ximerakis = pd.read_excel(
            sources.ximerakis,
            sheet_name=ximerakis_label,
            engine="openpyxl",
        )
        required_ximerakis = {"Gene", "padj", "logFC_Young_to_Old"}
        if missing := required_ximerakis.difference(ximerakis.columns):
            raise ValueError(
                f"Ximerakis sheet {ximerakis_label!r} is missing columns: "
                f"{sorted(missing)}"
            )
        ximerakis_effect = pd.to_numeric(
            ximerakis["logFC_Young_to_Old"], errors="coerce"
        )
        ximerakis_cutoff = empirical_effect_cutoff(ximerakis_effect)
        ximerakis_mask = (
            pd.to_numeric(ximerakis["padj"], errors="coerce")
            <= adjusted_p_threshold
        ) & (ximerakis_effect.abs() > ximerakis_cutoff)
        ximerakis_genes = _clean_gene_set(
            ximerakis.loc[ximerakis_mask, "Gene"]
        )

        jin_label = labels["jin_allen"]
        jin_genes = _clean_gene_set(
            jin.loc[jin["supertype_label"].astype(str).eq(jin_label), "gene"]
        )

        studies = {
            "Tabula Muris Senis": tms_genes,
            "Ximerakis et al.": ximerakis_genes,
            "Jin et al. (Allen)": jin_genes,
        }
        gene_sets[matched_cell_type] = studies

        for study, genes, cutoff, source_label in (
            ("Tabula Muris Senis", tms_genes, tms_cutoff, tms_label),
            ("Ximerakis et al.", ximerakis_genes, ximerakis_cutoff, ximerakis_label),
            ("Jin et al. (Allen)", jin_genes, np.nan, jin_label),
        ):
            summary_rows.append(
                {
                    "matched_cell_type": matched_cell_type,
                    "study": study,
                    "source_cell_type": source_label,
                    "adjusted_p_threshold": (
                        adjusted_p_threshold if np.isfinite(cutoff) else np.nan
                    ),
                    "empirical_effect_cutoff": cutoff,
                    "n_selected_genes": len(genes),
                }
            )

    return gene_sets, pd.DataFrame(summary_rows)


def plot_figure1d_venn(
    gene_sets: Mapping[str, Mapping[str, set[str]]],
    *,
    output_dir: str | Path | None = None,
    dpi: int = 600,
    remove_legend: bool = True,
):
    """Plot one three-study Venn diagram per matched broad cell type."""

    import matplotlib.pyplot as plt
    from venn import venn

    output_path = Path(output_dir) if output_dir is not None else None
    if output_path is not None:
        output_path.mkdir(parents=True, exist_ok=True)

    figures = {}
    for cell_type, studies in gene_sets.items():
        fig, ax = plt.subplots(figsize=(5, 5))
        venn(dict(studies), ax=ax)
        legend = ax.get_legend()
        if remove_legend and legend is not None:
            legend.remove()
        ax.set_title(cell_type)
        fig.tight_layout()
        if output_path is not None:
            filename = cell_type.lower().replace(" ", "_")
            fig.savefig(
                output_path / f"fig1d_aging_gene_overlap_{filename}.png",
                dpi=dpi,
                bbox_inches="tight",
            )
        figures[cell_type] = fig
    return figures
