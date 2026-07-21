# Figure 1D: published aging-gene inputs

Figure 1D compares published aging-associated gene sets from three mouse brain
single-cell transcriptomic studies. These source workbooks are supplementary
data from other publications and are **not redistributed in this repository**.

Download the workbooks from the original publications, place them in
`data/published_aging_gene_sets/`, and retain the filenames below. Alternatively,
set the `SCALE_AWARE_ST_AGING_GENE_DIR` environment variable to the folder that
contains them.

| Study | Expected filename | Data used by the notebook |
|---|---|---|
| Tabula Muris Senis (2020) | `TMS6_DEGs.xlsx` | Sheets `DGE_result.tissue_cell.bh_p` and `DGE_result.tissue_cell.age_coef` |
| Ximerakis et al. (2019) | `PMID_31551601_Differential gene expression data between young and old cell types.xlsx` | Sheets `OLG`, `EC`, and `MG` |
| Jin et al. (2024) / Allen Brain Institute | `41586_2024_8350_MOESM4_ESM.xlsx` | Sheet `supertype table` |

Publication pages and data:

- Tabula Muris Senis: [A single-cell transcriptomic atlas characterizes ageing tissues in the mouse](https://doi.org/10.1038/s41586-020-2496-1) and its [processed data record](https://doi.org/10.6084/m9.figshare.8273102.v2).
- Ximerakis et al.: [Single-cell transcriptomic profiling of the aging mouse brain](https://doi.org/10.1038/s41593-019-0491-3). Download the supplementary workbook titled “Differential gene expression data between young and old cell types.”
- Jin et al.: [Brain-wide cell-type-specific transcriptomic signatures of healthy ageing in mice](https://doi.org/10.1038/s41586-024-08350-8). Download the supplementary workbook containing the `supertype table` sheet.

The notebook validates the expected sheets and columns and stops with a clear
message if a file is absent or incompatible. For Tabula Muris Senis and
Ximerakis et al., genes are selected using adjusted p-value ≤ 0.01 and a
study/cell-type-specific empirical effect-size cutoff. The published Jin et al.
gene lists are used directly.
