# Data downloads and expected layout

This is the central data-setup guide for the numbered notebooks. It separates
three kinds of input:

1. **Zenodo data generated for this study** go under `data/`.
2. **Third-party published data** go under `external_data/` and are not
   redistributed by this project.
3. **Small version-controlled resources** are already present under
   `resources/`.

Paths below are relative to the repository root. Environment-variable
overrides are described at the end of this guide.

## 1. Study data from Zenodo

During peer review, download the five archives using the private Zenodo preview
link supplied with the manuscript. DOI `10.5281/zenodo.21420404` is reserved but
will not resolve until the record is published. Extract each archive directly
under `data/`. The result should be:

```text
data/
|-- figure1/
|   `-- panel_objects/
|-- primary_merfish/
|   |-- analysis_objects/
|   |-- ccf_metadata/
|   `-- spatial_sections/
|-- simulation/
|   `-- source_inputs/
|-- cosmx/
|   `-- cosmx_data_ct_final.h5ad
`-- supplementary_files/
    |-- Additional_File_1.xlsx
    |-- Additional_File_2.pdf
    |-- Additional_File_3.xlsx
    |-- Additional_File_4.xlsx
    |-- Additional_File_5.xlsx
    `-- Additional_File_6.xlsx
```

The archive-level README on Zenodo lists every deposited filename. With this
layout, no data-related environment variables are needed.

## 2. Third-party published inputs

These data remain subject to the original providers' terms and must be
downloaded from the cited sources.

### Sun et al. MERFISH ageing atlas

Used by notebooks 06, 07, and 08:

```text
external_data/sun_et_al/
|-- aging_coronal.h5ad
|-- 2023-12-22736D-TableS7_GeneClassificationTrajectory.xlsx
|-- 2023-12-22736D-TableS8_GOBPAging.xlsx
`-- 2023-12-22736D-TableS11_GOBPClockGenes.xlsx
```

- Download `aging_coronal.h5ad` from the authors'
  [processed MERFISH Zenodo record](https://zenodo.org/records/13883177).
- Download the Supplementary Tables archive from
  [Spatial transcriptomic clocks reveal cell proximity effects in brain
  ageing](https://doi.org/10.1038/s41586-024-08334-8), then extract Tables S7,
  S8, and S11 without changing their filenames.

The AnnData file is required even in `published` mode because these notebooks
also calculate expression summaries or rank-based comparisons directly from
the published atlas. The supplementary tables are external comparison inputs,
not substitutes for this project's Additional File 5.

### Figure 1D published ageing-gene workbooks

Used by notebook 01:

```text
external_data/published_aging_gene_sets/
|-- TMS6_DEGs.xlsx
|-- PMID_31551601_Differential gene expression data between young and old cell types.xlsx
`-- 41586_2024_8350_MOESM4_ESM.xlsx
```

The publication links, required sheets, and selection rules are documented in
`documentation_code/FIGURE1D_DATA_DOWNLOADS.md`.

### Public Xenium datasets used for Figure 1 panel comparisons

These inputs are only required when
`RUN_PUBLIC_PLATFORM_COMPARISONS = True` in notebook 01:

```text
external_data/xenium/
|-- Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs/
|-- Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs/
|-- Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs/
|-- Xenium_V1_FFPE_wildtype_2_5_months_outs/
`-- Xenium_Prime_Mouse_Brain_Coronal_FF_outs/
```

- The three v1 replicate outputs are available from the 10x Genomics
  [Fresh Frozen Mouse Brain Replicates dataset](https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard).
- The 5K output is available from the 10x Genomics
  [Xenium Prime Fresh Frozen Mouse Brain dataset](https://www.10xgenomics.com/datasets/xenium-prime-fresh-frozen-mouse-brain).
- Obtain the 2.5-month FFPE wild-type output from the corresponding 10x public
  Xenium dataset page and retain the output-directory name shown above.

Download complete output bundles rather than the small Explorer-only subsets;
the notebook reads the output directories with Squidpy.

### SenNet CosMx inputs

Notebook 01's optional public-platform comparison expects:

```text
external_data/sennet_cosmx/SNT638.MWRV.378/
```

The complete raw-CosMx reconstruction workflow uses four young-old pairs. The
required SenNet collection, files, and pair-specific directory organization are
documented in
`documentation_code/cosmx_python_scripts/cosmx_download_instructions.md`.
Those raw inputs are not required by notebooks 09 and 10 when using the
deposited `data/cosmx/cosmx_data_ct_final.h5ad` object.

## 3. Resources already supplied by GitHub

Do not download or substitute a different atlas workbook for notebook 04. The
exact uint16 ontology workbook used by the analysis is supplied at:

```text
resources/allen_ccf/CCFv3OntologyStructure_u16.xlsx
```

Marker lists, frozen functional annotations, and module dictionaries are also
already present under `resources/`.

## Notebook input matrix

| Notebook | Zenodo study data | Third-party downloads |
|---|---|---|
| 01 | Figure 1 panel objects and representative MERFISH sections | Three ageing-gene workbooks; Xenium and SenNet CosMx only for optional public-platform comparisons |
| 02 | Additional File 3 | None |
| 03 | Primary MERFISH AnnData and/or Additional File 4 | None |
| 04 | Primary MERFISH objects, sections, CCF metadata, and Additional File 4 | None; the exact ontology is in GitHub resources |
| 05 | Additional File 4 | None in the default frozen-annotation route |
| 06 | Primary MERFISH AnnData and Additional Files 4–5 | Sun `aging_coronal.h5ad` |
| 07 | Additional File 5 | Sun `aging_coronal.h5ad` and Table S7 |
| 08 | Additional File 5 | Sun `aging_coronal.h5ad` and Tables S7, S8, and S11 |
| 09 | CosMx and primary MERFISH AnnData plus Additional Files 4 and 6 | None when using the deposited CosMx object |
| 10 | CosMx AnnData | None |

## Storing data somewhere else

The defaults above are recommended, but large files do not have to live in the
Git checkout. Set paths before starting Jupyter:

```text
SCALE_AWARE_ST_DATA_DIR=/absolute/path/to/extracted_zenodo_data
SCALE_AWARE_ST_EXTERNAL_DATA_DIR=/absolute/path/to/third_party_data
```

`SCALE_AWARE_ST_ADDITIONAL_DATA_DIR` defaults to
`SCALE_AWARE_ST_DATA_DIR/supplementary_files`; set it separately only when the
Additional Files are stored elsewhere. `SCALE_AWARE_ST_RESOURCES_DIR` and
`SCALE_AWARE_ST_AGING_GENE_DIR` can likewise override the small-resource and
published-ageing-workbook locations.
