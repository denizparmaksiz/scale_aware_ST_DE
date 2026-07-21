# Scale-aware spatial transcriptomics analyses

This repository contains the numbered analysis notebooks and supporting code for the manuscript. The notebooks remain the primary, readable analysis record; `src/scale_aware_st` and `kimlabspatial` contain only shared configuration, data-loading, and analysis helpers used by those notebooks.

## Start here

1. Clone this repository and open it at the repository root.
2. Download the accompanying data from [Zenodo (DOI: 10.5281/zenodo.21420404)](https://doi.org/10.5281/zenodo.21420404).
3. Place or link the downloaded files using the directory layout documented in the deposit README and `PUBLICATION_UPLOAD_MANIFEST.md`.
4. Install the repository helpers with `pip install -e .` from the repository root.
5. Use the environment export matching the analysis in `environment_exports/`.
6. Run notebooks from `documentation_code/` in numeric order as needed.

Supplementary workbooks and large manuscript data are distributed through Zenodo rather than committed to GitHub. By default, notebooks look for them under `Supplementary tables/` and `data/`, respectively. These locations can be overridden without editing notebook code:

```text
SCALE_AWARE_ST_REPO_ROOT
SCALE_AWARE_ST_DATA_DIR
SCALE_AWARE_ST_ADDITIONAL_DATA_DIR
SCALE_AWARE_ST_EXTERNAL_DATA_DIR
SCALE_AWARE_ST_RESULTS_DIR
SCALE_AWARE_ST_RESOURCES_DIR
SCALE_AWARE_ST_AGING_GENE_DIR
```

## Published results or full ALDEx recomputation

Notebooks that consume ALDEx results expose this flag near the top:

```python
ALDEX_RESULT_MODE = "published"  # choose "published" or "recompute"
```

- `published` loads the pooled, frozen ALDEx results from Additional Files 4–6 and skips the corresponding model-fitting/postprocessing step.
- `recompute` reads the individual per-stratum ALDEx result workbooks produced by the documented R workflow and reproduces the downstream pooled/postprocessed outputs.

The pooled supplementary worksheets and individual recomputation files have different physical layouts. The shared loader splits each pooled worksheet by its exact stratum label so downstream notebook code receives the same dictionary-of-tables structure in either mode.

The historical BLMM sensitivity run is the one exception: complete per-stratum BLMM files are not redistributed. The published engine comparison is in Additional File 4 (`ALDEx_lme_vs_blmm_engine`), while full BLMM recomputation requires the pinned development implementation documented in `documentation_code/aldex_main_crossstudy_repo/README.md`.

## Notebook index

1. `01_fig1_panel_composition.ipynb` — panel composition, spatial examples, and cross-study aging-gene overlap
2. `02_simulation_benchmarking_plots.ipynb` — simulation benchmark summaries and plots from Additional File 3
3. `03_main_merfish_aldex_prepost.ipynb` — primary MERFISH ALDEx input preparation and result postprocessing
4. `04_main_MERFISH_bio_figs.ipynb` — primary biological figures and representative spatial sections
5. `05_primaryMERFISH_functional_annotations.ipynb` — functional annotations and module summaries
6. `06_Sun_main_binary_DE.ipynb` — Sun et al. binary-age validation
7. `07_sun_continuous_raw_validation.ipynb` — Sun et al. continuous-age validation
8. `08_sun_supplementary_analyses.ipynb` — Sun et al. supplementary comparisons
9. `09_cosmx_figs_DE.ipynb` — CosMx differential-expression figures
10. `10_cosmx_merfish_module_analysis.ipynb` — cross-platform CosMx/MERFISH module analysis

External source datasets are not redistributed. Their download and placement instructions are kept beside the relevant workflows, including `documentation_code/FIGURE1D_DATA_DOWNLOADS.md` and `documentation_code/cosmx_python_scripts/cosmx_download_instructions.md`.

## Environments

The three Python environments used for the manuscript are documented as both Conda YAML exports and `pip freeze` snapshots in `environment_exports/`:

- `squidpy` — general spatial analyses and most notebooks
- `mapmycells` — MapMyCells label transfer
- `scvi_abc` — scANVI/ABC Atlas workflow

The HPC R and simulation environments are documented under `documentation_code/aldex_main_crossstudy_repo/environments/` and `documentation_code/simulation_benchmarking/`.

## Data publication

`PUBLICATION_UPLOAD_MANIFEST.md` is the authoritative file-by-file decision record for GitHub, Zenodo, and externally downloaded data. The [Zenodo deposit](https://doi.org/10.5281/zenodo.21420404) includes its own `README.md`, which documents the archive layout, and `SHA256SUMS.txt`, which records checksums for the deposited source files.

The deposit provides five logical groups (`figure1`, `primary_merfish`, `simulation`, `cosmx`, and `supplementary_files`) as separate ZIP archives, together with the Zenodo README and checksum file. Separate archives make selective downloading practical while preserving the directory structure expected by the notebooks.

## Validation

The lightweight tests in `tests/` cover pooled ALDEx splitting and aging-gene cutoff behavior. Repository validation also checks notebook Python syntax, retired filename references, the staged Zenodo layout, and staged-file checksums. Full end-to-end notebook execution requires the matching environment and the large deposited or external inputs.

The completed runtime-test record, including branches that were intentionally not rerun, is documented in `TESTING.md`.
