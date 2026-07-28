# Scale-aware spatial transcriptomics analyses

This repository contains the numbered analysis notebooks and supporting code for the manuscript. The notebooks remain the primary, readable analysis record; `src/scale_aware_st` and `kimlabspatial` contain only shared configuration, data-loading, and analysis helpers used by those notebooks.

## Start here

1. Clone this repository and open it at the repository root.
2. Download the accompanying data from [Zenodo (DOI: 10.5281/zenodo.21420404)](https://doi.org/10.5281/zenodo.21420404).
3. Extract the five Zenodo archives under `data/`, preserving the directories inside each archive.
4. Create the tested core Python environment using the instructions below.
5. Run notebooks from `documentation_code/` using the registered `scale-aware-st-de` kernel.

Supplementary workbooks and large manuscript data are distributed through Zenodo rather than committed to GitHub. By default, notebooks use a single extracted Zenodo tree under `data/`; supplementary workbooks are therefore expected at `data/supplementary_files/`. Third-party inputs that cannot be redistributed belong under `external_data/`. See `DATA_DOWNLOADS.md` for the complete layout, download links, and notebook-by-notebook requirements.

All locations can be overridden without editing notebook code:

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

External source datasets are not redistributed. `DATA_DOWNLOADS.md` is the central download and placement guide. The Figure 1D and raw CosMx workflows also retain their more detailed guides under `documentation_code/`.

## Installation

The numbered notebooks were tested with **Python 3.10.14**. The core dependency
set deliberately pins the manuscript-era NumPy, Zarr, SpatialData, Scanpy, and
Squidpy versions; installing current unpinned releases is not supported. The
Squidpy build is pinned to commit
`df8e042264a99e489573273e449d585e3ff6143c`.

### Recommended Conda installation

Run these commands from the repository root:

```bash
conda env create -f environment_exports/environment-core.yml
conda activate scale-aware-st-de
```

The environment file creates Python 3.10.14 and installs this repository in
editable mode with its tested core dependencies. Editable installation means
changes to `src/scale_aware_st/` and `kimlabspatial/` are used immediately;
it does not copy the source code into the environment.

If the Conda environment already exists, update it from the repository root:

```bash
conda env update -f environment_exports/environment-core.yml --prune
```

### Manual installation into a fresh Python 3.10 environment

```bash
conda create -n scale-aware-st-de python=3.10.14 pip git
conda activate scale-aware-st-de
python -m pip install -e .
```

`pip install -e .` now installs the full core environment required by the
numbered notebooks, rather than only the small repository-helper dependency
set. The same portable pins are listed in
`environment_exports/requirements-core.txt`.

### Register the Jupyter kernel

```bash
python -m ipykernel install --user --name scale-aware-st-de --display-name "Python (scale-aware-st-de)"
```

Use `Python (scale-aware-st-de)` when opening notebooks under
`documentation_code/`. A Jupyter frontend may be installed separately if one
is not already available, for example with `python -m pip install jupyterlab`.

### Verify the installation

```bash
python -m pip check
python -c "import scanpy, squidpy, spatialdata, zarr; import kimlabspatial.preprocessing, scale_aware_st; print('Environment OK')"
```

To run the repository tests as well:

```bash
python -m pip install -e ".[test]"
python -m pytest
```

### Optional Scrublet support

Scrublet is only needed for raw-data preprocessing when doublet detection is
explicitly enabled, such as `pool_adatas(..., scrub=True)`. It is not required
for the numbered notebooks using deposited processed data.

```bash
python -m pip install -e ".[scrublet]"
```

Alternatively, install
`environment_exports/requirements-scrublet.txt`. On native Windows, Annoy may
require Microsoft Visual C++ Build Tools because a compatible prebuilt wheel
may not be available. Linux or WSL is therefore recommended for this optional
raw-preprocessing route. Without Scrublet, the package imports normally and
raises an informative error only if doublet detection is requested.

### Other manuscript environments

The historical exports remain under `environment_exports/` as provenance:

- `squidpy` was the main analysis environment and underlies the portable core
  environment above.
- `mapmycells` is needed only to recompute the MapMyCells stage of CosMx
  preprocessing.
- `scvi_abc` is needed only to recompute the scANVI/ABC Atlas stage of CosMx
  preprocessing.

Neither `mapmycells` nor `scvi_abc` is required to run the numbered notebooks
from the deposited processed CosMx object. The HPC R and simulation environments
are documented under `documentation_code/aldex_main_crossstudy_repo/environments/`
and `documentation_code/simulation_benchmarking/`.

## Data publication

`PUBLICATION_UPLOAD_MANIFEST.md` is the authoritative file-by-file decision record for GitHub, Zenodo, and externally downloaded data. The [Zenodo deposit](https://doi.org/10.5281/zenodo.21420404) includes its own `README.md`, which documents the archive layout, and `SHA256SUMS.txt`, which records checksums for the deposited source files.

The deposit provides five logical groups (`figure1`, `primary_merfish`, `simulation`, `cosmx`, and `supplementary_files`) as separate ZIP archives, together with the Zenodo README and checksum file. Separate archives make selective downloading practical while preserving the directory structure expected by the notebooks.

## Validation

The lightweight tests in `tests/` cover pooled ALDEx splitting and aging-gene cutoff behavior. Repository validation also checks notebook Python syntax, retired filename references, the staged Zenodo layout, and staged-file checksums. Full end-to-end notebook execution requires the matching environment and the large deposited or external inputs.

The completed runtime-test record, including branches that were intentionally not rerun, is documented in `TESTING.md`.
