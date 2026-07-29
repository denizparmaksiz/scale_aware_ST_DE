# Runtime testing record

The numbered public notebooks were smoke-tested sequentially in the exported
`squidpy` Conda environment (Python 3.10.14) using the staged Zenodo inputs and
the public validation datasets listed below. Each Python code cell was executed
in notebook order with a non-interactive Matplotlib backend.

## Completed tests

| Notebook | Tested route | Result |
|---|---|---|
| `01_fig1_panel_composition.ipynb` | Deposited Figure 1 panel objects, Figure 1D public aging-gene workbooks, and deposited MERFISH spatial sections | Passed |
| `02_simulation_benchmarking_plots.ipynb` | Published Additional File 3 route | Passed |
| `03_main_merfish_aldex_prepost.ipynb` | Published pooled ALDEx route | Passed |
| `04_main_MERFISH_bio_figs.ipynb` | Published ALDEx route plus full processed-AnnData and spatial-section workflow | Passed |
| `05_primaryMERFISH_functional_annotations.ipynb` | Published annotations and module resources | Passed |
| `06_Sun_main_binary_DE.ipynb` | Published pooled results and pseudobulk route | Passed |
| `07_sun_continuous_raw_validation.ipynb` | Published results plus the 1,453,144-cell Sun et al. object | Passed |
| `08_sun_supplementary_analyses.ipynb` | Published ALDEx/GO/pseudobulk route and Sun supplementary tables | Passed |
| `09_cosmx_figs_DE.ipynb` | Published results plus the 63,356-cell processed CosMx object | Passed |
| `10_cosmx_merfish_module_analysis.ipynb` | Frozen gene annotations and processed CosMx object | Passed |

Notebook 10 produced 85 module-score rows. All 85 keys and values matched the
`Shared_module_scores` sheet in Additional File 6; maximum numerical differences
were at floating-point rounding scale (at most approximately `9e-16`).

## Explicit limitations

- Notebook 1's optional public Xenium and raw SenNet CosMx platform-comparison
  blocks were not executed because those downloads were not staged. The notebook
  defaults to `RUN_PUBLIC_PLATFORM_COMPARISONS = False`; download the documented
  inputs and set the flag to `True` to run those blocks.
- The `ALDEX_RESULT_MODE = "recompute"` branches were not used to rerun all R
  ALDEx models during this smoke test. The pooled published-result loading and
  downstream Python analyses were tested. Full model recomputation requires the
  documented R/HPC environments and individual per-stratum outputs.
- Complete per-stratum BLMM sensitivity outputs are not redistributed. The
  frozen engine comparison is provided in Additional File 4.

Warnings stating that `FigureCanvasAgg` is non-interactive are expected during
headless smoke testing and do not indicate analysis failures.

## Mode-contract validation

After the published/recompute documentation and preflight checks were added,
all code cells in notebooks 03–09 were parsed again as Python and all notebooks
were parsed as valid notebook JSON. Static regression tests verify that every
affected notebook contains a mode preflight and that workflows crossing the
external R/ALDEx boundary reject an empty workbook collection with a specific
next-step message.

The analytical bodies of the notebooks were not changed. Nevertheless, every
notebook affected by this feedback block (03-09) received a full published-mode
regression run after the changes, and all seven passed:

- Notebook 03 loaded 34 cell type x region, 6 cell-type, and 33 anterior pooled
  strata.
- Notebook 04 completed the frozen pooled ALDEx, processed AnnData, ontology,
  and spatial-section workflow.
- Notebook 05 loaded all 34 cell type x region result sheets and completed the
  frozen annotation workflow.
- Notebook 06 loaded 35 pooled Sun binary-age strata, retained 21 ALDEx
  comparisons, and loaded 7 primary and 7 Sun pseudobulk cell types.
- Notebook 07 loaded 18 pooled continuous-age cell-type strata and the
  1,453,144-cell Sun et al. object.
- Notebook 08 loaded 120 pooled continuous-age subregion strata and completed
  the deposited pseudobulk and trajectory comparisons.
- Notebook 09 loaded the 63,356-cell CosMx object and completed the frozen ALDEx
  and deposited pseudobulk comparisons for five matched cell types.
