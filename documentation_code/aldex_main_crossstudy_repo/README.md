# ALDEx3 analyses: primary dataset and cross-study comparisons

This directory contains the R and SLURM workflows used for the manuscript's main-dataset and cross-study ALDEx3 analyses. Simulation benchmarking and Python postprocessing are intentionally outside this directory.

## Repository layout

- `R/aldex_helpers.R`: shared non-scientific utilities.
- `R/01_...R` through `R/08_...R`: run-specific analyses.
- `hpc/submit_aldex.sh`: shared SLURM launcher.
- `config/aldex_runs.csv`: manifest linking each manuscript analysis to its script, engine, model, inputs, and outputs.
- `docs/repo_planning_ALDEx.md`: remaining documentation and environment tasks.

## Input format

Each analysis expects one raw-count CSV and one metadata Excel file for each analysis stratum.

### Count matrix

- Rows are genes.
- Columns are cells/observations.
- Values are raw nonnegative transcript counts.
- The first CSV column contains gene names and is read as row names.
- Column names must match observation identifiers in the metadata.

### Metadata

- The first column contains observation identifiers and is read as row names.
- Required columns depend on the analysis; see `config/aldex_runs.csv`.
- The scripts align metadata rows to count-matrix columns before fitting.

The date suffixes in the original filenames are retained as configurable script arguments (`--input-date` and `--output-date`). Local desktop paths in the original run-tracking spreadsheet and historical GPFS paths in the submitted HPC scripts refer to equivalent logical input/output collections; neither is hard-coded in the cleaned scripts.

## Running on HPC

Example:

```bash
sbatch hpc/submit_aldex.sh \
  R/01_primary_merfish_ct_region_blmm.R \
  /gpfs/path/to/celltype_region \
  /gpfs/path/to/results \
  blmm
```

The launcher currently reproduces the historical cluster setup:

- `R/4.5.0`
- main ALDEx3 library: `$HOME/R/envs-4.5.0/ALDEx3`
- `blmm` library: `$HOME/R/envs-4.5.0/ALDEx3_blmm`

These are site-specific paths and should be changed for another cluster.

## ALDEx3 and the `blmm` engine

The manuscript analyses used two historical ALDEx3 installations:

- **Standard ALDEx3:** version `0.4.0`, commit `8c05ad40c41279dffa05dc808167ffcd53207740`.
- **BLMM ALDEx3:** version `1.1.0`, branch `fast-lmms`, commit `68004e7e8f89b825c8d8b11166ffb76696605b03`.

Several manuscript analyses use `method = "blmm"`. At present, that engine must be pulled separately from the `fast-lmms` development branch; installing the standard/main ALDEx3 release alone is not sufficient for those scripts. The `blmm` implementation is expected to be integrated into the main branch, but the manuscript workflow should remain pinned to the historical commit above until equivalence with a future integrated release has been verified.

Do not replace the historical `blmm` dependency with another mixed-model engine and assume identical results.

## Monte Carlo sampling

Final manuscript analyses use the values recorded in `config/aldex_runs.csv`, generally 2,000 Monte Carlo samples. The Sun continuous-age cell-type analysis was confirmed to use 1,000 samples.

The original scripts did not set a fixed random seed. The cleaned scripts preserve that behavior. A seed should not be added to the manuscript reproduction workflow without documenting that change.

## Remaining checks before publication

1. Export `sessionInfo()` and/or an `renv.lock` from each HPC R library.
2. Verify each cleaned script against at least one representative historical input/output pair.
3. Add a small desktop tutorial separately; it should use reduced MC sampling and clearly state that it will not reproduce the final manuscript results.
