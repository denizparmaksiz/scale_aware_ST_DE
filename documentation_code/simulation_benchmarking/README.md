# Simulation benchmarking scripts

These scripts no longer contain computer- or user-specific paths. By default,
the deposited source inputs are read from
`zenodo_data/simulation/source_inputs/` at the repository root and generated
files are written under `results/simulation/`.

Set `SCALE_AWARE_ST_DATA_DIR` to the extracted Zenodo data directory used by the
numbered notebooks. Alternatively, set `SIMULATION_SOURCE_DIR` directly to the
directory containing:

- `adata_nn_raw_counts_transposed.csv`
- `adata_nn_raw_counts_transposed_MD.xlsx`
- `updated_sim_design.rds`

Set `SIMULATION_WORK_DIR` when the generated matrices and model results should
be written elsewhere, which is recommended on an HPC scratch filesystem. The
following optional variables override individual locations:

- `SIMULATION_RESULTS_DIR`
- `SIMULATION_ALDEX_DIR`
- `SIMULATION_MAST_DIR`
- `SIMULATION_NB_DIR`
- `SIMULATION_SCANPY_DIR`
- `SIMULATION_BENCHMARK_DIR`

The Slurm submission scripts locate the adjacent `R/` and `python/` directories
automatically. Therefore, users do not need to edit a `cd` command before
submitting them. Cluster-specific module names, partitions, memory requests,
and environment locations may still need adjustment for the target HPC system.
