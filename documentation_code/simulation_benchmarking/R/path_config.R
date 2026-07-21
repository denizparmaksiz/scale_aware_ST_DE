# Shared portable paths for the simulation benchmarking scripts.
# Every location can be overridden without editing code.

env_path <- function(name, default) {
  normalizePath(Sys.getenv(name, unset = default), mustWork = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)[1]
if (is.na(script_file_arg)) stop("Run this analysis with Rscript so its location can be resolved.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg), mustWork = FALSE))
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)

data_root <- env_path("SCALE_AWARE_ST_DATA_DIR", file.path(repo_root, "zenodo_data"))
simulation_source_dir <- env_path(
  "SIMULATION_SOURCE_DIR",
  file.path(data_root, "simulation", "source_inputs")
)
simulation_work_dir <- env_path(
  "SIMULATION_WORK_DIR",
  file.path(repo_root, "results", "simulation")
)
simulation_generated_dir <- env_path(
  "SIMULATION_RESULTS_DIR",
  file.path(simulation_work_dir, "generated")
)
simulation_aldex_dir <- env_path(
  "SIMULATION_ALDEX_DIR",
  file.path(simulation_work_dir, "aldex")
)
simulation_mast_dir <- env_path(
  "SIMULATION_MAST_DIR",
  file.path(simulation_work_dir, "mast")
)
simulation_nb_dir <- env_path(
  "SIMULATION_NB_DIR",
  file.path(simulation_work_dir, "negative_binomial")
)
simulation_scanpy_dir <- env_path(
  "SIMULATION_SCANPY_DIR",
  file.path(simulation_work_dir, "scanpy")
)
simulation_benchmark_dir <- env_path(
  "SIMULATION_BENCHMARK_DIR",
  file.path(simulation_work_dir, "benchmarking")
)
