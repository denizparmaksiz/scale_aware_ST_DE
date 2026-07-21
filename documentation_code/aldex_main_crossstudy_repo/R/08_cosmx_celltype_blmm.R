#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ALDEx3)
  library(lme4)
  library(readxl)
  library(stringr)
  library(writexl)
})
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "R/08_cosmx_celltype_blmm.R"
source(file.path(dirname(normalizePath(script_path)), "aldex_helpers.R"))

config <- parse_cli_args(list(
  input_dir = Sys.getenv("ALDEX_INPUT_DIR", unset = "cosmx_data_ct_aldex_mem"),
  output_dir = Sys.getenv("ALDEX_OUTPUT_DIR", unset = "results/cosmx_celltype_blmm"),
  nsample = 2000,
  n_cores = 10,
  detection_threshold = 0.12,
  gamma_value = 0.80,
  scale_sd = 0.25,
  input_date = "05212026",
  output_date = "05212026"
))
config$output_dir <- ensure_output_dir(config$output_dir)
logger <- make_logger("cosmx_celltype_blmm")

files <- list.files(config$input_dir, pattern = "^cell_type_.*\\.csv$", full.names = FALSE)
file_info <- str_match(files, "^cell_type_(.*?)_(.*?)_\\d+\\.csv")
colnames(file_info) <- c("full", "region", "celltype")
file_info <- as.data.frame(file_info, stringsAsFactors = FALSE)
file_info <- file_info[!is.na(file_info$region), , drop = FALSE]
file_info$base_id <- paste(file_info$region, file_info$celltype, sep = "_")

for (id in file_info$base_id) {
  count_path <- file.path(config$input_dir, sprintf("cell_type_%s_%s.csv", id, config$input_date))
  metadata_path <- file.path(config$input_dir, sprintf("cell_type_%s_md_%s.xlsx", id, config$input_date))

  counts <- as.data.frame(read.csv(count_path, row.names = 1, check.names = FALSE))
  metadata <- as.data.frame(read_excel(metadata_path))
  rownames(metadata) <- metadata[[1]]
  metadata <- metadata[, -1, drop = FALSE]
  metadata$age_binary <- ifelse(metadata$age == "Yng", 0, 1)
  metadata$pair <- factor(metadata$pair)
  metadata$Area.um2 <- as.numeric(metadata$Area.um2)

  aligned <- align_counts_metadata(counts, metadata, id, logger$log)
  if (is.null(aligned)) next
  counts <- aligned$counts
  metadata <- aligned$metadata

  use_pair_re <- nlevels(metadata$pair) > 1
  formula <- if (use_pair_re) ~ age_binary + Area.um2 + (1 | pair) else ~ age_binary + Area.um2
  method <- if (use_pair_re) "blmm" else "lm"
  filtered <- filter_genes_by_detection(counts, config$detection_threshold)
  if (nrow(filtered) == 0) next

  fixed_design <- model.matrix(lme4::nobars(formula), data = metadata)
  if (ncol(fixed_design) > nrow(fixed_design)) {
    logger$log(sprintf("[%s] skipped: P > N", id))
    next
  }

  log_scale <- build_informed_log_scale(metadata$age_binary, config$nsample,
                                        config$gamma_value, config$scale_sd)
  return_parameters <- c("X", "estimate", "std.error", "p.val",
                         "p.val.adj", "logComp", "logScale")
  if (method == "blmm") return_parameters <- c(return_parameters, "random.effects")

  fit <- aldex(Y = filtered, X = formula, data = metadata, method = method,
               nsample = config$nsample, scale = log_scale,
               n.cores = config$n_cores, return.pars = return_parameters)
  results <- build_results_table(fit)
  output_path <- file.path(config$output_dir,
    sprintf("cosmx_data_ct_blmm_%s_results_%s.xlsx", id, config$output_date))
  write_xlsx(results, output_path)
  logger$log(sprintf("[%s] wrote: %s", id, output_path))
  rm(counts, metadata, aligned, filtered, fixed_design, log_scale, fit, results)
  gc()
}
