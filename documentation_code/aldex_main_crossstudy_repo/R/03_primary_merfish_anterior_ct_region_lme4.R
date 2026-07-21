#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ALDEx3)
  library(lme4)
  library(readxl)
  library(stringr)
  library(writexl)
})
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "R/03_primary_merfish_anterior_ct_region_lme4.R"
source(file.path(dirname(normalizePath(script_path)), "aldex_helpers.R"))

config <- parse_cli_args(list(
  input_dir = Sys.getenv("ALDEX_INPUT_DIR", unset = "celltype_region"),
  output_dir = Sys.getenv("ALDEX_OUTPUT_DIR", unset = "results/primary_merfish_anterior_ct_region_lme4"),
  nsample = 2000,
  n_cores = 10,
  detection_threshold = 0.20,
  gamma_value = 0.80,
  scale_sd = 0.25,
  input_date = "07142025",
  output_date = "07142025"
))
config$output_dir <- ensure_output_dir(config$output_dir)
logger <- make_logger("primary_merfish_anterior_ct_region_lme4")

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
  metadata$age_binary <- ifelse(metadata$Age == "Yng", 0, 1)
  metadata$batch <- factor(metadata$batch)
  metadata$volume <- as.numeric(metadata$volume)
  metadata <- metadata[metadata$batch %in% c("1", "2", "3"), , drop = FALSE]

  aligned <- align_counts_metadata(counts, metadata, id, logger$log)
  if (is.null(aligned)) next
  counts <- aligned$counts
  metadata <- aligned$metadata

  use_batch_re <- nlevels(metadata$batch) > 1
  formula <- if (use_batch_re) {
    ~ age_binary + volume + (1 | batch)
  } else {
    ~ age_binary + volume
  }
  method <- if (use_batch_re) "lme4" else "lm"

  filtered <- filter_genes_by_detection(counts, config$detection_threshold)
  if (nrow(filtered) == 0) {
    logger$log(sprintf("[%s] skipped: no genes met detection threshold", id))
    next
  }

  fixed_design <- model.matrix(lme4::nobars(formula), data = metadata)
  if (ncol(fixed_design) > nrow(fixed_design)) {
    logger$log(sprintf("[%s] skipped: P > N (P=%d, N=%d)", id, ncol(fixed_design), nrow(fixed_design)))
    next
  }

  log_scale <- build_informed_log_scale(metadata$age_binary, config$nsample,
                                        config$gamma_value, config$scale_sd)
  return_parameters <- c("X", "estimate", "std.error", "p.val",
                         "p.val.adj", "logComp", "logScale")
  if (method == "lme4") return_parameters <- c(return_parameters, "random.effects")

  logger$log(sprintf("[%s] running method=%s, nsample=%d", id, method, config$nsample))
  logger$memory(sprintf("%s | before ALDEx3", id))

  fit <- aldex(Y = filtered, X = formula, data = metadata, method = method,
               nsample = config$nsample, scale = log_scale,
               n.cores = config$n_cores, return.pars = return_parameters)
  results <- build_results_table(fit)

  output_path <- file.path(config$output_dir,
    sprintf("aldex_mem_cell_type_anterior%s_results_%s.xlsx", id, config$output_date))
  write_xlsx(results, output_path)
  logger$log(sprintf("[%s] wrote: %s", id, output_path))
  logger$memory(sprintf("%s | after write", id))

  rm(counts, metadata, aligned, filtered, fixed_design, log_scale, fit, results)
  gc()
}
