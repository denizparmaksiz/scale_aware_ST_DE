#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ALDEx3)
  library(lme4)
  library(readxl)
  library(stringr)
})
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "R/07_sun_continuous_subregion_blmm.R"
source(file.path(dirname(normalizePath(script_path)), "aldex_helpers.R"))

config <- parse_cli_args(list(
  input_dir = Sys.getenv("ALDEX_INPUT_DIR", unset = "cross_study_cont_subregion"),
  output_dir = Sys.getenv("ALDEX_OUTPUT_DIR", unset = "results/sun_continuous_subregion_blmm"),
  nsample = 2000,
  n_cores = 10,
  detection_threshold = 0.20,
  gamma_value = 0.80,
  scale_sd = 0.25,
  input_date = "05042026",
  output_date = "05042026"
))
config$output_dir <- ensure_output_dir(config$output_dir)
logger <- make_logger("sun_continuous_subregion_blmm")

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
  metadata$age_cont <- as.numeric(metadata$age)
  metadata$age_scaled <- (metadata$age_cont - min(metadata$age_cont)) /
                         (max(metadata$age_cont) - min(metadata$age_cont))
  metadata$batch <- factor(metadata$batch)
  metadata$slide_id <- factor(metadata$slide_id)
  metadata$volume <- as.numeric(metadata$volume)

  aligned <- align_counts_metadata(counts, metadata, id, logger$log)
  if (is.null(aligned)) next
  counts <- aligned$counts
  metadata <- aligned$metadata

  if (nlevels(metadata$slide_id) < 2) {
    logger$log(sprintf("[%s] skipped: only one slide_id", id))
    next
  }
  if (length(unique(metadata$age_scaled)) < 2) {
    logger$log(sprintf("[%s] skipped: no age variation", id))
    next
  }

  formula <- ~ age_scaled + volume + (1 | slide_id)
  filtered <- filter_genes_by_detection(counts, config$detection_threshold)
  if (nrow(filtered) == 0) {
    logger$log(sprintf("[%s] skipped: no genes met detection threshold", id))
    next
  }

  log_scale <- build_informed_log_scale(metadata$age_scaled, config$nsample,
                                        config$gamma_value, config$scale_sd)
  logger$log(sprintf("[%s] running blmm, nsample=%d", id, config$nsample))
  logger$memory(sprintf("%s | before ALDEx3", id))
  fit <- tryCatch(
    aldex(filtered, formula, data = metadata, method = "blmm",
          nsample = config$nsample, scale = log_scale,
          n.cores = config$n_cores,
          return.pars = c("estimate", "std.error", "p.val", "p.val.adj")),
    error = function(e) {
      logger$log(sprintf("[%s] ERROR: %s", id, e$message))
      NULL
    }
  )
  if (is.null(fit)) next

  results <- build_results_table(fit)
  output_path <- file.path(config$output_dir,
    sprintf("cross_study_ct_%s_cont_it_results_%s.rds", id, config$output_date))
  saveRDS(results, output_path)
  logger$log(sprintf("[%s] wrote: %s", id, output_path))
  logger$memory(sprintf("%s | after write", id))
  rm(counts, metadata, aligned, filtered, log_scale, fit, results)
  gc()
}
