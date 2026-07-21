library(dplyr)
library(Matrix)
library(glmmTMB)
library(stringr)
library(parallel)
library(writexl)

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), mustWork = FALSE)), "path_config.R"))

in_dir <- simulation_generated_dir
out_dir <- simulation_nb_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
target_sim <- suppressWarnings(as.integer(Sys.getenv("TARGET_SIM")))
if (is.na(target_sim)) stop("TARGET_SIM not set")
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))

fit_nb_gene <- function(counts, metadata, gene_name) {
  model_data <- metadata
  model_data$y <- as.numeric(counts)
  model_data$offset_log_depth <- log(model_data$observed_total)
  if (all(model_data$y == 0)) return(NULL)
  fit <- tryCatch(
    glmmTMB(y ~ age_binary + offset(offset_log_depth) + (1 | batch),
            family = nbinom2, data = model_data),
    error = function(error) NULL
  )
  if (is.null(fit)) return(NULL)
  coefficient_table <- summary(fit)$coefficients$cond
  if (!"age_binary" %in% rownames(coefficient_table)) return(NULL)
  coefficient <- coefficient_table["age_binary", , drop = FALSE]
  data.frame(
    gene = gene_name,
    age_coef = coefficient[1, "Estimate"],
    age_se = coefficient[1, "Std. Error"],
    age_z = coefficient[1, "z value"],
    age_pval = coefficient[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

manifest_files <- list.files(in_dir, pattern = "^manifest_.*_sim_\\d{3}_N_\\d{3}\\.rds$")
dataset_table <- bind_rows(lapply(manifest_files, function(filename) {
  match <- str_match(filename, "^manifest_(.+)_sim_(\\d{3})_N_(\\d{3})\\.rds$")
  if (any(is.na(match))) return(NULL)
  data.frame(ct_region = match[2], sim_id = as.integer(match[3]), N_batches = as.integer(match[4]))
}))
dataset_table <- dataset_table[dataset_table$sim_id == target_sim, , drop = FALSE]
dataset_table <- dataset_table[order(dataset_table$N_batches), , drop = FALSE]

for (i in seq_len(nrow(dataset_table))) {
  rep_id <- sprintf("%s_sim_%03d_N_%03d", dataset_table$ct_region[i],
                    dataset_table$sim_id[i], dataset_table$N_batches[i])
  message("[", i, "/", nrow(dataset_table), "] ", rep_id)
  obs <- read.csv(file.path(in_dir, paste0("obs_", rep_id, ".csv")))
  var <- read.csv(file.path(in_dir, paste0("var_", rep_id, ".csv")))
  counts <- Matrix::readMM(file.path(in_dir, paste0("counts_", rep_id, "_observed.mtx")))
  rownames(counts) <- var$gene
  colnames(counts) <- obs$cell_id
  metadata <- obs
  rownames(metadata) <- metadata$cell_id
  metadata$batch <- factor(metadata$batch)
  metadata$age_binary <- as.numeric(metadata$age_binary)
  filtered_counts <- counts[as.logical(var$keep_gene), , drop = FALSE]
  if (nrow(filtered_counts) == 0L) next
  result_list <- mclapply(seq_len(nrow(filtered_counts)), function(gene_index) {
    fit_nb_gene(filtered_counts[gene_index, ], metadata, rownames(filtered_counts)[gene_index])
  }, mc.cores = ncores)
  results <- bind_rows(result_list)
  if (nrow(results) == 0L) next
  results$age_padj <- p.adjust(results$age_pval, method = "BH")
  saveRDS(results, file.path(out_dir, paste0("nb_", rep_id, "_observed.rds")))
}


# ============================================================
# COMBINE CURRENTLY AVAILABLE NEGATIVE-BINOMIAL RESULTS
# ============================================================

result_files <- list.files(
  out_dir,
  pattern = "^nb_.*_observed\\.rds$",
  full.names = TRUE
)

if (length(result_files) > 0L) {
  combined_results <- bind_rows(lapply(result_files, function(path) {
    result <- readRDS(path)
    match <- str_match(
      basename(path),
      "^nb_(.+)_sim_(\\d+)_N_(\\d+)_observed\\.rds$"
    )
    if (any(is.na(match))) stop("Unexpected NB filename: ", basename(path))

    result$ct_region <- match[2]
    result$sim_id <- as.integer(match[3])
    result$N_batches <- as.integer(match[4])
    result$condition <- "observed"
    result$source_file <- basename(path)
    result
  }))

  write_xlsx(
    list(observed = combined_results),
    file.path(out_dir, "nb_observed_simulation_results.xlsx")
  )
}
