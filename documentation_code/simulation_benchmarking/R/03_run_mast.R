library(dplyr)
library(Matrix)
library(MAST)
library(lme4)
library(stringr)
library(writexl)

in_dir <- "/gpfs/Home/dpp5572/simulation_results/"
out_dir <- "/gpfs/Home/dpp5572/simulation_results_mast/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
target_sim <- suppressWarnings(as.integer(Sys.getenv("TARGET_SIM")))
if (is.na(target_sim)) stop("TARGET_SIM not set")

make_sca <- function(counts_gxc, metadata) {
  metadata <- metadata[colnames(counts_gxc), , drop = FALSE]
  stopifnot(identical(rownames(metadata), colnames(counts_gxc)))
  feature_data <- data.frame(primerid = rownames(counts_gxc), row.names = rownames(counts_gxc))
  MAST::FromMatrix(exprsArray = log2(as.matrix(counts_gxc) + 1),
                   cData = metadata, fData = feature_data)
}

extract_mast_results <- function(fit, coefficient = "age_binary") {
  table_all <- as.data.frame(summary(fit, doLRT = coefficient)$datatable)
  table_age <- table_all[table_all$contrast == coefficient, , drop = FALSE]
  hurdle <- table_age[table_age$component == "H", c("primerid", "Pr(>Chisq)"), drop = FALSE]
  colnames(hurdle) <- c("gene", "hurdle_pval")
  continuous <- table_age[table_age$component == "C", c("primerid", "coef", "ci.lo", "ci.hi"), drop = FALSE]
  colnames(continuous) <- c("gene", "age_cont_coef", "age_cont_ci_lo", "age_cont_ci_hi")
  discrete <- table_age[table_age$component == "D", c("primerid", "coef", "ci.lo", "ci.hi"), drop = FALSE]
  colnames(discrete) <- c("gene", "age_disc_coef", "age_disc_ci_lo", "age_disc_ci_hi")
  result <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE),
                   list(hurdle, continuous, discrete))
  result$hurdle_padj <- p.adjust(result$hurdle_pval, method = "BH")
  result
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
  sca <- make_sca(filtered_counts, metadata)
  fit <- tryCatch(
    zlm(~ age_binary + (1 | batch), sca, method = "glmer", ebayes = FALSE),
    error = function(error) { message("  fit failed: ", error$message); NULL }
  )
  if (is.null(fit)) next
  saveRDS(extract_mast_results(fit),
          file.path(out_dir, paste0("mast_", rep_id, "_observed.rds")))
}


# ============================================================
# COMBINE CURRENTLY AVAILABLE MAST RESULTS
# ============================================================

result_files <- list.files(
  out_dir,
  pattern = "^mast_.*_observed\\.rds$",
  full.names = TRUE
)

if (length(result_files) > 0L) {
  combined_results <- bind_rows(lapply(result_files, function(path) {
    result <- readRDS(path)
    match <- str_match(
      basename(path),
      "^mast_(.+)_sim_(\\d+)_N_(\\d+)_observed\\.rds$"
    )
    if (any(is.na(match))) stop("Unexpected MAST filename: ", basename(path))

    result$ct_region <- match[2]
    result$sim_id <- as.integer(match[3])
    result$N_batches <- as.integer(match[4])
    result$condition <- "observed"
    result$source_file <- basename(path)
    result
  }))

  write_xlsx(
    list(observed = combined_results),
    file.path(out_dir, "mast_observed_simulation_results.xlsx")
  )
}
