library(ALDEx3)
library(Matrix)
library(dplyr)
library(stringr)

in_dir <- "/gpfs/Home/dpp5572/simulation_results/"
out_dir <- "/gpfs/Home/dpp5572/simulation_results_aldex/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

S_draws <- 250
gamma_val <- 0.80
scale_sd <- 0.25
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
target_sim <- suppressWarnings(as.integer(Sys.getenv("TARGET_SIM")))
if (is.na(target_sim)) stop("TARGET_SIM not set")

tss_scale <- function(X, logComp, gamma = 0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  LambdaScale <- matrix(rnorm(P * nsample, 0, gamma), P, nsample)
  t(X) %*% LambdaScale
}

extract_aldex_results <- function(aldex_result) {
  estimate <- t(apply(aldex_result$estimate, c(1, 2), mean, na.rm = TRUE))
  std_error <- t(apply(aldex_result$std.error, c(1, 2), mean, na.rm = TRUE))
  p_value <- t(aldex_result$p.val)
  p_adjusted <- t(aldex_result$p.val.adj)
  results <- data.frame(gene = rownames(estimate), stringsAsFactors = FALSE)
  for (effect in colnames(estimate)) {
    results[[paste0(effect, ":est")]] <- estimate[, effect]
    results[[paste0(effect, ":estSE")]] <- std_error[, effect]
    results[[paste0(effect, ":pval")]] <- p_value[, effect]
    results[[paste0(effect, ":pval.adj")]] <- p_adjusted[, effect]
  }
  results
}

manifest_files <- list.files(in_dir, pattern = "^manifest_.*_sim_\\d{3}_N_\\d{3}\\.rds$", full.names = FALSE)
dataset_table <- bind_rows(lapply(manifest_files, function(filename) {
  match <- str_match(filename, "^manifest_(.+)_sim_(\\d{3})_N_(\\d{3})\\.rds$")
  if (any(is.na(match))) return(NULL)
  data.frame(ct_region = match[2], sim_id = as.integer(match[3]),
             N_batches = as.integer(match[4]), stringsAsFactors = FALSE)
}))
dataset_table <- dataset_table[dataset_table$sim_id == target_sim, , drop = FALSE]
dataset_table <- dataset_table[order(dataset_table$N_batches), , drop = FALSE]
if (nrow(dataset_table) == 0L) stop("No datasets found for sim_id = ", target_sim)

for (i in seq_len(nrow(dataset_table))) {
  rep_id <- sprintf("%s_sim_%03d_N_%03d", dataset_table$ct_region[i],
                    dataset_table$sim_id[i], dataset_table$N_batches[i])
  message("[", i, "/", nrow(dataset_table), "] ", rep_id)
  obs <- read.csv(file.path(in_dir, paste0("obs_", rep_id, ".csv")), stringsAsFactors = FALSE)
  var <- read.csv(file.path(in_dir, paste0("var_", rep_id, ".csv")), stringsAsFactors = FALSE)
  counts <- Matrix::readMM(file.path(in_dir, paste0("counts_", rep_id, "_observed.mtx")))
  rownames(counts) <- var$gene
  colnames(counts) <- obs$cell_id
  metadata <- obs
  rownames(metadata) <- metadata$cell_id
  metadata$batch <- factor(metadata$batch)
  metadata$age_binary <- as.numeric(metadata$age_binary)
  filtered_counts <- as.matrix(counts[as.logical(var$keep_gene), , drop = FALSE])
  storage.mode(filtered_counts) <- "integer"
  if (nrow(filtered_counts) == 0L) next

  model_formula <- ~ age_binary + (1 | batch)
  return_parameters <- c("X", "estimate", "std.error", "p.val", "p.val.adj",
                         "logComp", "logScale", "random.effects")
  result_tss <- aldex(
    Y = filtered_counts, X = model_formula, data = metadata, method = "blmm",
    nsample = S_draws, scale = tss_scale, n.cores = ncores,
    return.pars = return_parameters
  )
  saveRDS(extract_aldex_results(result_tss),
          file.path(out_dir, paste0("aldex_", rep_id, "_observed_tss.rds")))

  delta_age <- log2(gamma_val)
  theta_perp <- rnorm(S_draws, mean = delta_age, sd = scale_sd)
  informed_scale <- metadata$age_binary %*% t(theta_perp)
  result_informed <- aldex(
    Y = filtered_counts, X = model_formula, data = metadata, method = "blmm",
    nsample = S_draws, scale = informed_scale, n.cores = ncores,
    return.pars = return_parameters
  )
  saveRDS(extract_aldex_results(result_informed),
          file.path(out_dir, paste0("aldex_", rep_id, "_observed_it.rds")))
}


# ============================================================
# COMBINE CURRENTLY AVAILABLE ALDEx RESULTS
# ============================================================

result_files <- list.files(
  out_dir,
  pattern = "^aldex_.*_(observed_tss|observed_it)\\.rds$",
  full.names = TRUE
)

if (length(result_files) > 0L) {
  combined_results <- bind_rows(lapply(result_files, function(path) {
    result <- readRDS(path)
    match <- str_match(
      basename(path),
      "^aldex_(.+)_sim_(\\d+)_N_(\\d+)_(observed_tss|observed_it)\\.rds$"
    )
    if (any(is.na(match))) stop("Unexpected ALDEx filename: ", basename(path))

    result$ct_region <- match[2]
    result$sim_id <- as.integer(match[3])
    result$N_batches <- as.integer(match[4])
    result$condition <- match[5]
    result$source_file <- basename(path)
    result
  }))

  write_xlsx(
    list(
      observed_tss = filter(combined_results, condition == "observed_tss"),
      observed_it = filter(combined_results, condition == "observed_it")
    ),
    file.path(out_dir, "aldex_observed_simulation_results.xlsx")
  )
}
