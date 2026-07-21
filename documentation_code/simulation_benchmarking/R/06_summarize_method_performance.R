library(dplyr)
library(readr)
library(stringr)
library(writexl)

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), mustWork = FALSE)), "path_config.R"))

alpha <- 0.05
beta <- 0.5
aldex_dir <- simulation_aldex_dir
mast_dir <- simulation_mast_dir
nb_dir <- simulation_nb_dir
scanpy_dir <- simulation_scanpy_dir
manifest_dir <- simulation_generated_dir
out_dir <- simulation_benchmark_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_beta <- function(ppv, power, beta = 0.5) {
  ifelse(is.na(ppv) | is.na(power) | (ppv == 0 & power == 0), NA_real_,
         (1 + beta^2) * ppv * power / (beta^2 * ppv + power))
}

compute_confusion <- function(results, truth_genes, adjusted_p_column, alpha = 0.05) {
  genes <- as.character(results$gene)
  is_de <- genes %in% truth_genes
  significant <- !is.na(results[[adjusted_p_column]]) & results[[adjusted_p_column]] < alpha
  TP <- sum(significant & is_de)
  FP <- sum(significant & !is_de)
  FN <- sum(!significant & is_de)
  tibble(
    n_genes = length(genes),
    n_de_in_results = sum(is_de),
    TP = TP, FP = FP, FN = FN,
    power = ifelse(sum(is_de) > 0, TP / sum(is_de), NA_real_),
    emp_fdr = ifelse(TP + FP > 0, FP / (TP + FP), NA_real_),
    ppv = ifelse(TP + FP > 0, TP / (TP + FP), NA_real_)
  )
}

manifest_files <- list.files(manifest_dir, pattern = "^manifest_.*_sim_\\d{3}_N_\\d{3}\\.rds$", full.names = TRUE)
dataset_table <- bind_rows(lapply(manifest_files, function(path) {
  match <- str_match(basename(path), "^manifest_(.+)_sim_(\\d{3})_N_(\\d{3})\\.rds$")
  if (any(is.na(match))) return(NULL)
  data.frame(manifest_file = path, ct_region = match[2], sim_id = as.integer(match[3]),
             N_batches = as.integer(match[4]))
}))

summaries <- list()
for (i in seq_len(nrow(dataset_table))) {
  row <- dataset_table[i, ]
  rep_id <- sprintf("%s_sim_%03d_N_%03d", row$ct_region, row$sim_id, row$N_batches)
  manifest <- readRDS(row$manifest_file)
  truth_genes <- unique(as.character(manifest$de_table$gene))
  method_specs <- list(
    list(method = "ALDEx", condition = "observed_tss",
         file = file.path(aldex_dir, paste0("aldex_", rep_id, "_observed_tss.rds")),
         padj = "age_binary:pval.adj"),
    list(method = "ALDEx", condition = "observed_it",
         file = file.path(aldex_dir, paste0("aldex_", rep_id, "_observed_it.rds")),
         padj = "age_binary:pval.adj"),
    list(method = "MAST", condition = "observed",
         file = file.path(mast_dir, paste0("mast_", rep_id, "_observed.rds")),
         padj = "hurdle_padj"),
    list(method = "NB", condition = "observed",
         file = file.path(nb_dir, paste0("nb_", rep_id, "_observed.rds")),
         padj = "age_padj"),
    list(method = "Scanpy", condition = "observed",
         file = file.path(scanpy_dir, paste0("scanpy_", rep_id, "_observed.csv")),
         padj = "age_perm:pval.adj")
  )
  for (spec in method_specs) {
    if (!file.exists(spec$file)) next
    result <- if (grepl("\\.rds$", spec$file)) readRDS(spec$file) else read_csv(spec$file, show_col_types = FALSE)
    confusion <- compute_confusion(result, truth_genes, spec$padj, alpha)
    summaries[[length(summaries) + 1]] <- tibble(
      method = spec$method, condition = spec$condition, ct_region = row$ct_region,
      sim_id = row$sim_id, N_batches = row$N_batches
    ) %>% bind_cols(confusion) %>% mutate(f0_5 = f_beta(ppv, power, beta))
  }
}

plot_ready <- bind_rows(summaries)
average_by_N <- plot_ready %>%
  group_by(method, condition, N_batches) %>%
  summarise(
    n_sims = n_distinct(sim_id),
    power = mean(power, na.rm = TRUE),
    emp_fdr = mean(emp_fdr, na.rm = TRUE),
    ppv = mean(ppv, na.rm = TRUE),
    f0_5 = mean(f0_5, na.rm = TRUE),
    .groups = "drop"
  )
write_xlsx(list(simulation_level = plot_ready, average_by_N = average_by_N),
           file.path(out_dir, "simulation_method_performance.xlsx"))
