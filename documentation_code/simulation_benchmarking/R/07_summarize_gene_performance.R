library(dplyr)
library(readr)
library(stringr)
library(writexl)

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), mustWork = FALSE)), "path_config.R"))

aldex_dir <- simulation_aldex_dir
mast_dir <- simulation_mast_dir
nb_dir <- simulation_nb_dir
scanpy_dir <- simulation_scanpy_dir
manifest_dir <- simulation_generated_dir
out_file <- file.path(simulation_benchmark_dir, "simulation_gene_performance.xlsx")
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)

manifest_files <- list.files(manifest_dir, pattern = "^manifest_.*_sim_\\d{3}_N_\\d{3}\\.rds$", full.names = TRUE)
dataset_table <- bind_rows(lapply(manifest_files, function(path) {
  match <- str_match(basename(path), "^manifest_(.+)_sim_(\\d{3})_N_(\\d{3})\\.rds$")
  if (any(is.na(match))) return(NULL)
  data.frame(manifest_file = path, ct_region = match[2], sim_id = as.integer(match[3]),
             N_batches = as.integer(match[4]))
}))

all_results <- list()
for (i in seq_len(nrow(dataset_table))) {
  row <- dataset_table[i, ]
  rep_id <- sprintf("%s_sim_%03d_N_%03d", row$ct_region, row$sim_id, row$N_batches)
  manifest <- readRDS(row$manifest_file)
  truth <- manifest$de_table %>%
    transmute(gene = gene, injected = TRUE, true_effect = true_log2fc)

  append_result <- function(result, method, condition, coefficient, p_value, p_adjusted) {
    result %>%
      left_join(truth, by = "gene") %>%
      mutate(
        injected = ifelse(is.na(injected), FALSE, injected),
        true_effect = ifelse(is.na(true_effect), 0, true_effect)
      ) %>%
      transmute(
        gene, injected, true_effect,
        coef = .data[[coefficient]],
        pval = .data[[p_value]],
        padj = .data[[p_adjusted]],
        method = method, condition = condition,
        ct_region = row$ct_region, sim_id = row$sim_id, N_batches = row$N_batches
      )
  }

  specs <- list(
    list(file = file.path(aldex_dir, paste0("aldex_", rep_id, "_observed_tss.rds")),
         type = "rds", method = "ALDEx", condition = "observed_tss",
         coef = "age_binary:est", p = "age_binary:pval", padj = "age_binary:pval.adj"),
    list(file = file.path(aldex_dir, paste0("aldex_", rep_id, "_observed_it.rds")),
         type = "rds", method = "ALDEx", condition = "observed_it",
         coef = "age_binary:est", p = "age_binary:pval", padj = "age_binary:pval.adj"),
    list(file = file.path(mast_dir, paste0("mast_", rep_id, "_observed.rds")),
         type = "rds", method = "MAST", condition = "observed",
         coef = "age_cont_coef", p = "hurdle_pval", padj = "hurdle_padj"),
    list(file = file.path(nb_dir, paste0("nb_", rep_id, "_observed.rds")),
         type = "rds", method = "NB", condition = "observed",
         coef = "age_coef", p = "age_pval", padj = "age_padj"),
    list(file = file.path(scanpy_dir, paste0("scanpy_", rep_id, "_observed.csv")),
         type = "csv", method = "Scanpy", condition = "observed",
         coef = "age_perm:logfc", p = "age_perm:pval", padj = "age_perm:pval.adj")
  )

  for (spec in specs) {
    if (!file.exists(spec$file)) next
    result <- if (spec$type == "rds") readRDS(spec$file) else read_csv(spec$file, show_col_types = FALSE)
    all_results[[length(all_results) + 1]] <- append_result(
      result, spec$method, spec$condition, spec$coef, spec$p, spec$padj
    )
  }
}

write_xlsx(list(per_gene = bind_rows(all_results)), out_file)
