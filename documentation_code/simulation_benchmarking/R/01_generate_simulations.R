library(dplyr)
library(Matrix)
library(seqgendiff)
library(readxl)
library(data.table)

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), mustWork = FALSE)), "path_config.R"))

curr_dir <- simulation_source_dir
count_file <- file.path(curr_dir, "adata_nn_raw_counts_transposed.csv")
metadata_file <- file.path(curr_dir, "adata_nn_raw_counts_transposed_MD.xlsx")
design_file <- file.path(curr_dir, "updated_sim_design.rds")
out_dir <- simulation_generated_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

detection_threshold <- 0.20
prop_de <- 0.20
prop_down_old <- 0.80
log2fc_min <- 0.5
log2fc_max <- 4

sample_cells_by_batch <- function(meta_ct, batch_col = "batch", N_batches, seed) {
  set.seed(seed)
  batches <- unique(meta_ct[[batch_col]])
  if (length(batches) == 0L) stop("No batches found in the selected stratum.")
  sampled_batches <- sample(batches, size = N_batches, replace = N_batches > length(batches))
  idx_list <- lapply(sampled_batches, function(batch_now) which(meta_ct[[batch_col]] == batch_now))
  list(
    cell_idx = unlist(idx_list, use.names = FALSE),
    sampled_batches = sampled_batches,
    batch_draw_id = rep(seq_along(sampled_batches), vapply(idx_list, length, integer(1)))
  )
}

duplicate_subset <- function(counts_ct, meta_ct, cell_idx, batch_draw_id) {
  counts_sub <- counts_ct[, cell_idx, drop = FALSE]
  meta_sub <- meta_ct[cell_idx, , drop = FALSE]
  original_ids <- as.character(colnames(counts_sub))
  new_ids <- make.unique(paste0("draw", batch_draw_id, "__orig_", original_ids))
  colnames(counts_sub) <- new_ids
  rownames(meta_sub) <- new_ids
  meta_sub$cell_id <- new_ids
  meta_sub$cell_id_original <- original_ids
  meta_sub$batch_draw_id <- batch_draw_id
  list(counts = counts_sub, metadata = meta_sub)
}

permute_age_within_batch_draw <- function(metadata, age_col = "Age", draw_col = "batch_draw_id", seed) {
  set.seed(seed)
  permuted_age <- metadata[[age_col]]
  for (draw_now in unique(metadata[[draw_col]])) {
    idx <- which(metadata[[draw_col]] == draw_now)
    permuted_age[idx] <- sample(metadata[[age_col]][idx], size = length(idx), replace = FALSE)
  }
  permuted_age
}

make_effects <- function(counts_sub, prop_de, prop_down_old, log2fc_min, log2fc_max,
                         detection_threshold, seed) {
  set.seed(seed)
  genes <- rownames(counts_sub)
  detection_fraction <- rowMeans(counts_sub > 0, na.rm = TRUE)
  eligible_genes <- genes[detection_fraction >= detection_threshold]
  if (length(eligible_genes) == 0L) stop("No genes pass the expression-fraction threshold.")
  n_de <- max(1L, ceiling(prop_de * length(eligible_genes)))
  de_genes <- sample(eligible_genes, size = n_de, replace = FALSE)
  n_down <- round(prop_down_old * n_de)
  n_up <- n_de - n_down
  effect_signs <- sample(c(rep(-1, n_down), rep(1, n_up)))
  effect_sizes <- effect_signs * runif(n_de, min = log2fc_min, max = log2fc_max)
  coefficient_vector <- setNames(rep(0, length(genes)), genes)
  coefficient_vector[de_genes] <- effect_sizes
  coefficient_matrix <- matrix(coefficient_vector, ncol = 1, dimnames = list(genes, "age01"))
  list(
    coefficient_matrix = coefficient_matrix,
    de_table = data.frame(gene = de_genes, true_log2fc = effect_sizes, stringsAsFactors = FALSE)
  )
}

draw_depths_from_batch_pool <- function(original_counts, original_metadata, simulated_metadata,
                                        batch_col = "batch", seed) {
  set.seed(seed)
  original_totals <- colSums(original_counts)
  totals_by_batch <- split(original_totals, as.character(original_metadata[[batch_col]]))
  drawn_depths <- vapply(seq_len(nrow(simulated_metadata)), function(i) {
    batch_now <- as.character(simulated_metadata[[batch_col]][i])
    depth_pool <- totals_by_batch[[batch_now]]
    if (is.null(depth_pool) || length(depth_pool) == 0L) stop("No depth pool found for batch ", batch_now)
    sample(depth_pool, size = 1, replace = TRUE)
  }, numeric(1))
  names(drawn_depths) <- rownames(simulated_metadata)
  drawn_depths
}

multinomial_resample_counts <- function(latent_counts, depth_draws, seed) {
  set.seed(seed)
  if (!identical(colnames(latent_counts), names(depth_draws))) {
    stop("Depth-draw names do not match latent-count columns.")
  }
  observed_counts <- matrix(0L, nrow = nrow(latent_counts), ncol = ncol(latent_counts),
                            dimnames = dimnames(latent_counts))
  for (j in seq_len(ncol(latent_counts))) {
    latent_column <- latent_counts[, j]
    latent_total <- sum(latent_column)
    observed_depth <- as.integer(depth_draws[j])
    if (latent_total <= 0 || observed_depth <= 0) next
    observed_counts[, j] <- as.integer(
      rmultinom(1, size = observed_depth, prob = latent_column / latent_total)[, 1]
    )
  }
  observed_counts
}

simulate_one_dataset <- function(counts, metadata, design_row) {
  ct_region <- design_row$ct_region
  N_batches <- design_row$N_batches
  keep_cells <- metadata$ct_region == ct_region
  metadata_ct <- metadata[keep_cells, , drop = FALSE]
  counts_ct <- counts[, keep_cells, drop = FALSE]

  sampled <- sample_cells_by_batch(metadata_ct, N_batches = N_batches, seed = design_row$seed_resample)
  duplicated <- duplicate_subset(counts_ct, metadata_ct, sampled$cell_idx, sampled$batch_draw_id)
  counts_null <- duplicated$counts
  metadata_sim <- duplicated$metadata
  metadata_sim$age_perm <- permute_age_within_batch_draw(metadata_sim, seed = design_row$seed_permute)
  metadata_sim$age_binary <- ifelse(metadata_sim$age_perm == "Yng", 0L, 1L)

  design_fixed <- matrix(metadata_sim$age_binary, ncol = 1,
                         dimnames = list(rownames(metadata_sim), "age01"))
  effects <- make_effects(
    counts_null, prop_de, prop_down_old, log2fc_min, log2fc_max,
    detection_threshold, design_row$seed_effects
  )

  set.seed(design_row$seed_thin)
  thinning_result <- seqgendiff::thin_diff(
    mat = counts_null,
    design_fixed = design_fixed,
    coef_fixed = effects$coefficient_matrix,
    design_perm = NULL,
    coef_perm = NULL,
    target_cor = NULL,
    type = "thin",
    change_colnames = FALSE
  )
  counts_latent <- thinning_result$mat
  dimnames(counts_latent) <- dimnames(counts_null)

  depth_draws <- draw_depths_from_batch_pool(
    counts_ct, metadata_ct, metadata_sim, seed = design_row$seed_depth
  )
  counts_observed <- multinomial_resample_counts(
    counts_latent, depth_draws, seed = design_row$seed_multinom
  )

  detection_null <- rowMeans(counts_null > 0, na.rm = TRUE)
  detection_observed <- rowMeans(counts_observed > 0, na.rm = TRUE)
  keep_gene <- detection_null >= detection_threshold & detection_observed >= detection_threshold
  effect_map <- setNames(effects$de_table$true_log2fc, effects$de_table$gene)

  obs_df <- metadata_sim %>%
    mutate(
      cell_id = rownames(metadata_sim),
      depth_drawn = as.numeric(depth_draws[cell_id]),
      latent_total = as.numeric(colSums(counts_latent)[cell_id]),
      observed_total = as.numeric(colSums(counts_observed)[cell_id])
    ) %>%
    select(cell_id, cell_id_original, age_perm, age_binary, batch, batch_draw_id,
           depth_drawn, latent_total, observed_total)

  var_df <- data.frame(
    gene = rownames(counts_observed),
    keep_gene = keep_gene,
    is_de = rownames(counts_observed) %in% effects$de_table$gene,
    true_log2fc = unname(effect_map[rownames(counts_observed)]),
    stringsAsFactors = FALSE
  )
  var_df$true_log2fc[is.na(var_df$true_log2fc)] <- 0

  stopifnot(
    identical(obs_df$cell_id, colnames(counts_observed)),
    identical(var_df$gene, rownames(counts_observed)),
    all(counts_observed >= 0),
    all(counts_observed == round(counts_observed))
  )

  manifest <- list(
    ct_region = ct_region,
    N_batches = N_batches,
    sim_id = design_row$sim_id,
    seeds = design_row[c("seed_resample", "seed_permute", "seed_thin", "seed_effects",
                         "seed_depth", "seed_multinom")],
    sampled_batches = sampled$sampled_batches,
    de_table = effects$de_table,
    keep_genes = var_df$gene[var_df$keep_gene],
    parameters = list(
      detection_threshold = detection_threshold,
      prop_de = prop_de,
      prop_down_old = prop_down_old,
      log2fc_min = log2fc_min,
      log2fc_max = log2fc_max
    ),
    age_permutation_table = table(metadata_sim$age_perm)
  )

  list(counts_observed = counts_observed, obs = obs_df, var = var_df, manifest = manifest)
}

required_design_columns <- c(
  "ct_region", "N_batches", "sim_id", "seed_resample", "seed_permute",
  "seed_thin", "seed_effects", "seed_depth", "seed_multinom"
)
simulation_design <- as.data.frame(readRDS(design_file), stringsAsFactors = FALSE)
missing_design_columns <- setdiff(required_design_columns, colnames(simulation_design))
if (length(missing_design_columns) > 0) {
  stop("Simulation design is missing columns: ", paste(missing_design_columns, collapse = ", "))
}

count_df <- fread(count_file, data.table = FALSE, header = TRUE, integer64 = "double")
gene_names <- count_df[[1]]
count_df <- count_df[, -1, drop = FALSE]
counts <- as.matrix(count_df)
rownames(counts) <- gene_names

metadata <- read_excel(metadata_file) |> as.data.frame()
rownames(metadata) <- as.character(metadata[[1]])
metadata <- metadata[, -1, drop = FALSE]
metadata$batch <- factor(metadata$batch)
metadata$volume <- as.numeric(metadata$volume)

common_cells <- intersect(colnames(counts), rownames(metadata))
if (length(common_cells) == 0L) stop("No matching cells between counts and metadata.")
common_cells <- colnames(counts)[colnames(counts) %in% common_cells]
counts <- counts[, common_cells, drop = FALSE]
metadata <- metadata[common_cells, , drop = FALSE]
stopifnot(identical(colnames(counts), rownames(metadata)))

target_sim <- suppressWarnings(as.integer(Sys.getenv("TARGET_SIM")))
if (!is.na(target_sim)) {
  simulation_design <- simulation_design[simulation_design$sim_id == target_sim, , drop = FALSE]
}
if (nrow(simulation_design) == 0L) stop("No simulation-design rows selected.")

for (i in seq_len(nrow(simulation_design))) {
  design_row <- simulation_design[i, , drop = FALSE]
  rep_id <- sprintf("%s_sim_%03d_N_%03d", design_row$ct_region, design_row$sim_id,
                    design_row$N_batches)
  message("[", i, "/", nrow(simulation_design), "] ", rep_id)
  result <- simulate_one_dataset(counts, metadata, design_row)
  Matrix::writeMM(Matrix(result$counts_observed, sparse = TRUE),
                  file.path(out_dir, paste0("counts_", rep_id, "_observed.mtx")))
  write.csv(result$obs, file.path(out_dir, paste0("obs_", rep_id, ".csv")),
            row.names = FALSE, quote = TRUE)
  write.csv(result$var, file.path(out_dir, paste0("var_", rep_id, ".csv")), row.names = FALSE)
  saveRDS(result$manifest, file.path(out_dir, paste0("manifest_", rep_id, ".rds")))
}
