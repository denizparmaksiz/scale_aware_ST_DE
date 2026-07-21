# Shared utilities for manuscript ALDEx3 analyses.
# These functions only handle command-line arguments, logging, validation,
# result-table construction, and output writing. Analysis-specific model
# definitions remain in the individual run scripts.

parse_cli_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) return(defaults)
  if (length(args) %% 2 != 0) {
    stop("Arguments must be supplied as --name value pairs.")
  }
  out <- defaults
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[[i]])
    key <- gsub("-", "_", key)
    if (!key %in% names(out)) stop("Unknown argument: ", args[[i]])
    value <- args[[i + 1]]
    if (is.numeric(out[[key]])) value <- as.numeric(value)
    if (is.logical(out[[key]])) value <- tolower(value) %in% c("true", "1", "yes")
    out[[key]] <- value
  }
  out
}

make_logger <- function(prefix) {
  run_id <- sample.int(1e6, 1)
  log_file <- sprintf("%s_debug_%06d.log", prefix, run_id)
  mem_file <- sprintf("%s_memory_%06d.log", prefix, run_id)

  log_line <- function(message) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), message, "\n",
        file = log_file, append = TRUE)
  }

  track_memory <- function(stage = "") {
    mem_kb <- tryCatch(
      as.numeric(system("grep VmRSS /proc/$$/status | awk '{print $2}'", intern = TRUE)),
      error = function(e) NA_real_
    )
    mem_gb <- if (is.na(mem_kb)) NA_real_ else round(mem_kb / 1024^2, 3)
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", stage,
        " | Memory (GB): ", mem_gb, "\n", file = mem_file, append = TRUE)
  }

  list(log = log_line, memory = track_memory,
       log_file = log_file, memory_file = mem_file)
}

align_counts_metadata <- function(counts, metadata, id, log_line) {
  common <- intersect(colnames(counts), rownames(metadata))
  if (length(common) == 0) {
    log_line(sprintf("[%s] ERROR: no overlapping observation IDs", id))
    return(NULL)
  }
  list(
    counts = counts[, common, drop = FALSE],
    metadata = droplevels(metadata[common, , drop = FALSE])
  )
}

filter_genes_by_detection <- function(counts, threshold) {
  keep <- if (threshold == 0) {
    rep(TRUE, nrow(counts))
  } else {
    rowMeans(counts > 0, na.rm = TRUE) >= threshold
  }
  counts[keep, , drop = FALSE]
}

build_results_table <- function(aldex_result) {
  estimate <- apply(aldex_result$estimate, c(1, 2), mean, na.rm = TRUE)
  std_error <- apply(aldex_result$std.error, c(1, 2), mean, na.rm = TRUE)
  p_value <- aldex_result$p.val
  p_adjusted <- aldex_result$p.val.adj

  estimate <- t(estimate)
  std_error <- t(std_error)
  p_value <- t(p_value)
  p_adjusted <- t(p_adjusted)

  results <- data.frame(gene = rownames(estimate), stringsAsFactors = FALSE)
  for (effect in colnames(estimate)) {
    results[[paste0(effect, ":est")]] <- estimate[, effect]
    results[[paste0(effect, ":estSE")]] <- std_error[, effect]
    results[[paste0(effect, ":pval")]] <- p_value[, effect]
    results[[paste0(effect, ":pval.adj")]] <- p_adjusted[, effect]
  }
  results
}

build_informed_log_scale <- function(design_variable, nsample,
                                     gamma_value = 0.80,
                                     scale_sd = 0.25) {
  theta_perp <- rnorm(nsample, mean = log2(gamma_value), sd = scale_sd)
  as.matrix(design_variable) %*% t(theta_perp)
}

ensure_output_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}
