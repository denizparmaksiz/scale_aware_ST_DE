# GO biological-process enrichment of ALDEx results
#
# Adapted from the GO enrichment workflow accompanying Sun et al.
# Original source: add the Sun et al. repository URL and citation here.
#
# Modifications for the present analysis:
# - accepts an Excel workbook containing ALDEx results, with one analysis per sheet;
# - selects significant increasing and decreasing genes using ALDEx effect estimates
#   and Benjamini-Hochberg-adjusted p-values;
# - retains the original topGO "classic" Fisher-test workflow;
# - writes all enriched GO terms to a multi-sheet Excel workbook.

library(topGO)
library(org.Mm.eg.db)
library(readxl)
library(writexl)

# ============================================================
# SETTINGS
# ============================================================

config <- list(
  aldex_results_path = "path/to/aldex_results.xlsx",
  output_dir = "path/to/output_directory",
  output_filename = "go_aldex_BP.xlsx",
  ontology = "BP",
  effect_column = "age_scaled:est",
  adjusted_p_column = "age_scaled:pval.adj",
  adjusted_p_threshold = 0.05,
  node_size = 1
)

# ============================================================
# HELPERS
# ============================================================

validate_aldex_sheet <- function(data,
                                 sheet_name,
                                 effect_column,
                                 adjusted_p_column) {
  required_columns <- c("gene", effect_column, adjusted_p_column)
  missing_columns <- setdiff(required_columns, colnames(data))

  if (length(missing_columns) > 0L) {
    stop(
      "Sheet '", sheet_name, "' is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

make_result_name <- function(sheet_name, direction) {
  analysis_name <- sub("^[^_]*_", "", sheet_name)
  paste0(analysis_name, "_", direction)
}

run_direction_enrichment <- function(data,
                                     selected,
                                     result_name,
                                     ontology,
                                     go_to_genes,
                                     node_size) {
  if (!any(selected)) {
    message("No selected genes for ", result_name, "; skipping.")
    return(NULL)
  }

  if (all(selected)) {
    message("All genes selected for ", result_name, "; skipping.")
    return(NULL)
  }

  gene_indicator <- factor(
    as.integer(selected),
    levels = c(0, 1)
  )
  names(gene_indicator) <- data$gene

  selected_genes <- data$gene[selected]

  go_data <- methods::new(
    "topGOdata",
    ontology = ontology,
    allGenes = gene_indicator,
    annot = annFUN.GO2genes,
    GO2genes = go_to_genes,
    nodeSize = node_size
  )

  term_to_genes <- genesInTerm(go_data)
  fisher_result <- runTest(
    go_data,
    algorithm = "classic",
    statistic = "fisher"
  )

  n_terms <- nrow(topGO::termStat(go_data))
  enrichment <- GenTable(
    go_data,
    Fisher = fisher_result,
    orderBy = "Fisher",
    topNodes = n_terms,
    numChar = 1000
  )

  enrichment$Fisher <- as.numeric(enrichment$Fisher)
  enrichment <- enrichment[
    !is.na(enrichment$Fisher) & enrichment$Fisher < 0.05,
    ,
    drop = FALSE
  ]

  if (nrow(enrichment) == 0L) {
    message("No enriched GO terms for ", result_name, ".")
    return(NULL)
  }

  enrichment$genes <- vapply(
    enrichment$GO.ID,
    function(go_id) {
      paste(
        sort(intersect(selected_genes, term_to_genes[[go_id]])),
        collapse = " "
      )
    },
    character(1)
  )

  enrichment
}

# ============================================================
# ALDEX GO ENRICHMENT
# ============================================================

run_aldex_enrichment <- function(excel_path,
                                 ontology = "BP",
                                 output_path,
                                 effect_column = "age_scaled:est",
                                 adjusted_p_column = "age_scaled:pval.adj",
                                 adjusted_p_threshold = 0.05,
                                 node_size = 1) {
  if (!file.exists(excel_path)) {
    stop("ALDEx result workbook not found: ", excel_path)
  }

  output_parent <- dirname(output_path)
  dir.create(output_parent, showWarnings = FALSE, recursive = TRUE)

  sheet_names <- excel_sheets(excel_path)
  if (length(sheet_names) == 0L) {
    stop("The ALDEx workbook contains no worksheets.")
  }

  go_to_genes <- annFUN.org(
    whichOnto = ontology,
    feasibleGenes = NULL,
    mapping = "org.Mm.eg.db",
    ID = "symbol"
  )

  all_results <- list()

  for (sheet_name in sheet_names) {
    message("Processing sheet: ", sheet_name)

    data <- read_excel(excel_path, sheet = sheet_name)
    validate_aldex_sheet(
      data,
      sheet_name,
      effect_column,
      adjusted_p_column
    )

    data <- data[
      complete.cases(
        data[["gene"]],
        data[[effect_column]],
        data[[adjusted_p_column]]
      ),
      ,
      drop = FALSE
    ]

    if (nrow(data) == 0L) {
      message("No complete rows in ", sheet_name, "; skipping.")
      next
    }

    directions <- list(
      increasing = data[[effect_column]] > 0 &
        data[[adjusted_p_column]] < adjusted_p_threshold,
      decreasing = data[[effect_column]] < 0 &
        data[[adjusted_p_column]] < adjusted_p_threshold
    )

    for (direction in names(directions)) {
      result_name <- make_result_name(sheet_name, direction)

      enrichment <- run_direction_enrichment(
        data = data,
        selected = directions[[direction]],
        result_name = result_name,
        ontology = ontology,
        go_to_genes = go_to_genes,
        node_size = node_size
      )

      if (!is.null(enrichment)) {
        all_results[[result_name]] <- enrichment
      }
    }
  }

  if (length(all_results) == 0L) {
    warning("No enriched GO terms were identified; no workbook was written.")
    return(invisible(all_results))
  }

  write_xlsx(all_results, path = output_path)
  message("Saved GO enrichment results to: ", output_path)

  invisible(all_results)
}

# ============================================================
# RUN ANALYSIS
# ============================================================

run_enrichments <- function(config) {
  required_settings <- c("aldex_results_path", "output_dir", "output_filename")
  missing_settings <- required_settings[
    !vapply(config[required_settings], nzchar, logical(1))
  ]

  if (length(missing_settings) > 0L) {
    stop(
      "Missing configuration values: ",
      paste(missing_settings, collapse = ", ")
    )
  }

  output_path <- file.path(
    config$output_dir,
    config$output_filename
  )

  message("Running ALDEx GO enrichment (", config$ontology, ")")

  run_aldex_enrichment(
    excel_path = config$aldex_results_path,
    ontology = config$ontology,
    output_path = output_path,
    effect_column = config$effect_column,
    adjusted_p_column = config$adjusted_p_column,
    adjusted_p_threshold = config$adjusted_p_threshold,
    node_size = config$node_size
  )
}

if (sys.nframe() == 0L) {
  run_enrichments(config)
}
