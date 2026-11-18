# Pre-binning Kraken report data

read_pre_binning_kraken_report <- function(path) {
  sample <- sub("[.]tsv$", "", basename(path))
  report <- readr::read_tsv(
    path,
    col_names = c(
      "report_percent", "clade_reads", "direct_reads",
      "rank", "tax_id", "taxon"
    ),
    col_types = readr::cols(
      report_percent = readr::col_double(),
      clade_reads = readr::col_double(),
      direct_reads = readr::col_double(),
      rank = readr::col_character(),
      tax_id = readr::col_double(),
      taxon = readr::col_character()
    ),
    show_col_types = FALSE,
    trim_ws = TRUE
  ) %>%
    dplyr::mutate(taxon = stringr::str_trim(taxon))

  get_taxon_value <- function(tax_id, column) {
    value <- report[report$tax_id == tax_id, column, drop = TRUE]
    if (length(value) != 1) {
      stop("Expected one Kraken entry for taxon ID ", tax_id, " in ", path)
    }
    value
  }

  unclassified_reads <- get_taxon_value(0, "clade_reads")
  classified_reads <- get_taxon_value(1, "clade_reads")
  bacteria_reads <- get_taxon_value(2, "clade_reads")
  neisseria_reads <- get_taxon_value(482, "clade_reads")
  ng_reads <- get_taxon_value(485, "clade_reads")

  tibble::tibble(
    sample = sample,
    total_reads = unclassified_reads + classified_reads,
    unclassified_reads = unclassified_reads,
    unclassified_percent = get_taxon_value(0, "report_percent"),
    bacteria_reads = bacteria_reads,
    bacteria_percent = get_taxon_value(2, "report_percent"),
    neisseria_reads = neisseria_reads,
    neisseria_percent = get_taxon_value(482, "report_percent"),
    ng_reads = ng_reads,
    ng_percent = get_taxon_value(485, "report_percent"),
    bacteria_as_neisseria_percent = 100 * neisseria_reads / bacteria_reads,
    neisseria_as_ng_percent = 100 * ng_reads / neisseria_reads
  )
}

load_pre_binning_kraken_data <- function(sample_pool_data = NULL) {
  kraken_files <- list.files(
    "./data/NG_probe/kraken2/pre-binning",
    pattern = "[.]tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  pre_binning_kraken_summary <- purrr::map_dfr(
    kraken_files,
    read_pre_binning_kraken_report
  ) %>%
    dplyr::arrange(sample)

  if (!is.null(sample_pool_data)) {
    if (
      anyDuplicated(sample_pool_data$sample) ||
        !setequal(pre_binning_kraken_summary$sample, sample_pool_data$sample)
    ) {
      stop("Pre-binning Kraken reports do not match the pool metadata.")
    }

    pre_binning_kraken_summary <- pre_binning_kraken_summary %>%
      dplyr::left_join(
        sample_pool_data %>% dplyr::select(sample, pool_level),
        by = "sample"
      )
  }

  percent_columns <- grep(
    "percent$",
    names(pre_binning_kraken_summary),
    value = TRUE
  )

  if (
    length(kraken_files) != 12 ||
      nrow(pre_binning_kraken_summary) != 12 ||
      anyDuplicated(pre_binning_kraken_summary$sample) ||
      anyNA(pre_binning_kraken_summary) ||
      any(pre_binning_kraken_summary$total_reads <= 0) ||
      any(pre_binning_kraken_summary[, percent_columns] < 0) ||
      any(pre_binning_kraken_summary[, percent_columns] > 100)
  ) {
    stop("Pre-binning Kraken summaries could not be validated.")
  }

  list(pre_binning_kraken_summary = pre_binning_kraken_summary)
}
