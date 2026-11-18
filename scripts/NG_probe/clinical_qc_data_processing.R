# Clinical isolate QC and probe pool data

pilot_report_dir <- paste0(
  "./data/NG_probe/",
  "2026-10-20_pilot-strains_NG-probe_Report"
)

clean_clinical_qc_summary <- function(qcd_data) {
  qcd_data %>%
    dplyr::rename(
      initial_total_reads = initial_reads,
      num_contigs_ge_1000_bp = `# contigs (>= 1000 bp)`,
      num_contigs = `# contigs`,
      largest_contig = `Largest contig`,
      total_length = `Total length`,
      gc_percent = `GC (%)`,
      n50 = N50,
      complete_buscos = `Complete BUSCOs`,
      single_copy_buscos = `Single-copy BUSCOs`,
      duplicated_buscos = `Duplicated BUSCOs`,
      fragmented_buscos = `Fragmented BUSCOs`,
      missing_buscos = `Missing BUSCOs`,
      total_busco_groups_searched = `Total BUSCO groups searched`
    )
}

clean_probe_data <- function(probe_data) {
  probe_data %>%
    dplyr::rename(
      pool_id = Pool,
      sample_id = `Sample ID`,
      original_percent = `Original%`,
      post_capture_percent = `Post-capture%`,
      enrichment = Enrichment,
      raw_mapping = `Raw mapping`,
      raw_total = `Raw total`
    ) %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      original_percent = readr::parse_number(original_percent),
      post_capture_percent = readr::parse_number(post_capture_percent),
      pool_level = dplyr::case_when(
        pool_id == "A" ~ "High",
        pool_id == "B" ~ "Medium",
        pool_id == "C" ~ "Low",
        TRUE ~ NA_character_
      ),
      pool_level = factor(pool_level, levels = c("Low", "Medium", "High"))
    )
}

load_reference_coverage_summary <- function() {
  readr::read_tsv(
    file.path(pilot_report_dir, "bowtie2_coverage_summary.tsv"),
    show_col_types = FALSE
  ) %>%
    dplyr::rename(
      sample = Sample,
      coverage_1x = `>=1x`,
      coverage_5x = `>=5x`,
      coverage_10x = `>=10x`,
      coverage_25x = `>=25x`,
      coverage_50x = `>=50x`
    )
}

load_clinical_qc_data <- function(include_pool_data = TRUE) {
  qcd_data <- readr::read_csv(
    file.path(
      pilot_report_dir,
      "2026-10-20_pilot-strains_NG-probe_QC_summary.csv"
    ),
    show_col_types = FALSE
  ) %>%
    clean_clinical_qc_summary()

  qcd_data$sample_id <- stringr::str_split(
    qcd_data$sample,
    "_",
    simplify = TRUE
  )[, 5]

  if (anyDuplicated(qcd_data$sample) ||
      anyDuplicated(qcd_data$sample_id) ||
      anyNA(qcd_data$sample_id) ||
      any(!nzchar(qcd_data$sample_id))) {
    stop("Clinical QC data require one sample ID per row.")
  }

  if (include_pool_data) {
    probe_data <- readr::read_csv(
      "./data/NG_probe/metadata/neisseria_enrichment.csv",
      show_col_types = FALSE,
      col_types = readr::cols(`Sample ID` = readr::col_character())
    ) %>%
      clean_probe_data()

    if (anyDuplicated(probe_data$sample_id) ||
        anyNA(probe_data$pool_id) ||
        anyNA(probe_data$pool_level) ||
        anyNA(probe_data$original_percent) ||
        anyNA(probe_data$post_capture_percent) ||
        any(!dplyr::between(probe_data$original_percent, 0, 100)) ||
        any(!dplyr::between(probe_data$post_capture_percent, 0, 100)) ||
        !setequal(qcd_data$sample_id, probe_data$sample_id)) {
      stop("Probe enrichment data are incomplete or do not match the clinical QC samples.")
    }

    qcd_data <- dplyr::left_join(qcd_data, probe_data, by = "sample_id")
  }

  qcd_data <- qcd_data %>%
    dplyr::mutate(
      qc_passing = dplyr::between(
        total_length,
        2200000 * 0.85,
        2200000 * 1.15
      ),
      ng_species = post_map_top_species_1 == "Neisseria gonorrhoeae"
    )

  qcd_data_passing <- qcd_data %>%
    dplyr::filter(qc_passing, ng_species)

  if (anyNA(qcd_data$qc_passing) ||
      anyNA(qcd_data$ng_species) ||
      any(!qcd_data_passing$qc_passing) ||
      any(!qcd_data_passing$ng_species)) {
    stop("Clinical QC flags could not be assigned consistently.")
  }

  clinical_qc_data <- list(
    qcd_data = qcd_data,
    qcd_data_passing = qcd_data_passing
  )

  if (include_pool_data) {
    sample_pool_data <- qcd_data %>%
      dplyr::select(sample, sample_id, pool_id, pool_level) %>%
      dplyr::distinct()

    reference_coverage_summary <- load_reference_coverage_summary()

    coverage_columns <- c(
      "coverage_1x", "coverage_5x", "coverage_10x",
      "coverage_25x", "coverage_50x"
    )
    coverage_values <- as.matrix(
      reference_coverage_summary[, coverage_columns]
    )

    if (anyDuplicated(reference_coverage_summary$sample) ||
        !setequal(reference_coverage_summary$sample, sample_pool_data$sample) ||
        anyNA(coverage_values) ||
        any(!is.finite(coverage_values)) ||
        any(coverage_values < 0 | coverage_values > 100) ||
        any(apply(coverage_values, 1, function(x) any(diff(x) > 0)))) {
      stop("Reference coverage data are invalid or do not match the clinical QC samples.")
    }

    reference_coverage_summary <- reference_coverage_summary %>%
      dplyr::left_join(sample_pool_data, by = "sample")

    clinical_qc_data$probe_data <- probe_data
    clinical_qc_data$sample_pool_data <- sample_pool_data
    clinical_qc_data$reference_coverage_summary <- reference_coverage_summary
  }

  clinical_qc_data
}
