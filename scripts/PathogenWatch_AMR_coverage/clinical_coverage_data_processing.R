# Clinical supplement data

clinical_coverage_input_paths <- c(
  chromosome = paste0(
    "./data/reference_genome_coverage/coverage/curated_datasets/",
    "FA1090_genome_coverage_summary_statistics.RDS"
  ),
  position = paste0(
    "./data/PathogenWatch_AMR_coverage/coverage/curated_datasets/",
    "Pathogenwatch_gene_coverage_by_position_chromosome_normalized.RDS"
  ),
  locus = paste0(
    "./data/PathogenWatch_AMR_coverage/coverage/curated_datasets/",
    "Pathogenwatch_gene_coverage_summary_statistics_chromosome_normalized.RDS"
  ),
  breadth = paste0(
    "./data/PathogenWatch_AMR_coverage/coverage/curated_datasets/",
    "PathogenWatch_coverage_master_summary_file.RDS"
  )
)

load_clinical_coverage_sources <- function(all_isolates) {
  chromosome <- readRDS(clinical_coverage_input_paths[["chromosome"]])
  position <- readRDS(clinical_coverage_input_paths[["position"]])
  locus <- readRDS(clinical_coverage_input_paths[["locus"]])
  breadth <- readRDS(clinical_coverage_input_paths[["breadth"]])

  source_isolates <- list(
    chromosome = chromosome$isolate,
    position = position$isolate,
    locus = locus$isolate,
    breadth = breadth$isolate
  )
  source_sets_match <- vapply(
    source_isolates,
    function(x) setequal(unique(x), all_isolates),
    logical(1)
  )

  expected_breadth <- position %>%
    dplyr::group_by(isolate, locus = genotype) %>%
    dplyr::summarise(
      expected_1 = mean(Total_Depth >= 1) * 100,
      expected_5 = mean(Total_Depth >= 5) * 100,
      expected_10 = mean(Total_Depth >= 10) * 100,
      expected_15 = mean(Total_Depth >= 15) * 100,
      expected_20 = mean(Total_Depth >= 20) * 100,
      expected_25 = mean(Total_Depth >= 25) * 100,
      expected_50 = mean(Total_Depth >= 50) * 100,
      .groups = "drop"
    )

  breadth_check <- breadth %>%
    dplyr::left_join(
      expected_breadth,
      by = c("isolate", "locus"),
      relationship = "one-to-one"
    )
  thresholds <- c(1L, 5L, 10L, 15L, 20L, 25L, 50L)
  breadth_matches <- vapply(
    thresholds,
    function(depth_threshold) {
      observed <- breadth_check[[paste0("percent_coverage_", depth_threshold)]]
      expected <- breadth_check[[paste0("expected_", depth_threshold)]]
      isTRUE(all.equal(observed, expected, tolerance = 1e-12))
    },
    logical(1)
  )

  chromosome_denominator <- stats::setNames(
    chromosome$median_non_zero,
    chromosome$isolate
  )
  position_denominator_matches <- isTRUE(all.equal(
    position$chromosome_median_non_zero_depth,
    unname(chromosome_denominator[position$isolate]),
    tolerance = 1e-12
  ))
  locus_denominator_matches <- isTRUE(all.equal(
    locus$chromosome_median_non_zero_depth,
    unname(chromosome_denominator[locus$isolate]),
    tolerance = 1e-12
  ))
  position_normalization_matches <- isTRUE(all.equal(
    position$positional_chromosome_median_non_zero_normalized_depth,
    position$Total_Depth / position$chromosome_median_non_zero_depth,
    tolerance = 1e-12
  ))
  locus_normalization_matches <- isTRUE(all.equal(
    locus$loci_chromosome_median_non_zero_normalized_depth,
    locus$loci_median_non_zero_depth /
      locus$chromosome_median_non_zero_depth,
    tolerance = 1e-12
  ))

  breadth_columns <- paste0("percent_coverage_", thresholds)
  breadth_matrix <- as.matrix(breadth[, breadth_columns])

  if (
    length(all_isolates) != 12 ||
      !all(source_sets_match) ||
      nrow(chromosome) != length(all_isolates) ||
      nrow(locus) != length(all_isolates) * 35 ||
      nrow(breadth) != length(all_isolates) * 35 ||
      dplyr::n_distinct(position$genotype) != 35 ||
      anyDuplicated(chromosome$isolate) ||
      anyDuplicated(position[c("isolate", "genotype", "position")]) ||
      anyDuplicated(locus[c("isolate", "genotype")]) ||
      anyDuplicated(breadth[c("isolate", "locus")]) ||
      anyNA(position[c(
        "Total_Depth",
        "positional_chromosome_median_non_zero_normalized_depth"
      )]) ||
      anyNA(locus[c(
        "loci_median_non_zero_depth",
        "loci_chromosome_median_non_zero_normalized_depth"
      )]) ||
      anyNA(breadth_matrix) ||
      any(!is.finite(breadth_matrix)) ||
      any(breadth_matrix < 0 | breadth_matrix > 100) ||
      any(apply(breadth_matrix, 1, function(x) any(diff(x) > 0))) ||
      !all(breadth_matches) ||
      !position_denominator_matches ||
      !locus_denominator_matches ||
      !position_normalization_matches ||
      !locus_normalization_matches
  ) {
    stop("Curated clinical coverage sources could not be validated.")
  }

  list(
    chromosome = chromosome,
    position = position,
    locus = locus,
    breadth = breadth
  )
}

prepare_clinical_pathogenwatch_coverage <- function(
  clinical_coverage_sources,
  isolates
) {
  chromosome_coverage <- clinical_coverage_sources$chromosome %>%
    dplyr::filter(isolate %in% isolates) %>%
    dplyr::transmute(
      isolate,
      chromosome_median_depth = median,
      chromosome_median_nonzero_depth = median_non_zero
    ) %>%
    dplyr::arrange(match(isolate, isolates))

  coverage_qc <- clinical_coverage_sources$position %>%
    dplyr::filter(isolate %in% isolates) %>%
    dplyr::transmute(
      isolate,
      locus = genotype,
      locus_position = as.integer(position),
      depth = Total_Depth,
      normalized_depth =
        positional_chromosome_median_non_zero_normalized_depth
    )

  clinical_locus_coverage <- clinical_coverage_sources$locus %>%
    dplyr::filter(isolate %in% isolates) %>%
    dplyr::transmute(
      isolate,
      locus = genotype,
      median_nonzero_depth = loci_median_non_zero_depth,
      normalized_depth =
        loci_chromosome_median_non_zero_normalized_depth
    )

  if (nrow(coverage_qc) == 0 ||
      anyNA(coverage_qc$locus_position) ||
      anyDuplicated(coverage_qc[c("isolate", "locus", "locus_position")]) ||
      any(!is.finite(coverage_qc$depth)) ||
      any(coverage_qc$depth < 0) ||
      any(!is.finite(coverage_qc$normalized_depth)) ||
      any(coverage_qc$normalized_depth < 0) ||
      !setequal(coverage_qc$isolate, isolates) ||
      anyDuplicated(clinical_locus_coverage[c("isolate", "locus")]) ||
      any(!is.finite(clinical_locus_coverage$median_nonzero_depth)) ||
      any(clinical_locus_coverage$median_nonzero_depth < 0) ||
      any(!is.finite(clinical_locus_coverage$normalized_depth)) ||
      any(clinical_locus_coverage$normalized_depth < 0) ||
      !setequal(clinical_locus_coverage$isolate, isolates) ||
      !setequal(clinical_locus_coverage$locus, coverage_qc$locus)) {
    stop("Clinical Pathogenwatch coverage data could not be validated.")
  }

  list(
    chromosome_coverage = chromosome_coverage,
    coverage_qc = coverage_qc,
    clinical_locus_coverage = clinical_locus_coverage
  )
}

load_coverage_breadth_data <- function(clinical_coverage_sources, isolates) {
  coverage_breadth_data <- clinical_coverage_sources$breadth %>%
    dplyr::filter(isolate %in% isolates) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("percent_coverage_"),
      names_to = "depth_threshold",
      names_prefix = "percent_coverage_",
      values_to = "percent_coverage"
    ) %>%
    dplyr::transmute(
      isolate,
      locus,
      depth_threshold = as.integer(depth_threshold),
      percent_coverage
    )

  expected_depth_thresholds <- c(1L, 5L, 10L, 15L, 20L, 25L, 50L)
  threshold_counts <- coverage_breadth_data %>%
    dplyr::count(isolate, locus, name = "num_depth_thresholds")
  if (nrow(coverage_breadth_data) == 0 ||
      anyDuplicated(coverage_breadth_data[c("isolate", "locus", "depth_threshold")]) ||
      !setequal(coverage_breadth_data$isolate, isolates) ||
      !setequal(coverage_breadth_data$depth_threshold, expected_depth_thresholds) ||
      any(threshold_counts$num_depth_thresholds != length(expected_depth_thresholds)) ||
      anyNA(coverage_breadth_data$percent_coverage) ||
      any(coverage_breadth_data$percent_coverage < 0) ||
      any(coverage_breadth_data$percent_coverage > 100)) {
    stop("Clinical coverage breadth data could not be validated.")
  }

  coverage_breadth_data
}

prepare_isolate_metadata <- function(qcd_data_passing, isolates) {
  isolate_metadata <- qcd_data_passing %>%
    dplyr::filter(sample %in% isolates) %>%
    dplyr::transmute(
      isolate = sample,
      sample_id,
      pool_id,
      pool_level
    ) %>%
    dplyr::arrange(match(isolate, isolates))

  if (anyDuplicated(isolate_metadata$isolate) ||
      !identical(isolate_metadata$isolate, isolates) ||
      anyNA(isolate_metadata)) {
    stop("Clinical isolate metadata could not be validated.")
  }

  isolate_metadata
}

prepare_locus_metadata <- function(coverage_qc) {
  locus_metadata <- coverage_qc %>%
    dplyr::group_by(locus) %>%
    dplyr::summarise(locus_length = max(locus_position), .groups = "drop")

  observed_lengths <- coverage_qc %>%
    dplyr::count(isolate, locus, name = "observed_length") %>%
    dplyr::left_join(locus_metadata, by = "locus", relationship = "many-to-one")

  if (anyDuplicated(locus_metadata$locus) ||
      anyNA(locus_metadata$locus_length) ||
      any(locus_metadata$locus_length <= 0) ||
      any(observed_lengths$observed_length != observed_lengths$locus_length)) {
    stop("Clinical Pathogenwatch locus metadata could not be validated.")
  }

  locus_metadata
}

prepare_clinical_supplement_data <- function(
  qcd_data_passing,
  isolates,
  coverage_qc,
  coverage_breadth_data
) {
  list(
    coverage_breadth_data = coverage_breadth_data,
    isolate_metadata = prepare_isolate_metadata(qcd_data_passing, isolates),
    locus_metadata = prepare_locus_metadata(coverage_qc)
  )
}
