# Clinical AMRSearch data

clean_amrsearch_gene_data <- function(amrsearch_gene_data) {
  amrsearch_gene_data %>%
    dplyr::rename(
      ref_id = refId,
      ref_start = refStart,
      contig_id = queryId,
      query_start = queryStart,
      percent_identity = pid,
      e_value = evalue,
      ref_stop = refStop,
      query_stop = queryStop,
      ref_length = refLength,
      variant_name = name
    ) %>%
    dplyr::mutate(
      contig_start = pmin(query_start, query_stop),
      contig_stop = pmax(query_start, query_stop)
    )
}

prepare_amrsearch_shared_data <- function(
  amrsearch_gene_data,
  amrsearch_matrix,
  qcd_data,
  qcd_data_passing,
  assembly_statistics,
  mlst_data
) {
  isolates <- sort(qcd_data_passing$sample)

  amrsearch_gene_data <- clean_amrsearch_gene_data(amrsearch_gene_data)

  amrsearch_matrix <- amrsearch_matrix %>%
    dplyr::left_join(
      qcd_data %>% dplyr::select(isolate = sample, pool_level),
      by = "isolate",
      relationship = "many-to-one"
    ) %>%
    dplyr::left_join(
      mlst_data %>% dplyr::select(isolate = sample, ST),
      by = "isolate",
      relationship = "many-to-one"
    ) %>%
    dplyr::filter(isolate %in% isolates)

  amrsearch_contig_data <- amrsearch_gene_data %>%
    dplyr::filter(isolate %in% isolates) %>%
    dplyr::group_by(ref_id) %>%
    dplyr::filter(entry == "gene" | !any(entry == "gene", na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      assembly_statistics,
      by = "contig_id",
      relationship = "many-to-one"
    ) %>%
    dplyr::mutate(
      dist_from_start = contig_start - 1,
      dist_from_end = length - contig_stop,
      min_dist_to_edge = pmin(dist_from_start, dist_from_end)
    ) %>%
    dplyr::group_by(genotype_name) %>%
    dplyr::mutate(
      label = paste0(gsub("_", " ", genotype_name), " (n=", dplyr::n(), ")")
    ) %>%
    dplyr::ungroup()

  amrsearch_genotype_data <- amrsearch_gene_data %>%
    dplyr::select(genotype_name, entry) %>%
    dplyr::distinct()

  gene_calls <- amrsearch_gene_data %>%
    dplyr::filter(entry == "gene", isolate %in% isolates) %>%
    dplyr::distinct(isolate, ref_id)

  genotype_columns <- intersect(
    names(amrsearch_matrix),
    amrsearch_genotype_data$genotype_name
  )

  contig_sample_ids <- stringr::str_split(
    amrsearch_contig_data$isolate,
    "_",
    simplify = TRUE
  )[, 5]

  if (anyDuplicated(assembly_statistics$contig_id) ||
      anyDuplicated(amrsearch_matrix$isolate) ||
      !setequal(amrsearch_matrix$isolate, isolates) ||
      anyNA(amrsearch_matrix$pool_level) ||
      anyNA(amrsearch_matrix$ST) ||
      anyNA(amrsearch_gene_data$genotype_name) ||
      any(!amrsearch_gene_data$entry %in% c("gene", "variant")) ||
      anyDuplicated(amrsearch_genotype_data$genotype_name) ||
      length(genotype_columns) == 0 ||
      anyNA(amrsearch_matrix[, genotype_columns]) ||
      any(!unlist(amrsearch_matrix[, genotype_columns]) %in% c(0, 1)) ||
      nrow(amrsearch_contig_data) == 0 ||
      anyNA(amrsearch_contig_data$length) ||
      any(amrsearch_contig_data$length <= 0) ||
      any(amrsearch_contig_data$contig_start < 1) ||
      any(amrsearch_contig_data$contig_stop > amrsearch_contig_data$length) ||
      any(as.character(amrsearch_contig_data$sample_id) != contig_sample_ids)) {
    stop("Clinical AMRSearch data could not be matched or validated.")
  }

  list(
    amrsearch_matrix = amrsearch_matrix,
    amrsearch_gene_data = amrsearch_gene_data,
    amrsearch_contig_data = amrsearch_contig_data,
    amrsearch_genotype_data = amrsearch_genotype_data,
    gene_calls = gene_calls
  )
}
