# Build processed AMRSearch data

suppressPackageStartupMessages(library(tidyverse))

source("./lib/processed_data_functions.R")
source("./scripts/amrsearch_10_20_2026_probes/amrsearch_parsing_functions.R")
source("./scripts/amrsearch_10_20_2026_probes/amrsearch_data_processing.R")

amrsearch_files <- list.files(
  "./data/amrsearch_10_20_2026_probes",
  pattern = "[.]jsn$",
  full.names = TRUE
)
results <- lapply(amrsearch_files, get_amrsearch_data)

resistance_profiles <- results %>%
  lapply(function(x) x$resistance_profile) %>%
  dplyr::bind_rows()
amrsearch_gene_data_raw <- results %>%
  lapply(function(x) x$genotype_profile) %>%
  dplyr::bind_rows()

gene_profile <- amrsearch_gene_data_raw %>%
  dplyr::filter(entry == "gene")
mutation_profile <- amrsearch_gene_data_raw %>%
  dplyr::filter(entry == "variant")

resistance_wide <- resistance_profiles %>%
  dplyr::select(isolate, antibiotic, state) %>%
  tidyr::pivot_wider(names_from = antibiotic, values_from = state) %>%
  dplyr::select(isolate, sort(colnames(.)[-1]))
gene_wide <- gene_profile %>%
  dplyr::select(isolate, refId) %>%
  dplyr::distinct() %>%
  dplyr::mutate(present = 1) %>%
  tidyr::pivot_wider(
    names_from = refId,
    values_from = present,
    values_fill = 0
  ) %>%
  dplyr::select(isolate, sort(colnames(.)[-1]))
mutation_wide <- mutation_profile %>%
  dplyr::select(isolate, genotype_name) %>%
  dplyr::distinct() %>%
  dplyr::mutate(present = 1) %>%
  tidyr::pivot_wider(
    names_from = genotype_name,
    values_from = present,
    values_fill = 0
  ) %>%
  dplyr::select(isolate, sort(colnames(.)[-1]))
amrsearch_matrix_raw <- resistance_wide %>%
  dplyr::left_join(gene_wide, by = "isolate") %>%
  dplyr::left_join(mutation_wide, by = "isolate")

if (
  length(amrsearch_files) != 12 ||
    !identical(dim(amrsearch_gene_data_raw), c(226L, 15L)) ||
    !identical(dim(amrsearch_matrix_raw), c(12L, 47L)) ||
    anyDuplicated(amrsearch_matrix_raw$isolate)
) {
  stop("The 2026-10-20 AMRSearch reports do not match the expected data.")
}

pilot_qc_data <- readRDS("./data/NG_probe/curated/pilot_qc_data.RDS")
pilot_assembly_data <- readRDS(
  "./data/NG_probe/curated/pilot_assembly_data.RDS"
)
validate_processed_data(pilot_qc_data)
validate_processed_data(pilot_assembly_data)

amrsearch_data <- prepare_amrsearch_shared_data(
  amrsearch_gene_data = amrsearch_gene_data_raw,
  amrsearch_matrix = amrsearch_matrix_raw,
  qcd_data = pilot_qc_data$qcd_data,
  qcd_data_passing = pilot_qc_data$qcd_data_passing,
  assembly_statistics = pilot_assembly_data$assembly_statistics,
  mlst_data = pilot_assembly_data$mlst_data
)
amrsearch_data$gene_calls <- amrsearch_data$gene_calls %>%
  dplyr::rename(locus = ref_id)
amrsearch_data$provenance <- create_data_provenance(
  builder = "scripts/amrsearch_10_20_2026_probes/build_amrsearch_data.R",
  input_files = c(
    "data/amrsearch_10_20_2026_probes/*.jsn",
    "data/NG_probe/curated/pilot_qc_data.RDS",
    "data/NG_probe/curated/pilot_assembly_data.RDS"
  ),
  pilot_run_id = "2026-10-20_pilot-strains_NG-probe"
)

save_processed_data(
  amrsearch_data,
  "./data/amrsearch_10_20_2026_probes/curated/amrsearch_data.RDS"
)

message(
  "AMRSearch processed data created in ",
  "data/amrsearch_10_20_2026_probes/curated"
)
