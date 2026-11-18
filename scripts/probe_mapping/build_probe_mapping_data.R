# Build shared probe mapping data

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Biostrings))

source("./lib/processed_data_functions.R")
source("./scripts/probe_mapping/probe_mapping_data_processing.R")

reference_blast_data <- readRDS(
  "./data/reference_genome/curated/reference_blast_data.RDS"
)
validate_processed_data(reference_blast_data)
probe_mapping_data <- load_probe_mapping_data(
  reference_blast_data$pathogenwatch_blast
)
probe_mapping_data$provenance <- create_data_provenance(
  builder = "scripts/probe_mapping/build_probe_mapping_data.R",
  input_files = c(
    unname(probe_mapping_input_paths),
    "data/reference_genome/curated/reference_blast_data.RDS"
  ),
  coverage_run_id =
    "240808-903-Stefan_Green-D10891NGNRH-Ngonorrhoeae_probe_mapping"
)

save_processed_data(
  probe_mapping_data,
  "./data/probe_mapping/curated/probe_mapping_data.RDS"
)

message("Probe mapping data created in data/probe_mapping/curated")
