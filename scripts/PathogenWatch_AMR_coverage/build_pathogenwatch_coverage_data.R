# Build processed Pathogenwatch coverage data

suppressPackageStartupMessages(library(tidyverse))

source("./lib/processed_data_functions.R")
source("./scripts/PathogenWatch_AMR_coverage/clinical_mutation_functions.R")
source("./scripts/PathogenWatch_AMR_coverage/clinical_coverage_data_processing.R")

pilot_run_id <- "2026-10-20_pilot-strains_NG-probe"
coverage_run_id <- "2026-10-20_pilot-strains_NG-probe_coverage"

pilot_qc_data <- readRDS("./data/NG_probe/curated/pilot_qc_data.RDS")
validate_processed_data(pilot_qc_data)
all_isolates <- sort(pilot_qc_data$qcd_data$sample)
isolates <- sort(pilot_qc_data$qcd_data_passing$sample)

clinical_coverage_sources <- load_clinical_coverage_sources(all_isolates)
clinical_coverage_data <- prepare_clinical_pathogenwatch_coverage(
  clinical_coverage_sources,
  isolates
)
coverage_breadth_data <- load_coverage_breadth_data(
  clinical_coverage_sources,
  isolates
)
clinical_supplement_data <- prepare_clinical_supplement_data(
  pilot_qc_data$qcd_data_passing,
  isolates,
  clinical_coverage_data$coverage_qc,
  coverage_breadth_data
)

if (
  length(isolates) != 9 ||
    !setequal(clinical_coverage_data$coverage_qc$isolate, isolates) ||
    !setequal(clinical_supplement_data$coverage_breadth_data$isolate, isolates)
) {
  stop("Clinical coverage does not match the active QC-passing isolates.")
}

pathogenwatch_coverage_data <- c(
  list(isolates = isolates),
  clinical_coverage_data,
  clinical_supplement_data,
  list(
    provenance = create_data_provenance(
      builder = paste0(
        "scripts/PathogenWatch_AMR_coverage/",
        "build_pathogenwatch_coverage_data.R"
      ),
      input_files = c(
        unname(clinical_coverage_input_paths),
        "data/NG_probe/curated/pilot_qc_data.RDS"
      ),
      pilot_run_id = pilot_run_id,
      coverage_run_id = coverage_run_id,
      coverage_status = "current"
    )
  )
)
save_processed_data(
  pathogenwatch_coverage_data,
  paste0(
    "./data/PathogenWatch_AMR_coverage/curated/",
    "pathogenwatch_coverage_data.RDS"
  )
)

config_file <- "./data/PathogenWatch_AMR_coverage/database/485.toml"
config <- toml::read_toml(config_file)
pathogenwatch_mutation_data <- prepare_pathogenwatch_mutation_data(config)
pathogenwatch_mutation_data$provenance <- create_data_provenance(
  builder = paste0(
    "scripts/PathogenWatch_AMR_coverage/",
    "build_pathogenwatch_coverage_data.R"
  ),
  input_files = config_file,
  coverage_run_id = "Pathogenwatch_485_database"
)
save_processed_data(
  pathogenwatch_mutation_data,
  paste0(
    "./data/PathogenWatch_AMR_coverage/curated/",
    "pathogenwatch_mutation_data.RDS"
  )
)

message(
  "Pathogenwatch processed data created in ",
  "data/PathogenWatch_AMR_coverage/curated"
)
