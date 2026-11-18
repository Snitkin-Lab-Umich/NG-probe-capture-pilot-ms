# Build processed NG probe pilot data

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))

source("./lib/processed_data_functions.R")
source("./scripts/NG_probe/clinical_qc_data_processing.R")
source("./scripts/NG_probe/assembly_data_processing.R")
source("./scripts/NG_probe/figure_3_species_data_processing.R")
source("./scripts/NG_probe/kraken_report_data_processing.R")

pilot_run_id <- "2026-10-20_pilot-strains_NG-probe"
pilot_report <- file.path(
  pilot_report_dir,
  "2026-10-20_pilot-strains_NG-probe_QC_summary.csv"
)
pilot_coverage <- file.path(pilot_report_dir, "bowtie2_coverage_summary.tsv")
probe_metadata <- "./data/NG_probe/metadata/neisseria_enrichment.csv"

pilot_qc_data <- load_clinical_qc_data()
if (
  nrow(pilot_qc_data$qcd_data) != 12 ||
    dplyr::n_distinct(pilot_qc_data$qcd_data$sample) != 12 ||
    nrow(pilot_qc_data$qcd_data_passing) != 9
) {
  stop("The active pilot must contain 12 samples and 9 QC-passing isolates.")
}
pilot_qc_data$provenance <- create_data_provenance(
  builder = "scripts/NG_probe/build_ng_probe_data.R",
  input_files = c(pilot_report, pilot_coverage, probe_metadata),
  pilot_run_id = pilot_run_id
)
save_processed_data(
  pilot_qc_data,
  "./data/NG_probe/curated/pilot_qc_data.RDS"
)

pilot_assembly_data <- load_pilot_assembly_data()
if (
  !setequal(
    pilot_assembly_data$mlst_data$sample,
    pilot_qc_data$qcd_data$sample
  ) ||
    !setequal(
      pilot_assembly_data$assembly_statistics$sample,
      pilot_qc_data$qcd_data$sample
    )
) {
  stop("Pilot assembly and MLST samples do not match the pilot QC report.")
}
pilot_assembly_data$provenance <- create_data_provenance(
  builder = "scripts/NG_probe/build_ng_probe_data.R",
  input_files = c(
    "data/NG_probe/spades/bowtie2/*/*_contigs_l1000.fasta",
    "data/NG_probe/mlst/bowtie2/*/report.tsv"
  ),
  pilot_run_id = pilot_run_id
)
save_processed_data(
  pilot_assembly_data,
  "./data/NG_probe/curated/pilot_assembly_data.RDS"
)

figure_3_species_data <- prepare_figure_3_species_data(
  pilot_qc_data$qcd_data,
  pilot_qc_data$sample_pool_data
)
figure_3_species_data$provenance <- create_data_provenance(
  builder = "scripts/NG_probe/build_ng_probe_data.R",
  input_files = pilot_report,
  pilot_run_id = pilot_run_id
)
save_processed_data(
  figure_3_species_data,
  "./data/NG_probe/curated/figure_3_species_data.RDS"
)

pre_binning_kraken_data <- load_pre_binning_kraken_data(
  pilot_qc_data$sample_pool_data
)
if (!setequal(
  pre_binning_kraken_data$pre_binning_kraken_summary$sample,
  pilot_qc_data$qcd_data$sample
)) {
  stop("Pre-binning Kraken reports do not match the active pilot samples.")
}
pre_binning_kraken_data$provenance <- create_data_provenance(
  builder = "scripts/NG_probe/build_ng_probe_data.R",
  input_files = "data/NG_probe/kraken2/pre-binning/*/*.tsv",
  pilot_run_id = pilot_run_id
)
save_processed_data(
  pre_binning_kraken_data,
  "./data/NG_probe/curated/pre_binning_kraken_data.RDS"
)

message("NG probe processed data created in data/NG_probe/curated")
