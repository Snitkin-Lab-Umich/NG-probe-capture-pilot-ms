# Build shared Pathogenwatch data

# Packages
packages <- c("tidyverse")
invisible(lapply(packages, library, character.only = TRUE, quietly = TRUE))

# Data processing
source("./lib/processed_data_functions.R")
source("./scripts/pathogenwatch/01__pathogenwatch_summary_functions.R")

profiles_350 <- readRDS("./data/representative_isolates/source/selections/profiles_350.RDS")

source("./scripts/pathogenwatch/02__pathogenwatch_data_processing.R")
source("./scripts/pathogenwatch/03__chosen_genomes_data_processing.R")

pathogenwatch_provenance <- create_data_provenance(
  builder = "scripts/pathogenwatch/04__build_pathogenwatch_data.R",
  input_files = c(
    "../../data/pathogenwatch_metadata/Neisseria_gonorrhoeae_amrsearch.csv",
    "../../data/ng_pathogenwatch_files/ng_snps.RDS",
    "../../data/pathogenwatch_metadata/Neisseria_gonorrhoeae_mlst.csv",
    "../../data/pathogenwatch_metadata/Neisseria_gonorrhoeae_metadata.csv",
    "data/phylogeny/source/cognac_nj.tre",
    "./data/representative_isolates/source/selections/profiles_350.RDS"
  )
)

# Shared figure data
pathogenwatch_data <- list(
  ng_amrsearch = ng_amrsearch,
  ng_STs = ng_STs,
  ng_sample_metadata = ng_sample_metadata,
  ng_meta = ng_meta,
  collapsed_snps_binary = collapsed_snps_binary,
  ng_meta_gene_grouped = ng_meta_gene_grouped,
  chosen_genomes = chosen_genomes,
  pathogenwatch_isolate_data = pathogenwatch_isolate_data,
  resistance_variables = resistance_variables,
  mutation_variables = mutation_variables,
  chosen_genome_data = chosen_genome_data,
  chosen_genome_ids = chosen_genome_ids,
  provenance = pathogenwatch_provenance
)

save_processed_data(
  pathogenwatch_data,
  "./data/pathogenwatch/curated/pathogenwatch_data.RDS"
)
