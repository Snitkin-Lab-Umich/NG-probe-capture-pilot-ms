# Build shared reference BLAST data

# Packages
packages <- c("tidyverse")
invisible(lapply(packages, library, character.only = TRUE, quietly = TRUE))

# Processed data functions
source("./lib/processed_data_functions.R")

# Pathogenwatch reference-gene hits against the ASM684v1 reference genome
pathogenwatch_blast <- read_table(
  paste0(
    "./data/reference_genome/PathogenWatch_gene_detection/",
    "GCF_000006845.1_ASM684v1_genomic_blast_clean.tsv"
  ),
  show_col_types = FALSE,
  progress = FALSE
)

reference_blast_data <- list(
  pathogenwatch_blast = pathogenwatch_blast,
  provenance = create_data_provenance(
    builder = "scripts/reference_genome/build_reference_blast_data.R",
    input_files = paste0(
      "data/reference_genome/PathogenWatch_gene_detection/",
      "GCF_000006845.1_ASM684v1_genomic_blast_clean.tsv"
    )
  )
)

save_processed_data(
  reference_blast_data,
  "./data/reference_genome/curated/reference_blast_data.RDS"
)
