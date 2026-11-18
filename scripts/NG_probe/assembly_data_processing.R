# Clinical assembly and MLST data

load_pilot_assembly_statistics <- function() {
  assembly_files <- list.files(
    "./data/NG_probe/spades/bowtie2",
    pattern = "_contigs_l1000[.]fasta$",
    recursive = TRUE,
    full.names = TRUE
  )

  assembly_statistics <- purrr::map_dfr(assembly_files, function(path) {
    sample <- sub("_contigs_l1000[.]fasta$", "", basename(path))
    sample_id <- stringr::str_split(sample, "_", simplify = TRUE)[, 5]
    assembly <- Biostrings::readDNAStringSet(path)

    tibble::tibble(
      contig_id = names(assembly),
      length = BiocGenerics::width(assembly),
      sample = sample,
      sample_id = sample_id
    )
  })

  if (
    length(assembly_files) != 12 ||
      dplyr::n_distinct(assembly_statistics$sample) != 12 ||
      anyDuplicated(assembly_statistics$contig_id) ||
      anyNA(assembly_statistics) ||
      any(assembly_statistics$length < 1000)
  ) {
    stop("Pilot assembly statistics could not be validated.")
  }

  assembly_statistics
}

load_pilot_mlst_data <- function() {
  mlst_files <- list.files(
    "./data/NG_probe/mlst/bowtie2",
    pattern = "report[.]tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  mlst_data <- purrr::map_dfr(mlst_files, function(path) {
    readr::read_tsv(
      path,
      col_names = c(
        "filename", "scheme", "ST", "abcZ", "adK", "aroE",
        "fumC", "gdh", "pdhC", "pgm"
      ),
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    )
  }) %>%
    dplyr::mutate(
      sample = sub("_contigs_l1000[.]fasta$", "", basename(filename)),
      sample_id = stringr::str_split(sample, "_", simplify = TRUE)[, 5]
    )

  if (
    length(mlst_files) != 12 ||
      nrow(mlst_data) != 12 ||
      anyDuplicated(mlst_data$sample) ||
      anyNA(mlst_data$sample) ||
      anyNA(mlst_data$sample_id)
  ) {
    stop("Pilot MLST data could not be validated.")
  }

  mlst_data
}

load_pilot_assembly_data <- function() {
  list(
    assembly_statistics = load_pilot_assembly_statistics(),
    mlst_data = load_pilot_mlst_data()
  )
}
