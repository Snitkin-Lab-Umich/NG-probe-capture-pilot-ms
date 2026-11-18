# Chosen genomes

if (!exists("profiles_350", inherits = FALSE)) {
  stop("profiles_350 must be loaded from data/representative_isolates/source/selections/profiles_350.RDS")
}

if (
  !is.data.frame(profiles_350) ||
    !all(
      c("cl_prof_id", "cluster", "isolates", "iso_by_qual") %in%
        colnames(profiles_350)
    ) ||
    anyNA(profiles_350$iso_by_qual) ||
    anyDuplicated(profiles_350$iso_by_qual)
) {
  stop("profiles_350 must contain unique selected isolates and profile assignments")
}

chosen_genomes <- profiles_350$iso_by_qual %>% unname()

chosen_genome_data <- profiles_350 %>%
  dplyr::transmute(
    profile_id = as.character(cl_prof_id),
    cluster_id = cluster,
    genome_id = as.character(iso_by_qual)
  ) %>%
  dplyr::left_join(
    pathogenwatch_isolate_data %>%
      dplyr::select(pathogenwatch_genome_id, genome_id),
    by = "genome_id"
  ) %>%
  dplyr::select(
    genome_id,
    pathogenwatch_genome_id,
    profile_id,
    cluster_id
  )

if (
  anyNA(chosen_genome_data$pathogenwatch_genome_id) ||
  anyDuplicated(chosen_genome_data$genome_id) ||
    !setequal(chosen_genome_data$genome_id, chosen_genomes)
) {
  stop("Chosen genomes must map once to the canonical Pathogenwatch isolate data")
}

chosen_genome_ids <- chosen_genome_data$genome_id

pathogenwatch_isolate_data <- pathogenwatch_isolate_data %>%
  dplyr::mutate(is_chosen = genome_id %in% chosen_genome_ids)
