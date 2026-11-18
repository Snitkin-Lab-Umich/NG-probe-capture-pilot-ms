# Figure 3 species data

species_long_tbl <- function(data, stage_prefix, stage_label) {
  species_cols <- paste0(stage_prefix, "_top_species_", 1:3)
  abund_cols <- paste0(stage_prefix, "_relative_abundance_", 1:3)

  data %>%
    select(sample, all_of(species_cols), all_of(abund_cols)) %>%
    pivot_longer(
      cols = c(all_of(species_cols), all_of(abund_cols)),
      names_to = c(".value", "rank"),
      names_pattern = paste0(stage_prefix, "_(top_species|relative_abundance)_([123])")
    ) %>%
    dplyr::rename(species = top_species, rel_abund = relative_abundance) %>%
    mutate(
      rank = as.integer(rank),
      rel_abund = as.numeric(rel_abund),
      species = as.character(species),
      stage = stage_label
    ) %>%
    select(sample, stage, rank, species, rel_abund)
}

prepare_figure_3_species_data <- function(qcd_data, sample_pool_data) {
  top3_all <- bind_rows(
    species_long_tbl(qcd_data, "pre_map", "Initial"),
    species_long_tbl(qcd_data, "post_map", "Post read mapping")
  ) %>%
    left_join(
      sample_pool_data %>% select(sample, pool_level),
      by = "sample"
    )

  ng_abund_initial <- top3_all %>%
    filter(
      stage == "Initial",
      species == "Neisseria gonorrhoeae"
    ) %>%
    group_by(sample) %>%
    arrange(rank, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      sample,
      stage,
      relative_abundance = rel_abund,
      top_id = rank,
      pool_level
    )

  list(
    top3_all = top3_all,
    ng_abund_initial = ng_abund_initial
  )
}
