# Build shared phylogeny and profile data

# Packages
packages <- c("tidyverse")

# Load packages
invisible(lapply(packages, library, character.only = TRUE, quietly = TRUE))

# Processed data functions
source("./lib/processed_data_functions.R")

# Representative isolate data
profiles_350 <- readRDS("./data/representative_isolates/source/selections/profiles_350.RDS")

# Phylogeny and 350-AA profile data
source("./scripts/phylogeny/03__phylogeny_profile_data_processing.R")

phylogeny_profile_data <- list(
  tree = tree,
  profiles_350 = profiles_350,
  ng_STs = ng_STs,
  sequence_types = sequence_types,
  tree_tip_data = tree_tip_data,
  profile_assignments = profile_assignments,
  genome_groups = genome_groups,
  tree_genome_groups = tree_genome_groups,
  major_sts = major_sts,
  selected_350 = selected_350,
  selected_genome_ids = selected_genome_ids,
  tree_350 = tree_350,
  phylogeny_track_data = phylogeny_track_data,
  phylogeny_tracks = phylogeny_tracks,
  cluster_track_levels = cluster_track_levels,
  profile_track_levels = profile_track_levels
)

# Profile assignments across thresholds
source("./scripts/phylogeny/04__threshold_profile_data_processing.R")

cluster_str <- readRDS(
  "./data/phylogeny/source/distinct_profiles/cluster_str.RDS"
)
cluster_stats <- readRDS(
  "./data/phylogeny/source/distinct_profiles/cluster_stats.RDS"
)
cluster_results_list <- readRDS(
  "./data/phylogeny/source/distinct_profiles/cluster_results_list.RDS"
)

profile_stats <- purrr::map_dfr(names(cluster_results_list), function(threshold) {
  collapsed_profiles <- cluster_results_list[[threshold]]$cluster_collapse
  tibble(
    threshold = as.numeric(threshold),
    profile_groups = nrow(collapsed_profiles),
    profile_clusters = sum(collapsed_profiles$count > 1),
    profile_singletons = sum(collapsed_profiles$count == 1)
  )
})

threshold_group_counts <- cluster_stats %>%
  left_join(profile_stats, by = "threshold") %>%
  pivot_longer(
    cols = c(
      total_cluster_singleton_ct, no_cluster, not_in_cluster_ct,
      profile_groups, profile_clusters, profile_singletons
    ),
    names_to = "metric",
    values_to = "count"
  ) %>%
  mutate(
    category = case_when(
      metric %in% c("total_cluster_singleton_ct", "profile_groups") ~ "Groups",
      metric %in% c("no_cluster", "profile_clusters") ~ "Clusters",
      TRUE ~ "Singletons"
    ),
    type = if_else(
      metric %in% c(
        "total_cluster_singleton_ct", "no_cluster", "not_in_cluster_ct"
      ),
      "AA-based clustering",
      "AA+Trait-based clustering"
    )
  )

profile_medians <- purrr::map_dfr(names(cluster_results_list), function(threshold) {
  collapsed_profiles <- cluster_results_list[[threshold]]$cluster_collapse
  tibble(
    threshold = as.numeric(threshold),
    profile_size_med = median(collapsed_profiles$count[collapsed_profiles$count > 1])
  )
})

threshold_group_medians <- cluster_stats %>%
  left_join(profile_medians, by = "threshold") %>%
  select(threshold, cluster_size_med, profile_size_med) %>%
  pivot_longer(
    cols = c(cluster_size_med, profile_size_med),
    names_to = "metric",
    values_to = "size"
  )

make_group_size_data <- function(data, type) {
  sizes <- apply(data, 2, function(x) as.vector(table(x)))
  purrr::imap_dfr(sizes, function(size, threshold) {
    tibble(
      threshold = as.numeric(threshold),
      size = size[size != 1],
      type = type
    )
  })
}

threshold_group_sizes <- bind_rows(
  make_group_size_data(profile_id, "AA+trait-based clustering"),
  make_group_size_data(cluster_str, "AA-based clustering")
)

threshold_group_size_medians <- threshold_group_sizes %>%
  group_by(threshold, type) %>%
  summarise(median_size = median(size, na.rm = TRUE), .groups = "drop")

threshold_provenance <- create_data_provenance(
  builder = "scripts/phylogeny/05__build_phylogeny_data.R",
  input_files = c(
    "data/phylogeny/source/distinct_profiles/profile_id.RDS",
    "data/phylogeny/source/distinct_profiles/cluster_str.RDS",
    "data/phylogeny/source/distinct_profiles/cluster_stats.RDS",
    "data/phylogeny/source/distinct_profiles/cluster_results_list.RDS"
  )
)

threshold_profile_data <- list(
  profile_id = profile_id,
  threshold_profiles = threshold_profiles,
  threshold_350_profiles = threshold_350_profiles,
  threshold_350_cluster_levels = threshold_350_cluster_levels,
  threshold_350_profile_levels = threshold_350_profile_levels,
  selected_threshold_350_profiles = selected_threshold_350_profiles,
  threshold_group_counts = threshold_group_counts,
  threshold_group_medians = threshold_group_medians,
  threshold_group_sizes = threshold_group_sizes,
  threshold_group_size_medians = threshold_group_size_medians,
  provenance = threshold_provenance
)

save_processed_data(
  threshold_profile_data,
  "./data/phylogeny/curated/threshold_profile_data.RDS"
)

# Ordered Pathogenwatch annotations for the selected-isolate tree
load_processed_data(
  "./data/pathogenwatch/curated/pathogenwatch_data.RDS",
  c(
    "pathogenwatch_isolate_data",
    "resistance_variables",
    "mutation_variables"
  )
)

if (
  anyNA(pathogenwatch_isolate_data$genome_id) ||
    anyDuplicated(pathogenwatch_isolate_data$genome_id) ||
    !all(
      c(resistance_variables, mutation_variables) %in%
        colnames(pathogenwatch_isolate_data)
    ) ||
    !all(tree_350$tip.label %in% pathogenwatch_isolate_data$genome_id)
) {
  stop("Pathogenwatch annotations must have unique genome IDs and expected columns")
}

selected_phylogeny_annotations <- tibble(
  genome_id = tree_350$tip.label,
  selected_tree_order = seq_along(tree_350$tip.label)
) %>%
  left_join(
    tree_genome_groups %>% select(genome_id, sequence_type, st_track),
    by = "genome_id"
  ) %>%
  left_join(
    pathogenwatch_isolate_data %>%
      select(
        genome_id,
        all_of(resistance_variables),
        all_of(mutation_variables)
      ),
    by = "genome_id"
  )

if (
  anyNA(selected_phylogeny_annotations$sequence_type) ||
    nrow(selected_phylogeny_annotations) != length(tree_350$tip.label)
) {
  stop("Selected phylogeny annotations do not match the selected tree")
}

phylogeny_profile_data$selected_phylogeny_annotations <-
  selected_phylogeny_annotations
phylogeny_profile_data$resistance_variables <- resistance_variables
phylogeny_profile_data$mutation_variables <- mutation_variables
phylogeny_provenance <- create_data_provenance(
  builder = "scripts/phylogeny/05__build_phylogeny_data.R",
  input_files = c(
    "data/phylogeny/source/cognac_nj.tre",
    "./data/representative_isolates/source/selections/profiles_350.RDS",
    "../../data/pathogenwatch_metadata/Neisseria_gonorrhoeae_mlst.csv",
    "data/pathogenwatch/curated/pathogenwatch_data.RDS"
  )
)
phylogeny_profile_data$provenance <- phylogeny_provenance

save_processed_data(
  phylogeny_profile_data,
  "./data/phylogeny/curated/phylogeny_profile_data.RDS"
)
