# Profile assignments across thresholds

profile_id <- readRDS(
  "./data/phylogeny/source/distinct_profiles/profile_id.RDS"
)

# Check and clean profile assignments
if (
  !is.data.frame(profile_id) ||
    !"350" %in% colnames(profile_id) ||
    is.null(rownames(profile_id)) ||
    anyNA(rownames(profile_id)) ||
    any(rownames(profile_id) == "") ||
    anyDuplicated(rownames(profile_id))
) {
  stop("profile_id must contain threshold 350 and unique genome row names")
}

threshold_profiles <- profile_id %>%
  rownames_to_column("genome_id")

threshold_350_profiles <- threshold_profiles %>%
  transmute(
    genome_id,
    profile = as.character(`350`),
    aa_cluster = sub("\\..*$", "", profile)
  )

if (anyNA(threshold_350_profiles$profile)) {
  stop("Threshold-350 profiles cannot contain missing assignments")
}

# Use the shared tree order when phylogeny data have already been loaded
if (exists("tree_tip_data", inherits = FALSE)) {
  if (!setequal(tree_tip_data$genome_id, threshold_350_profiles$genome_id)) {
    stop("Tree tips and threshold-350 profile assignments do not match")
  }

  threshold_350_profiles <- tree_tip_data %>%
    select(genome_id) %>%
    left_join(threshold_350_profiles, by = "genome_id")
}

threshold_350_cluster_levels <- unique(
  threshold_350_profiles$aa_cluster
)
threshold_350_profile_levels <- unique(
  threshold_350_profiles$profile
)

if (exists("tree_350", inherits = FALSE)) {
  selected_threshold_350_profiles <- tibble(
    genome_id = tree_350$tip.label,
    selected_tree_order = seq_along(tree_350$tip.label)
  ) %>%
    left_join(threshold_350_profiles, by = "genome_id")

  if (anyNA(selected_threshold_350_profiles$profile)) {
    stop("Every selected tree tip must have a threshold-350 profile")
  }
}
