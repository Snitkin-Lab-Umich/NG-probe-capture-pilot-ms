# Phylogeny and 350-AA profile data

# Load data
if (!exists("tree", inherits = FALSE)) {
  tree <- ape::read.tree(
    "./data/phylogeny/source/cognac_nj.tre"
  )
}

if (!exists("profiles_350", inherits = FALSE)) {
  stop("profiles_350 must be loaded from data/representative_isolates/source/selections/profiles_350.RDS")
}

if (!exists("ng_STs", inherits = FALSE)) {
  ng_STs <- read.csv(
    "../../data/pathogenwatch_metadata/Neisseria_gonorrhoeae_mlst.csv"
  )
}

# Check shared inputs
if (
  !inherits(tree, "phylo") ||
    length(tree$tip.label) == 0 ||
    anyNA(tree$tip.label) ||
    anyDuplicated(tree$tip.label)
) {
  stop("tree must be a phylo object with unique, non-missing tip labels")
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

if (
  !is.data.frame(ng_STs) ||
    !all(c("Genome.Name", "ST") %in% colnames(ng_STs)) ||
    anyNA(ng_STs$Genome.Name) ||
    anyDuplicated(ng_STs$Genome.Name)
) {
  stop("ng_STs must contain one row per non-missing Genome.Name and an ST column")
}

# Generate canonical sequence type data
sequence_types <- ng_STs %>%
  transmute(
    genome_id = as.character(Genome.Name),
    sequence_type = as.character(ST)
  )

# Generate genome-level profile data
genome_groups <- profiles_350 %>%
  select(
    profile = cl_prof_id,
    aa_cluster = cluster,
    isolates,
    selected_genome = iso_by_qual
  ) %>%
  separate_rows(isolates, sep = ",") %>%
  transmute(
    genome_id = str_trim(isolates),
    profile = as.character(profile),
    aa_cluster = as.character(aa_cluster),
    selected = genome_id == selected_genome
  ) %>%
  left_join(sequence_types, by = "genome_id")

# Check data
if (
  nrow(genome_groups) != length(tree$tip.label) ||
    n_distinct(genome_groups$genome_id) != length(tree$tip.label) ||
    !setequal(tree$tip.label, genome_groups$genome_id)
) {
  stop("Tree tips and 350-AA profile assignments do not match")
}
if (sum(genome_groups$selected) != nrow(profiles_350)) {
  stop("Selected-isolate membership does not match profiles_350")
}

# Identify major sequence types
major_sts <- genome_groups %>%
  count(sequence_type, name = "original_genomes") %>%
  arrange(desc(original_genomes), sequence_type) %>%
  filter(original_genomes >= 100) %>%
  pull(sequence_type)

genome_groups <- genome_groups %>%
  mutate(
    st_class = if_else(
      sequence_type %in% major_sts,
      "Major ST",
      "Minor ST"
    ),
    st_track = if_else(
      sequence_type %in% major_sts,
      sequence_type,
      "Other ST"
    )
  )

# Retain the established Figure 1 data and create canonical ordered data
profile_assignments <- genome_groups
genome_groups <- profile_assignments %>%
  transmute(
    genome_id,
    profile,
    aa_cluster,
    selected,
    ST = sequence_type,
    st_class,
    ST_track = st_track
  )

tree_tip_data <- tibble(
  genome_id = tree$tip.label,
  tree_order = seq_along(tree$tip.label)
)

tree_genome_groups <- tree_tip_data %>%
  left_join(profile_assignments, by = "genome_id")

if (anyNA(tree_genome_groups$profile)) {
  stop("Every tree tip must have a 350-AA profile assignment")
}

# Generate selected isolate data and tree
selected_350 <- profiles_350 %>%
  transmute(
    genome_id = as.character(iso_by_qual),
    profile = as.character(cl_prof_id),
    aa_cluster = as.character(cluster)
  )

selected_genome_ids <- selected_350$genome_id
tree_350 <- ape::keep.tip(tree, selected_genome_ids)

# Generate consistently numbered annotation tracks
cluster_number_lookup <- tree_genome_groups %>%
  filter(sequence_type %in% major_sts) %>%
  distinct(sequence_type, aa_cluster) %>%
  group_by(sequence_type) %>%
  mutate(cluster_number = row_number()) %>%
  ungroup()

profile_number_lookup <- tree_genome_groups %>%
  filter(sequence_type %in% major_sts) %>%
  distinct(sequence_type, aa_cluster, profile) %>%
  group_by(sequence_type, aa_cluster) %>%
  mutate(profile_number = row_number()) %>%
  ungroup()

cluster_track_levels <- paste0(
  "Cluster ",
  seq_len(max(cluster_number_lookup$cluster_number))
)
profile_track_levels <- paste0(
  "Profile ",
  seq_len(max(profile_number_lookup$profile_number))
)

phylogeny_track_data <- tree_genome_groups %>%
  left_join(
    cluster_number_lookup,
    by = c("sequence_type", "aa_cluster")
  ) %>%
  left_join(
    profile_number_lookup,
    by = c("sequence_type", "aa_cluster", "profile")
  ) %>%
  mutate(
    cluster_track = if_else(
      sequence_type %in% major_sts,
      paste0("Cluster ", cluster_number),
      "Cluster gap"
    ),
    profile_track = if_else(
      sequence_type %in% major_sts,
      paste0("Profile ", profile_number),
      "Profile gap"
    )
  ) %>%
  select(
    genome_id,
    tree_order,
    profile,
    aa_cluster,
    selected,
    sequence_type,
    st_class,
    st_track,
    cluster_number,
    profile_number,
    cluster_track,
    profile_track
  )

phylogeny_tracks <- phylogeny_track_data %>%
  select(genome_id, st_track, cluster_track, profile_track) %>%
  column_to_rownames("genome_id") %>%
  as.data.frame()

colnames(phylogeny_tracks) <- c(
  "Major ST",
  "350-AA cluster",
  "Profile"
)
