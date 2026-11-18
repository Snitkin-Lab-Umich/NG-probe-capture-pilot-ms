# Pathogenwatch data

# Tree
if (!exists("tree", inherits = FALSE)) {
  tree <- ape::read.tree(
    "./data/phylogeny/source/cognac_nj.tre"
  )
}

if (
  !inherits(tree, "phylo") ||
    length(tree$tip.label) == 0 ||
    anyNA(tree$tip.label) ||
    anyDuplicated(tree$tip.label)
) {
  stop("tree must be a phylo object with unique, non-missing tip labels")
}

# Relevant metadata
ng_amrsearch <- readr::read_csv(
  "../../data/pathogenwatch_metadata/Neisseria_gonorrhoeae_amrsearch.csv",
  show_col_types = FALSE
) %>%
  subset(`Genome Name` %in% tree$tip.label)

ng_snps <- readRDS("../../data/ng_pathogenwatch_files/ng_snps.RDS") %>%
  subset(`Genome.Name` %in% tree$tip.label) %>%
  dplyr::mutate(`Genome Name` = Genome.Name)

if (!exists("ng_STs", inherits = FALSE)) {
  ng_STs <- read.csv(
    "../../data/pathogenwatch_metadata/Neisseria_gonorrhoeae_mlst.csv"
  )
}

if (
  !is.data.frame(ng_STs) ||
    !all(c("Genome.Name", "ST") %in% colnames(ng_STs)) ||
    anyNA(ng_STs$Genome.Name) ||
    anyDuplicated(ng_STs$Genome.Name)
) {
  stop("ng_STs must contain one row per non-missing Genome.Name and an ST column")
}

ng_sample_metadata <- readr::read_csv(
  "../../data/pathogenwatch_metadata/Neisseria_gonorrhoeae_metadata.csv",
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
) %>%
  subset(displayname %in% tree$tip.label)

if (
  anyNA(ng_sample_metadata$id) ||
    anyNA(ng_sample_metadata$displayname) ||
    anyDuplicated(ng_sample_metadata$id) ||
    anyDuplicated(ng_sample_metadata$displayname)
) {
  stop("Pathogenwatch metadata must contain unique, non-missing isolate IDs and names")
}

# Generate grouped data
ng_meta <- dplyr::left_join(
  ng_amrsearch[, c(-1, -3, -4)],
  ng_snps[, c(-1, -3, -4)],
  by = dplyr::join_by(`Genome Name` == `Genome.Name`)
) %>%
  as.data.frame()

colnames(ng_meta)[1] <- "Isolate"
rownames(ng_meta) <- ng_meta$Isolate

# NG SNP gene grouping
snp_data <- ng_snps[, 5:66]
prefixes <- unique(stringr::str_extract(colnames(snp_data), "^[^_]+"))

collapsed_snps <- as.data.frame(sapply(prefixes, function(prefix) {
  matched_cols <- grep(paste0("^", prefix, "_"), colnames(snp_data), value = TRUE)

  if (length(matched_cols) > 0) {
    rowSums(snp_data[, matched_cols, drop = FALSE], na.rm = TRUE)
  } else {
    snp_data[, prefix, drop = FALSE]
  }
}))

rownames(collapsed_snps) <- rownames(ng_snps)

collapsed_snps_binary <- collapsed_snps %>%
  dplyr::mutate_all(~ as.factor(ifelse(. > 0, "1", "0")))

colnames(collapsed_snps_binary)[colnames(collapsed_snps_binary) == "X16S"] <- "16S_rDNA"
colnames(collapsed_snps_binary)[colnames(collapsed_snps_binary) == "X23S"] <- "23S_rDNA"

ng_meta_gene_grouped <- dplyr::left_join(
  ng_meta[, 1:9],
  collapsed_snps_binary %>%
    as.data.frame() %>%
    dplyr::mutate(Isolate = rownames(.))
) %>%
  as.data.frame()

rownames(ng_meta_gene_grouped) <- ng_meta_gene_grouped$Isolate

# Canonical isolate data
pathogenwatch_collection_data <- clean_pathogenwatch_collection_metadata(
  ng_sample_metadata
)

pathogenwatch_resistance_data <- ng_amrsearch %>%
  dplyr::transmute(
    pathogenwatch_genome_id = as.character(`Genome ID`),
    genome_id = as.character(`Genome Name`),
    dplyr::across(dplyr::all_of(resistance_variables))
  )

pathogenwatch_sequence_type_data <- ng_STs %>%
  dplyr::transmute(
    pathogenwatch_genome_id = as.character(Genome.ID),
    genome_id = as.character(Genome.Name),
    sequence_type = ST
  ) %>%
  dplyr::filter(genome_id %in% tree$tip.label)

mutation_variables <- colnames(collapsed_snps_binary)

pathogenwatch_isolate_mutation_data <- dplyr::bind_cols(
  ng_snps %>%
    dplyr::transmute(
      pathogenwatch_genome_id = as.character(Genome.ID),
      genome_id = as.character(Genome.Name)
    ),
  collapsed_snps_binary %>%
    tibble::as_tibble(.name_repair = "minimal")
)

pathogenwatch_isolate_data <- pathogenwatch_collection_data %>%
  dplyr::left_join(
    pathogenwatch_resistance_data,
    by = c("pathogenwatch_genome_id", "genome_id")
  ) %>%
  dplyr::left_join(
    pathogenwatch_sequence_type_data,
    by = c("pathogenwatch_genome_id", "genome_id")
  ) %>%
  dplyr::left_join(
    pathogenwatch_isolate_mutation_data,
    by = c("pathogenwatch_genome_id", "genome_id")
  )

if (
  nrow(pathogenwatch_isolate_data) != nrow(ng_sample_metadata) ||
    anyNA(pathogenwatch_isolate_data$pathogenwatch_genome_id) ||
    anyNA(pathogenwatch_isolate_data$genome_id) ||
    anyDuplicated(pathogenwatch_isolate_data$pathogenwatch_genome_id) ||
    anyDuplicated(pathogenwatch_isolate_data$genome_id) ||
    !setequal(pathogenwatch_collection_data$genome_id, pathogenwatch_resistance_data$genome_id) ||
    !setequal(pathogenwatch_collection_data$genome_id, pathogenwatch_sequence_type_data$genome_id) ||
    !setequal(
      pathogenwatch_collection_data$genome_id,
      pathogenwatch_isolate_mutation_data$genome_id
    )
) {
  stop("Canonical Pathogenwatch data must contain one complete row per genome_id")
}

if (
  anyNA(pathogenwatch_isolate_data[, resistance_variables]) ||
    anyNA(pathogenwatch_isolate_data[, mutation_variables])
) {
  stop("Canonical Pathogenwatch resistance and mutation data must be complete")
}
