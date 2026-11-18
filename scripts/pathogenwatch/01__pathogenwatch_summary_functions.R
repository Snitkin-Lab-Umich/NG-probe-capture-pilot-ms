# Pathogenwatch summary functions

# Resistance variables
resistance_variables <- c(
  "Azithromycin", "Ceftriaxone", "Cefixime", "Ciprofloxacin",
  "Penicillin", "Sulfonamides", "Spectinomycin", "Tetracycline"
)

# Legacy alias for Figure S2
resistance_vars <- resistance_variables

# Collection metadata cleaning
get_pathogenwatch_collection_year <- function(metadata) {
  ifelse(
    is.na(metadata$year),
    metadata$`Collection date`,
    ifelse(
      is.na(metadata$`Collection date`) & is.na(metadata$year),
      "Unknown",
      metadata$year
    )
  )
}

standardize_pathogenwatch_country <- function(country) {
  dplyr::recode(
    country,
    "Brasil" = "Brazil",
    "HongKong" = "Hong Kong",
    "The Netherlands" = "Netherlands",
    "UK" = "United Kingdom",
    "USA" = "United States",
    "Viet Nam" = "Vietnam",
    .default = country,
    .missing = "Unknown"
  )
}

clean_pathogenwatch_collection_metadata <- function(metadata) {
  collection_year_recode <- get_pathogenwatch_collection_year(metadata)
  collection_country <- standardize_pathogenwatch_country(metadata$Country)
  collection_country <- dplyr::recode(
    collection_country,
    "Scotland" = "United Kingdom",
    "Caribbean" = "Dominican Republic",
    .default = collection_country
  )

  metadata %>%
    dplyr::transmute(
      pathogenwatch_genome_id = as.character(id),
      genome_id = as.character(displayname),
      collection_year = suppressWarnings(as.numeric(collection_year_recode)),
      collection_decade = floor(collection_year / 10) * 10,
      collection_country = collection_country,
      collection_continent = countrycode::countrycode(
        collection_country,
        origin = "country.name",
        destination = "continent",
        warn = FALSE
      )
    )
}

# Collection year
recode_collection_year <- function(metadata) {
  metadata$year_recode <- get_pathogenwatch_collection_year(metadata)

  metadata
}

# Collection country
recode_collection_country <- function(metadata) {
  metadata %>%
    dplyr::mutate(
      Country_recode = standardize_pathogenwatch_country(Country)
    )
}

# Sequence type distribution
summarise_sequence_types <- function(
  metadata,
  n_top = 51,
  sequence_type_column = "ST"
) {
  STs <- metadata %>%
    dplyr::filter(!is.na(.data[[sequence_type_column]])) %>%
    dplyr::count(ST = .data[[sequence_type_column]], name = "n") %>%
    dplyr::mutate(prop = n / sum(n), perc = prop * 100) %>%
    dplyr::arrange(dplyr::desc(n))

  STs_sorted <- STs %>%
    dplyr::arrange(dplyr::desc(n))

  top_STs_vec <- STs_sorted$ST[1:n_top]

  STs_summary <- STs %>%
    dplyr::mutate(ST_group = ifelse(ST %in% top_STs_vec, as.character(ST), "Other STs")) %>%
    dplyr::group_by(ST_group) %>%
    dplyr::summarise(n = sum(n), perc = 100 * n / sum(STs$n), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(perc))

  STs_summary$ST_group <- factor(STs_summary$ST_group, levels = STs_summary$ST_group)

  STs_summary
}

# Resistance distribution
summarise_resistance <- function(
  metadata,
  isolate_column,
  resistance_levels,
  resistance_columns = resistance_variables
) {
  reshaped_resistance <- reshape2::melt(
    metadata %>%
      dplyr::select(dplyr::all_of(c(isolate_column, resistance_columns))),
    id.vars = isolate_column
  )

  reshaped_resistance$value <- factor(
    reshaped_resistance$value,
    levels = resistance_levels
  )

  reshaped_resistance$variable <- factor(
    reshaped_resistance$variable,
    levels = reshaped_resistance %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(n = sum(value %in% c("INTERMEDIATE", "RESISTANT")), .groups = "drop") %>%
      dplyr::arrange(n) %>%
      dplyr::pull(variable)
  )

  reshaped_resistance %>%
    dplyr::group_by(variable, value) %>%
    dplyr::summarise(
      n = dplyr::n(),
      perc = 100 * n / length(unique(reshaped_resistance[[isolate_column]])),
      .groups = "drop"
    )
}

# Mutational distribution by gene
summarise_mutation_genes <- function(
  metadata,
  isolates = NULL,
  isolate_column = "Isolate",
  mutation_columns = NULL
) {
  if (is.null(mutation_columns)) {
    mutation_columns <- setdiff(
      colnames(metadata),
      c(
        "Isolate", resistance_variables,
        "Genome Name.y", "ST", "ST_recode"
      )
    )
  }

  mutation_metadata <- metadata %>%
    dplyr::transmute(
      Isolate = .data[[isolate_column]],
      dplyr::across(dplyr::all_of(mutation_columns))
    )

  if (!is.null(isolates)) {
    mutation_metadata <- mutation_metadata %>%
      subset(Isolate %in% isolates)
  }

  reshaped_mutation <- reshape2::melt(
    mutation_metadata,
    id.vars = "Isolate"
  )

  reshaped_mutation$value <- factor(reshaped_mutation$value, levels = c(0, 1))

  reshaped_mutation$variable <- factor(
    reshaped_mutation$variable,
    levels = reshaped_mutation %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(n = sum(value %in% 1), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::pull(variable)
  )

  reshaped_mutation_gene_collapsed <- reshaped_mutation %>%
    dplyr::mutate(
      value = as.numeric(as.character(value)),
      gene = sub("_.*", "", variable)
    ) %>%
    dplyr::group_by(Isolate, gene) %>%
    dplyr::summarise(
      present = as.numeric(any(value == 1)),
      .groups = "drop"
    ) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n = sum(present),
      perc = 100 * n / dplyr::n_distinct(reshaped_mutation$Isolate),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n))

  reshaped_mutation_gene_collapsed$gene <- factor(
    reshaped_mutation_gene_collapsed$gene,
    levels = reshaped_mutation_gene_collapsed$gene
  )

  reshaped_mutation_gene_collapsed
}
