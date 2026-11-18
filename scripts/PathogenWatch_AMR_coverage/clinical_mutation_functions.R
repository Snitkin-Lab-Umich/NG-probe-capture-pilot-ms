# Shared Pathogenwatch mutation parsing functions

parse_pathogenwatch_variant <- function(v, detect_aa_insertion = FALSE) {
  if (is.na(v) || v == "") {
    return(tibble(
      position = NA_integer_, ref = NA_character_,
      alt = NA_character_, mut_type = "presence_absence"
    ))
  }
  if (v %in% c("truncated", "disrupted", "frameshift")) {
    return(tibble(
      position = NA_integer_, ref = NA_character_,
      alt = v, mut_type = "structural"
    ))
  }
  if (grepl("^ins\\d+", v)) {
    pos <- as.integer(regmatches(v, regexpr("\\d+", v)))
    alt <- gsub("^ins\\d+", "", v)
    mut_type <- if (detect_aa_insertion && grepl("^[A-Z]+$", alt)) {
      "AA_insertion"
    } else {
      "insertion"
    }
    return(tibble(
      position = pos, ref = NA_character_,
      alt = alt, mut_type = mut_type
    ))
  }
  if (grepl("^-{2,}\\d+[A-Za-z]", v)) {
    pos <- as.integer(regmatches(v, regexpr("\\d+", v)))
    alt <- gsub("^-+\\d+", "", v)
    mut_type <- if (grepl("^[A-Z]+$", alt)) "AA_insertion" else "insertion"
    return(tibble(position = pos, ref = NA_character_, alt = alt, mut_type = mut_type))
  }
  if (grepl("^[A-Za-z]-\\d+", v)) {
    ref <- substr(v, 1, 1)
    pos <- -as.integer(regmatches(v, regexpr("\\d+", v)))
    alt <- gsub("^[A-Za-z]-\\d+", "", v)
    mut_type <- if (alt %in% c("del", "-")) "promoter_deletion" else "promoter_substitution"
    return(tibble(position = pos, ref = ref, alt = alt, mut_type = mut_type))
  }
  if (grepl("^[A-Z]\\d+(-|del)$", v)) {
    return(tibble(
      position = as.integer(regmatches(v, regexpr("\\d+", v))),
      ref = substr(v, 1, 1), alt = "-", mut_type = "AA_deletion"
    ))
  }
  if (grepl("^[a-z]\\d+(-|del)$", v)) {
    return(tibble(
      position = as.integer(regmatches(v, regexpr("\\d+", v))),
      ref = substr(v, 1, 1), alt = "-", mut_type = "nucleotide_deletion"
    ))
  }
  if (grepl("^[A-Z]\\d+[A-Z]+$", v)) {
    return(tibble(
      position = as.integer(regmatches(v, regexpr("\\d+", v))),
      ref = substr(v, 1, 1),
      alt = gsub("^[A-Z]\\d+", "", v),
      mut_type = "AA_substitution"
    ))
  }
  if (grepl("^[a-z]\\d+[a-z]$", v)) {
    return(tibble(
      position = as.integer(regmatches(v, regexpr("\\d+", v))),
      ref = substr(v, 1, 1),
      alt = substr(v, nchar(v), nchar(v)),
      mut_type = "nucleotide_substitution"
    ))
  }
  tibble(
    position = NA_integer_, ref = NA_character_,
    alt = NA_character_, mut_type = "unknown"
  )
}

add_pathogenwatch_gene_positions <- function(amr_df, config) {
  gene_info <- purrr::map_dfr(config$genes, ~ tibble(
    gene = .x$name,
    gene_length = nchar(.x$sequence)
  ))

  amr_df %>%
    left_join(gene_info, by = "gene") %>%
    mutate(position_gene_aware = case_when(
      is.na(position) ~ NA_integer_,
      position < 0 ~ gene_length + position,
      mut_type == "AA_substitution" ~ (position - 1L) * 3L + 1L,
      mut_type == "AA_deletion" ~ (position - 1L) * 3L + 1L,
      mut_type == "AA_insertion" ~ (position - 1L) * 3L + 1L,
      TRUE ~ position
    ))
}

format_pathogenwatch_loci <- function(locus) {
  namedf <- unique(locus) %>%
    as.character() %>%
    stringr::str_split("_", simplify = TRUE)

  namedf[, 1] <- sapply(namedf[, 1], function(x) {
    paste0("<i>", x, "</i>", collapse = "")
  })
  namedf[, 2] <- sapply(namedf[, 2], function(x) {
    ifelse(x == "rDNA", paste0("<i>", x, "</i>", collapse = ""), x)
  })

  renamed <- apply(namedf, 1, paste0, collapse = " ") %>%
    trimws() %>%
    `names<-`(unique(as.character(locus))) %>%
    gsub("<i>blaTEM</i>", "<i>bla</i><sub>TEM</sub>", .)

  list(
    renamed = renamed,
    locus_levels = sort(unique(renamed)) %>% rev()
  )
}

add_contig_edge_distances <- function(data) {
  data %>%
    mutate(
      dist_from_start = contig_start - 1,
      dist_from_end = length - contig_stop,
      min_dist_to_edge = pmin(dist_from_start, dist_from_end)
    )
}

prepare_pathogenwatch_mutation_data <- function(config) {
  mutation_entries <- function(entries) {
    purrr::map_dfr(entries, function(m) {
      if (!is.list(m) || is.null(m$variants)) return(NULL)
      tibble(gene = m$gene, variant = unlist(m$variants))
    })
  }

  amr_df <- purrr::map_dfr(config$mechanisms, function(entry) {
    bind_rows(
      mutation_entries(entry$members),
      purrr::map_dfr(entry$phenotypes, ~ mutation_entries(.x$modifiers))
    )
  }) %>%
    distinct(gene, variant) %>%
    mutate(parsed = purrr::map(
      variant,
      ~ parse_pathogenwatch_variant(.x, detect_aa_insertion = TRUE)
    )) %>%
    tidyr::unnest(parsed)

  amr_df <- add_pathogenwatch_gene_positions(amr_df, config)

  mutation_positions <- amr_df %>%
    filter(!is.na(position_gene_aware)) %>%
    distinct(locus = gene, locus_position = position_gene_aware)

  legacy_mutation_positions <- purrr::map_dfr(
    config$mechanisms,
    function(entry) mutation_entries(entry$members)
  ) %>%
    distinct(gene, variant) %>%
    mutate(parsed = purrr::map(
      variant,
      ~ parse_pathogenwatch_variant(.x, detect_aa_insertion = FALSE)
    )) %>%
    tidyr::unnest(parsed) %>%
    add_pathogenwatch_gene_positions(config) %>%
    filter(!is.na(position_gene_aware)) %>%
    distinct(gene, position = position_gene_aware)

  list(
    amr_df = amr_df,
    mutation_positions = mutation_positions,
    legacy_mutation_positions = legacy_mutation_positions
  )
}
