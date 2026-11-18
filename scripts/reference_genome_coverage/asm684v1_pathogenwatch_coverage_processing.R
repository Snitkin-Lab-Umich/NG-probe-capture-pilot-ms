# ASM684v1 Pathogenwatch reference coverage

median_depth_value <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

median_nonzero_value <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  if (!any(x > 0, na.rm = TRUE)) return(0)
  median(x[x > 0 & !is.na(x)])
}

save_asm684v1_coverage_data <- function(x, name, tsv = FALSE) {
  output_dir <- "./data/reference_genome_coverage/curated_datasets"
  output_prefix <- file.path(output_dir, "ASM684v1_PathogenWatch_")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, paste0(output_prefix, name, ".RDS"))
  if (tsv) readr::write_tsv(x, paste0(output_prefix, name, ".tsv"))
}

prepare_gene_vs_chromosome_coverage <- function(
  isolates,
  coverage_qc,
  clinical_locus_coverage,
  gene_calls,
  locus_coverage
) {
  comparison_loci <- coverage_qc %>%
    dplyr::semi_join(gene_calls, by = c("isolate", "locus")) %>%
    dplyr::distinct(locus) %>%
    dplyr::pull(locus)

  gene_based_coverage <- clinical_locus_coverage %>%
    dplyr::filter(locus %in% comparison_loci) %>%
    dplyr::transmute(
      isolate,
      locus,
      normalized_depth,
      coverage_method = "Gene-based coverage"
    )

  chromosome_based_coverage <- locus_coverage %>%
    dplyr::filter(qc_passing, isolate %in% isolates) %>%
    dplyr::transmute(
      isolate,
      locus,
      normalized_depth = median_nonzero_normalized_depth,
      coverage_method = "Reference chromosome-based coverage"
    )

  shared_loci <- sort(intersect(
    gene_based_coverage$locus,
    chromosome_based_coverage$locus
  ))
  gene_only_loci <- sort(setdiff(
    gene_based_coverage$locus,
    chromosome_based_coverage$locus
  ))
  chromosome_only_loci <- sort(setdiff(
    chromosome_based_coverage$locus,
    gene_based_coverage$locus
  ))
  locus_order <- c(shared_loci, gene_only_loci, chromosome_only_loci)

  gene_vs_chromosome_coverage <- dplyr::bind_rows(
    gene_based_coverage,
    chromosome_based_coverage
  ) %>%
    dplyr::mutate(
      locus = factor(locus, levels = locus_order),
      coverage_method = factor(
        coverage_method,
        levels = c(
          "Gene-based coverage",
          "Reference chromosome-based coverage"
        )
      )
    )

  method_counts <- gene_vs_chromosome_coverage %>%
    dplyr::count(coverage_method, isolate, name = "num_loci")

  if (length(locus_order) == 0 ||
      anyDuplicated(gene_vs_chromosome_coverage[c(
        "coverage_method", "isolate", "locus"
      )]) ||
      anyNA(gene_vs_chromosome_coverage) ||
      any(!is.finite(gene_vs_chromosome_coverage$normalized_depth)) ||
      any(gene_vs_chromosome_coverage$normalized_depth < 0) ||
      !setequal(gene_vs_chromosome_coverage$isolate, isolates) ||
      any(method_counts$num_loci <= 0)) {
    stop("Gene and chromosome coverage data could not be validated.")
  }

  gene_vs_chromosome_coverage
}

prepare_asm684v1_pathogenwatch_coverage_data <- function(
  qcd_data_passing,
  amr_df,
  pathogenwatch_blast,
  pathogenwatch_locus_metadata
) {
  blast_hits <- pathogenwatch_blast %>%
    mutate(
      source_row = row_number(),
      gene = sseqid,
      chromosome = qseqid,
      start = pmin(qstart, qend),
      end = pmax(qstart, qend),
      strand = if_else(
        sign(qend - qstart) == sign(send - sstart),
        "+",
        "-"
      ),
      interval_length = end - start + 1L,
      gene_span = abs(send - sstart) + 1L,
      partial_hit = gene_span < slen,
      hit_id = paste(
        gene, chromosome, qstart, qend, sstart, send, source_row,
        sep = ":"
      )
    ) %>%
    group_by(gene) %>%
    mutate(hit_number = match(start, sort(unique(start)))) %>%
    ungroup()

  chromosome_coverage <- readRDS(
    "./data/reference_genome_coverage/FA1090_genome_coverage_summary_statistics.RDS"
  )

  isolate_data <- chromosome_coverage %>%
    transmute(
      isolate,
      chromosome_median = median,
      chromosome_median_non_zero = median_non_zero,
      qc_passing = isolate %in% qcd_data_passing$sample,
      depth_file = paste0(
        "./data/reference_genome_coverage/", isolate, "_summary"
      )
    )

  if (anyDuplicated(isolate_data$isolate) ||
      any(!is.finite(isolate_data$chromosome_median_non_zero) |
            isolate_data$chromosome_median_non_zero <= 0) ||
      !setequal(
        isolate_data$isolate[isolate_data$qc_passing],
        qcd_data_passing$sample
      ) ||
      any(!file.exists(isolate_data$depth_file))) {
    stop("Chromosome normalization data could not be validated.")
  }

  reference <- Biostrings::readDNAStringSet(
    "./data/reference_genome/GCF_000006845.1_ASM684v1_genomic.fasta"
  )
  names(reference) <- sub(" .*", "", names(reference))

  genes <- Biostrings::readDNAStringSet(
    paste0(
      "./data/PathogenWatch_AMR_coverage/database/concatenated_genes/",
      "PathogenWatch_reference_genes.fasta"
    )
  )
  names(genes) <- sub(" .*", "", names(genes))

  alignment_matrix <- pwalign::nucleotideSubstitutionMatrix(
    match = 2,
    mismatch = -3,
    baseOnly = TRUE
  )

  map_hit <- function(i) {
    hit <- blast_hits[i, ]
    query <- Biostrings::subseq(
      reference[[hit$chromosome]],
      hit$start,
      hit$end
    )
    subject <- Biostrings::subseq(
      genes[[hit$gene]],
      min(hit$sstart, hit$send),
      max(hit$sstart, hit$send)
    )

    if (hit$qstart > hit$qend) query <- Biostrings::reverseComplement(query)
    if (hit$sstart > hit$send) subject <- Biostrings::reverseComplement(subject)

    if (hit$gaps > 0) {
      alignment <- pwalign::pairwiseAlignment(
        query,
        subject,
        type = "global",
        substitutionMatrix = alignment_matrix,
        gapOpening = 4,
        gapExtension = 2
      )
      query_aligned <- as.character(pwalign::alignedPattern(alignment))
      subject_aligned <- as.character(pwalign::alignedSubject(alignment))
    } else {
      query_aligned <- as.character(query)
      subject_aligned <- as.character(subject)
    }

    query_bases <- strsplit(query_aligned, "")[[1]]
    subject_bases <- strsplit(subject_aligned, "")[[1]]
    gap <- query_bases == "-" | subject_bases == "-"
    actual <- c(
      length(query_bases),
      sum(query_bases == subject_bases & !gap),
      sum(query_bases != subject_bases & !gap),
      sum(gap),
      sum(diff(c(FALSE, query_bases == "-")) == 1) +
        sum(diff(c(FALSE, subject_bases == "-")) == 1)
    )
    expected <- as.numeric(
      hit[1, c("length", "nident", "mismatch", "gaps", "gapopen")]
    )

    if (!identical(actual + 0, expected)) {
      stop("Alignment does not match BLAST hit: ", hit$hit_id)
    }

    tibble(
      hit_id = hit$hit_id,
      gene = hit$gene,
      chromosome = hit$chromosome,
      alignment_column = seq_along(query_bases),
      chromosome_position = if_else(
        query_bases != "-",
        hit$qstart +
          (cumsum(query_bases != "-") - 1L) * sign(hit$qend - hit$qstart),
        NA_real_
      ),
      position = if_else(
        subject_bases != "-",
        hit$sstart +
          (cumsum(subject_bases != "-") - 1L) * sign(hit$send - hit$sstart),
        NA_real_
      ),
      query_base = query_bases,
      subject_base = subject_bases
    )
  }

  coordinate_map <- map_dfr(seq_len(nrow(blast_hits)), map_hit)
  alignment_provenance <- tibble(
    method = "endpoint-constrained global alignment",
    match = 2,
    mismatch = -3,
    gap_opening = 4,
    gap_extension = 2,
    Biostrings = as.character(packageVersion("Biostrings")),
    pwalign = as.character(packageVersion("pwalign"))
  )

  save_asm684v1_coverage_data(blast_hits, "hit_catalogue", TRUE)
  save_asm684v1_coverage_data(coordinate_map, "coordinate_map")
  save_asm684v1_coverage_data(
    alignment_provenance,
    "alignment_provenance",
    TRUE
  )
  save_asm684v1_coverage_data(isolate_data, "isolate_data", TRUE)

  hit_loci <- coordinate_map %>%
    filter(!is.na(chromosome_position)) %>%
    transmute(Locus = paste0(chromosome, ":", chromosome_position)) %>%
    distinct()

  coverage_by_hit_position <- map_dfr(
    seq_len(nrow(isolate_data)),
    function(i) {
      sample <- isolate_data[i, ]
      depth <- data.table::fread(
        sample$depth_file,
        select = c("Locus", "Total_Depth")
      ) %>%
        as_tibble() %>%
        semi_join(hit_loci, by = "Locus") %>%
        transmute(
          chromosome = sub(":.*$", "", Locus),
          chromosome_position = as.integer(sub("^.*:", "", Locus)),
          depth = Total_Depth
        )

      coordinate_map %>%
        left_join(
          depth,
          by = c("chromosome", "chromosome_position"),
          relationship = "many-to-one"
        ) %>%
        mutate(
          isolate = sample$isolate,
          qc_passing = sample$qc_passing,
          chromosome_median_non_zero = sample$chromosome_median_non_zero,
          normalized_depth = depth / chromosome_median_non_zero
        )
    }
  )

  hit_coverage <- coverage_by_hit_position %>%
    filter(!is.na(chromosome_position)) %>%
    group_by(
      isolate,
      qc_passing,
      gene,
      hit_id,
      chromosome_median_non_zero
    ) %>%
    summarise(
      observed_bases = sum(!is.na(depth)),
      median_depth = median_depth_value(depth),
      median_non_zero_depth = median_nonzero_value(depth),
      zero_depth_fraction = if (all(is.na(depth))) {
        NA_real_
      } else {
        mean(depth == 0, na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    mutate(
      median_normalized_depth = median_depth / chromosome_median_non_zero,
      median_non_zero_normalized_depth =
        median_non_zero_depth / chromosome_median_non_zero
    ) %>%
    left_join(
      blast_hits %>%
        select(hit_id, hit_number, interval_length, strand, partial_hit),
      by = "hit_id"
    )

  save_asm684v1_coverage_data(
    coverage_by_hit_position,
    "hit_positional_coverage"
  )
  save_asm684v1_coverage_data(hit_coverage, "hit_coverage_summary", TRUE)

  gene_coverage <- hit_coverage %>%
    group_by(isolate, qc_passing, gene) %>%
    summarise(
      n_hits = n(),
      n_observed_hits = sum(!is.na(median_non_zero_depth)),
      median_depth = median_depth_value(median_depth),
      median_non_zero_depth = median_depth_value(median_non_zero_depth),
      median_normalized_depth = median_depth_value(median_normalized_depth),
      median_non_zero_normalized_depth =
        median_depth_value(median_non_zero_normalized_depth),
      .groups = "drop"
    )

  gene_positional_coverage <- coverage_by_hit_position %>%
    filter(!is.na(position)) %>%
    group_by(isolate, qc_passing, gene, position) %>%
    summarise(
      n_mapped_hits = sum(!is.na(depth)),
      n_positive_hits = sum(depth > 0, na.rm = TRUE),
      depth = median_nonzero_value(depth),
      normalized_depth = median_nonzero_value(normalized_depth),
      .groups = "drop"
    )

  stats <- gene_positional_coverage %>%
    filter(qc_passing) %>%
    group_by(gene, position) %>%
    summarise(
      n_isolates = sum(!is.na(normalized_depth)),
      n_positive_isolates = sum(normalized_depth > 0, na.rm = TRUE),
      median_normalized_depth = median_nonzero_value(normalized_depth),
      .groups = "drop"
    )

  save_asm684v1_coverage_data(gene_coverage, "gene_coverage_summary", TRUE)
  save_asm684v1_coverage_data(
    gene_positional_coverage,
    "gene_positional_coverage"
  )
  save_asm684v1_coverage_data(stats, "heatmap_statistics", TRUE)

  mutation_positions <- amr_df %>%
    filter(
      gene %in% blast_hits$gene,
      !is.na(position_gene_aware),
      position_gene_aware >= 1,
      position_gene_aware <= gene_length
    ) %>%
    distinct(gene, position = position_gene_aware)

  genes_without_reference_hits <- tibble(gene = names(genes)) %>%
    filter(!gene %in% blast_hits$gene)

  save_asm684v1_coverage_data(amr_df, "database_mutations", TRUE)
  save_asm684v1_coverage_data(
    genes_without_reference_hits,
    "genes_without_reference_hits",
    TRUE
  )

  combo_genes <- blast_hits %>%
    distinct(gene, slen) %>%
    arrange(gene)
  combo_positions <- seq_len(max(combo_genes$slen))
  combo_matrix <- matrix(
    NA_real_,
    nrow(combo_genes),
    length(combo_positions),
    dimnames = list(combo_genes$gene, combo_positions)
  )
  combo_matrix[cbind(
    match(stats$gene, combo_genes$gene),
    stats$position
  )] <- stats$median_normalized_depth
  combo_boxes <- lapply(combo_genes$gene, function(gene_name) {
    gene_coverage %>%
      filter(qc_passing, gene == gene_name) %>%
      pull(median_non_zero_normalized_depth)
  })

  save_asm684v1_coverage_data(
    list(
      genes = combo_genes,
      matrix = combo_matrix,
      boxplots = combo_boxes,
      mutation_positions = mutation_positions
    ),
    "figure_data"
  )

  blast_hit_data <- blast_hits %>%
    transmute(
      hit_id,
      hit_number,
      source_row,
      reference_id = isolate_no,
      chromosome,
      chromosome_start = qstart,
      chromosome_end = qend,
      interval_start = start,
      interval_end = end,
      locus = gene,
      locus_length = slen,
      locus_start = sstart,
      locus_end = send,
      alignment_length = length,
      percent_identity = pident,
      identical_bases = nident,
      mismatches = mismatch,
      gap_openings = gapopen,
      gaps,
      e_value = evalue,
      percent_coverage = coverage,
      strand,
      interval_length,
      locus_span = gene_span,
      partial_hit
    )

  coordinate_data <- coordinate_map %>%
    transmute(
      hit_id,
      locus = gene,
      chromosome,
      alignment_column,
      chromosome_position = as.integer(chromosome_position),
      locus_position = as.integer(position),
      chromosome_base = query_base,
      locus_base = subject_base
    )

  isolate_coverage_data <- isolate_data %>%
    dplyr::rename(chromosome_median_nonzero = chromosome_median_non_zero)

  hit_positional_coverage <- coverage_by_hit_position %>%
    dplyr::rename(
      locus = gene,
      locus_position = position,
      chromosome_base = query_base,
      locus_base = subject_base,
      chromosome_median_nonzero = chromosome_median_non_zero
    )

  hit_coverage_data <- hit_coverage %>%
    dplyr::rename(
      locus = gene,
      chromosome_median_nonzero = chromosome_median_non_zero,
      median_nonzero_depth = median_non_zero_depth,
      median_nonzero_normalized_depth =
        median_non_zero_normalized_depth
    )

  locus_coverage_data <- gene_coverage %>%
    dplyr::rename(
      locus = gene,
      median_nonzero_depth = median_non_zero_depth,
      median_nonzero_normalized_depth =
        median_non_zero_normalized_depth
    )

  locus_positional_coverage <- gene_positional_coverage %>%
    dplyr::rename(locus = gene, locus_position = position)

  heatmap_statistics <- stats %>%
    dplyr::rename(locus = gene, locus_position = position)

  mutation_position_data <- mutation_positions %>%
    dplyr::rename(locus = gene, locus_position = position)

  reference_locus_metadata <- combo_genes %>%
    transmute(locus = gene, locus_length = as.integer(slen))

  reference_locus_lengths <- pathogenwatch_locus_metadata %>%
    dplyr::semi_join(reference_locus_metadata, by = "locus") %>%
    dplyr::arrange(match(locus, reference_locus_metadata$locus))

  loci_without_reference_hits <- genes_without_reference_hits %>%
    dplyr::rename(locus = gene)

  if (anyDuplicated(blast_hit_data$hit_id) ||
      anyDuplicated(hit_coverage_data[c("isolate", "hit_id")]) ||
      anyDuplicated(
        locus_positional_coverage[c("isolate", "locus", "locus_position")]
      ) ||
      anyDuplicated(heatmap_statistics[c("locus", "locus_position")]) ||
      !setequal(
        hit_coverage_data$isolate[hit_coverage_data$qc_passing],
        qcd_data_passing$sample
      ) ||
      anyDuplicated(pathogenwatch_locus_metadata$locus) ||
      !identical(
        reference_locus_lengths$locus,
        reference_locus_metadata$locus
      ) ||
      !identical(
        as.integer(reference_locus_lengths$locus_length),
        reference_locus_metadata$locus_length
      ) ||
      !all(mutation_position_data$locus %in% reference_locus_metadata$locus) ||
      any(
        mutation_position_data$locus_position >
          reference_locus_metadata$locus_length[
            match(
              mutation_position_data$locus,
              reference_locus_metadata$locus
            )
          ]
      )) {
    stop("ASM684v1 Pathogenwatch coverage data could not be validated.")
  }

  list(
    blast_hit_data = blast_hit_data,
    coordinate_data = coordinate_data,
    alignment_provenance = alignment_provenance,
    isolate_coverage_data = isolate_coverage_data,
    hit_positional_coverage = hit_positional_coverage,
    hit_coverage = hit_coverage_data,
    hit_coverage_qc = hit_coverage_data %>% filter(qc_passing),
    locus_coverage = locus_coverage_data,
    locus_coverage_qc = locus_coverage_data %>% filter(qc_passing),
    locus_positional_coverage = locus_positional_coverage,
    heatmap_statistics = heatmap_statistics,
    mutation_positions = mutation_position_data,
    locus_metadata = pathogenwatch_locus_metadata,
    reference_loci = reference_locus_metadata$locus,
    loci_without_reference_hits = loci_without_reference_hits
  )
}
