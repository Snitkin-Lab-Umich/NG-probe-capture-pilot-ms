# Probe mapping data

probe_mapping_input_paths <- c(
  probe_fasta = paste0(
    "./data/probe_mapping/",
    "240808-903-Stefan_Green-D10891NGNRH-Ngonorrhoeae.fas"
  ),
  bait_bed = "./data/probe_mapping/ref-139Kbaits.sorted.merge.bed",
  reference_gff = paste0(
    "./data/probe_mapping/",
    "GCF_000006845.1_ASM684v1_genomic.gff"
  ),
  reference_fasta = paste0(
    "./data/probe_mapping/",
    "GCF_000006845.1_ASM684v1_genomic.fna"
  )
)

read_reference_bed <- function(path) {
  readr::read_tsv(
    path,
    col_names = c("seqid", "start0", "end0"),
    col_types = readr::cols(
      seqid = readr::col_character(),
      start0 = readr::col_integer(),
      end0 = readr::col_integer()
    ),
    progress = FALSE
  )
}

reference_bed_to_granges <- function(data) {
  GenomicRanges::GRanges(
    seqnames = data$seqid,
    ranges = IRanges::IRanges(start = data$start0 + 1L, end = data$end0)
  )
}

collapse_reference_metadata <- function(x) {
  if (inherits(x, "List") || is.list(x)) {
    return(vapply(
      as.list(x),
      function(value) {
        value <- as.character(value)
        if (length(value) == 0L) NA_character_ else paste(value, collapse = ",")
      },
      character(1)
    ))
  }
  as.character(x)
}

reference_overlap_width <- function(query, subject) {
  subject <- GenomicRanges::reduce(subject, ignore.strand = TRUE)
  hits <- GenomicRanges::findOverlaps(query, subject, ignore.strand = TRUE)
  answer <- integer(length(query))
  if (length(hits) == 0L) return(answer)

  overlap_widths <- BiocGenerics::width(
    GenomicRanges::pintersect(
      query[S4Vectors::queryHits(hits)],
      subject[S4Vectors::subjectHits(hits)],
      ignore.strand = TRUE
    )
  )
  sums <- rowsum(overlap_widths, S4Vectors::queryHits(hits), reorder = FALSE)
  answer[as.integer(rownames(sums))] <- as.integer(sums[, 1])
  answer
}

make_relative_gene_bins <- function(gene_range, number_of_bins = 10L) {
  gene_width <- BiocGenerics::width(gene_range)
  boundaries <- floor(seq(0, gene_width, length.out = number_of_bins + 1L))
  relative_start <- boundaries[seq_len(number_of_bins)] + 1L
  relative_end <- boundaries[seq_len(number_of_bins) + 1L]

  if (as.character(BiocGenerics::strand(gene_range)) == "-") {
    bin_start <- BiocGenerics::end(gene_range) - relative_end + 1L
    bin_end <- BiocGenerics::end(gene_range) - relative_start + 1L
  } else {
    bin_start <- BiocGenerics::start(gene_range) + relative_start - 1L
    bin_end <- BiocGenerics::start(gene_range) + relative_end - 1L
  }

  GenomicRanges::GRanges(
    seqnames = as.character(GenomicRanges::seqnames(gene_range)),
    ranges = IRanges::IRanges(start = bin_start, end = bin_end),
    strand = BiocGenerics::strand(gene_range)
  )
}

clean_reference_locus <- function(x) {
  x <- gsub("(16S|23S)_rDNA", "\\1_rRNA_gene", x)
  gsub("_", " ", x)
}

load_probe_mapping_data <- function(pathogenwatch_blast) {
  bait_file <- probe_mapping_input_paths[["bait_bed"]]
  gff_file <- probe_mapping_input_paths[["reference_gff"]]
  fasta_file <- probe_mapping_input_paths[["reference_fasta"]]

  baits <- read_reference_bed(bait_file) %>%
    reference_bed_to_granges() %>%
    GenomicRanges::reduce(ignore.strand = TRUE)
  annotations <- rtracklayer::import(gff_file, format = "GFF3")
  genome <- Biostrings::readDNAStringSet(fasta_file, use.names = TRUE)
  names(genome) <- sub(" .*", "", names(genome))
  genome_lengths <- setNames(as.integer(BiocGenerics::width(genome)), names(genome))
  genes <- annotations[annotations$type %in% c("gene", "pseudogene")]
  gene_meta <- S4Vectors::mcols(genes)
  gene_covered_bp <- reference_overlap_width(genes, baits)

  gene_coverage <- tibble::tibble(
    gene_index = seq_along(genes),
    seqid = as.character(GenomicRanges::seqnames(genes)),
    start = BiocGenerics::start(genes),
    end = BiocGenerics::end(genes),
    strand = as.character(BiocGenerics::strand(genes)),
    feature_type = as.character(gene_meta$type),
    gene_id = collapse_reference_metadata(gene_meta$ID),
    locus_tag = collapse_reference_metadata(gene_meta$locus_tag),
    gene_symbol = collapse_reference_metadata(gene_meta$gene),
    name = collapse_reference_metadata(gene_meta$Name),
    gene_biotype = collapse_reference_metadata(gene_meta$gene_biotype),
    gene_length_bp = BiocGenerics::width(genes),
    covered_bp = gene_covered_bp
  ) %>%
    dplyr::mutate(
      coverage_percent = 100 * covered_bp / gene_length_bp,
      coverage_status = dplyr::case_when(
        covered_bp == gene_length_bp ~ "Fully covered",
        covered_bp == 0L ~ "Not covered",
        TRUE ~ "Partially covered"
      ),
      coverage_status = factor(
        coverage_status,
        levels = c("Fully covered", "Partially covered", "Not covered")
      )
    )

  reference_seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames = names(genome_lengths),
    seqlengths = genome_lengths
  )
  genome_bins <- GenomicRanges::tileGenome(
    reference_seqinfo,
    tilewidth = 10000L,
    cut.last.tile.in.chrom = TRUE
  )
  bin_covered_bp <- reference_overlap_width(genome_bins, baits)
  bin_coverage <- tibble::tibble(
    start = BiocGenerics::start(genome_bins),
    end = BiocGenerics::end(genome_bins),
    midpoint_mb = (start + end) / 2e6,
    coverage_percent = 100 * bin_covered_bp / BiocGenerics::width(genome_bins)
  )

  gene_track_data <- gene_coverage %>%
    dplyr::mutate(start_mb = start / 1e6, end_mb = end / 1e6)
  bait_track_data <- as.data.frame(baits) %>%
    dplyr::transmute(start_mb = start / 1e6, end_mb = end / 1e6)
  partial_genes <- gene_coverage %>%
    dplyr::filter(coverage_status == "Partially covered")

  gene_bin_ranges <- GenomicRanges::GRangesList(
    lapply(seq_along(genes), function(index) make_relative_gene_bins(genes[index]))
  ) %>%
    unlist(use.names = FALSE)
  S4Vectors::mcols(gene_bin_ranges)$gene_index <- rep(
    seq_along(genes),
    each = 10L
  )
  S4Vectors::mcols(gene_bin_ranges)$relative_bin <- rep(
    seq_len(10L),
    times = length(genes)
  )
  gene_bin_covered_bp <- reference_overlap_width(gene_bin_ranges, baits)
  relative_gene_coverage <- tibble::tibble(
    gene_index = S4Vectors::mcols(gene_bin_ranges)$gene_index,
    relative_bin = S4Vectors::mcols(gene_bin_ranges)$relative_bin,
    coverage_percent = 100 * gene_bin_covered_bp /
      BiocGenerics::width(gene_bin_ranges)
  ) %>%
    dplyr::filter(gene_index %in% partial_genes$gene_index) %>%
    dplyr::mutate(
      relative_position = factor(
        relative_bin,
        levels = seq_len(10L),
        labels = paste0(seq(0, 90, 10), "-", seq(10, 100, 10), "%")
      )
    )

  blast_annot <- pathogenwatch_blast %>%
    dplyr::filter(sseqid != "mtrR_promoter_mosaic_1") %>%
    dplyr::mutate(q_lo = pmin(qstart, qend), q_hi = pmax(qstart, qend)) %>%
    dplyr::left_join(
      gene_coverage %>%
        dplyr::mutate(g_lo = pmin(start, end), g_hi = pmax(start, end)),
      by = dplyr::join_by(overlaps(q_lo, q_hi, g_lo, g_hi)),
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(
      overlap_bp = pmax(0, pmin(q_hi, g_hi) - pmax(q_lo, g_lo) + 1),
      frac_of_hit = overlap_bp / (q_hi - q_lo + 1),
      frac_of_gene = overlap_bp / (g_hi - g_lo + 1)
    )

  amr_gene_summary <- blast_annot %>%
    dplyr::group_by(qseqid, qstart, qend) %>%
    dplyr::slice_max(overlap_bp, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(isolate_no, seqid, sseqid) %>%
    dplyr::arrange(q_lo, .by_group = TRUE) %>%
    dplyr::mutate(
      locus = dplyr::if_else(
        sseqid %in% c("16S_rDNA", "23S_rDNA"),
        paste0(sseqid, "_", dplyr::row_number()),
        sseqid
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      !locus %in% c("mtrR_promoter", "macAB_promoter", "norM_promoter")
    )

  reference_genome_summary <- tibble::tibble(
    genome_bp = sum(genome_lengths),
    covered_bp = sum(BiocGenerics::width(baits))
  ) %>%
    dplyr::mutate(
      uncovered_bp = genome_bp - covered_bp,
      genome_covered_percent = 100 * covered_bp / genome_bp
    )

  reference_gene_summary <- gene_coverage %>%
    dplyr::summarise(
      total_genes = dplyr::n(),
      covered_genes = sum(covered_bp > 0),
      covered_genes_percent = 100 * mean(covered_bp > 0),
      fully_covered_genes = sum(coverage_status == "Fully covered"),
      fully_covered_genes_percent = 100 * mean(coverage_status == "Fully covered"),
      partially_covered_genes = sum(coverage_status == "Partially covered"),
      uncovered_genes = sum(coverage_status == "Not covered")
    )

  partial_gene_summary <- partial_genes %>%
    dplyr::summarise(
      partially_covered_genes = dplyr::n(),
      partial_genes_at_least_80_percent = sum(coverage_percent >= 80),
      partial_genes_at_least_80_percent_value =
        100 * mean(coverage_percent >= 80)
    )

  amr_gene_coverage <- amr_gene_summary %>%
    dplyr::transmute(
      locus = as.character(locus),
      reference_locus = sseqid,
      locus_tag,
      gene_symbol,
      gene_length_bp,
      covered_bp,
      coverage_percent,
      coverage_status
    )

  amr_reference_summary <- amr_gene_coverage %>%
    dplyr::filter(!reference_locus %in% c("16S_rDNA", "23S_rDNA")) %>%
    dplyr::summarise(
      resistance_loci = dplyr::n_distinct(reference_locus),
      loci_with_probe_coverage = dplyr::n_distinct(
        reference_locus[covered_bp > 0]
      ),
      fully_covered_loci = dplyr::n_distinct(
        reference_locus[coverage_status == "Fully covered"]
      )
    )

  amr_gene_segments <- purrr::map_dfr(seq_len(nrow(amr_gene_summary)), function(index) {
    gene_row <- amr_gene_summary[index, ]
    gene_range <- genes[gene_row$gene_index]
    covered_ranges <- GenomicRanges::intersect(
      baits,
      gene_range,
      ignore.strand = TRUE
    )
    if (length(covered_ranges) == 0L) {
      return(tibble::tibble(
        locus = character(),
        segment_start = numeric(),
        segment_end = numeric()
      ))
    }

    if (gene_row$strand == "-") {
      segment_start <- 100 * (
        BiocGenerics::end(gene_range) - BiocGenerics::end(covered_ranges)
      ) / BiocGenerics::width(gene_range)
      segment_end <- 100 * (
        BiocGenerics::end(gene_range) - BiocGenerics::start(covered_ranges) + 1L
      ) / BiocGenerics::width(gene_range)
    } else {
      segment_start <- 100 * (
        BiocGenerics::start(covered_ranges) - BiocGenerics::start(gene_range)
      ) / BiocGenerics::width(gene_range)
      segment_end <- 100 * (
        BiocGenerics::end(covered_ranges) - BiocGenerics::start(gene_range) + 1L
      ) / BiocGenerics::width(gene_range)
    }

    tibble::tibble(
      locus = gene_row$locus,
      segment_start = segment_start,
      segment_end = segment_end
    )
  })

  amr_locus_levels <- amr_gene_summary %>%
    dplyr::arrange(dplyr::desc(coverage_percent), locus) %>%
    dplyr::pull(locus) %>%
    rev() %>%
    clean_reference_locus()
  amr_gene_segments <- amr_gene_segments %>%
    dplyr::mutate(
      locus = factor(
        clean_reference_locus(as.character(locus)),
        levels = amr_locus_levels
      ),
      segment_midpoint = (segment_start + segment_end) / 2,
      segment_width = segment_end - segment_start
    )
  amr_gene_summary <- amr_gene_summary %>%
    dplyr::mutate(
      locus = factor(
        clean_reference_locus(as.character(locus)),
        levels = amr_locus_levels
      )
    )

  if (
    nrow(reference_genome_summary) != 1 ||
      reference_genome_summary$genome_bp != 2153922 ||
      reference_genome_summary$covered_bp != 2063614 ||
      reference_genome_summary$uncovered_bp < 0 ||
      nrow(gene_coverage) != 2250 ||
      nrow(bin_coverage) != 216 ||
      nrow(bait_track_data) != 857 ||
      nrow(amr_gene_coverage) != 24 ||
      amr_reference_summary$resistance_loci != 16 ||
      anyNA(gene_coverage$coverage_percent) ||
      any(!dplyr::between(gene_coverage$coverage_percent, 0, 100)) ||
      any(reference_gene_summary$covered_genes +
          reference_gene_summary$uncovered_genes !=
          reference_gene_summary$total_genes) ||
      partial_gene_summary$partially_covered_genes !=
          reference_gene_summary$partially_covered_genes
  ) {
    stop("Probe mapping data do not match the current inputs")
  }

  list(
    reference_genome_summary = reference_genome_summary,
    gene_coverage = gene_coverage,
    reference_gene_summary = reference_gene_summary,
    partial_gene_summary = partial_gene_summary,
    amr_gene_coverage = amr_gene_coverage,
    amr_reference_summary = amr_reference_summary,
    bin_coverage = bin_coverage,
    gene_track_data = gene_track_data,
    bait_track_data = bait_track_data,
    partial_genes = partial_genes,
    relative_gene_coverage = relative_gene_coverage,
    amr_gene_summary = amr_gene_summary,
    amr_gene_segments = amr_gene_segments,
    reference_limit_mb = sum(genome_lengths) / 1e6
  )
}
