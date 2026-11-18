# Build all processed data

builders <- c(
  "./scripts/pathogenwatch/04__build_pathogenwatch_data.R",
  "./scripts/phylogeny/05__build_phylogeny_data.R",
  "./scripts/NG_probe/build_ng_probe_data.R",
  "./scripts/amrsearch_10_20_2026_probes/build_amrsearch_data.R",
  paste0(
    "./scripts/PathogenWatch_AMR_coverage/",
    "build_pathogenwatch_coverage_data.R"
  ),
  "./scripts/reference_genome/build_reference_blast_data.R",
  "./scripts/probe_mapping/build_probe_mapping_data.R"
)

for (builder in builders) {
  message("Running ", builder)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    builder
  )

  if (!identical(status, 0L)) {
    stop("Processed data builder failed: ", builder)
  }
}

message("All processed data were created successfully.")
