# Processed data functions

validate_processed_data <- function(data) {
  if (!is.list(data) || is.null(names(data)) || any(names(data) == "")) {
    stop("Processed data must be a named list.")
  }

  if (!"provenance" %in% names(data) || !is.list(data$provenance)) {
    stop("Processed data must contain a provenance list.")
  }

  required_provenance <- c(
    "builder", "input_files", "generated_at", "pilot_run_id",
    "coverage_run_id", "coverage_status"
  )
  missing_provenance <- setdiff(required_provenance, names(data$provenance))

  if (length(missing_provenance) > 0) {
    stop(
      "Processed data provenance is missing: ",
      paste(missing_provenance, collapse = ", ")
    )
  }

  invisible(data)
}

create_data_provenance <- function(
  builder,
  input_files,
  pilot_run_id = NA_character_,
  coverage_run_id = NA_character_,
  coverage_status = "current"
) {
  list(
    builder = builder,
    input_files = unname(input_files),
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    pilot_run_id = pilot_run_id,
    coverage_run_id = coverage_run_id,
    coverage_status = coverage_status
  )
}

save_processed_data <- function(data, path, overwrite = TRUE) {
  validate_processed_data(data)

  output_dir <- dirname(path)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(path) && !overwrite) {
    stop("Processed data already exist: ", path)
  }

  temporary_file <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = output_dir
  )
  on.exit(unlink(temporary_file), add = TRUE)

  saveRDS(data, temporary_file)
  if (!file.rename(temporary_file, path)) {
    stop("Processed data could not be moved into place: ", path)
  }

  invisible(path)
}

load_processed_data <- function(
  path,
  objects = NULL,
  envir = parent.frame(),
  include_provenance = FALSE
) {
  if (!file.exists(path)) {
    stop(
      "Processed data not found: ", path, "\n",
      "Run scripts/build_all_processed_data.R first."
    )
  }

  processed_data <- readRDS(path)
  validate_processed_data(processed_data)

  available_objects <- setdiff(names(processed_data), "provenance")
  if (is.null(objects)) {
    objects <- available_objects
  }

  missing_objects <- setdiff(objects, available_objects)
  if (length(missing_objects) > 0) {
    stop(
      "Processed data are missing: ",
      paste(missing_objects, collapse = ", ")
    )
  }

  selected_data <- processed_data[objects]
  if (include_provenance) {
    selected_data$provenance <- processed_data$provenance
  }

  existing_objects <- intersect(
    names(selected_data),
    ls(envir = envir, all.names = TRUE)
  )
  changed_objects <- existing_objects[!vapply(
    existing_objects,
    function(object) {
      identical(get(object, envir = envir), selected_data[[object]])
    },
    logical(1)
  )]

  if (length(changed_objects) > 0) {
    stop(
      "Processed data would overwrite different object(s): ",
      paste(changed_objects, collapse = ", ")
    )
  }

  list2env(
    selected_data[setdiff(names(selected_data), existing_objects)],
    envir = envir
  )
  invisible(processed_data)
}
