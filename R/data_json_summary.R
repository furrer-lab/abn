#' Compute summary metadata for data JSON
#'
#' Computes the derived `metadata.summary` block for a column-oriented data JSON
#' object. The summary contains row and column counts, total missing values, and
#' per-column missing and unique-value counts. Missing values are JSON `null` or
#' scalar R `NA` values after parsing.
#'
#' @param x A JSON string, path to a JSON file, or parsed list containing
#'   top-level `metadata` and `data` fields.
#' @return A list containing derived summary metadata.
#' @keywords internal
compute_data_json_summary <- function(x) {
  x <- parse_data_json_input(x)
  validate_data_json_shape(x)

  column_order <- data_json_as_character_vector(x$metadata$column_order)
  data <- x$data
  column_lengths <- vapply(data[column_order], length, integer(1))
  if (length(unique(column_lengths)) != 1) {
    stop("data columns must all have the same length", call. = FALSE)
  }

  columns <- list()
  missing_total <- 0L
  for (column_name in column_order) {
    values <- data[[column_name]]
    missing <- data_json_missing_vector(values)
    missing_count <- sum(missing)
    missing_total <- missing_total + missing_count
    non_missing_values <- values[!missing]
    columns[[column_name]] <- list(
      missing_count = missing_count,
      unique_count = length(unique(vapply(non_missing_values, data_json_value_key, character(1))))
    )
  }

  list(
    row_count = unname(column_lengths[[1]]),
    column_count = length(column_order),
    missing_count = missing_total,
    columns = columns
  )
}

#' Validate or repair data JSON summary metadata
#'
#' `metadata.summary` is optional derived metadata for data JSON objects. When it
#' is present, this function recomputes the summary from `data` and verifies that
#' the stored values match. Stale summaries are treated as errors by default
#' because they can mislead importers and downstream tools.
#'
#' Set `on_error = "return"` to receive diagnostics without signalling an error,
#' or `on_error = "warning"` to warn and still return diagnostics. Set
#' `repair = TRUE` to replace stale summary metadata with recomputed values. If
#' `metadata.summary` is absent, validation succeeds; with `repair = TRUE`, the
#' missing summary block is added.
#'
#' @param x A JSON string, path to a JSON file, or parsed list containing
#'   top-level `metadata` and `data` fields.
#' @param on_error One of `"error"`, `"warning"`, or `"return"`.
#' @param repair If `TRUE`, recompute and replace or add `metadata.summary`.
#' @return A diagnostics list containing validity, summaries, and optionally a
#'   repaired object. The list contains `valid`, `repaired`, `errors`,
#'   `warnings`, `summary_expected`, `summary_found`, and `object`.
#' @examples
#' data_json <- list(
#'   metadata = list(
#'     schema_version = "abn-data-v1",
#'     orientation = "columnar",
#'     column_order = c("G1", "B1"),
#'     columns = list(
#'       list(name = "G1", semantic_type = "continuous"),
#'       list(name = "B1", semantic_type = "binary")
#'     )
#'   ),
#'   data = list(
#'     G1 = list(0.1, -0.4, 1.2),
#'     B1 = list("0", NULL, "1")
#'   )
#' )
#'
#' diagnostics <- validate_data_json_summary(data_json, on_error = "return")
#' diagnostics$valid
#'
#' repaired <- repair_data_json_summary(data_json)
#' repaired$metadata$summary$row_count
#' @export
validate_data_json_summary <- function(x, on_error = c("error", "warning", "return"),
                                       repair = FALSE) {
  on_error <- match.arg(on_error)
  object <- parse_data_json_input(x)
  expected <- compute_data_json_summary(object)
  found <- object$metadata$summary

  errors <- character()
  if (is.null(found)) {
    # Summary metadata is optional. When absent, validation succeeds and repair
    # can materialize the derived summary for callers that want it persisted.
    errors <- character()
  } else {
    errors <- compare_data_json_summary(found, expected)
  }

  valid <- length(errors) == 0L
  repaired <- FALSE
  if (repair && (is.null(found) || !valid)) {
    object$metadata$summary <- expected
    repaired <- TRUE
  }

  warnings <- character()
  if (!valid) {
    if (on_error == "error" && !repair) {
      stop(paste(c("data JSON summary is inconsistent", errors), collapse = "; "),
           call. = FALSE)
    }
    if (on_error == "warning") {
      warnings <- errors
      warning(paste(c("data JSON summary is inconsistent", errors), collapse = "; "),
              call. = FALSE)
    }
  }

  list(
    valid = valid,
    repaired = repaired,
    errors = errors,
    warnings = warnings,
    summary_expected = expected,
    summary_found = found,
    object = object
  )
}

#' Repair data JSON summary metadata
#'
#' Recomputes summary metadata from `data` and stores it at `metadata.summary`.
#' Existing stale summary values are replaced. If the summary is absent, it is
#' added.
#'
#' @param x A JSON string, path to a JSON file, or parsed list.
#' @return A data JSON object with recomputed `metadata.summary`.
#' @seealso [validate_data_json_summary()]
#' @export
repair_data_json_summary <- function(x) {
  validate_data_json_summary(x, on_error = "return", repair = TRUE)$object
}

parse_data_json_input <- function(x) {
  if (is.list(x)) return(x)
  if (!is.character(x) || length(x) != 1L) {
    stop("x must be a parsed list, JSON string, or path to a JSON file", call. = FALSE)
  }
  if (file.exists(x)) {
    x <- paste(readLines(x, warn = FALSE), collapse = "\n")
  }
  jsonlite::fromJSON(x, simplifyVector = FALSE)
}

validate_data_json_shape <- function(x) {
  if (is.null(x$metadata) || !is.list(x$metadata)) {
    stop("data JSON must contain a metadata object", call. = FALSE)
  }
  if (is.null(x$data) || !is.list(x$data)) {
    stop("data JSON must contain a data object", call. = FALSE)
  }
  if (is.null(x$metadata$column_order)) {
    stop("metadata.column_order is required", call. = FALSE)
  }
  column_order <- data_json_as_character_vector(x$metadata$column_order)
  if (anyDuplicated(column_order)) {
    stop("metadata.column_order must not contain duplicate columns", call. = FALSE)
  }
  data_names <- names(x$data)
  if (!setequal(column_order, data_names)) {
    stop("metadata.column_order and data column names must match", call. = FALSE)
  }
  invisible(TRUE)
}

data_json_missing_vector <- function(values) {
  vapply(values, function(value) {
    is.null(value) || (length(value) == 1L && is.atomic(value) && is.na(value))
  }, logical(1))
}

data_json_value_key <- function(value) {
  if (is.null(value)) return("<NULL>")
  if (length(value) == 1L && is.atomic(value) && is.na(value)) return("<NA>")
  paste0(typeof(value), ":", as.character(value))
}

compare_data_json_summary <- function(found, expected) {
  errors <- character()
  for (field in c("row_count", "column_count", "missing_count")) {
    if (!identical(as.integer(found[[field]]), as.integer(expected[[field]]))) {
      errors <- c(errors, paste0("summary.", field, " mismatch"))
    }
  }

  found_columns <- found$columns
  expected_columns <- expected$columns
  if (is.null(found_columns) || !setequal(names(found_columns), names(expected_columns))) {
    errors <- c(errors, "summary.columns mismatch")
    return(errors)
  }

  for (column_name in names(expected_columns)) {
    for (field in c("missing_count", "unique_count")) {
      if (!identical(as.integer(found_columns[[column_name]][[field]]),
                     as.integer(expected_columns[[column_name]][[field]]))) {
        errors <- c(errors, paste0("summary.columns.", column_name, ".", field,
                                   " mismatch"))
      }
    }
  }
  errors
}

data_json_as_character_vector <- function(x) {
  as.character(unlist(x, use.names = FALSE))
}
