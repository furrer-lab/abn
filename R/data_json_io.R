#' Export raw ABN data to JSON
#'
#' Serializes a data frame plus ABN distribution metadata into the data JSON
#' structure. The payload is column-oriented: `metadata` describes how to
#' interpret the values, while `data` contains only scalar values.
#'
#' @param data.df A data frame containing node columns and, optionally, a group
#'   column.
#' @param data.dists A named list or vector giving ABN distributions for node
#'   columns. Names must match non-group node columns in `data.df` and in the
#'   same order.
#' @param group.var Optional name of the grouping column in `data.df`.
#' @param file Optional path to write the JSON document. If `NULL`, the JSON
#'   string is returned.
#' @param pretty If `TRUE`, pretty-print the JSON output.
#' @param include_summary If `TRUE`, include recomputed `metadata.summary`.
#' @return A JSON string, or invisibly the output file path when `file` is used.
#' @details
#' The columnar `bn-data` format is specified by
#' `inst/schemas/bn-data.schema.json` and described in
#' `vignette("bn-data-json-specification", package = "abn")`.
#' Zero-row datasets are supported; zero-column datasets are rejected.
#' Raw categorical samples may have one observed category. Declared factor
#' levels, including unused levels, are preserved; unknown category domains in
#' empty or all-missing label columns are represented by an empty array.
#' Collection fields remain arrays even when they contain zero or one element.
#'
#' Missing observations are JSON `null`; numeric `NaN`, `Inf` and `-Inf` use
#' the strings `"NaN"`, `"Infinity"` and `"-Infinity"`. The abn adapter reserves
#' these strings and rejects them as ordinary category labels.
#' R classes, storage types and row names are recorded in `metadata.adapters.r`.
#' Row labels are always strings; `row_names_type` explicitly distinguishes
#' character, integer and automatic R row names without guessing from labels.
#' @export
export_abnData <- function(data.df, data.dists, group.var = NULL, file = NULL,
                           pretty = TRUE, include_summary = TRUE) {
  validate_abn_data_inputs(data.df, data.dists, group.var)

  data.dists <- as.list(data.dists)
  column_order <- colnames(data.df)
  metadata <- list(
    schema_version = "bn-data",
    issuer = "abn::export_abnData",
    orientation = "columnar",
    column_order = as.list(column_order),
    missing_values = list(
      missing = NULL,
      not_a_number = "NaN",
      positive_infinity = "Infinity",
      negative_infinity = "-Infinity"
    ),
    columns = lapply(column_order, function(column_name) {
      data_json_column_descriptor(
        column_name = column_name,
        values = data.df[[column_name]],
        distribution = data.dists[[column_name]],
        is_group = identical(column_name, group.var)
      )
    }),
    data_dists = data.dists,
    group_var = group.var,
    adapters = list(
      r = list(
        row_names = as.list(row.names(data.df)),
        row_names_type = if (.row_names_info(data.df, type = 1L) < 0L) {
          "automatic"
        } else {
          typeof(attr(data.df, "row.names"))
        },
        columns = data_json_r_adapter_columns(data.df)
      )
    )
  )

  object <- list(
    metadata = metadata,
    data = data_json_data_payload(data.df, data.dists, group.var)
  )
  validate_abn_data_json_object(object)
  if (isTRUE(include_summary)) {
    object$metadata$summary <- compute_data_json_summary(object)
  }

  json <- jsonlite::toJSON(object, auto_unbox = TRUE, pretty = pretty,
                           null = "null", digits = NA)
  if (!is.null(file)) {
    writeLines(json, con = file)
    return(invisible(file))
  }
  json
}

#' Import raw ABN data from JSON
#'
#' Reconstructs a data frame and associated ABN metadata from a data JSON
#' document produced by [export_abnData()].
#'
#' @param file Optional path to a JSON file.
#' @param json Optional JSON string. If both `file` and `json` are supplied,
#'   `json` is used.
#' @param validate If `TRUE`, check document structure, column descriptors,
#'   distribution/encoding consistency, category membership, values, nullability,
#'   and optional summary metadata before reconstruction. Runtime checks do not
#'   require the optional `jsonvalidate` package.
#' @return A list with `data.df`, `data.dists`, `group.var`, and `metadata`.
#' @details
#' Numeric, integer, boolean and categorical columns are reconstructed using
#' generic column descriptors. Ordered and unused categories are retained.
#' An optional `metadata.adapters.r` block restores R-specific column storage
#' and row-name types. Without it, label columns become factors using their
#' declared categories, and custom row names are unavailable, so R generates
#' default row names.
#'
#' Supplied row labels remain strings unless `metadata.adapters.r.row_names_type`
#' explicitly requests integer or automatic R row names. For legacy adapters
#' without this marker, labels such as `"001"` are preserved as strings.
#' The abn adapter supports the documented `NaN` and infinity tokens for number
#' columns and reserves them against use as ordinary labels.
#' Zero-row and all-missing categorical datasets are supported; zero-column
#' datasets are rejected. See
#' `vignette("bn-data-json-specification", package = "abn")` for the format.
#' @export
import_abnData <- function(file = NULL, json = NULL, validate = TRUE) {
  if (is.null(file) && is.null(json)) {
    stop("Either 'file' or 'json' must be provided", call. = FALSE)
  }
  if (!is.null(json)) {
    object <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  } else {
    if (!file.exists(file)) {
      stop(sprintf("File '%s' does not exist", file), call. = FALSE)
    }
    json <- paste(readLines(file, warn = FALSE), collapse = "\n")
    object <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  }

  if (isTRUE(validate)) {
    validate_abn_data_json_object(object)
    validate_data_json_summary(object, on_error = "error")
  }

  column_order <- data_json_as_character_vector(object$metadata$column_order)
  columns_by_name <- data_json_columns_by_name(object$metadata$columns)
  r_columns <- object$metadata$adapters$r$columns %||% list()
  data_list <- stats::setNames(vector("list", length(column_order)), column_order)
  for (column_name in column_order) {
    data_list[[column_name]] <- data_json_reconstruct_column(
      values = object$data[[column_name]],
      column = columns_by_name[[column_name]],
      r_column = r_columns[[column_name]]
    )
  }

  data.df <- as.data.frame(data_list, stringsAsFactors = FALSE,
                           optional = TRUE, check.names = FALSE)
  data.df <- data_json_restore_row_names(data.df, object$metadata$adapters$r)

  list(
    data.df = data.df,
    data.dists = object$metadata$data_dists,
    group.var = object$metadata$group_var,
    metadata = object$metadata
  )
}

data_json_restore_row_names <- function(data.df, r_adapter) {
  row_names <- r_adapter$row_names
  if (is.null(row_names)) return(data.df)
  # Labels are always strings. Old adapters without an explicit type marker
  # retain character row names, even if every label looks like a number.
  labels <- vapply(as.list(row_names), function(label) {
    if (!is.character(label) || length(label) != 1L || is.na(label)) {
      stop("metadata.adapters.r.row_names must contain string labels", call. = FALSE)
    }
    label
  }, character(1), USE.NAMES = FALSE)
  if (length(labels) != nrow(data.df) || anyDuplicated(labels)) {
    stop("metadata.adapters.r.row_names must be unique and match the row count",
         call. = FALSE)
  }
  type <- r_adapter$row_names_type %||% "character"
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("automatic", "integer", "character")) {
    stop("metadata.adapters.r.row_names_type must be automatic, integer or character",
         call. = FALSE)
  }
  if (identical(type, "automatic")) {
    if (!identical(labels, as.character(seq_len(nrow(data.df))))) {
      stop("Automatic row names must match the default row sequence", call. = FALSE)
    }
    row.names(data.df) <- NULL
  } else if (identical(type, "integer")) {
    integers <- suppressWarnings(as.integer(labels))
    if (anyNA(integers) || !identical(as.character(integers), labels)) {
      stop("Integer row names must contain exact, representable integer labels",
           call. = FALSE)
    }
    row.names(data.df) <- integers
  } else {
    row.names(data.df) <- labels
  }
  data.df
}

validate_abn_data_inputs <- function(data.df, data.dists, group.var) {
  if (!is.data.frame(data.df)) {
    stop("data.df must be a data frame", call. = FALSE)
  }
  if (ncol(data.df) == 0L) {
    stop("data.df must contain at least one column; zero-row data are supported",
         call. = FALSE)
  }
  if (is.null(data.dists) || is.null(names(data.dists))) {
    stop("data.dists must be a named list or vector", call. = FALSE)
  }
  data.dists <- as.list(data.dists)
  valid_distributions <- c("binomial", "multinomial", "poisson", "gaussian")
  if (any(!unlist(data.dists) %in% valid_distributions)) {
    stop("data.dists contains unsupported distributions", call. = FALSE)
  }
  if (!is.null(group.var)) {
    if (!is.character(group.var) || length(group.var) != 1L || is.na(group.var)) {
      stop("group.var must be a single column name", call. = FALSE)
    }
    if (!group.var %in% colnames(data.df)) {
      stop("group.var must be a column in data.df", call. = FALSE)
    }
  }

  node_columns <- setdiff(colnames(data.df), group.var)
  if (!identical(names(data.dists), node_columns)) {
    stop("data.dists names must match non-group data.df columns in order",
         call. = FALSE)
  }
  invisible(TRUE)
}

validate_abn_data_json_object <- function(object) {
  abn_json_assert_object(object, "data document", c("metadata", "data"))
  abn_json_assert_object(object$metadata, "metadata",
                          c("schema_version", "orientation", "column_order", "columns", "data_dists"))
  abn_json_assert_object(object$data, "data")
  abn_json_assert_array(object$metadata$column_order, "metadata.column_order")
  for (name in object$metadata$column_order) {
    abn_json_assert_string(name, "metadata.column_order entry", nonempty = TRUE)
  }
  abn_json_assert_array(object$metadata$columns, "metadata.columns")
  for (column in object$metadata$columns) {
    abn_json_assert_object(column, "column descriptor",
                            c("name", "role", "semantic_type", "value_encoding", "nullable"))
    for (field in c("name", "role", "semantic_type", "value_encoding")) {
      abn_json_assert_string(column[[field]], paste0("column descriptor ", field), nonempty = TRUE)
    }
    abn_json_assert_boolean(column$nullable, "column descriptor nullable")
    if (!is.null(column$ordered)) abn_json_assert_boolean(column$ordered, "column descriptor ordered")
    if (!is.null(column$model_distribution)) {
      abn_json_assert_string(column$model_distribution, "column descriptor model_distribution")
    }
    if (!is.null(column$categories)) {
      abn_json_assert_array(column$categories, "column descriptor categories")
      for (category in column$categories) abn_json_assert_string(category, "category")
    }
  }
  abn_json_assert_object(object$metadata$data_dists, "metadata.data_dists")
  for (distribution in object$metadata$data_dists) {
    abn_json_assert_string(distribution, "metadata.data_dists entry", nonempty = TRUE)
  }
  for (column in object$data) abn_json_assert_array(column, "data column")
  if (!is.null(object$metadata$adapters)) {
    abn_json_assert_object(object$metadata$adapters, "metadata.adapters")
    if (!is.null(object$metadata$adapters$r)) {
      abn_json_assert_object(object$metadata$adapters$r, "metadata.adapters.r")
      r_columns <- object$metadata$adapters$r$columns
      if (!is.null(r_columns)) {
        abn_json_assert_object(r_columns, "metadata.adapters.r.columns")
        for (r_column in r_columns) abn_json_assert_object(r_column, "R column adapter")
      }
    }
  }
  validate_data_json_shape(object)
  if (!identical(object$metadata$schema_version, "bn-data")) {
    stop("unsupported data JSON schema version", call. = FALSE)
  }
  if (!identical(object$metadata$orientation, "columnar")) {
    stop("data JSON orientation must be 'columnar'", call. = FALSE)
  }
  validate_data_json_document_semantics(object)
  invisible(TRUE)
}

validate_data_json_type <- function(value, what, kind) {
  ok <- switch(kind,
               list = is.list(value),
               character = is.character(value),
               logical = is.logical(value),
               numeric = is.numeric(value),
               FALSE)
  if (!isTRUE(ok)) stop(what, " must be JSON ", kind, call. = FALSE)
  invisible(TRUE)
}

validate_data_json_document_semantics <- function(object) {
  validate_data_json_type(object$metadata, "metadata", "list")
  if (!(is.list(object$metadata$column_order) ||
        is.character(object$metadata$column_order))) {
    stop("metadata.column_order must be a JSON array or character vector", call. = FALSE)
  }
  column_order <- data_json_as_character_vector(object$metadata$column_order)
  if (any(!nzchar(column_order))) {
    stop("metadata.column_order must not contain empty column names", call. = FALSE)
  }

  columns_by_name <- data_json_columns_by_name(object$metadata$columns)
  column_order <- data_json_as_character_vector(object$metadata$column_order)
  if (!identical(names(columns_by_name), column_order)) {
    stop("metadata.columns must match metadata.column_order", call. = FALSE)
  }
  special_tokens <- c("NaN", "Infinity", "-Infinity")

  for (column_name in column_order) {
    column <- columns_by_name[[column_name]]
    if (!is.list(column)) {
      stop("metadata.columns entries must be JSON objects", call. = FALSE)
    }
    if (!identical(as.character(column$name), column_name)) {
      stop("metadata column name does not match its column_order entry", call. = FALSE)
    }
    if (!column$role %in% c("node", "group")) {
      stop("metadata.columns entries need role 'node' or 'group'", call. = FALSE)
    }
    if (!column$semantic_type %in% c("continuous", "count", "binary", "categorical", "unknown")) {
      stop("metadata column semantic_type is unsupported", call. = FALSE)
    }
    if (!column$value_encoding %in% c("labels", "integer", "number", "boolean", "string")) {
      stop("metadata column value_encoding is unsupported", call. = FALSE)
    }
    if (!is.logical(column$nullable) || length(column$nullable) != 1L) {
      stop("metadata column nullable must be a JSON boolean", call. = FALSE)
    }
    if (!is.null(column$ordered) &&
        (!is.logical(column$ordered) || length(column$ordered) != 1L)) {
      stop("metadata column ordered must be a JSON boolean", call. = FALSE)
    }
    if (identical(column$value_encoding, "labels") && is.null(column$categories)) {
      stop("label-encoded columns must declare categories", call. = FALSE)
    }
    if (!is.null(column$categories)) {
      categories <- as.character(unlist(column$categories, use.names = FALSE))
      if (anyNA(categories) || anyDuplicated(categories)) {
        stop("metadata column categories must contain unique string labels", call. = FALSE)
      }
      if (any(categories %in% special_tokens)) {
        stop("metadata column categories must not use the special numeric tokens",
             call. = FALSE)
      }
    }
  }

  missing_values <- object$metadata$missing_values
  if (!is.null(missing_values)) {
    if (!is.list(missing_values) || !is.null(missing_values$missing)) {
      stop("metadata.missing_values.missing must be JSON null", call. = FALSE)
    }
    tokens <- vapply(c("not_a_number", "positive_infinity", "negative_infinity"),
                     function(token) {
      value <- missing_values[[token]]
      if (!is.character(value) || length(value) != 1L) {
        stop("metadata.missing_values.", token, " must be a scalar string",
             call. = FALSE)
      }
      value
    }, character(1))
    if (!identical(unname(tokens), special_tokens) || anyDuplicated(tokens)) {
      stop("metadata.missing_values must declare the documented NaN and infinity tokens",
           call. = FALSE)
    }
  }

  group_var <- object$metadata$group_var
  data.dists <- object$metadata$data_dists
  if (!is.list(data.dists)) {
    stop("metadata.data_dists must be a JSON object", call. = FALSE)
  }
  if (!is.null(group_var) &&
      (!is.character(group_var) || length(group_var) != 1L || !nzchar(group_var))) {
    stop("metadata.group_var must be a single column name or null", call. = FALSE)
  }
  node_columns <- setdiff(column_order, group_var)
  if (!identical(names(data.dists), node_columns)) {
    stop("metadata.data_dists names must match node columns in order", call. = FALSE)
  }
  dists <- unlist(data.dists, use.names = FALSE)
  if (any(!dists %in% c("binomial", "multinomial", "poisson", "gaussian"))) {
    stop("metadata.data_dists contains unsupported distributions", call. = FALSE)
  }
  role_columns <- vapply(columns_by_name, function(column) column$role, character(1))
  if (!is.null(group_var)) {
    if (!group_var %in% column_order) {
      stop("metadata.group_var must name a column", call. = FALSE)
    }
    if (!identical(role_columns[[group_var]], "group")) {
      stop("metadata.group_var must have role 'group'", call. = FALSE)
    }
    if (any(role_columns[setdiff(names(role_columns), group_var)] == "group")) {
      stop("metadata.columns contains a group column not named by metadata.group_var",
           call. = FALSE)
    }
  } else if (any(role_columns == "group")) {
    stop("group columns require metadata.group_var", call. = FALSE)
  }

  expected_semantics <- c(gaussian = "continuous", poisson = "count",
                          binomial = "binary", multinomial = "categorical")
  expected_encodings <- c(gaussian = "number", poisson = "integer",
                          binomial = "labels", multinomial = "labels")
  for (column_name in column_order) {
    column <- columns_by_name[[column_name]]
    if (identical(column$role, "group")) {
      if (!identical(column$semantic_type, "categorical") ||
          !identical(column$value_encoding, "labels")) {
        stop("group columns must be categorical label columns", call. = FALSE)
      }
      if (!is.null(column$model_distribution)) {
        stop("group columns must not declare model_distribution", call. = FALSE)
      }
      next
    }
    distribution <- data.dists[[column_name]]
    if (distribution %in% c("binomial", "multinomial")) {
      n_categories <- length(as.character(column$categories))
      if (distribution == "binomial" && n_categories > 2L) {
        stop("binomial columns must declare at most two categories", call. = FALSE)
      }
    }
    if (!identical(column$model_distribution, distribution)) {
      stop("metadata.columns model_distribution for '", column_name,
           "' must match metadata.data_dists", call. = FALSE)
    }
    if (!identical(column$semantic_type, expected_semantics[[distribution]])) {
      stop("metadata.columns semantic_type for '", column_name,
           "' is incompatible with model_distribution", call. = FALSE)
    }
    if (!identical(column$value_encoding, expected_encodings[[distribution]])) {
      stop("metadata.columns value_encoding for '", column_name,
           "' is incompatible with model_distribution", call. = FALSE)
    }
    if (isTRUE(column$ordered) && !identical(column$value_encoding, "labels")) {
      stop("ordered is only supported for label-encoded columns", call. = FALSE)
    }
  }

  for (column_name in column_order) {
    values <- object$data[[column_name]]
    missing <- data_json_missing_vector(values)
    if (any(missing) && !isTRUE(columns_by_name[[column_name]]$nullable)) {
      stop("column '", column_name, "' contains missing values but is not nullable",
           call. = FALSE)
    }
    validate_data_json_values(values, columns_by_name[[column_name]], column_name)
  }

  adapters <- object$metadata$adapters
  if (!is.null(adapters)) {
    if (!is.list(adapters)) {
      stop("metadata.adapters must be a JSON object", call. = FALSE)
    }
    r_columns <- adapters$r$columns
    if (!is.null(r_columns)) {
      if (!is.list(r_columns)) {
        stop("metadata.adapters.r.columns must be a JSON object", call. = FALSE)
      }
      for (column_name in names(r_columns)) {
        if (!column_name %in% column_order) {
          stop("metadata.adapters.r.columns references unknown column '", column_name,
               "'", call. = FALSE)
        }
        r_column <- r_columns[[column_name]]
        column <- columns_by_name[[column_name]]
        r_levels <- as.character(r_column$levels %||% character())
        if (length(r_levels) > 0L && !is.null(column$categories) &&
            !identical(r_levels, as.character(column$categories))) {
          stop("metadata.adapters.r.columns levels disagree with generic categories for '",
               column_name, "'", call. = FALSE)
        }
        if (!is.null(r_column$ordered) && !is.null(column$ordered) &&
            !identical(isTRUE(r_column$ordered), isTRUE(column$ordered))) {
          stop("metadata.adapters.r.columns ordered disagrees with generic metadata for '",
               column_name, "'", call. = FALSE)
        }
        r_class <- as.character(r_column$class %||% character())
        r_type <- as.character(r_column$typeof %||% "")
        if (identical(column$value_encoding, "labels") &&
            !("factor" %in% r_class || "ordered" %in% r_class)) {
          if (r_type %in% c("integer", "double")) {
            values <- suppressWarnings(as.numeric(as.character(column$categories)))
            if (anyNA(values) || any(!is.finite(values))) {
              stop("metadata.adapters.r.columns labels are not numeric for '",
                   column_name, "'", call. = FALSE)
            }
            if (identical(r_type, "integer") &&
                any(abs(values - round(values)) > .Machine$double.eps^0.5 |
                    values < -2147483648 |
                    values > 2147483647)) {
              stop("metadata.adapters.r.columns labels are not integer-compatible for '",
                   column_name, "'", call. = FALSE)
            }
          } else if (identical(r_type, "logical") &&
                     !identical(sort(unique(as.character(column$categories))), c("FALSE", "TRUE")) &&
                     !all(as.character(column$categories) %in% c("FALSE", "TRUE"))) {
            stop("metadata.adapters.r.columns labels are not logical-compatible for '",
                 column_name, "'", call. = FALSE)
          }
        }
      }
    }
  }
  invisible(TRUE)
}

validate_data_json_values <- function(values, column, column_name) {
  if (!is.list(values)) {
    values <- as.list(values)
  }
  encoding <- column$value_encoding
  categories <- as.character(column$categories %||% character())
  for (raw_value in values) {
    value <- raw_value
    if (is.null(value)) next
    if (!is.atomic(value) || length(value) != 1L || is.na(value)) {
      stop("column '", column_name,
           "' must contain scalar, non-NA values or JSON null", call. = FALSE)
    }
    if (encoding == "number") {
      if (is.numeric(value)) {
        if (!is.finite(value)) {
          stop("column '", column_name,
               "' uses a raw non-finite number; use NaN or Infinity tokens instead",
               call. = FALSE)
        }
      } else if (!is.character(value) || !value %in% c("NaN", "Infinity", "-Infinity")) {
        stop("column '", column_name,
             "' must contain numbers or the documented special tokens", call. = FALSE)
      }
    } else if (encoding == "integer") {
      if (!is.numeric(value) || !is.finite(value) ||
          abs(value - round(value)) > .Machine$double.eps^0.5 ||
          value < 0 || value > .Machine$integer.max) {
        stop("column '", column_name, "' must contain non-negative integers", call. = FALSE)
      }
    } else if (encoding == "boolean") {
      if (!is.logical(value)) {
        stop("column '", column_name, "' must contain booleans", call. = FALSE)
      }
    } else if (encoding == "labels") {
      if (!is.character(value)) {
        stop("column '", column_name, "' must contain label strings", call. = FALSE)
      }
      if (!value %in% categories) {
        stop("column '", column_name, "' contains undeclared category '", value,
             "'", call. = FALSE)
      }
    } else if (encoding == "string") {
      if (!is.character(value)) {
        stop("column '", column_name, "' must contain strings", call. = FALSE)
      }
      if (value %in% c("NaN", "Infinity", "-Infinity")) {
        stop("column '", column_name,
             "' uses a reserved special token as an ordinary string", call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}


data_json_column_descriptor <- function(column_name, values, distribution, is_group) {
  descriptor <- list(
    name = column_name,
    role = if (is_group) "group" else "node",
    semantic_type = data_json_semantic_type(distribution, values, is_group),
    model_distribution = if (is_group) NULL else distribution,
    value_encoding = data_json_value_encoding(distribution, values, is_group),
    nullable = any(is.na(values))
  )
  encoding <- data_json_value_encoding(distribution, values, is_group)
  if (identical(encoding, "labels")) {
    labels <- if (is.factor(values)) levels(values) else if (is.logical(values)) {
      c("FALSE", "TRUE")
    } else {
      unique(as.character(values[!is.na(values)]))
    }
    descriptor$categories <- as.list(labels)
    descriptor$ordered <- isTRUE(is.ordered(values))
  } else if (is.factor(values)) {
    descriptor$categories <- as.list(levels(values))
    descriptor$ordered <- is.ordered(values)
  }
  descriptor
}

data_json_semantic_type <- function(distribution, values, is_group) {
  if (is_group) return("categorical")
  switch(distribution,
         gaussian = "continuous",
         poisson = "count",
         binomial = "binary",
         multinomial = "categorical",
         if (is.factor(values)) "categorical" else "unknown")
}

data_json_value_encoding <- function(distribution, values, is_group) {
  if (is.factor(values) || is_group || distribution %in% c("binomial", "multinomial")) {
    return("labels")
  }
  if (distribution == "poisson" || is.integer(values)) return("integer")
  if (is.logical(values)) return("boolean")
  if (is.numeric(values)) return("number")
  "string"
}

data_json_r_adapter_columns <- function(data.df) {
  columns <- list()
  for (column_name in colnames(data.df)) {
    values <- data.df[[column_name]]
    columns[[column_name]] <- list(
      class = as.list(class(values)),
      typeof = typeof(values),
      levels = if (is.factor(values)) as.list(levels(values)) else NULL,
      ordered = if (is.factor(values)) is.ordered(values) else NULL
    )
  }
  columns
}

data_json_data_payload <- function(data.df, data.dists = NULL, group.var = NULL) {
  payload <- list()
  data.dists <- as.list(data.dists %||% stats::setNames(vector("list", ncol(data.df)), colnames(data.df)))
  for (column_name in colnames(data.df)) {
    values <- data.df[[column_name]]
    is_group <- identical(column_name, group.var)
    encoding <- data_json_value_encoding(data.dists[[column_name]], values, is_group)
    if (identical(encoding, "labels")) values <- as.character(values)
    payload[[column_name]] <- lapply(values, data_json_scalar_value)
  }
  payload
}

data_json_scalar_value <- function(value) {
  if (is.numeric(value) && is.nan(value)) return("NaN")
  if (is.numeric(value) && identical(value, Inf)) return("Infinity")
  if (is.numeric(value) && identical(value, -Inf)) return("-Infinity")
  if (length(value) == 0L || is.na(value)) return(NULL)
  unname(value)
}

data_json_columns_by_name <- function(columns) {
  if (is.null(columns) || !is.list(columns)) {
    stop("metadata.columns must be a list", call. = FALSE)
  }
  names <- vapply(columns, function(column) as.character(column$name), character(1))
  if (anyDuplicated(names)) {
    stop("metadata.columns must not contain duplicate names", call. = FALSE)
  }
  stats::setNames(columns, names)
}

data_json_decode_scalar <- function(value, encoding) {
  if (is.null(value)) {
    return(switch(encoding,
                  number = NA_real_,
                  integer = NA_integer_,
                  boolean = NA,
                  NA_character_))
  }
  switch(encoding,
         number = if (identical(value, "NaN")) NaN else
           if (identical(value, "Infinity")) Inf else
             if (identical(value, "-Infinity")) -Inf else as.numeric(value),
         integer = as.integer(value),
         boolean = as.logical(value),
         labels = as.character(value),
         string = as.character(value),
         as.character(value))
}

data_json_reconstruct_column <- function(values, column, r_column) {
  values <- unname(values)
  r_class <- as.character(r_column$class %||% character())
  r_type <- as.character(r_column$typeof %||% "character")
  encoding <- column$value_encoding

  if (encoding == "labels") {
    labels <- vapply(values, data_json_decode_scalar, character(1), encoding = encoding)
    levels <- as.character(column$categories)
    if ("ordered" %in% r_class || isTRUE(column$ordered) || isTRUE(r_column$ordered)) {
      return(ordered(labels, levels = levels))
    }
    if ("factor" %in% r_class) {
      return(factor(labels, levels = levels))
    }
    if (length(r_class) > 0L && identical(r_class, "character")) {
      return(unname(labels))
    }
    if (identical(r_type, "integer")) {
      out <- suppressWarnings(as.integer(labels))
      out[is.na(labels)] <- NA_integer_
      return(out)
    }
    if (identical(r_type, "double")) {
      out <- suppressWarnings(as.numeric(labels))
      out[is.na(labels)] <- NA_real_
      return(out)
    }
    if (identical(r_type, "logical")) {
      out <- ifelse(is.na(labels), NA, labels == "TRUE")
      return(out)
    }
    return(factor(labels, levels = levels))
  }
  if (encoding == "number") {
    decoded <- vapply(values, data_json_decode_scalar, numeric(1), encoding = encoding)
    return(decoded)
  }
  if (encoding == "integer") {
    out <- vapply(values, data_json_decode_scalar, integer(1), encoding = encoding)
    if (identical(r_type, "double")) return(as.numeric(out))
    return(out)
  }
  if (encoding == "boolean") {
    return(vapply(values, data_json_decode_scalar, logical(1), encoding = encoding))
  }
  if (identical(r_type, "integer")) {
    return(as.integer(vapply(values, data_json_decode_scalar, integer(1),
                             encoding = "integer")))
  }
  if (identical(r_type, "double")) {
    return(as.numeric(vapply(values, data_json_decode_scalar, numeric(1),
                             encoding = "number")))
  }
  if (identical(r_type, "logical")) {
    return(as.logical(vapply(values, data_json_decode_scalar, logical(1),
                             encoding = "boolean")))
  }
  vapply(values, data_json_decode_scalar, character(1), encoding = encoding)
}

data_json_null_to_na <- function(value) {
  if (is.null(value)) return(NA)
  value
}
