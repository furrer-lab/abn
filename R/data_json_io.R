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
#' @export
export_abnData <- function(data.df, data.dists, group.var = NULL, file = NULL,
                           pretty = TRUE, include_summary = TRUE) {
  validate_abn_data_inputs(data.df, data.dists, group.var)

  data.dists <- as.list(data.dists)
  column_order <- colnames(data.df)
  metadata <- list(
    schema_version = "bn-data-v1",
    issuer = "abn::export_abnData",
    orientation = "columnar",
    column_order = column_order,
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
        row_names = row.names(data.df),
        columns = data_json_r_adapter_columns(data.df)
      )
    )
  )

  object <- list(
    metadata = metadata,
    data = data_json_data_payload(data.df)
  )
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
#' @param validate If `TRUE`, validate schema shape and optional summary
#'   metadata before reconstruction.
#' @return A list with `data.df`, `data.dists`, `group.var`, and `metadata`.
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
  row_names <- object$metadata$adapters$r$row_names
  if (!is.null(row_names) && length(row_names) == nrow(data.df)) {
    row.names(data.df) <- as.character(row_names)
  }

  list(
    data.df = data.df,
    data.dists = object$metadata$data_dists,
    group.var = object$metadata$group_var,
    metadata = object$metadata
  )
}

validate_abn_data_inputs <- function(data.df, data.dists, group.var) {
  if (!is.data.frame(data.df)) {
    stop("data.df must be a data frame", call. = FALSE)
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
  validate_data_json_shape(object)
  if (!identical(object$metadata$schema_version, "bn-data-v1")) {
    stop("unsupported data JSON schema version", call. = FALSE)
  }
  if (!identical(object$metadata$orientation, "columnar")) {
    stop("data JSON orientation must be 'columnar'", call. = FALSE)
  }
  columns_by_name <- data_json_columns_by_name(object$metadata$columns)
  column_order <- data_json_as_character_vector(object$metadata$column_order)
  if (!identical(names(columns_by_name), column_order)) {
    stop("metadata.columns must match metadata.column_order", call. = FALSE)
  }
  group.var <- object$metadata$group_var
  data.dists <- object$metadata$data_dists
  node_columns <- setdiff(column_order, group.var)
  if (!identical(names(data.dists), node_columns)) {
    stop("metadata.data_dists names must match node columns in order",
         call. = FALSE)
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
  if (is.factor(values)) {
    descriptor$categories <- levels(values)
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
      class = class(values),
      typeof = typeof(values),
      levels = if (is.factor(values)) levels(values) else NULL,
      ordered = if (is.factor(values)) is.ordered(values) else NULL
    )
  }
  columns
}

data_json_data_payload <- function(data.df) {
  payload <- list()
  for (column_name in colnames(data.df)) {
    values <- data.df[[column_name]]
    if (is.factor(values)) values <- as.character(values)
    payload[[column_name]] <- lapply(values, data_json_scalar_value)
  }
  payload
}

data_json_scalar_value <- function(value) {
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

data_json_reconstruct_column <- function(values, column, r_column) {
  values <- lapply(values, data_json_null_to_na)
  r_class <- as.character(r_column$class %||% character())
  r_type <- as.character(r_column$typeof %||% "character")

  if ("ordered" %in% r_class) {
    return(ordered(unlist(values, use.names = FALSE),
                   levels = as.character(r_column$levels)))
  }
  if ("factor" %in% r_class) {
    levels <- r_column$levels %||% column$categories
    return(factor(unlist(values, use.names = FALSE),
                  levels = as.character(levels)))
  }
  if (identical(r_type, "integer")) {
    return(as.integer(unlist(values, use.names = FALSE)))
  }
  if (identical(r_type, "double")) {
    return(as.numeric(unlist(values, use.names = FALSE)))
  }
  if (identical(r_type, "logical")) {
    return(as.logical(unlist(values, use.names = FALSE)))
  }
  if (length(r_class) == 0L && identical(column$value_encoding, "labels")) {
    return(factor(unlist(values, use.names = FALSE),
                  levels = as.character(column$categories)))
  }
  as.character(unlist(values, use.names = FALSE))
}

data_json_null_to_na <- function(value) {
  if (is.null(value)) return(NA)
  value
}
