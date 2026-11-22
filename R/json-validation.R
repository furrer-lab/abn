# Small structural checks shared by the JSON adapters. Never coerce first:
# coercion can turn a malformed scalar or object into apparently valid data.
abn_json_assert_object <- function(value, where, required = character(0)) {
  if (!is.list(value) || (length(value) > 0L && is.null(names(value)))) {
    stop(where, " must be a JSON object", call. = FALSE)
  }
  if (anyDuplicated(names(value))) {
    stop(where, " contains duplicate keys", call. = FALSE)
  }
  missing <- setdiff(required, names(value))
  if (length(missing)) {
    stop(where, " is missing required field '", missing[1L], "'", call. = FALSE)
  }
  invisible(TRUE)
}

abn_json_assert_array <- function(value, where) {
  if (!is.list(value) || (length(value) > 0L && !is.null(names(value)))) {
    stop(where, " must be a JSON array", call. = FALSE)
  }
  invisible(TRUE)
}

abn_json_assert_string <- function(value, where, nonempty = FALSE) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      (nonempty && !nzchar(value))) {
    stop(where, " must be a scalar string", call. = FALSE)
  }
  invisible(TRUE)
}

abn_json_assert_boolean <- function(value, where) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(where, " must be a scalar boolean", call. = FALSE)
  }
  invisible(TRUE)
}
