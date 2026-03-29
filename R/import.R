#' Import abnFit object from structured JSON format
#'
#' @description
#' Imports a fitted Additive Bayesian Network (ABN) model from a structured JSON
#' format, reconstructing an abnFit object suitable for further analysis.
#'
#' @param file Character string specifying a file path containing the JSON
#'    representation of the model. Alternatively, a JSON string can be provided
#'    directly via the `json` parameter.
#' @param json Optional character string containing the JSON representation
#'    of the model. If provided, `file` is ignored.
#' @param validate Logical, whether to validate the imported object against
#'    the abnFit class requirements. Default is TRUE.
#' @param ... Additional import options (currently unused).
#'
#' @return An object of class \code{abnFit} representing the imported model.
#'
#' @importFrom jsonlite fromJSON
#' @export
import_abnFit <- function(file = NULL, json = NULL, validate = TRUE, ...) {
  # Input validation
  if (is.null(file) && is.null(json)) {
    stop("Either 'file' or 'json' must be provided")
  }

  if (!is.null(file) && !is.null(json)) {
    warning("Both 'file' and 'json' provided; using 'json' parameter")
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for JSON import.")
  }

  # Load JSON
  if (!is.null(json)) {
    json_list <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  } else {
    if (!file.exists(file)) {
      stop(sprintf("File '%s' does not exist", file))
    }
    # Use warn = FALSE to suppress "incomplete final line" warnings from readLines
    json_content <- paste(readLines(file, warn = FALSE), collapse = "\n")
    json_list <- jsonlite::fromJSON(json_content, simplifyVector = FALSE)
  }

  # Validate JSON structure
  validate_json_structure(json_list)

  # Determine method (default to mle if not specified)
  method <- json_list$method %||% "mle"

  # Reconstruct abnFit object based on method
  if (method == "mle") {
    abn_fit <- reconstruct_abnfit_mle(json_list)
  } else if (method == "bayes") {
    abn_fit <- reconstruct_abnfit_bayes(json_list)
  } else {
    stop(sprintf("Unsupported method '%s'. Supported methods are 'mle' and 'bayes'.", method))
  }

  # Optional validation
  if (validate) {
    validate_abnfit_object(abn_fit)
  }

  return(abn_fit)
}

#' Validate JSON structure for abnFit import
#' @keywords internal
validate_json_structure <- function(json_list) {
  required_keys <- c("variables", "parameters", "arcs")
  if (!is.list(json_list) || !all(required_keys %in% names(json_list))) {
    stop("Invalid JSON structure: must contain 'variables', 'parameters', and 'arcs' components")
  }
  if (!is.null(json_list$linkFunctions)) {
    if (!is.list(json_list$linkFunctions)) {
      stop("Invalid JSON structure: 'linkFunctions' must be an array")
    }
  }
  has_lf_id <- any(vapply(json_list$parameters, function(p) {
    !is.null(p$link_function_id)
  }, logical(1)))
  if (has_lf_id && is.null(json_list$linkFunctions)) {
    stop("Invalid JSON structure: 'linkFunctions' is required when 'link_function_id' is used in parameters")
  }
  invisible(TRUE)
}

#' Reconstruct abnFit object for MLE method from JSON
#' @keywords internal
reconstruct_abnfit_mle <- function(json_list) {
  variables <- json_list$variables

  id_to_name <- stats::setNames(
    vapply(variables, function(x) as.character(x$attribute_name), character(1)),
    vapply(variables, function(x) as.character(x$variable_id), character(1))
  )
  variable_names <- unname(id_to_name)
  model_types <- vapply(variables, function(x) as.character(x$model_type), character(1))
  n <- length(variable_names)

  lf_lookup <- list()
  if (!is.null(json_list$linkFunctions) && is.list(json_list$linkFunctions)) {
    for (lf in json_list$linkFunctions) {
      lf_lookup[[as.character(lf$id)]] <- as.character(lf$name)
    }
  }

  data_dists <- stats::setNames(model_types, variable_names)

  data_list <- lapply(seq_along(variables), function(i) {
    dist <- model_types[i]
    if (dist == "gaussian") {
      numeric(0)
    } else if (dist == "binomial") {
      factor(levels = c("0", "1"))
    } else if (dist == "poisson") {
      integer(0)
    } else if (dist == "multinomial") {
      states <- variables[[i]]$states
      if (is.null(states) || length(states) == 0) {
        state_names <- c("1", "2")
      } else {
        state_names <- vapply(states, function(x) as.character(x$value_name), character(1))
      }
      factor(levels = state_names)
    } else {
      numeric(0)
    }
  })
  data_df <- as.data.frame(stats::setNames(data_list, variable_names))

  dag_matrix <- matrix(0, nrow = n, ncol = n,
                       dimnames = list(variable_names, variable_names))

  if (!is.null(json_list$arcs) && length(json_list$arcs) > 0) {
    for (arc in json_list$arcs) {
      src_raw <- as.character(arc$source_variable_id)
      tgt_raw <- as.character(arc$target_variable_id)
      src <- id_to_name[src_raw]
      if (length(src) == 0 || is.na(src)) src <- src_raw
      tgt <- id_to_name[tgt_raw]
      if (length(tgt) == 0 || is.na(tgt)) tgt <- tgt_raw
      if (src %in% variable_names && tgt %in% variable_names) {
        dag_matrix[src, tgt] <- 1
      }
    }
  }

  abnDag <- list(
    dag = dag_matrix,
    data.df = data_df,
    data.dists = data_dists
  )
  class(abnDag) <- "abnDag"

  multinomial.states <- list()
  for (i in seq_along(variables)) {
    if (model_types[i] == "multinomial") {
      states <- variables[[i]]$states
      if (!is.null(states) && length(states) > 0) {
        var_name <- variable_names[i]
        multinomial.states[[var_name]] <- lapply(states, function(s) {
          list(
            state_id = as.character(s$state_id),
            value_name = as.character(s$value_name),
            is_baseline = isTRUE(s$is_baseline)
          )
        })
      }
    }
  }

  coef_list <- stats::setNames(vector("list", n), variable_names)
  stderror_list <- stats::setNames(vector("list", n), variable_names)

  if (!is.null(json_list$parameters)) {
    for (param in json_list$parameters) {
      raw_target_id <- as.character(param$source$variable_id)
      target_name <- id_to_name[raw_target_id]
      if (length(target_name) == 0 || is.na(target_name)) target_name <- raw_target_id
      if (!(target_name %in% variable_names)) next

      if (!is.null(param$link_function_name)) {
        lf_name <- as.character(param$link_function_name)
      } else if (!is.null(param$link_function_id)) {
        lf_name <- lf_lookup[[as.character(param$link_function_id)]]
        if (is.null(lf_name)) lf_name <- "unknown"
      } else {
        lf_name <- "unknown"
      }

      coeffs_input <- if (is.list(param$coefficients)) param$coefficients else list()
      vals <- numeric()
      ses <- numeric()
      c_names <- character()

      for (coeff in coeffs_input) {
        type <- as.character(coeff$condition_type)
        if (type == "intercept") {
          c_name <- paste0(target_name, "|intercept")
        } else if (type == "linear_term") {
          parent_raw <- if (!is.null(coeff$conditions) && length(coeff$conditions) > 0) {
            as.character(coeff$conditions[[1]]$parent_variable_id)
          } else {
            "unknown"
          }
          parent_name <- id_to_name[parent_raw]
          if (length(parent_name) == 0 || is.na(parent_name)) parent_name <- parent_raw
          c_name <- paste0(target_name, "|", parent_name)
        } else {
          c_name <- paste0(target_name, "|", type)
        }

        vals <- c(vals, as.numeric(coeff$value))
        ses <- c(ses, as.numeric(coeff$stderr %||% NA))
        c_names <- c(c_names, c_name)
      }

      coef_list[[target_name]] <- matrix(vals, nrow = 1, dimnames = list(NULL, c_names))
      stderror_list[[target_name]] <- matrix(ses, nrow = 1, dimnames = list(NULL, c_names))
    }
  }

  for (var_id in variable_names) {
    if (is.null(coef_list[[var_id]])) {
      coef_list[[var_id]] <- matrix(numeric(0), nrow = 0, ncol = 0)
      stderror_list[[var_id]] <- matrix(numeric(0), nrow = 0, ncol = 0)
    }
  }

  abn_fit <- list(
    abnDag = abnDag,
    coef = coef_list,
    Stderror = stderror_list,
    method = "mle",
    multinomial.states = multinomial.states,
    scenario_id = json_list$scenario_id %||% json_list$scenarioId %||% NULL,
    label = json_list$label %||% NULL,
    call = match.call()
  )

  class(abn_fit) <- "abnFit"
  return(abn_fit)
}

#' Reconstruct abnFit object for Bayesian method from JSON
#' @keywords internal
reconstruct_abnfit_bayes <- function(json_list) {
  warning("Bayesian model import is not fully implemented yet. Returning MLE-structured object.")
  abn_fit <- reconstruct_abnfit_mle(json_list)
  abn_fit$method <- "bayes"
  return(abn_fit)
}

#' Validate abnFit object meets class requirements
#' @keywords internal
validate_abnfit_object <- function(object) {
  if (!inherits(object, "abnFit")) {
    stop("Imported object is not of class 'abnFit'")
  }
  if (!inherits(object$abnDag, "abnDag")) {
    stop("abnDag component is not of class 'abnDag'")
  }
  invisible(TRUE)
}

#' Default null operator
#' @keywords internal
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}
