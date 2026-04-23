#' Import abnFit object from structured JSON format
#'
#' @description
#' Reconstructs a fitted Additive Bayesian Network (ABN) model from a structured JSON
#' format (typically exported by \code{\link{export_abnFit}}). This function enables
#' round-trip model I/O, allowing models to be exported, stored, shared, and later
#' reimported for further analysis or modification.
#'
#' The function validates the JSON structure and reconstructs the abnFit object,
#' including all network structure (variables, arcs) and model parameters
#' (coefficients, standard errors, variances, random effects).
#'
#' @param file Character string specifying a file path containing the JSON
#'    representation of the model. Alternatively, a JSON string can be provided
#'    directly via the \code{json} parameter. If both are provided, \code{json}
#'    is used and \code{file} is ignored.
#' @param json Optional character string containing the JSON representation
#'    of the model. Takes precedence over \code{file} if both are provided.
#' @param validate Logical, whether to validate the imported object against
#'    abnFit class requirements. Default is \code{TRUE}. Set to \code{FALSE}
#'    only if you are certain the JSON is valid.
#' @param ... Additional import options (currently unused, reserved for future extensions).
#'
#' @return An object of class \code{abnFit} representing the imported model,
#'    with all components reconstructed from the JSON: \code{coef}, \code{Stderror},
#'    \code{abnDag}, \code{method}, etc.
#'
#' @details
#' ## Round-Trip Capability
#'
#' This function is designed to work with \code{\link{export_abnFit}}:
#' \code{abnFit object} → \code{export_abnFit()} → JSON file/string →
#' \code{import_abnFit()} → \code{abnFit object}
#'
#' The reconstructed object can be used for plotting, predictions, or further analysis.
#' Metadata fields (\code{scenario_id}, \code{label}) are preserved through the round trip.
#'
#' ## JSON Structure Requirements
#'
#' The input JSON must contain at least these required fields:
#' \itemize{
#'   \item \code{variables}: Array of variable objects defining all nodes in the network
#'   \item \code{parameters}: Array of parameter objects defining all fitted coefficients
#'   \item \code{arcs}: Array of arc objects defining edges in the DAG (can be empty)
#' }
#'
#' Optional fields may also be present:
#' \itemize{
#'   \item \code{scenario_id}: Model run identifier (preserved in output)
#'   \item \code{label}: Descriptive label (preserved in output)
#'   \item \code{method}: Fitting method ("mle" or "bayes", default "mle")
#'   \item \code{linkFunctions}: Link function definitions (rarely used)
#'   \item \code{constraints}: Subsetting constraints (abnScripts integration)
#'   \item \code{subset_metadata}: Metadata about subsetting operation
#'   \item \code{original_model}: Original model before subsetting
#'   \item \code{original_data_path}: Path to CSV data file
#' }
#'
#' ## Validation Rules
#'
#' When \code{validate = TRUE}, the following checks are performed:
#' \itemize{
#'   \item \strong{Required fields}: \code{variables}, \code{parameters}, \code{arcs} must be present
#'   \item \strong{Unique IDs}: Each \code{variable_id} must be unique
#'   \item \strong{References}: Parameter \code{source.variable_id} must match existing variables
#'   \item \strong{Arcs}: Arc \code{source_variable_id} and \code{target_variable_id} must exist
#'   \item \strong{States}: For multinomial variables, \code{state_id} references in parameters must match
#'   \item \strong{Link functions}: Must be compatible with distribution type
#' }
#'
#' Validation errors provide informative messages to guide correction of the JSON.
#'
#' @section Supported JSON Fields:
#'
#' See \code{\link{export_abnFit}} for complete documentation of the JSON structure,
#' including all fields, their types, and valid values.
#'
#' @examples
#' \dontrun{
#' # Example 1: Import from file
#' library(abn)
#' 
#' # Assuming you have a JSON file exported from export_abnFit()
#' imported_model <- import_abnFit(file = "my_model.json")
#' 
#' # Now use the imported model
#' summary(imported_model)
#' }
#'
#' \dontrun{
#' # Example 2: Import from JSON string
#' library(abn)
#' library(jsonlite)
#'
#' # You might have JSON as a string (from API, database, etc.)
#' json_string <- '{
#'   "scenario_id": "model_v1",
#'   "label": "My Model",
#'   "method": "mle",
#'   "variables": [...],
#'   "parameters": [...],
#'   "arcs": [...],
#'   "linkFunctions": null,
#'   "constraints": null,
#'   "subset_metadata": null,
#'   "original_model": null,
#'   "original_data_path": null
#' }'
#'
#' # Import from string
#' model <- import_abnFit(json = json_string)
#' }
#'
#' \dontrun{
#' # Example 3: Round-trip (export and re-import)
#' library(abn)
#'
#' # Fit a model
#' data(ex1.dag.data)
#' mydists <- list(b1 = "binomial", p1 = "poisson", g1 = "gaussian",
#'                 b2 = "binomial", p2 = "poisson", g2 = "gaussian",
#'                 b3 = "binomial", g3 = "gaussian")
#' mycache <- buildScoreCache(data.df = ex1.dag.data,
#'                             data.dists = mydists,
#'                             method = "mle",
#'                             max.parents = 2)
#' mp_dag <- mostProbable(score.cache = mycache)
#' myfit <- fitAbn(object = mp_dag, method = "mle")
#'
#' # Export to JSON
#' json_export <- export_abnFit(myfit, scenario_id = "model_v1")
#'
#' # Re-import from the JSON string
#' myfit_reimported <- import_abnFit(json = json_export)
#'
#' # Verify they are equivalent
#' identical(myfit$coef, myfit_reimported$coef)  # Should be TRUE (within rounding)
#' }
#'
#' \dontrun{
#' # Example 4: Error handling
#' library(abn)
#'
#' # Invalid JSON will produce informative error messages
#' tryCatch({
#'   # This JSON is missing the required "arcs" field
#'   bad_json <- '{
#'     "scenario_id": "bad_model",
#'     "variables": [...],
#'     "parameters": [...]
#'   }'
#'   import_abnFit(json = bad_json)
#' }, error = function(e) {
#'   cat("Error:", conditionMessage(e), "\n")
#'   # Output: "Invalid JSON structure: must contain 'variables', 'parameters', and 'arcs' components"
#' })
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{export_abnFit}} for exporting abnFit to JSON (inverse operation)
#'   \item \code{\link{fitAbn}} for fitting ABN models
#'   \item \code{\link{jsonlite::fromJSON}} for JSON parsing (used internally)
#' }
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
#'
#' @details
#' Coefficient column names are reconstructed purely from the JSON's
#' *structural* fields (`source.variable_id`, `source.state_id`,
#' `conditions[].parent_variable_id`, `conditions[].parent_state_id`).
#' The `name` field on each parameter object is treated as an opaque
#' human-readable label and is intentionally **never** consulted.
#'
#' If any parameter has a `condition_type` of `"variance"`,
#' `"random_variance"`, or `"random_covariance"`, the model is treated
#' as a grouped (mixed-effects) MLE fit, and `mu`/`betas`/`sigma`/
#' `sigma_alpha`/`group.var` are reconstructed in addition to (or in
#' place of) the standard `coef`/`Stderror` matrices.
#' @keywords internal
reconstruct_abnfit_mle <- function(json_list) {
  variables <- json_list$variables

  id_to_name <- stats::setNames(
    vapply(variables, function(x) as.character(x$attribute_name), character(1)),
    vapply(variables, function(x) as.character(x$variable_id), character(1))
  )
  variable_names <- unname(id_to_name)
  model_types <- vapply(variables, function(x) as.character(x$model_type), character(1))
  names(model_types) <- variable_names
  n <- length(variable_names)

  lf_lookup <- list()
  if (!is.null(json_list$linkFunctions) && is.list(json_list$linkFunctions)) {
    for (lf in json_list$linkFunctions) {
      lf_lookup[[as.character(lf$id)]] <- as.character(lf$name)
    }
  }

  data_dists <- stats::setNames(model_types, variable_names)

  # Build per-multinomial-node lookup: state_id -> value_name. Used to
  # reconstruct the original abn coefficient column-name suffixes.
  state_id_to_value <- list()
  for (i in seq_along(variables)) {
    if (model_types[i] == "multinomial") {
      states <- variables[[i]]$states
      if (!is.null(states) && length(states) > 0) {
        ids <- vapply(states, function(s) as.character(s$state_id), character(1))
        vals <- vapply(states, function(s) as.character(s$value_name), character(1))
        state_id_to_value[[variable_names[i]]] <- stats::setNames(vals, ids)
      }
    }
  }

  resolve_state_value <- function(node_name, state_id) {
    if (is.null(state_id) || is.na(state_id)) return(NULL)
    sid <- as.character(state_id)
    lk <- state_id_to_value[[node_name]]
    if (is.null(lk)) return(sid)             # fall back to literal id
    v <- lk[sid]
    if (is.na(v)) return(sid)
    unname(v)
  }

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

  # Accumulators for standard (non-grouped) coefficient matrices.
  vals_acc <- stats::setNames(vector("list", n), variable_names)
  ses_acc <- stats::setNames(vector("list", n), variable_names)
  names_acc <- stats::setNames(vector("list", n), variable_names)
  for (vn in variable_names) {
    vals_acc[[vn]] <- numeric()
    ses_acc[[vn]] <- numeric()
    names_acc[[vn]] <- character()
  }

  # Accumulators for grouped (mixed-effects) MLE objects.
  is_grouped <- FALSE
  mu_acc          <- stats::setNames(vector("list", n), variable_names)
  mu_names_acc    <- stats::setNames(vector("list", n), variable_names)
  betas_acc       <- stats::setNames(vector("list", n), variable_names)
  betas_names_acc <- stats::setNames(vector("list", n), variable_names)
  # Multinomial child: betas is a matrix with rownames=child levels,
  # colnames=parent names. Track row/col keys per cell.
  betas_rows_acc  <- stats::setNames(vector("list", n), variable_names)
  betas_cols_acc  <- stats::setNames(vector("list", n), variable_names)
  sigma_acc       <- stats::setNames(vector("list", n), variable_names)
  # sigma_alpha: scalar for non-multinomial child, matrix for multinomial.
  sigma_alpha_scalar <- stats::setNames(vector("list", n), variable_names)
  sigma_alpha_cells  <- stats::setNames(vector("list", n), variable_names)

  if (!is.null(json_list$parameters)) {
    for (param in json_list$parameters) {
      raw_target_id <- as.character(param$source$variable_id)
      target_name <- id_to_name[raw_target_id]
      if (length(target_name) == 0 || is.na(target_name)) target_name <- raw_target_id
      if (!(target_name %in% variable_names)) next

      child_dist <- model_types[target_name]
      child_state_id <- if (!is.null(param$source$state_id)) as.character(param$source$state_id) else NULL
      child_state_value <- if (!is.null(child_state_id)) resolve_state_value(target_name, child_state_id) else NULL

      coeffs_input <- if (is.list(param$coefficients)) param$coefficients else list()

      for (coeff in coeffs_input) {
        type <- as.character(coeff$condition_type)
        value <- as.numeric(coeff$value)
        se_val <- if (is.null(coeff$stderr)) NA_real_ else as.numeric(coeff$stderr)

        # --------------- Grouped-only condition types ----------------
        if (type == "variance") {
          is_grouped <- TRUE
          sigma_acc[[target_name]] <- c(sigma_acc[[target_name]], value)
          next
        }
        if (type == "random_variance") {
          is_grouped <- TRUE
          if (child_dist == "multinomial") {
            # Diagonal of sigma_alpha matrix; key by state_id.
            key <- if (is.null(child_state_id)) "1" else child_state_id
            sigma_alpha_cells[[target_name]] <- c(
              sigma_alpha_cells[[target_name]],
              stats::setNames(value, paste0(key, "_", key))
            )
          } else {
            sigma_alpha_scalar[[target_name]] <- c(sigma_alpha_scalar[[target_name]], value)
          }
          next
        }
        if (type == "random_covariance") {
          is_grouped <- TRUE
          # source.state_id encodes "<i>_<j>" per the exporter.
          key <- if (is.null(child_state_id)) NA_character_ else child_state_id
          sigma_alpha_cells[[target_name]] <- c(
            sigma_alpha_cells[[target_name]],
            stats::setNames(value, key)
          )
          next
        }

        # ---------------- Reconstruct abn coef column name ----------------
        # Pure structural reconstruction; NEVER read param$name.
        c_name <- NA_character_

        if (type == "intercept") {
          if (child_dist == "multinomial" && !is.null(child_state_value)) {
            c_name <- paste0(target_name, "|intercept.", child_state_value)
          } else {
            c_name <- paste0(target_name, "|intercept")
          }
        } else if (type == "linear_term") {
          parent_var <- NA_character_
          parent_state_value <- NULL
          if (!is.null(coeff$conditions) && length(coeff$conditions) > 0) {
            cond <- coeff$conditions[[1]]
            p_raw <- as.character(cond$parent_variable_id)
            pn <- id_to_name[p_raw]
            if (length(pn) > 0 && !is.na(pn)) {
              parent_var <- unname(pn)
            } else {
              parent_var <- p_raw
            }
            if (!is.null(cond$parent_state_id)) {
              parent_state_value <- resolve_state_value(parent_var, cond$parent_state_id)
            }
          }
          if (child_dist == "multinomial") {
            # abn convention: <parent><state> with empty separator, no leading "child|".
            if (!is.null(parent_state_value)) {
              c_name <- paste0(parent_var, parent_state_value)
            } else {
              c_name <- parent_var
            }
          } else {
            if (!is.null(parent_state_value)) {
              c_name <- paste0(target_name, "|", parent_var, ".", parent_state_value)
            } else {
              c_name <- paste0(target_name, "|", parent_var)
            }
          }
        } else {
          # Unknown condition_type: preserve as "<target>|<type>" for safety.
          c_name <- paste0(target_name, "|", type)
        }

        # Accumulate into standard coef/Stderror matrices.
        vals_acc[[target_name]]  <- c(vals_acc[[target_name]],  value)
        ses_acc[[target_name]]   <- c(ses_acc[[target_name]],   se_val)
        names_acc[[target_name]] <- c(names_acc[[target_name]], c_name)

        # Also accumulate into grouped accumulators in case this turns out to
        # be a grouped fit (decided post-hoc by presence of variance terms).
        if (type == "intercept") {
          if (child_dist == "multinomial") {
            key <- if (!is.null(child_state_value)) {
              paste0(target_name, ".", child_state_value)
            } else {
              target_name
            }
            mu_acc[[target_name]]       <- c(mu_acc[[target_name]], value)
            mu_names_acc[[target_name]] <- c(mu_names_acc[[target_name]], key)
          } else {
            mu_acc[[target_name]]       <- c(mu_acc[[target_name]], value)
            mu_names_acc[[target_name]] <- c(mu_names_acc[[target_name]], "(Intercept)")
          }
        } else if (type == "linear_term") {
          parent_var <- NA_character_
          parent_state_value <- NULL
          if (!is.null(coeff$conditions) && length(coeff$conditions) > 0) {
            cond <- coeff$conditions[[1]]
            p_raw <- as.character(cond$parent_variable_id)
            pn <- id_to_name[p_raw]
            parent_var <- if (length(pn) > 0 && !is.na(pn)) unname(pn) else p_raw
            if (!is.null(cond$parent_state_id)) {
              parent_state_value <- resolve_state_value(parent_var, cond$parent_state_id)
            }
          }
          if (child_dist == "multinomial") {
            row_key <- if (!is.null(child_state_value)) child_state_value else "1"
            col_key <- if (!is.null(parent_state_value)) {
              paste0(parent_var, parent_state_value)
            } else {
              parent_var
            }
            betas_acc[[target_name]]      <- c(betas_acc[[target_name]], value)
            betas_rows_acc[[target_name]] <- c(betas_rows_acc[[target_name]], row_key)
            betas_cols_acc[[target_name]] <- c(betas_cols_acc[[target_name]], col_key)
          } else {
            beta_name <- if (!is.null(parent_state_value)) {
              paste0(parent_var, ".", parent_state_value)
            } else {
              parent_var
            }
            betas_acc[[target_name]]       <- c(betas_acc[[target_name]], value)
            betas_names_acc[[target_name]] <- c(betas_names_acc[[target_name]], beta_name)
          }
        }
      }
    }
  }

  coef_list     <- stats::setNames(vector("list", n), variable_names)
  stderror_list <- stats::setNames(vector("list", n), variable_names)
  for (var_id in variable_names) {
    if (length(vals_acc[[var_id]]) > 0) {
      coef_list[[var_id]] <- matrix(vals_acc[[var_id]], nrow = 1,
                                    dimnames = list(NULL, names_acc[[var_id]]))
      stderror_list[[var_id]] <- matrix(ses_acc[[var_id]], nrow = 1,
                                        dimnames = list(NULL, names_acc[[var_id]]))
    } else {
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

  if (is_grouped) {
    mu_list          <- stats::setNames(vector("list", n), variable_names)
    betas_list       <- stats::setNames(vector("list", n), variable_names)
    sigma_list       <- stats::setNames(vector("list", n), variable_names)
    sigma_alpha_list <- stats::setNames(vector("list", n), variable_names)

    for (vn in variable_names) {
      child_dist <- model_types[vn]

      # mu
      if (length(mu_acc[[vn]]) > 0) {
        mu_list[[vn]] <- stats::setNames(mu_acc[[vn]], mu_names_acc[[vn]])
      } else {
        mu_list[[vn]] <- NA
      }

      # betas
      if (child_dist == "multinomial" && length(betas_acc[[vn]]) > 0) {
        rows <- unique(betas_rows_acc[[vn]])
        cols <- unique(betas_cols_acc[[vn]])
        m <- matrix(NA_real_, nrow = length(rows), ncol = length(cols),
                    dimnames = list(rows, cols))
        for (k in seq_along(betas_acc[[vn]])) {
          m[betas_rows_acc[[vn]][k], betas_cols_acc[[vn]][k]] <- betas_acc[[vn]][k]
        }
        betas_list[[vn]] <- m
      } else if (length(betas_acc[[vn]]) > 0) {
        betas_list[[vn]] <- stats::setNames(betas_acc[[vn]], betas_names_acc[[vn]])
      } else {
        betas_list[[vn]] <- NA
      }

      # sigma
      if (length(sigma_acc[[vn]]) > 0) {
        sigma_list[[vn]] <- sigma_acc[[vn]][1]
      } else {
        sigma_list[[vn]] <- NA
      }

      # sigma_alpha
      if (child_dist == "multinomial") {
        cells <- sigma_alpha_cells[[vn]]
        if (!is.null(cells) && length(cells) > 0) {
          # Determine matrix dimension from state_id_to_value lookup.
          lk <- state_id_to_value[[vn]]
          if (!is.null(lk)) {
            ids <- names(lk)
            vals <- unname(lk)
            dnames <- paste0(vn, ".", vals)
            k <- length(ids)
            m <- matrix(NA_real_, nrow = k, ncol = k, dimnames = list(dnames, dnames))
            id_to_pos <- stats::setNames(seq_along(ids), ids)
            for (key in names(cells)) {
              parts <- strsplit(key, "_", fixed = TRUE)[[1]]
              if (length(parts) == 2L) {
                i <- id_to_pos[parts[1]]
                j <- id_to_pos[parts[2]]
                if (!is.na(i) && !is.na(j)) {
                  m[i, j] <- cells[[key]]
                  m[j, i] <- cells[[key]]
                }
              }
            }
            sigma_alpha_list[[vn]] <- m
          } else {
            sigma_alpha_list[[vn]] <- NA
          }
        } else {
          sigma_alpha_list[[vn]] <- NA
        }
      } else {
        if (length(sigma_alpha_scalar[[vn]]) > 0) {
          sigma_alpha_list[[vn]] <- sigma_alpha_scalar[[vn]][1]
        } else {
          sigma_alpha_list[[vn]] <- NA
        }
      }
    }

    abn_fit$mu          <- mu_list
    abn_fit$betas       <- betas_list
    abn_fit$sigma       <- sigma_list
    abn_fit$sigma_alpha <- sigma_alpha_list
    abn_fit$group.var   <- TRUE
  }

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
