#' Import abnFit object from structured JSON format
#'
#' @description
#' Imports a fitted Additive Bayesian Network (ABN) model from a structured JSON
#' format, reconstructing an abnFit object suitable for further analysis.
#'
#' @param file Character string specifying a file path containing the JSON
#'   representation of the model. Alternatively, a JSON string can be provided
#'   directly via the `json` parameter.
#' @param json Optional character string containing the JSON representation
#'   of the model. If provided, `file` is ignored.
#' @param validate Logical, whether to validate the imported object against
#'   the abnFit class requirements. Default is TRUE.
#' @param ... Additional import options (currently unused, reserved for future extensions).
#'
#' @return An object of class \code{abnFit} representing the imported model.
#'
#' @details
#' This function provides a standardized way to import fitted ABN models from JSON,
#' enabling model sharing, archiving, and integration with external tools or
#' databases. The import function expects the JSON to follow the same structure
#' produced by \code{\link{export_abnFit}}.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{export_abnFit}} for exporting ABN models to JSON
#'   \item \code{\link{fitAbn}} for fitting ABN models
#' }
#'
#' @examples
#' \dontrun{
#' # Import from file
#' imported_model <- import_abnFit(file = "path/to/model.json")
#'
#' # Import from JSON string
#' json_string <- '{"scenario_id": "test", "label": "Test Model", ...}'
#' imported_model <- import_abnFit(json = json_string)
#'
#' # Validate the imported model
#' if (inherits(imported_model, "abnFit")) {
#'   print(imported_model)
#' }
#' }
#'
#' @export
import_abnFit <- function(file = NULL, json = NULL, validate = TRUE, ...) {
  # Input validation
  if (is.null(file) && is.null(json)) {
    stop("Either 'file' or 'json' must be provided")
  }
  
  if (!is.null(file) && !is.null(json)) {
    warning("Both 'file' and 'json' provided; using 'json' parameter")
  }
  
  # Load JSON
  if (!is.null(json)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Package 'jsonlite' is required for JSON import. Please install it.")
    }
    json_list <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  } else {
    if (!file.exists(file)) {
      stop(sprintf("File '%s' does not exist", file))
    }
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Package 'jsonlite' is required for JSON import. Please install it.")
    }
    json_list <- jsonlite::fromJSON(paste(readLines(file), collapse = "\n"), simplifyVector = FALSE)
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
#'
#' @param json_list List object parsed from JSON
#' @keywords internal
validate_json_structure <- function(json_list) {
  required_keys <- c("variables", "parameters", "arcs")
  if (!is.list(json_list) || !all(required_keys %in% names(json_list))) {
    stop("Invalid JSON structure: must contain 'variables', 'parameters', and 'arcs' components")
  }
  
  # Handle label/name field - prefer label, fallback to name, warn if both present
  if (!is.null(json_list$label) && !is.null(json_list$name)) {
    warning("Both 'label' and 'name' fields present; using 'label' and ignoring 'name'")
  }
  
  # Validate variables structure
  if (!is.null(json_list$variables)) {
    for (i in seq_along(json_list$variables)) {
      var <- json_list$variables[[i]]
      if (!all(c("variable_id", "attribute_name", "model_type") %in% names(var))) {
        stop(sprintf("Variable %d missing required fields", i))
      }
      # Ignore transform_params if present
      if (!is.null(var$states) && !is.list(var$states)) {
        stop(sprintf("Variable %d states must be a list", i))
      }
    }
  }
  
  # Handle linkFunctions lookup table if present
  link_functions <- json_list$linkFunctions %||% list()
  
  # Validate parameters structure (accept both embedded and normalized formats)
  if (!is.null(json_list$parameters)) {
    for (i in seq_along(json_list$parameters)) {
      param <- json_list$parameters[[i]]
      
      # Handle link function resolution
      link_func_name <- NULL
      if (!is.null(param$link_function_name)) {
        # Direct link function name
        link_func_name <- param$link_function_name
      } else if (!is.null(param$link_function_id) && 
                 length(link_functions) > 0 && 
                 !is.null(link_functions[[param$link_function_id]])) {
        # Look up from linkFunctions table
        lf <- link_functions[[param$link_function_id]]
        link_func_name <- lf$label %||% lf$name %||% NULL
        if (is.null(link_func_name)) {
          warning(sprintf("Could not resolve link function for ID %s", param$link_function_id))
        }
      }
      
      if (is.null(link_func_name)) {
        stop(sprintf("Parameter %d missing link function information", i))
      }
      
      # Check for required fields (accepting both formats)
      required_param_fields <- c("parameter_id", "name", "source", "coefficients")
      if (!all(required_param_fields %in% names(param))) {
        stop(sprintf("Parameter %d missing required fields: %s", i, 
                     paste(setdiff(required_param_fields, names(param)), collapse=", ")))
      }
      
      if (!is.list(param$source) || !"variable_id" %in% names(param$source)) {
        stop(sprintf("Parameter %d source missing variable_id", i))
      }
      
      # Handle both embedded and normalized coefficient formats
      if (is.list(param$coefficients)) {
        # Embedded format (current export format)
        for (j in seq_along(param$coefficients)) {
          coeff <- param$coefficients[[j]]
          if (!all(c("value", "condition_type", "conditions") %in% names(coeff))) {
            stop(sprintf("Parameter %d coefficient %d missing required fields", i, j))
          }
          if (!is.list(coeff$conditions)) {
            stop(sprintf("Parameter %d coefficient %d conditions must be a list", i, j))
          }
        }
      } else if (is.data.frame(param$coefficients)) {
        # Normalized format (from abstract schema) - convert to embedded
        # We'll handle this in the reconstruction phase
      } else {
        stop(sprintf("Parameter %d coefficients must be a list or data frame", i))
      }
    }
  }
  
  # Validate arcs structure
  if (!is.null(json_list$arcs)) {
    for (i in seq_along(json_list$arcs)) {
      arc <- json_list$arcs[[i]]
      if (!all(c("source_variable_id", "target_variable_id") %in% names(arc))) {
        stop(sprintf("Arc %d missing required fields", i))
      }
    }
  }
  
  invisible(TRUE)
}

#' Reconstruct abnFit object for MLE method from JSON
#'
#' @param json_list List object parsed from JSON
#' @keywords internal
reconstruct_abnfit_mle <- function(json_list) {
  # Extract basic information
  scenario_id <- json_list$scenario_id %||% NULL
  # Handle label/name field - prefer label, fallback to name
  label <- json_list$label %||% (json_list$name %||% NULL)
  
  # Handle label/name conflict
  if (!is.null(json_list$label) && !is.null(json_list$name)) {
    warning("Both 'label' and 'name' fields present; using 'label' and ignoring 'name'")
  }
  
  # Build variables information
  variables <- json_list$variables
  variable_ids <- sapply(variables, function(x) x$variable_id)
  attribute_names <- sapply(variables, function(x) x$attribute_name)
  model_types <- sapply(variables, function(x) x$model_type)
  
  # Create data distribution list
  data_dists <- set_names(model_types, variable_ids)
  
  # Create empty data frame with appropriate column types
  # We'll create a minimal dataset just for structure - actual values aren't needed for the object structure
  n_vars <- length(variable_ids)
  data_df <- as.data.frame(set_names(vector("list", n_vars), variable_ids))
  
  # Initialize each column with appropriate type based on distribution
  for (i in seq_along(variable_ids)) {
    var_id <- variable_ids[i]
    dist <- model_types[i]
    if (dist == "gaussian") {
      data_df[[var_id]] <- numeric(0)
    } else if (dist == "binomial") {
      data_df[[var_id]] <- factor(0, levels = c(0, 1))
    } else if (dist == "poisson") {
      data_df[[var_id]] <- integer(0)
    } else if (dist == "multinomial") {
      # Get states if available
      states <- variables[[i]]$states %||% list()
      state_names <- if (length(states) > 0) {
        sapply(states, function(x) x$value_name)
      } else {
        c("state1", "state2")  # fallback
      }
      data_df[[var_id]] <- factor(character(0), levels = state_names)
    }
  }
  
  # Build DAG from arcs
  n <- length(variable_ids)
  dag_matrix <- matrix(0, nrow = n, ncol = n, 
                       dimnames = list(variable_ids, variable_ids))
  
  if (!is.null(json_list$arcs) && length(json_list$arcs) > 0) {
    for (arc in json_list$arcs) {
      source <- arc$source_variable_id
      target <- arc$target_variable_id
      if (source %in% variable_ids && target %in% variable_ids) {
        dag_matrix[source, target] <- 1
      }
    }
  }
  
  # Create abnDag structure (simplified)
  abnDag <- list(
    dag = dag_matrix,
    data.df = data_df,
    data.dists = data_dists
  )
  class(abnDag) <- "abnDag"
  
  # Build coefficients and standard errors from parameters
  # Initialize empty lists for each node
  coef_list <- set_names(vector("list", n), variable_ids)
  stderror_list <- set_names(vector("list", n), variable_ids)
  
  # Process each parameter
  if (!is.null(json_list$parameters)) {
    for (param in json_list$parameters) {
      target_var <- param$source$variable_id
      if (!(target_var %in% variable_ids)) {
        next  # Skip if variable not found
      }
      
      # Initialize coefficient matrix for this variable if needed
      if (is.null(coef_list[[target_var]])) {
        coef_list[[target_var]] <- matrix(NA, nrow = 0, ncol = 0)
        stderror_list[[target_var]] <- matrix(NA, nrow = 0, ncol = 0)
      }
      
      # Process each coefficient - handle both embedded and normalized formats
      if (is.list(param$coefficients)) {
        // Embedded format (current export format)
        for (coeff in param$coefficients) {
          // Create coefficient name
          if (coeff$condition_type == "intercept") {
            coeff_name <- paste0(target_var, "|intercept")
          } else if (coeff$condition_type == "linear_term") {
            // Find parent variable from conditions
            parent_var <- NULL
            if (!is.null(coeff$conditions) && length(coeff$conditions) > 0) {
              parent_var <- coeff$conditions[[1]]$parent_variable_id
            }
            if (is.null(parent_var)) {
              parent_var <- "unknown"
            }
            coeff_name <- paste0(target_var, "|", parent_var)
          } else {
            // For other types, use a generic name
            coeff_name <- paste0(target_var, "|", coeff$condition_type, "_", 
                                 length(coef_list[[target_var]]) + 1)
          }
          
          // Resize matrices to accommodate new coefficient
          current_coef <- coef_list[[target_var]]
          current_se <- stderror_list[[target_var]]
          
          new_ncol <- ncol(current_coef) + 1
          coef_list[[target_var]] <- cbind(current_coef, 
                                           c(rep(NA, nrow(current_coef)), 
                                             coeff$value))
          stderror_list[[target_var]] <- cbind(current_se, 
                                               c(rep(NA, nrow(current_se)), 
                                                 coeff$stderr %||% NULL))
          
          // Set column names
          colnames(coef_list[[target_var]]) <- 
            c(if (ncol(current_coef) > 0) colnames(current_coef) else character(0), 
              coeff_name)
          colnames(stderror_list[[target_var]]) <- 
            c(if (ncol(current_se) > 0) colnames(current_se) else character(0), 
              coeff_name)
        }
      } else if (is.data.frame(param$coefficients)) {
        // Normalized format (from abstract schema) - convert to embedded
        coeff_df <- param$coefficients
        for (j in seq_len(nrow(coeff_df))) {
          coeff_row <- coeff_df[j, ]
          
          // Extract coefficient information
          coeff_value <- coeff_row$value
          coeff_stderr <- coeff_row$stderr %||% NULL
          coeff_condition_type <- coeff_row$condition_type
          coeff_conditions <- coeff_row$conditions
          
          // Create coefficient name based on condition type
          if (coeff_condition_type == "intercept") {
            coeff_name <- paste0(target_var, "|intercept")
          } else if (coeff_condition_type == "linear_term") {
            // Find parent variable from conditions
            parent_var <- NULL
            if (!is.null(coeff_conditions) && length(coeff_conditions) > 0) {
              // Handle conditions as list or data frame
              if (is.list(coeff_conditions)) {
                parent_var <- coeff_conditions$parent_variable_id
              } else if (is.data.frame(coeff_conditions) && nrow(coeff_conditions) > 0) {
                parent_var <- coeff_conditions$parent_variable_id[1]
              }
            }
            if (is.null(parent_var)) {
              parent_var <- "unknown"
            }
            coeff_name <- paste0(target_var, "|", parent_var)
          } else {
            // For other types, use a generic name
            coeff_name <- paste0(target_var, "|", coeff_condition_type, "_", 
                                 length(coef_list[[target_var]]) + 1)
          }
          
          // Resize matrices to accommodate new coefficient
          current_coef <- coef_list[[target_var]]
          current_se <- stderror_list[[target_var]]
          
          new_ncol <- ncol(current_coef) + 1
          coef_list[[target_var]] <- cbind(current_coef, 
                                           c(rep(NA, nrow(current_coef)), 
                                             coeff_value))
          stderror_list[[target_var]] <- cbind(current_se, 
                                               c(rep(NA, nrow(current_se)), 
                                                 coeff_stderr))
          
          // Set column names
          colnames(coef_list[[target_var]]) <- 
            c(if (ncol(current_coef) > 0) colnames(current_coef) else character(0), 
              coeff_name)
          colnames(stderror_list[[target_var]]) <- 
            c(if (ncol(current_se) > 0) colnames(current_se) else character(0), 
              coeff_name)
        }
      }
          if (is.null(parent_var)) {
            parent_var <- "unknown"
          }
          coeff_name <- paste0(target_var, "|", parent_var)
        } else {
          # For other types, use a generic name
          coeff_name <- paste0(target_var, "|", coeff$condition_type, "_", 
                               length(coef_list[[target_var]]) + 1)
        }
        
        # Resize matrices to accommodate new coefficient
        current_coef <- coef_list[[target_var]]
        current_se <- stderror_list[[target_var]]
        
        new_ncol <- ncol(current_coef) + 1
        coef_list[[target_var]] <- cbind(current_coef, 
                                         c(rep(NA, nrow(current_coef)), 
                                           coeff$value))
        stderror_list[[target_var]] <- cbind(current_se, 
                                             c(rep(NA, nrow(current_se)), 
                                               coeff$stderr %||% NULL))
        
        # Set column names
        colnames(coef_list[[target_var]]) <- 
          c(if (ncol(current_coef) > 0) colnames(current_coef) else character(0), 
            coeff_name)
        colnames(stderror_list[[target_var]]) <- 
          c(if (ncol(current_se) > 0) colnames(current_se) else character(0), 
            coeff_name)
      }
    }
  }
  
  # Ensure all variables have coefficient matrices (even if empty)
  for (var_id in variable_ids) {
    if (is.null(coef_list[[var_id]])) {
      coef_list[[var_id]] <- matrix(NA, nrow = 0, ncol = 0)
      stderror_list[[var_id]] <- matrix(NA, nrow = 0, ncol = 0)
    }
    
    # Set dimension names if needed
    if (nrow(coef_list[[var_id]]) == 0 && ncol(coef_list[[var_id]]) == 0) {
      dimnames(coef_list[[var_id]]) <- list(NULL, character(0))
      dimnames(stderror_list[[var_id]]) <- list(NULL, character(0))
    }
  }
  
  # Create abnFit object
  abn_fit <- list(
    abnDag = abnDag,
    coef = coef_list,
    Stderror = stderror_list,
    method = "mle",
    fit.stats = list(),  # Empty fit stats
    model.dic = NA,      # Placeholder
    call = match.call()
  )
  
  class(abn_fit) <- "abnFit"
  return(abn_fit)
}

#' Reconstruct abnFit object for Bayesian method from JSON
#'
#' @param json_list List object parsed from JSON
#' @keywords internal
reconstruct_abnfit_bayes <- function(json_list) {
  # For now, we'll create a basic structure and issue a warning
  # that full Bayesian import is not yet implemented
  warning("Bayesian model import is not fully implemented yet. Creating placeholder structure.")
  
  # Reuse the MLE reconstruction for basic structure
  abn_fit <- reconstruct_abnfit_mle(json_list)
  abn_fit$method <- "bayes"
  
  return(abn_fit)
}

#' Helper function to determine link function from distribution
#'
#' @keywords internal
get_link_function <- function(distribution) {
  link_functions <- list(
    "gaussian" = "identity",
    "binomial" = "logit",
    "poisson" = "log",
    "multinomial" = "logit"
  )
  
  return(if (is.null(link_functions[[distribution]])) "identity" else link_functions[[distribution]])
}

#' Validate abnFit object meets class requirements
#'
#' @param object Object to validate
#' @keywords internal
validate_abnfit_object <- function(object) {
  if (!inherits(object, "abnFit")) {
    stop("Imported object is not of class 'abnFit'")
  }
  
  # Additional validation of required components
  required_components <- c("abnDag", "method")
  missing <- setdiff(required_components, names(object))
  if (length(missing) > 0) {
    stop(sprintf("abnFit object missing required components: %s", 
                 paste(missing, collapse = ", ")))
  }
  
  # Validate abnDag structure
  if (!inherits(object$abnDag, "abnDag")) {
    stop("abnDag component is not of class 'abnDag'")
  }
  
  required_dag_components <- c("dag", "data.df", "data.dists")
  missing_dag <- setdiff(required_dag_components, names(object$abnDag))
  if (length(missing_dag) > 0) {
    stop(sprintf("abnDag object missing required components: %s", 
                 paste(missing_dag, collapse = ", ")))
  }
  
  # Validate method
  if (!object$method %in% c("mle", "bayes")) {
    stop(sprintf("Invalid method '%s'. Must be 'mle' or 'bayes'.", object$method))
  }
  
  invisible(TRUE)
}