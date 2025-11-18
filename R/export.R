#' Export abnFit object to structured format
#'
#' @param object An object of class abnFit
#' @param format The export format, currently only "json" is supported.
#' @param include_network Whether to include network structure
#' @param file Optional file path to save JSON. If NULL, returns JSON string
#' @param pretty Logical, whether to format JSON nicely
#' @param scenario_id Optional identifier for the model run/scenario
#' @param label Optional name/label for the scenario
#' @param ... Additional export options
#' @return A JSON string or writes to file if 'file' is provided
export_abnFit <- function(object, format = "json", include_network = TRUE,
                          file = NULL, pretty = TRUE, scenario_id = NULL,
                          label = NULL, ...) {
  if (!inherits(object, "abnFit")) {
    stop("Input object must be of class 'abnFit'")
  }

  # Dispatch based on method
  if (object$method == "mle") {
    export_list <- export_abnFit_mle(object, format, include_network,
                                     scenario_id, label, ...)
  } else if (object$method == "bayes") {
    export_list <- export_abnFit_bayes(object, format, include_network,
                                       scenario_id, label, ...)
  } else {
    stop("Unsupported method in abnFit object. Supported methods are 'mle' and 'bayes'.")
  }

  # Convert to desired format
  if (format == "json") {
    return(export_to_json(export_list, format, file, pretty))
  } else {
    stop("Currently, only 'json' format is supported.")
  }
}

#' Helper function to convert export list to JSON
#' @param export_list The list to convert to JSON. Must contain variables, parameters, and arcs components, see details.
#' @inheritParams export_abnFit
#' @details The export_list must be a named list with the following components:
#' \itemize{
#' \item scenario_id: Optional identifier for the model run/scenario.
#' \item label: Optional name/label for the scenario.
#' \item variables: An array where each element represents a variable/node with its metadata, distribution type, and states (for categorical variables).
#' \item parameters: An array where each element represents a parameter with its link function, source variable, and coefficients.
#' \item arcs: An array with arc details, each containing source_variable_id and target_variable_id.
#' }
#' @keywords internal
export_to_json <- function(export_list, format, file = NULL, pretty = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for JSON export. Please install it.")
  }

  # Validate export list structure - scenario_id and label are optional
  required_keys <- c("variables", "parameters", "arcs")
  if (!is.list(export_list) || !all(required_keys %in% names(export_list))) {
    stop("export_list must be a named list with components: variables, parameters, arcs")
  }

  # Convert to JSON
  json_str <- jsonlite::toJSON(export_list, auto_unbox = TRUE, pretty = pretty, null = "null")

  # Write to file or return string
  if (!is.null(file)) {
    writeLines(json_str, con = file)
    return(invisible(file))
  } else {
    return(json_str)
  }
}

#' Export abnFit object fitted with MLE (non-mixed effects)
#' @inheritParams export_abnFit
#' @details This function handles abnFit objects fitted using Maximum Likelihood Estimation (MLE)
#' without mixed-effects. It extracts variables metadata, arc details, and node parameterisations.
#' @return A named list with components: scenario_id, label, variables, parameters, arcs.
#' @keywords internal
export_abnFit_mle <- function(object, format, include_network, scenario_id = NULL,
                              label = NULL, ...) {
  # Extract arc details first
  arcs_details <- export_abnFit_mle_arcs(object)

  # Extract variable and parameter details based on grouping
  if (!is.null(object$group.var)) {
    # With grouping (mixed-effects)
    result <- export_abnFit_mle_grouped_nodes(object, format, include_network, ...)
    variables_list <- result$variables
    parameters_list <- result$parameters
  } else {
    # Without grouping
    result <- export_abnFit_mle_nodes(object, format, include_network, ...)
    variables_list <- result$variables
    parameters_list <- result$parameters
  }

  # Create export list with scenario_id and label at the top
  export_structure <- list()

  # Add scenario_id and label first (will be null if not provided)
  export_structure$scenario_id <- scenario_id
  export_structure$label <- label

  # Add the main components
  export_structure$variables <- variables_list
  export_structure$parameters <- parameters_list
  export_structure$arcs <- arcs_details

  return(export_structure)
}

#' Export node information from abnFit objects fitted with MLE (non-mixed effects)
#'
#' @param object An object of class abnFit fitted with method = "mle"
#' @param ... Additional arguments (currently unused)
#'
#' @details This function extracts node parameterisation information from abnFit objects
#' that were fitted using the Maximum Likelihood Estimation (MLE) approach without
#' mixed-effects (i.e., no group.var specified). The function processes the coefficients
#' and standard errors stored in the abnFit object.
#'
#' The \code{coef} component contains the estimated regression coefficients for each node,
#' stored as a matrix where column names indicate the parameter names (e.g., "g2",
#' "m11", "b1|intercept"). These represent the linear model coefficients from the
#' generalized linear model fitted to each node given its parents in the DAG.
#'
#' The \code{Stderror} component contains the corresponding standard errors for each
#' coefficient, providing a measure of uncertainty in the parameter estimates. The
#' structure mirrors that of the \code{coef} component.
#'
#' For different distribution types:
#' \itemize{
#' \item Gaussian nodes: Include intercept and slope coefficients
#' \item Binomial/Poisson nodes: Include intercept and slope coefficients on logit/log scale
#' \item Multinomial nodes: Include category-specific intercepts (reference level omitted)
#' and coefficients, following standard multinomial logistic regression conventions
#' }
#'
#' @returns A named list with two components: variables and parameters.
#' Variables is an array where each element represents a variable with its metadata.
#' Parameters is an array where each element represents a parameter with its coefficients.
#'
#' @keywords internal
export_abnFit_mle_nodes <- function(object, ...) {
  # Input validation
  if (!inherits(object, "abnFit")) {
    stop("Object must be of class 'abnFit'", call. = FALSE)
  }

  if (object$method != "mle") {
    stop("This function only handles abnFit objects fitted with method = 'mle'", call. = FALSE)
  }

  if (is.null(object$coef) || is.null(object$Stderror)) {
    stop("abnFit object must contain 'coef' and 'Stderror' components", call. = FALSE)
  }

  # Initialize output
  variables_list <- list()
  parameters_list <- list()
  parameter_counter <- 1

  node_names <- names(object$coef)
  node_dists <- object$abnDag$data.dists

  # Get parent information from DAG
  dag_matrix <- as.matrix(object$abnDag)

  # Process each node
  for (node_id in node_names) {
    # Extract coefficients and standard errors for this node
    coef_mat <- object$coef[[node_id]]
    se_mat <- object$Stderror[[node_id]]

    # Convert matrices to named vectors
    coef_vec <- as.numeric(coef_mat)
    names(coef_vec) <- colnames(coef_mat)
    se_vec <- as.numeric(se_mat)
    names(se_vec) <- colnames(se_mat)

    # Get distribution type for this node
    distribution <- node_dists[[node_id]]

    # Determine link function based on distribution
    link_function <- get_link_function(distribution)

    # Get parent nodes for this child
    parent_nodes <- names(dag_matrix[, node_id])[dag_matrix[, node_id] == 1]

    # Create variable entry
    variable_entry <- list(
      variable_id = node_id,
      attribute_name = node_id,
      model_type = distribution
    )

    # Add states for categorical variables (multinomial)
    if (distribution == "multinomial") {
      variable_entry$states <- extract_states_from_data(object, node_id)
    } else {
      variable_entry$states <- NULL
    }

    variables_list[[length(variables_list) + 1]] <- variable_entry

    # Extract parameters based on distribution type
    param_result <- extract_parameters_by_distribution(
      coef_vec, se_vec, distribution, node_id,
      parent_nodes, parameter_counter, link_function
    )

    # Add parameters to list
    parameters_list <- c(parameters_list, param_result$parameters)
    parameter_counter <- param_result$next_counter
  }

  return(list(
    variables = variables_list,
    parameters = parameters_list
  ))
}

#' Helper function to determine link function from distribution
#' @keywords internal
get_link_function <- function(distribution) {
  link_functions <- list(
    "gaussian" = "identity",
    "binomial" = "logit",
    "poisson" = "log",
    "multinomial" = "logit"
  )

  return(link_functions[[distribution]] %||% "identity")
}

#' Helper function for NULL default
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Extract states for categorical variables from data
#' @keywords internal
extract_states_from_data <- function(object, node_id) {
  # Get the data for this variable
  data_col <- object$abnDag$data.df[[node_id]]

  if (is.factor(data_col) || is.character(data_col)) {
    unique_vals <- sort(unique(as.character(data_col)))
    states <- lapply(seq_along(unique_vals), function(i) {
      list(
        state_id = as.character(i),
        value_name = unique_vals[i],
        is_baseline = (i == 1)  # First level is baseline
      )
    })
    return(states)
  }

  return(NULL)
}

#' Helper function to extract parameters based on distribution type
#' @keywords internal
extract_parameters_by_distribution <- function(coef_vec, se_vec, distribution, node_id,
                                               parent_nodes, start_counter, link_function) {
  param_names <- names(coef_vec)
  parameters <- list()
  counter <- start_counter

  if (distribution %in% c("gaussian", "binomial", "poisson")) {
    # Find intercept parameter
    intercept_pattern <- paste0(node_id, "\\|intercept")
    intercept_idx <- grep(intercept_pattern, param_names, ignore.case = TRUE)

    if (length(intercept_idx) > 0) {
      intercept_name <- param_names[intercept_idx[1]]

      # Create parameter for intercept
      param_entry <- list(
        parameter_id = as.character(counter),
        name = "intercept",
        link_function_name = link_function,
        source = list(
          variable_id = node_id
        ),
        coefficients = list(
          list(
            value = unname(coef_vec[intercept_name]),
            stderr = unname(se_vec[intercept_name]),
            condition_type = "intercept",
            conditions = list()
          )
        )
      )

      parameters[[length(parameters) + 1]] <- param_entry
      counter <- counter + 1

      # Remove intercept from remaining parameters
      param_names <- param_names[-intercept_idx[1]]
      coef_remaining <- coef_vec[param_names]
      se_remaining <- se_vec[param_names]
    } else {
      coef_remaining <- coef_vec
      se_remaining <- se_vec
    }

    # Process parent coefficients as linear terms
    if (length(coef_remaining) > 0) {
      for (i in seq_along(coef_remaining)) {
        param_name <- names(coef_remaining)[i]
        parent_var <- parent_nodes[i]

        param_entry <- list(
          parameter_id = as.character(counter),
          name = param_name,
          link_function_name = link_function,
          source = list(
            variable_id = node_id
          ),
          coefficients = list(
            list(
              value = unname(coef_remaining[i]),
              stderr = unname(se_remaining[i]),
              condition_type = "linear_term",
              conditions = list(
                list(
                  parent_variable_id = parent_var,
                  parent_state_id = NULL
                )
              )
            )
          )
        )

        parameters[[length(parameters) + 1]] <- param_entry
        counter <- counter + 1
      }
    }

  } else if (distribution == "multinomial") {
    # For multinomial, parameters are category-specific
    intercept_pattern <- paste0(node_id, "\\|intercept\\.")
    intercept_idx <- grep(intercept_pattern, param_names)

    if (length(intercept_idx) > 0) {
      # Extract category numbers from intercept names
      intercept_names <- param_names[intercept_idx]
      categories <- gsub(paste0(".*", node_id, "\\|intercept\\."), "", intercept_names)

      # Process each category
      for (cat in unique(categories)) {
        # Intercept for this category
        cat_intercept_pattern <- paste0(node_id, "\\|intercept\\.", cat, "$")
        cat_intercept_idx <- grep(cat_intercept_pattern, names(coef_vec))

        if (length(cat_intercept_idx) > 0) {
          param_entry <- list(
            parameter_id = as.character(counter),
            name = paste0("prob_", cat),
            link_function_name = link_function,
            source = list(
              variable_id = node_id,
              state_id = cat
            ),
            coefficients = list(
              list(
                value = unname(coef_vec[cat_intercept_idx[1]]),
                stderr = unname(se_vec[cat_intercept_idx[1]]),
                condition_type = "intercept",
                conditions = list()
              )
            )
          )

          parameters[[length(parameters) + 1]] <- param_entry
          counter <- counter + 1
        }

        # Parent coefficients for this category
        cat_pattern <- paste0("\\.", cat, "$")
        cat_coef_idx <- grep(cat_pattern, param_names)
        cat_coef_idx <- setdiff(cat_coef_idx, cat_intercept_idx)

        if (length(cat_coef_idx) > 0) {
          for (idx in cat_coef_idx) {
            param_name <- param_names[idx]
            # Try to match with parent nodes
            parent_var <- NULL
            for (p in parent_nodes) {
              if (grepl(p, param_name)) {
                parent_var <- p
                break
              }
            }

            if (!is.null(parent_var)) {
              param_entry <- list(
                parameter_id = as.character(counter),
                name = param_name,
                link_function_name = link_function,
                source = list(
                  variable_id = node_id,
                  state_id = cat
                ),
                coefficients = list(
                  list(
                    value = unname(coef_vec[idx]),
                    stderr = unname(se_vec[idx]),
                    condition_type = "linear_term",
                    conditions = list(
                      list(
                        parent_variable_id = parent_var,
                        parent_state_id = NULL
                      )
                    )
                  )
                )
              )

              parameters[[length(parameters) + 1]] <- param_entry
              counter <- counter + 1
            }
          }
        }
      }
    }
  }

  return(list(
    parameters = parameters,
    next_counter = counter
  ))
}

#' Export node information from abnFit objects fitted with MLE (mixed effects)
#'
#' @param object An object of class abnFit fitted with method = "mle" and group.var specified.
#' @param ... Additional arguments (currently unused)
#'
#' @details This function extracts node parameterisation information from abnFit objects
#' that were fitted using the Maximum Likelihood Estimation (MLE) approach WITH
#' mixed-effects (i.e., group.var was specified).
#'
#' For mixed-effects models, the structure includes:
#' \itemize{
#' \item Fixed effects: Population-level intercepts and coefficients
#' \item Random effects: Group-level variance components (sigma, sigma_alpha)
#' }
#'
#' The export format follows the same variables/parameters structure, but parameters
#' will include both fixed-effect coefficients and random-effect variance components.
#'
#' @returns A named list with two components: variables and parameters.
#' Variables is an array where each element represents a variable with its metadata.
#' Parameters is an array where each element represents a parameter, including both
#' fixed-effect coefficients and random-effect variance components.
#'
#' @keywords internal
export_abnFit_mle_grouped_nodes <- function(object, ...) {
  # Input validation
  if (!inherits(object, "abnFit")) {
    stop("Object must be of class 'abnFit'", call. = FALSE)
  }

  if (object$method != "mle") {
    stop("This function only handles abnFit objects fitted with method = 'mle'", call. = FALSE)
  }

  if (is.null(object$group.var)) {
    stop("This function only handles grouped (mixed-effects) models", call. = FALSE)
  }

  if (is.null(object$mu) || is.null(object$betas) ||
      is.null(object$sigma) || is.null(object$sigma_alpha)) {
    stop("abnFit object must contain 'mu', 'betas', 'sigma', and 'sigma_alpha' components", call. = FALSE)
  }

  # Initialize output
  variables_list <- list()
  parameters_list <- list()
  parameter_counter <- 1

  node_names <- names(object$mu)
  node_dists <- object$abnDag$data.dists

  # Get parent information from DAG
  dag_matrix <- as.matrix(object$abnDag)

  # Process each node
  for (node_id in node_names) {
    # Extract parameters for this node
    mu_val <- object$mu[[node_id]]
    betas_val <- object$betas[[node_id]]
    sigma_val <- object$sigma[[node_id]]
    sigma_alpha_val <- object$sigma_alpha[[node_id]]

    # Get distribution type for this node
    distribution <- node_dists[[node_id]]

    # Determine link function based on distribution
    link_function <- get_link_function(distribution)

    # Get parent nodes for this child
    parent_nodes <- names(dag_matrix[, node_id])[dag_matrix[, node_id] == 1]

    # Create variable entry
    variable_entry <- list(
      variable_id = node_id,
      attribute_name = node_id,
      model_type = distribution
    )

    # Add states for categorical variables (multinomial)
    if (distribution == "multinomial") {
      variable_entry$states <- extract_states_from_data(object, node_id)
    } else {
      variable_entry$states <- NULL
    }

    variables_list[[length(variables_list) + 1]] <- variable_entry

    # Extract parameters for mixed-effects models
    param_result <- extract_parameters_mixed_effects(
      mu_val, betas_val, sigma_val, sigma_alpha_val,
      distribution, node_id, parent_nodes, parameter_counter, link_function
    )

    # Add parameters to list
    parameters_list <- c(parameters_list, param_result$parameters)
    parameter_counter <- param_result$next_counter
  }

  return(list(
    variables = variables_list,
    parameters = parameters_list
  ))
}

#' Helper function to extract parameters based on distribution type with grouping
#' @param mu Fixed-effect intercept(s) from mu component
#' @param betas Fixed-effect coefficients from betas component
#' @param sigma Residual variance from sigma component
#' @param sigma_alpha Random-effect variance/covariance from sigma_alpha component
#' @param distribution Node distribution type
#' @param node_id Node identifier
#' @keywords internal
extract_parameters_by_distribution_grouping <- function(mu, betas, sigma, sigma_alpha, distribution, node_id) {

  param_list <- list()

  # Handle fixed effects (intercept and coefficients)
  param_list$fixed_effects <- list()

  if (distribution %in% c("gaussian", "binomial", "poisson")) {

    # Single intercept for these distributions
    param_list$fixed_effects$intercept <- mu

    # Handle parent coefficients
    if (!is.null(betas) && !is.logical(betas) && length(betas) > 0) {
      # Convert to list format for JSON export
      if (is.matrix(betas)) {
        coeff_list <- list()
        for (i in seq_len(ncol(betas))) {
          coeff_name <- colnames(betas)[i]
          coeff_list[[coeff_name]] <- betas[, i]
        }
        param_list$fixed_effects$coefficients <- coeff_list
      } else if (is.vector(betas) && length(betas) > 0) {
        coeff_list <- as.list(betas)
        param_list$fixed_effects$coefficients <- coeff_list
      }
    } else {
      param_list$fixed_effects$coefficients <- list()
    }

  } else if (distribution == "multinomial") {

    # Handle category-specific intercepts for multinomial
    if (length(mu) > 1) {
      # Multiple categories
      categories <- list()
      category_names <- names(mu)

      for (i in seq_along(mu)) {
        cat_name <- category_names[i]
        # Extract category number from name (e.g., "m1.2" -> "2")
        cat_num <- gsub(".*\\.", "", cat_name)
        categories[[paste0("category_", cat_num)]] <- list(
          intercept = mu[i]
        )
      }
      param_list$fixed_effects$categories <- categories
    } else {
      param_list$fixed_effects$intercept <- mu
    }

    # Handle parent coefficients for multinomial
    if (!is.null(betas) && !is.logical(betas) && length(betas) > 0) {
      if (is.matrix(betas)) {
        # Matrix format for multinomial with parents
        coeff_list <- list()
        for (i in seq_len(nrow(betas))) {
          for (j in seq_len(ncol(betas))) {
            row_name <- rownames(betas)[i]
            col_name <- colnames(betas)[j]
            coeff_list[[paste0(row_name, "_", col_name)]] <- betas[i, j]
          }
        }
        param_list$fixed_effects$coefficients <- coeff_list
      } else {
        param_list$fixed_effects$coefficients <- as.list(betas)
      }
    }
  }

  # Handle random effects (variance components)
  param_list$random_effects <- list()

  # Residual variance (sigma)
  if (!is.null(sigma) && !is.logical(sigma) && length(sigma) > 0) {
    param_list$random_effects$sigma <- sigma
  }

  # Random intercept variance/covariance (sigma_alpha)
  if (!is.null(sigma_alpha) && !is.logical(sigma_alpha)) {
    if (is.matrix(sigma_alpha)) {
      # Convert matrix to list format for JSON export
      sigma_alpha_list <- list()
      if (nrow(sigma_alpha) == ncol(sigma_alpha)) {
        # Symmetric variance-covariance matrix
        for (i in seq_len(nrow(sigma_alpha))) {
          for (j in seq_len(ncol(sigma_alpha))) {
            row_name <- rownames(sigma_alpha)[i]
            col_name <- colnames(sigma_alpha)[j]
            if (i <= j) {  # Only store upper triangle (symmetric)
              element_name <- if (i == j) paste0("var_", row_name) else paste0("cov_", row_name, "_", col_name)
              sigma_alpha_list[[element_name]] <- sigma_alpha[i, j]
            }
          }
        }
      }
      param_list$random_effects$sigma_alpha <- sigma_alpha_list
    } else {
      # Scalar variance
      param_list$random_effects$sigma_alpha <- sigma_alpha
    }
  }

  return(param_list)
}

#' Export arc information from abnFit objects fitted with MLE
#' @inheritParams export_abnFit
#' @details This function extracts arc information from abnFit objects fitted using MLE.
#' It retrieves the source and target nodes for each arc in the directed acyclic graph (DAG).
#' Currently, frequency, significance, and constraint information are not included in the export.
#' @return An array containing arc details: source_variable_id and target_variable_id for each arc.
#' @keywords internal
export_abnFit_mle_arcs <- function(object, ...) {
  edgelist <- as.data.frame(object$abnDag)

  if (nrow(edgelist) == 0) {
    return(list())
  }

  # Create array of arc objects
  arcs <- lapply(seq_len(nrow(edgelist)), function(i) {
    list(
      source_variable_id = as.character(edgelist$from[i]),
      target_variable_id = as.character(edgelist$to[i])
    )
  })

  return(arcs)
}

#' Export abnFit object fitted with Bayesian methods
#'
#' @inheritParams export_abnFit
#'
#' @details This function handles abnFit objects fitted using Bayesian methods.
#' It will extract the posterior distributions and other Bayesian-specific information.
#'
#' The structure will follow the same variables/parameters/arcs format, but parameters
#' will include posterior summaries (mean, median, credible intervals) instead of
#' point estimates and standard errors.
#'
#' TODO: Implement the full extraction logic for Bayesian models, including:
#' - Posterior mean/median for parameters
#' - Credible intervals
#' - Convergence diagnostics (Rhat, ESS)
#' - Prior specifications
#'
#' @return A named list with components: scenario_id, label, variables, parameters, arcs.
#' @keywords internal
export_abnFit_bayes <- function(object, format, include_network,
                                scenario_id = NULL, label = NULL, ...) {
  # Input validation
  if (!inherits(object, "abnFit")) {
    stop("Object must be of class 'abnFit'", call. = FALSE)
  }

  if (object$method != "bayes") {
    stop("This function only handles abnFit objects fitted with method = 'bayes'", call. = FALSE)
  }

  # TODO: Extract variables, parameters, and arcs from Bayesian fit
  # This should:
  # 1. Extract posterior distributions for parameters
  # 2. Compute summary statistics (mean, median, quantiles)
  # 3. Include convergence diagnostics
  # 4. Format according to variables/parameters/arcs structure

  warning("Bayesian model export is not fully implemented yet. Returning placeholder structure.")

  # Placeholder structure
  export_structure <- list()
  export_structure$scenario_id <- scenario_id
  export_structure$label <- label
  export_structure$variables <- list()
  export_structure$parameters <- list()
  export_structure$arcs <- list()

  # TODO: Populate with actual Bayesian model information

  return(export_structure)
}
