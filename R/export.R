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
#'   that were fitted using the Maximum Likelihood Estimation (MLE) approach without
#'   mixed-effects (i.e., no group.var specified). The function processes the coefficients
#'   and standard errors stored in the abnFit object.
#'
#'   The \code{coef} component contains the estimated regression coefficients for each node,
#'   stored as a matrix where column names indicate the parameter names (e.g., "g2",
#'   "m11", "b1|intercept"). These represent the linear model coefficients from the
#'   generalized linear model fitted to each node given its parents in the DAG.
#'
#'   The \code{Stderror} component contains the corresponding standard errors for each
#'   coefficient, providing a measure of uncertainty in the parameter estimates. The
#'   structure mirrors that of the \code{coef} component.
#'
#'   For different distribution types:
#'   \itemize{
#'     \item Gaussian nodes: Include intercept and slope coefficients
#'     \item Binomial/Poisson nodes: Include intercept and slope coefficients on logit/log scale
#'     \item Multinomial nodes: Include category-specific intercepts (reference level omitted)
#'           and coefficients, following standard multinomial logistic regression conventions
#'   }
#'
#' @returns A named list where each element represents a node in the network.
#'   Each node is identified by its node ID and contains the following components:
#'   \item{nodeid}{A list representing a single node with components:}
#'   \item{label}{Character string. The display name/label of the node.}
#'   \item{distribution}{Character string. The statistical distribution type
#'     (e.g., "gaussian", "binomial", "poisson", "multinomial").}
#'   \item{df}{A named vector of degrees of freedom for the node's model parameters, if available.}
#'   \item{mse}{Mean Squared Error of the node's model, if available.}
#'   \item{sse}{Sum of Squared Errors of the node's model, if available.}
#'   \item{parameterisation}{A named list containing the estimated parameters
#'     for this node. The structure depends on the distribution type and fitting method:
#'     \itemize{
#'       \item For Gaussian nodes: \code{intercept}, \code{coefficients}, \code{stderr}
#'       \item For Binomial/Poisson nodes: \code{intercept}, \code{coefficients}, \code{stderr}
#'       \item For Multinomial nodes: category-specific parameter lists with \code{stderr}
#'     }}
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
  nodes_list <- list()
  node_names <- names(object$coef)
  node_dists <- object$abnDag$data.dists

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

    # Extract parameters based on distribution type
    param_list <- extract_parameters_by_distribution(coef_vec, se_vec, distribution, node_id)

    # Get degree of freedom if available
    if (!is.null(object$df)) {
      df_int <- object$df[[node_id]]
    }

    # Get mse
    if (!is.null(object$mse)) {
      mse_val <- object$mse[[node_id]]
    }

    # Get sse
    if (!is.null(object$sse)) {
      sse_val <- object$sse[[node_id]]
    }

    # Create node entry
    nodes_list[[node_id]] <- list(
      label = node_id,
      distribution = distribution,
      df = if (exists("df_int")) df_int else NULL,
      mse = if (exists("mse_val")) mse_val else NULL,
      sse = if (exists("sse_val")) sse_val else NULL,
      parameterisation = param_list
    )
  }

  return(nodes_list)
}

#' Helper function to extract parameters based on distribution type without grouping
#' @keywords internal
extract_parameters_by_distribution <- function(coef_vec, se_vec, distribution, node_id) {

  param_names <- names(coef_vec)

  # Initialize parameter list
  param_list <- list()

  if (distribution %in% c("gaussian", "binomial", "poisson")) {

    # Find intercept parameter
    intercept_pattern <- paste0(node_id, "\\|intercept")
    intercept_idx <- grep(intercept_pattern, param_names, ignore.case = TRUE)

    if (length(intercept_idx) > 0) {
      intercept_name <- param_names[intercept_idx[1]]
      param_list$intercept <- list(
        estimate = coef_vec[intercept_name],
        stderr = se_vec[intercept_name]
      )
      # Remove intercept from remaining parameters
      param_names <- param_names[-intercept_idx[1]]
      coef_remaining <- coef_vec[param_names]
      se_remaining <- se_vec[param_names]
    } else {
      coef_remaining <- coef_vec
      se_remaining <- se_vec
    }

    # Store remaining coefficients (parent effects)
    if (length(coef_remaining) > 0) {
      param_list$coefficients <- list()
      for (i in seq_along(coef_remaining)) {
        param_name <- names(coef_remaining)[i]
        param_list$coefficients[[param_name]] <- list(
          estimate = coef_remaining[i],
          stderr = se_remaining[i]
        )
      }
    }

  } else if (distribution == "multinomial") {

    # For multinomial, parameters are category-specific
    # Extract intercepts and coefficients by category
    intercept_pattern <- paste0(node_id, "\\|intercept\\.")
    intercept_idx <- grep(intercept_pattern, param_names)

    # Group parameters by category
    param_list$categories <- list()

    if (length(intercept_idx) > 0) {
      # Extract category numbers from intercept names
      intercept_names <- param_names[intercept_idx]
      categories <- gsub(paste0(".*", node_id, "\\|intercept\\."), "", intercept_names)

      for (cat in unique(categories)) {
        cat_pattern <- paste0("\\.", cat, "$")
        cat_params <- grep(cat_pattern, param_names)

        if (length(cat_params) > 0) {
          param_list$categories[[paste0("category_", cat)]] <- list()

          for (idx in cat_params) {
            param_name <- param_names[idx]
            clean_name <- gsub(paste0(node_id, "\\|"), "", param_name)
            param_list$categories[[paste0("category_", cat)]][[clean_name]] <- list(
              estimate = coef_vec[idx],
              stderr = se_vec[idx]
            )
          }
        }
      }
    }

    # Handle non-intercept multinomial parameters (parent effects)
    non_intercept_idx <- setdiff(seq_along(param_names), intercept_idx)
    if (length(non_intercept_idx) > 0) {
      remaining_params <- param_names[non_intercept_idx]

      for (idx in non_intercept_idx) {
        param_name <- param_names[idx]
        # Extract category information if present
        if (grepl("\\d+$", param_name)) {
          category <- gsub(".*([0-9]+)$", "\\1", param_name)
          clean_name <- gsub("([0-9]+)$", "", param_name)

          if (is.null(param_list$categories[[paste0("category_", category)]])) {
            param_list$categories[[paste0("category_", category)]] <- list()
          }

          param_list$categories[[paste0("category_", category)]][[clean_name]] <- list(
            estimate = coef_vec[idx],
            stderr = se_vec[idx]
          )
        }
      }
    }
  }

  return(param_list)
}

#' Export node information from abnFit objects fitted with MLE (mixed effects)
#' @param object An object of class abnFit fitted with method = "mle" and group.var specified.
#' @param ... Additional arguments (currently unused)
#' @details This function is a placeholder for exporting node information from abnFit objects.
#'  Currently, it raises an error indicating that export for grouped MLE models is not implemented.
#' @return This function does not return a value. It raises an error.
#' @keywords internal
export_abnFit_mle_grouped_nodes <- function(object, ...) {
  # Input validation
  if (is.null(object$mu) || is.null(object$betas) ||
      is.null(object$sigma) || is.null(object$sigma_alpha)) {
    stop("abnFit object must contain 'mu', 'betas', 'sigma', and 'sigma_alpha' components", call. = FALSE)
  }

  # Initialize output
  node_list <- list()
  node_names <- names(object$mu)
  node_dists <- object$abnDag$data.dists

  # Process each node
  for (node_id in node_names) {
    # Extract parameters for this node
    mu_val <- object$mu[[node_id]]
    betas_val <- object$betas[[node_id]]
    sigma_val <- object$sigma[[node_id]]
    sigma_alpha_val <- object$sigma_alpha[[node_id]]

    # Get distribution type for this node
    distribution <- node_dists[[node_id]]

    # Extract parameters based on distribution type
    param_list <- extract_parameters_by_distribution_grouping(
      mu_val, betas_val, sigma_val, sigma_alpha_val, distribution, node_id
    )
    # Get degree of freedom if available
    if (!is.null(object$df)) {
      df_int <- object$df[[node_id]]
    }

    # Get mse
    if (!is.null(object$mse)) {
      mse_val <- object$mse[[node_id]]
    }

    # Get sse
    if (!is.null(object$sse)) {
      sse_val <- object$sse[[node_id]]
    }

    # Create node entry
    node_list[[node_id]] <- list(
      label = node_id,
      distribution = node_dists[[node_id]],
      df = if (exists("df_int")) df_int else NULL,
      mse = if (exists("mse_val")) mse_val else NULL,
      sse = if (exists("sse_val")) sse_val else NULL,
      parameterisation = param_list
    )
  }

  return(node_list)
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

#' Export graph metadata from abnFit objects fitted with MLE
#' @inheritParams export_abnFit
#' @details This function extracts metadata about the graph from abnFit objects.
#' It retrieves information such as the fitting method, number of nodes, number of observations,
#' and model scores (AIC, BIC, log-likelihood). If a grouping variable was used during fitting,
#' it is also included.
#' @return A named list containing graph metadata.
#' @keywords internal
export_abnFit_mle_graph <- function(object, ...) {
  return(
    list(
      method = object$method,
      n_nodes = ncol(object$abnDag$data.df),
      n_observations = nrow(object$abnDag$data.df),
      scores = list(
        aic = object$aic,
        bic = object$bic,
        mlik = object$mlik,
        mdl = object$mdl
      ),
      groupVar = if(!is.null(object$group.var)) object$group.var else NULL,
      groupedVariables = if(!is.null(object$grouped.vars)) names(object$abnDag$data.dists)[object$grouped.vars] else NULL
    )
  )
}

#' Export arc information from abnFit objects fitted with MLE
#' @inheritParams export_abnFit
#' @details This function extracts arc information from abnFit objects fitted using MLE.
#' It retrieves the source and target nodes for each arc in the directed acyclic graph (DAG).
#' Currently, frequency, significance, and constraint information are not included in the export.
#' @return A named list containing arc details: source nodes, target nodes, frequency, significance, and constraint.
#' @keywords internal
export_abnFit_mle_arcs <- function(object, ...) {
  edgelist <- as.data.frame(object$abnDag)
  warning("Arc export currently does not include frequency, significance, or constraint information.")
  arcs <- list(
    source = edgelist$from,
    target = edgelist$to,
    frequency = NULL,
    significance = NULL,
    constraint = NULL
  )
  return(arcs)
}

#' Export abnFit object fitted with Bayesian methods
#' @inheritParams export_abnFit
#' @details This function is a placeholder for exporting abnFit objects fitted using Bayesian methods.
#' Currently, it raises an error indicating that export for Bayesian models is not implemented.
#' @return This function does not return a value. It raises an error.
#' @keywords internal
export_abnFit_bayes <- function(object, format, include_network, ...) {
  # Extract graph metadata

  # Extract arc details

  # Extract node details
  stop("Export for Bayesian models is not implemented yet.")
}
