#' Export abnFit object to structured format
#'
#' @param object An object of class abnFit
#' @param format The export format, currently only "json" is supported.
#' @param include_network Whether to include network structure
#' @param file Optional file path to save JSON. If NULL, returns JSON string
#' @param pretty Logical, whether to format JSON nicely
#' @param ... Additional export options
#' @return A JSON string or writes to file if 'file' is provided
export_abnFit <- function(object, format = "json", include_network = TRUE, file = NULL, pretty = TRUE, ...) {
  if (!inherits(object, "abnFit")) {
    stop("Input object must be of class 'abnFit'")
  }

  # Dispatch based on method
  if (object$method == "mle") {
    export_list <- export_abnFit_mle(object, format, include_network, ...)
  } else if (object$method == "bayes") {
    export_list <- export_abnFit_bayes(object, format, include_network, ...)
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
#' @param export_list The list to convert to JSON. Must contain graph, nodes, and arcs components, see details.
#' @inheritParams export_abnFit
#' @details The export_list must be a named list with the following components:
#' \itemize{
#' \item graph: A named list with metadata about the graph (e.g., method, n_nodes, n_observations, groupVar)
#' \item nodes: A named list where each element represents a node with its parameters and distribution. Each node should have:
#'  \itemize{
#'  \item label: Character string. The display name/label of the node.
#'  \item distribution: Character string. The statistical distribution type (e.g., "gaussian", "binomial", "poisson", "multinomial").
#'  \item parameterisation: A named list containing the estimated parameters for this node. The structure depends on the distribution type and fitting method.
#'  }
#'  \item arcs: A named list with arc details, including:
#'  \itemize{
#'  \item source: Character vector of source node IDs.
#'  \item target: Character vector of target node IDs.
#'  \item frequency: Numeric vector of arc frequencies (if available).
#'  \item significance: Numeric vector of arc significance values (if available).
#'  \item constraint: Character vector indicating if the arc was constrained (if available).
#'  }
#'  }
#' @keywords internal
export_to_json <- function(export_list, format, file = NULL, pretty = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for JSON export. Please install it.")
  }

  # Validate export list structure
  if (!is.list(export_list) || !all(c("graph", "nodes", "arcs") %in% names(export_list))) {
    stop("export_list must be a named list with components: graph, nodes, arcs")
  }
  if (!all(c("source", "target") %in% names(export_list$arcs))) {
    stop("export_list$arcs must contain at least 'source' and 'target' components")
  }
  if (!all(c("label", "distribution", "parameterisation") %in% names(export_list$nodes[[1]]))) {
    stop("Each node in export_list$nodes must contain 'label', 'distribution', and 'parameterisation' components")
  }
  if (!all(c("method", "n_nodes", "n_observations", "scores", "groupVar") %in% names(export_list$graph))) {
    stop("export_list$graph must contain 'method', 'n_nodes', 'n_observations', 'scores', and 'groupVar' components")
  }

  # Convert to JSON
  json_str <- jsonlite::toJSON(export_list, auto_unbox = TRUE, pretty = pretty)

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
#' without mixed-effects. It extracts graph metadata, arc details, and node parameterisations.
#' @return A named list with components: graph, nodes, arcs.
#' @keywords internal
export_abnFit_mle <- function(object, format, include_network, ...) {
  # Extract graph metadata
  graph_meta <- export_abnFit_mle_graph(object)

  # Extract arc details
  arcs_details <- export_abnFit_mle_arcs(object)

  # Extract node details
  if (!is.null(object$group.var)) {
    # With grouping
    node_details <- export_abnFit_mle_grouped_nodes(object, format, include_network, ...)
  } else if (is.null(object$group.var)) {
    # Without grouping
    node_details <- export_abnFit_mle_nodes(object, format, include_network, ...)
  } else {
    stop("Unexpected state for group.var in abnFit object.")
  }

  # Create export list
  return(
    list(
      graph = graph_meta,
      nodes = node_details,
      arcs = arcs_details
    )
  )
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

#' Helper function to extract parameters based on distribution type
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
  # Initialize output
  node_list <- list()
  node_names <- names(object$mu)
  node_dists <- object$abnDag$data.dists

  # Process each node
  for (node_id in node_names) {
    # Extract parameters (implementation later)

    # Get distribution type for this node
    distribution <- node_dists[[node_id]]

    # Extract parameters based on distribution type (implementation later)
    param_list <- list()  # Placeholder, implement extraction logic later

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
      parameterisation = NULL  # Placeholder, implement extraction logic later
    )
  }

  return(node_list)
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
      group_var = if(!is.null(object$group.var)) object$group.var else NULL
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
