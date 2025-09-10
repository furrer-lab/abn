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
#' \item graph: A named list with metadata about the graph (e.g., method, n_nodes, n_observations, group_var)
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
  if (!all(c("method", "n_nodes", "n_observations", "scores", "group_var") %in% names(export_list$graph))) {
    stop("export_list$graph must contain 'method', 'n_nodes', 'n_observations', 'scores', and 'group_var' components")
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

  return(json_str)
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
