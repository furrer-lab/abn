#' Export abnFit object to structured JSON format
#'
#' @description
#' Exports a fitted Additive Bayesian Network (ABN) model to a structured JSON
#' format suitable for storage, sharing, and interoperability with other analysis
#' tools. The export includes network structure (variables and arcs) and model
#' parameters (coefficients, variances, and their associated metadata).
#'
#' @param object An object of class \code{abnFit}, typically created by \code{\link{fitAbn}}.
#' @param format Character string specifying the export format. Currently, only
#'   \code{"json"} is supported.
#' @param include_network Logical, whether to include network structure (variables
#'   and arcs). Default is \code{TRUE}.
#' @param file Optional character string specifying a file path to save the JSON
#'   output. If \code{NULL} (default), the JSON string is returned.
#' @param pretty Logical, whether to format the JSON output with indentation for
#'   readability. Default is \code{TRUE}. Set to \code{FALSE} for more compact output.
#' @param scenario_id Optional character string or numeric identifier for the model
#'   run or scenario. Useful for tracking multiple model versions or experiments.
#'   Default is \code{NULL}.
#' @param label Optional character string providing a descriptive name or label
#'   for the scenario. Default is \code{NULL}.
#' @param ... Additional export options (currently unused, reserved for future extensions).
#'
#' @return
#' If \code{file} is \code{NULL}, returns a character string containing the JSON
#' representation of the model. If \code{file} is provided, writes the JSON to
#' the specified file and invisibly returns the file path.
#'
#' @details
#' This function provides a standardized way to export fitted ABN models to JSON,
#' facilitating model sharing, archiving, and integration with external tools or
#' databases. The JSON structure is designed to be both human-readable and
#' machine-parseable, following a flat architecture to avoid deep nesting.
#'
#' ## Supported Model Types
#'
#' The function handles different model fitting methods:
#' \itemize{
#'   \item \strong{MLE without grouping}: Standard maximum likelihood estimation
#'     for all supported distributions (Gaussian, Binomial, Poisson, Multinomial).
#'     Exports fixed-effect parameters with standard errors.
#'   \item \strong{MLE with grouping}: Generalized Linear Mixed Models (GLMM)
#'     with group-level random effects. Exports both fixed effects (mu, betas)
#'     and random effects (sigma, sigma_alpha).
#'   \item \strong{Bayesian}: Placeholder for future implementation of Bayesian
#'     model exports including posterior distributions.
#' }
#'
#' ## JSON Structure Overview
#'
#' The exported JSON follows a three-component structure:
#' \itemize{
#'   \item \strong{variables}: An array where each element represents a node/variable
#'     in the network with metadata including identifier, attribute name, distribution
#'     type, and states (for categorical variables).
#'   \item \strong{parameters}: An array where each element represents a model
#'     parameter (intercepts, coefficients, variances) with associated values,
#'     standard errors, link functions, and parent variable conditions.
#'   \item \strong{arcs}: An array where each element represents a directed edge
#'     in the network, specifying source and target variable identifiers.
#' }
#'
#' Additionally, optional top-level fields \code{scenario_id} and \code{label}
#' can be used to identify and describe the model.
#'
#' @section JSON Schema:
#'
#' ### Top-Level Fields
#'
#' \describe{
#'   \item{\code{scenario_id}}{Optional string or numeric identifier for the model
#'     run. Can be \code{null}.}
#'   \item{\code{label}}{Optional descriptive name for the model. Can be \code{null}.}
#'   \item{\code{variables}}{Array of variable objects (see Variables section).}
#'   \item{\code{parameters}}{Array of parameter objects (see Parameters section).}
#'   \item{\code{arcs}}{Array of arc objects (see Arcs section).}
#' }
#'
#' ### Variables Array
#'
#' Each variable object contains:
#' \describe{
#'   \item{\code{variable_id}}{Unique identifier for the variable (string). This ID
#'     is used throughout the JSON to reference this variable in parameters' \code{source}
#'     fields and in arcs' \code{source_variable_id}/\code{target_variable_id} fields.}
#'   \item{\code{attribute_name}}{Original attribute name from the data (string).}
#'   \item{\code{model_type}}{Distribution type: \code{"gaussian"}, \code{"binomial"},
#'     \code{"poisson"}, or \code{"multinomial"}.}
#'   \item{\code{states}}{Array of state objects for multinomial variables only.
#'     Each state has \code{state_id} (used to reference specific categories in
#'     parameters), \code{value_name} (the category label), and \code{is_baseline}
#'     (whether this is the reference category). \code{NULL} for continuous variables.}
#' }
#'
#' ### Parameters Array
#'
#' Each parameter object contains:
#' \describe{
#'   \item{\code{parameter_id}}{Unique identifier for the parameter (string).}
#'   \item{\code{name}}{Parameter name (e.g., \code{"intercept"}, \code{"prob_2"},
#'     coefficient name, \code{"sigma"}, \code{"sigma_alpha"}).}
#'   \item{\code{link_function_name}}{Link function: \code{"identity"} (Gaussian),
#'     \code{"logit"} (Binomial, Multinomial), or \code{"log"} (Poisson).}
#'   \item{\code{source}}{Object identifying which variable and state this parameter
#'     belongs to. Contains \code{variable_id} (required, references a variable from
#'     the variables array) and optional \code{state_id} (references a specific state
#'     for category-specific parameters in multinomial models).}
#'   \item{\code{coefficients}}{Array of coefficient objects (typically length 1),
#'     each with \code{value}, \code{stderr} (or \code{NULL} for mixed models),
#'     \code{condition_type}, and \code{conditions} array.}
#' }
#'
#' #### Coefficient Condition Types
#'
#' \itemize{
#'   \item \code{"intercept"}: Baseline parameter with no parent dependencies
#'   \item \code{"linear_term"}: Effect of a parent variable
#'   \item \code{"CPT_combination"}: Conditional probability table entry (future use)
#'   \item \code{"variance"}: Residual variance (Gaussian/Poisson only)
#'   \item \code{"random_variance"}: Random effect variance (mixed models)
#'   \item \code{"random_covariance"}: Random effect covariance (multinomial mixed models)
#' }
#'
#' ### Arcs Array
#'
#' Each arc object contains:
#' \describe{
#'   \item{\code{source_variable_id}}{Identifier of the parent/source node.}
#'   \item{\code{target_variable_id}}{Identifier of the child/target node.}
#' }
#'
#' @section Design Rationale:
#'
#' The JSON structure uses a flat architecture with three parallel arrays rather
#' than deeply nested objects. This design offers several advantages:
#' \itemize{
#'   \item \strong{Database compatibility}: Easy to store in relational or document
#'     databases with minimal transformation.
#'   \item \strong{Extensibility}: New parameter types or metadata can be added
#'     without restructuring existing fields.
#'   \item \strong{Parsability}: Simpler to query and transform programmatically.
#'   \item \strong{Flexibility}: Supports both CPT-style and GLM(M)-style models
#'     through the polymorphic \code{source} and \code{conditions} structure.
#' }
#'
#' Parameters are linked to variables through the \code{source.variable_id} field,
#' with optional \code{source.state_id} for category-specific parameters in
#' multinomial models. Parent dependencies are encoded in the \code{conditions}
#' array within each coefficient.
#'
#' @examples
#' \dontrun{
#' # Load example data and fit a model
#' library(abn)
#' data(ex1.dag.data)
#'
#' # Define distributions
#' mydists <- list(b1 = "binomial", p1 = "poisson", g1 = "gaussian",
#'                 b2 = "binomial", p2 = "poisson", g2 = "gaussian",
#'                 b3 = "binomial", g3 = "gaussian")
#'
#' # Build score cache
#' mycache <- buildScoreCache(data.df = ex1.dag.data,
#'                             data.dists = mydists,
#'                             method = "mle",
#'                             max.parents = 2)
#'
#' # Find most probable DAG
#' mp_dag <- mostProbable(score.cache = mycache)
#'
#' # Fit the model
#' myfit <- fitAbn(object = mp_dag, method = "mle")
#'
#' # Export to JSON string with metadata
#' json_export <- export_abnFit(myfit,
#'                              scenario_id = "example_model_v1",
#'                              label = "Example ABN Model")
#'
#' # View the structure
#' library(jsonlite)
#' parsed <- fromJSON(json_export)
#' str(parsed, max.level = 2)
#'
#' # Export to file
#' export_abnFit(myfit,
#'               file = "my_abn_model.json",
#'               scenario_id = "example_model_v1",
#'               label = "Example ABN Model",
#'               pretty = TRUE)
#'
#' # Export with compact formatting
#' compact_json <- export_abnFit(myfit, pretty = FALSE)
#'
#' # ---
#' # Mixed-effects model example
#' # (Requires data with grouping structure)
#'
#' # Add grouping variable
#' ex1.dag.data$group <- rep(1:5, length.out = nrow(ex1.dag.data))
#'
#' # Build cache with grouping
#' mycache_grouped <- buildScoreCache(data.df = ex1.dag.data,
#'                                     data.dists = mydists,
#'                                     method = "mle",
#'                                     group.var = "group",
#'                                     max.parents = 2)
#'
#' # Fit grouped model
#' myfit_grouped <- fitAbn(object = mp_dag,
#'                         method = "mle",
#'                         group.var = "group")
#'
#' # Export grouped model (includes random effects)
#' json_grouped <- export_abnFit(myfit_grouped,
#'                               scenario_id = "grouped_model_v1",
#'                               label = "Mixed Effects ABN")
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{fitAbn}} for fitting ABN models
#'   \item \code{\link{buildScoreCache}} for structure learning
#'   \item \code{\link{mostProbable}} for finding the most probable network structure
#' }
#'
#' @export
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
  dag_matrix <- as.matrix(object$abnDag$dag)

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

    # Get parent nodes for this child (child is node_id, parents are columns with 1s)
    node_idx <- which(colnames(dag_matrix) == node_id)
    parent_nodes <- names(dag_matrix[node_idx, ])[dag_matrix[node_idx, ] == 1]

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

  return(if (is.null(link_functions[[distribution]])) "identity" else link_functions[[distribution]])
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
  dag_matrix <- as.matrix(object$abnDag$dag)

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
    node_idx <- which(colnames(dag_matrix) == node_id)
    parent_nodes <- names(dag_matrix[node_idx, ])[dag_matrix[node_idx, ] == 1]

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

#' Helper function to extract parameters for mixed-effects models
#'
#' @param mu Fixed-effect intercept(s) from mu component
#' @param betas Fixed-effect coefficients from betas component
#' @param sigma Residual variance from sigma component
#' @param sigma_alpha Random-effect variance/covariance from sigma_alpha component
#' @param distribution Node distribution type
#' @param node_id Node identifier
#' @param parent_nodes Parent node identifiers
#' @param start_counter Starting parameter counter
#' @param link_function Link function name
#'
#' @details Extracts parameters from mixed-effects models following the new JSON structure.
#'
#' For each node, creates parameters for:
#' - Fixed-effect intercept (from mu)
#' - Fixed-effect coefficients for parents (from betas)
#' - Residual variance (from sigma, for Gaussian/Poisson)
#' - Random-effect variance/covariance (from sigma_alpha)
#'
#' @keywords internal
extract_parameters_mixed_effects <- function(mu, betas, sigma, sigma_alpha,
                                             distribution, node_id, parent_nodes,
                                             start_counter, link_function) {
  parameters <- list()
  counter <- start_counter

  # Handle different distribution types
  if (distribution %in% c("gaussian", "binomial", "poisson")) {
    # 1. Fixed-effect intercept (mu)
    if (!is.null(mu) && !is.na(mu) && length(mu) > 0) {
      param_entry <- list(
        parameter_id = as.character(counter),
        name = "intercept",
        link_function_name = link_function,
        source = list(variable_id = node_id),
        coefficients = list(
          list(
            value = as.numeric(mu)[1],
            stderr = NULL,
            condition_type = "intercept",
            conditions = list()
          )
        )
      )
      parameters[[length(parameters) + 1]] <- param_entry
      counter <- counter + 1
    }

    # 2. Fixed-effect coefficients (betas)
    if (!is.null(betas) && !is.logical(betas) && length(betas) > 0 && !all(is.na(betas))) {
      beta_names <- names(betas)
      for (i in seq_along(betas)) {
        beta_name <- beta_names[i]
        beta_value <- betas[i]

        # Match with parent nodes
        parent_var <- NULL
        parent_state <- NULL
        for (p in parent_nodes) {
          if (grepl(paste0("^", p), beta_name) || grepl(paste0(p, "[0-9]"), beta_name)) {
            parent_var <- p
            state_match <- regmatches(beta_name, regexpr("[0-9]+$", beta_name))
            parent_state <- if (length(state_match) > 0) state_match else NULL
            break
          }
        }

        if (!is.null(parent_var)) {
          param_entry <- list(
            parameter_id = as.character(counter),
            name = beta_name,
            link_function_name = link_function,
            source = list(variable_id = node_id),
            coefficients = list(
              list(
                value = as.numeric(beta_value),
                stderr = NULL,
                condition_type = "linear_term",
                conditions = list(
                  list(
                    parent_variable_id = parent_var,
                    parent_state_id = parent_state
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

    # 3. Residual variance (sigma) - only for Gaussian and Poisson
    if (distribution %in% c("gaussian", "poisson")) {
      if (!is.null(sigma) && !is.logical(sigma) && length(sigma) > 0 && !all(is.na(sigma))) {
        param_entry <- list(
          parameter_id = as.character(counter),
          name = "sigma",
          link_function_name = "identity",
          source = list(variable_id = node_id),
          coefficients = list(
            list(
              value = as.numeric(sigma)[1],
              stderr = NULL,
              condition_type = "variance",
              conditions = list()
            )
          )
        )
        parameters[[length(parameters) + 1]] <- param_entry
        counter <- counter + 1
      }
    }

    # 4. Random-effect variance (sigma_alpha)
    if (!is.null(sigma_alpha) && !all(is.na(sigma_alpha))) {
      param_entry <- list(
        parameter_id = as.character(counter),
        name = "sigma_alpha",
        link_function_name = "identity",
        source = list(variable_id = node_id),
        coefficients = list(
          list(
            value = if (is.matrix(sigma_alpha)) as.numeric(sigma_alpha[1,1]) else as.numeric(sigma_alpha)[1],
            stderr = NULL,
            condition_type = "random_variance",
            conditions = list()
          )
        )
      )
      parameters[[length(parameters) + 1]] <- param_entry
      counter <- counter + 1
    }

  } else if (distribution == "multinomial") {
    # For multinomial, mu contains category-specific intercepts
    if (!is.null(mu) && length(mu) > 0) {
      mu_names <- names(mu)
      categories <- gsub(paste0(".*", node_id, "\\."), "", mu_names)

      for (i in seq_along(mu)) {
        cat <- categories[i]
        param_entry <- list(
          parameter_id = as.character(counter),
          name = paste0("intercept_", cat),
          link_function_name = link_function,
          source = list(variable_id = node_id, state_id = cat),
          coefficients = list(
            list(
              value = as.numeric(mu[i]),
              stderr = NULL,
              condition_type = "intercept",
              conditions = list()
            )
          )
        )
        parameters[[length(parameters) + 1]] <- param_entry
        counter <- counter + 1
      }
    }

    # Fixed-effect coefficients (betas) - matrix format
    if (!is.null(betas) && is.matrix(betas) && !all(is.na(betas))) {
      categories <- rownames(betas)
      parent_names <- colnames(betas)

      for (i in seq_len(nrow(betas))) {
        cat <- categories[i]
        for (j in seq_len(ncol(betas))) {
          parent_name <- parent_names[j]

          parent_var <- NULL
          parent_state <- NULL
          for (p in parent_nodes) {
            if (grepl(paste0("^", p), parent_name)) {
              parent_var <- p
              state_match <- regmatches(parent_name, regexpr("[0-9]+$", parent_name))
              parent_state <- if (length(state_match) > 0) state_match else NULL
              break
            }
          }

          if (!is.null(parent_var)) {
            param_entry <- list(
              parameter_id = as.character(counter),
              name = paste0(parent_name, "_cat", cat),
              link_function_name = link_function,
              source = list(variable_id = node_id, state_id = cat),
              coefficients = list(
                list(
                  value = as.numeric(betas[i, j]),
                  stderr = NULL,
                  condition_type = "linear_term",
                  conditions = list(
                    list(
                      parent_variable_id = parent_var,
                      parent_state_id = parent_state
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

    # Random-effect variance-covariance matrix
    if (!is.null(sigma_alpha) && is.matrix(sigma_alpha)) {
      categories <- rownames(sigma_alpha)

      for (i in seq_len(nrow(sigma_alpha))) {
        for (j in i:ncol(sigma_alpha)) {
          cat_i <- gsub(".*~", "", categories[i])
          cat_j <- gsub(".*~", "", categories[j])

          param_entry <- list(
            parameter_id = as.character(counter),
            name = if (i == j) paste0("sigma_alpha_", cat_i) else paste0("sigma_alpha_", cat_i, "_", cat_j),
            link_function_name = "identity",
            source = list(
              variable_id = node_id,
              state_id = if (i == j) cat_i else paste0(cat_i, "_", cat_j)
            ),
            coefficients = list(
              list(
                value = as.numeric(sigma_alpha[i, j]),
                stderr = NULL,
                condition_type = if (i == j) "random_variance" else "random_covariance",
                conditions = list()
              )
            )
          )
          parameters[[length(parameters) + 1]] <- param_entry
          counter <- counter + 1
        }
      }
    }
  }

  return(list(parameters = parameters, next_counter = counter))
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
