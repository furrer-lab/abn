#' Compare two abnFit objects for equivalence
#'
#' Compares abnFit objects ignoring parameter and arc ordering.
#' This is useful for testing round-trip export/import operations where
#' the order of parameters or arcs may change but the model structure remains
#' equivalent.
#'
#' @param obj1 First abnFit object to compare
#' @param obj2 Second abnFit object to compare
#' @param tolerance Numerical tolerance for comparing coefficient values
#'   (default: 1e-6)
#'
#' @return Logical TRUE if objects are equivalent, FALSE otherwise
#'
#' @keywords internal
abnfit_objects_equal <- function(obj1, obj2, tolerance = 1e-6) {
  # Check class
  if (!inherits(obj1, "abnFit") || !inherits(obj2, "abnFit")) {
    message("Class check failed: obj1=", class(obj1), " obj2=", class(obj2))
    return(FALSE)
  }

  # Check method
  if (obj1$method != obj2$method) {
    message("Method mismatch: obj1$method=", obj1$method, " obj2$method=", obj2$method)
    return(FALSE)
  }

  # Check abnDag structure
  if (!abndag_equal(obj1$abnDag, obj2$abnDag)) {
    message("abnDag structure mismatch")
    return(FALSE)
  }

  # Check number of nodes
  if (length(obj1$coef) != length(obj2$coef)) {
    message("Node count mismatch: obj1 has ", length(obj1$coef), 
            " nodes, obj2 has ", length(obj2$coef), " nodes")
    return(FALSE)
  }

  # Check coefficients for each node (order-independent)
  node_names <- names(obj1$coef)
  if (!setequal(node_names, names(obj2$coef))) {
    message("Node names mismatch: obj1 nodes=", paste(node_names, collapse=","),
            " obj2 nodes=", paste(names(obj2$coef), collapse=","))
    return(FALSE)
  }

  for (node in node_names) {
    coef1 <- obj1$coef[[node]]
    coef2 <- obj2$coef[[node]]

    # Compare coefficient matrices
    if (!compare_coefficient_matrices(coef1, coef2, tolerance = tolerance)) {
      message("Coefficient mismatch for node ", node)
      message("  obj1 coef dim: ", paste(dim(coef1), collapse="x"),
              " obj2 coef dim: ", paste(dim(coef2), collapse="x"))
      return(FALSE)
    }

    # Compare standard errors
    se1 <- obj1$Stderror[[node]]
    se2 <- obj2$Stderror[[node]]

    if (!compare_coefficient_matrices(se1, se2, tolerance = tolerance)) {
      message("Standard error mismatch for node ", node)
      message("  obj1 SE dim: ", paste(dim(se1), collapse="x"),
              " obj2 SE dim: ", paste(dim(se2), collapse="x"))
      return(FALSE)
    }
  }

  # Check scenario_id
  if (!identical(obj1$scenario_id, obj2$scenario_id)) {
    return(FALSE)
  }

  # Check label
  if (!identical(obj1$label, obj2$label)) {
    return(FALSE)
  }

  # If both have multinomial.states, check them
  if (!is.null(obj1$multinomial.states) || !is.null(obj2$multinomial.states)) {
    if (!identical(obj1$multinomial.states, obj2$multinomial.states)) {
      return(FALSE)
    }
  }

  return(TRUE)
}

#' Compare two abnDag objects for equivalence
#'
#' Compares the structure of abnDag objects.
#'
#' @param dag1 First abnDag object
#' @param dag2 Second abnDag object
#'
#' @return Logical TRUE if DAGs are equivalent, FALSE otherwise
#'
#' @keywords internal
abndag_equal <- function(dag1, dag2) {
  if (!inherits(dag1, "abnDag") || !inherits(dag2, "abnDag")) {
    return(FALSE)
  }

  # Check variable names (order-independent)
  if (!setequal(colnames(dag1$dag), colnames(dag2$dag))) {
    return(FALSE)
  }

  # Check data.dists for all variables
  var_names <- colnames(dag1$dag)
  for (var in var_names) {
    if (!var %in% colnames(dag2$dag)) {
      return(FALSE)
    }
    if (dag1$data.dists[var] != dag2$data.dists[var]) {
      return(FALSE)
    }
  }

  # Check DAG structure (reorder rows/cols to match, then compare)
  dag1_reorder <- dag1$dag[var_names, var_names]
  dag2_reorder <- dag2$dag[var_names, var_names]

  if (!identical(dag1_reorder, dag2_reorder)) {
    return(FALSE)
  }

  return(TRUE)
}

#' Compare two coefficient/standard error matrices for equivalence
#'
#' Compares matrices allowing for numerical tolerance and column name matching.
#' Columns are matched by name, so order doesn't matter.
#'
#' @param mat1 First matrix
#' @param mat2 Second matrix
#' @param tolerance Numerical tolerance for comparison
#'
#' @return Logical TRUE if matrices are equivalent, FALSE otherwise
#'
#' @keywords internal
compare_coefficient_matrices <- function(mat1, mat2, tolerance = 1e-6) {
  # Both empty
  if (is.null(dim(mat1)) && is.null(dim(mat2))) {
    return(TRUE)
  }

  # Convert to matrix if needed
  if (!is.matrix(mat1)) {
    if (length(mat1) == 0) {
      mat1 <- matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
      mat1 <- as.matrix(mat1)
    }
  }

  if (!is.matrix(mat2)) {
    if (length(mat2) == 0) {
      mat2 <- matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
      mat2 <- as.matrix(mat2)
    }
  }

  # Check dimensions
  if (nrow(mat1) != nrow(mat2) || ncol(mat1) != ncol(mat2)) {
    return(FALSE)
  }

  # If both are empty
  if (nrow(mat1) == 0 && ncol(mat1) == 0) {
    return(TRUE)
  }

  # Check column names match (order-independent)
  col1 <- colnames(mat1)
  col2 <- colnames(mat2)

  if (is.null(col1) || is.null(col2)) {
    # If one has colnames and the other doesn't
    if (is.null(col1) != is.null(col2)) {
      return(FALSE)
    }
    # Both are NULL, compare as-is
    return(isTRUE(all.equal(mat1, mat2, tolerance = tolerance)))
  }

  if (!setequal(col1, col2)) {
    return(FALSE)
  }

  # Reorder mat2 to match mat1's column order
  mat2_reorder <- mat2[, col1, drop = FALSE]

  # Compare numerically
  return(isTRUE(all.equal(mat1, mat2_reorder, tolerance = tolerance)))
}
