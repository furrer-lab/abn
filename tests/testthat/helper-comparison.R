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
#'   (default: \code{.Machine$double.eps^0.5}, ~1.5e-8). Tightened from the
#'   previous 1e-6 to catch precision regressions in JSON round-trip.
#'
#' @return Logical TRUE if objects are equivalent, FALSE otherwise
#'
#' @keywords internal
abnfit_objects_equal <- function(obj1, obj2, tolerance = .Machine$double.eps^0.5) {
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

  is_grouped1 <- !is.null(obj1$group.var) || !is.null(obj1$mu)
  is_grouped2 <- !is.null(obj2$group.var) || !is.null(obj2$mu)

  if (!is_grouped1 && !is_grouped2) {
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
  }

  # Check scenario_id
  if (!identical(obj1$scenario_id, obj2$scenario_id)) {
    return(FALSE)
  }

  # Check label
  if (!identical(obj1$label, obj2$label)) {
    return(FALSE)
  }

  # If both have multinomial.states, check them. Older direct fits may omit this
  # derived helper even when multinomial states can be reconstructed from data.
  if (!is.null(obj1$multinomial.states) && !is.null(obj2$multinomial.states)) {
    if (!identical(obj1$multinomial.states, obj2$multinomial.states)) {
      return(FALSE)
    }
  }

  # Grouped (mixed-effects) MLE objects expose extra components. Only compare
  # them when at least one side advertises grouping; otherwise pre-grouped
  # tests (which only populate $coef / $Stderror) keep working unchanged.
  if (is_grouped1 || is_grouped2) {
    for (fld in c("mu", "betas", "sigma", "sigma_alpha")) {
      v1 <- obj1[[fld]]
      v2 <- obj2[[fld]]
      if (is.null(v1) && is.null(v2)) next
      if (xor(is.null(v1), is.null(v2))) {
        message("Grouped field '", fld, "' presence mismatch")
        return(FALSE)
      }
      if (length(v1) != length(v2) || !setequal(names(v1), names(v2))) {
        message("Grouped field '", fld, "' node-key mismatch")
        return(FALSE)
      }
      for (nm in names(v1)) {
        a <- v1[[nm]]
        b <- v2[[nm]]
        if (is.null(a) && is.null(b)) next
        if (absent_grouped_value(a) && absent_grouped_value(b)) {
          next
        }
        if (is.matrix(a) || is.matrix(b)) {
          if (identical(fld, "betas") && is.matrix(a) && is.matrix(b) &&
              !is.null(rownames(a)) && !is.null(rownames(b)) &&
              !setequal(rownames(a), rownames(b))) {
            normalize_multinomial_beta_rows <- function(rows, node) {
              rows <- sub(paste0("^", node, "\\."), "", rows)
              rows <- sub(paste0("^", node), "", rows)
              rows
            }
            rownames(a) <- normalize_multinomial_beta_rows(rownames(a), nm)
            rownames(b) <- normalize_multinomial_beta_rows(rownames(b), nm)
          }
          if (!compare_coefficient_matrices(a, b, tolerance = tolerance)) {
            message("Grouped field '", fld, "' for node '", nm, "': matrix mismatch")
            return(FALSE)
          }
        } else {
          if (!compare_grouped_numeric(a, b, tolerance = tolerance)) {
            message("Grouped field '", fld, "' for node '", nm, "': value mismatch")
            return(FALSE)
          }
        }
      }
    }
  }

  return(TRUE)
}

absent_grouped_value <- function(x) {
  length(x) == 0 || (length(x) == 1 && is.na(x))
}

compare_grouped_numeric <- function(a, b, tolerance = .Machine$double.eps^0.5) {
  if (length(a) <= 1 && length(b) <= 1) {
    return(isTRUE(all.equal(unname(as.numeric(a)), unname(as.numeric(b)),
                            tolerance = tolerance)))
  }
  if (!is.null(names(a)) && !is.null(names(b))) {
    if (!setequal(names(a), names(b))) return(FALSE)
    b <- b[names(a)]
  } else if (xor(is.null(names(a)), is.null(names(b)))) {
    return(FALSE)
  }
  isTRUE(all.equal(unname(as.numeric(a)), unname(as.numeric(b)),
                   tolerance = tolerance))
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
    if (!identical(unname(unlist(dag1$data.dists[var])),
                   unname(unlist(dag2$data.dists[var])))) {
      return(FALSE)
    }
  }

  # Check DAG structure (reorder rows/cols to match, then compare)
  dag1_reorder <- dag1$dag[var_names, var_names]
  dag2_reorder <- dag2$dag[var_names, var_names]

  if (!identical(dimnames(dag1_reorder), dimnames(dag2_reorder))) {
    return(FALSE)
  }
  if (!isTRUE(all.equal(unname(as.numeric(dag1_reorder)),
                        unname(as.numeric(dag2_reorder)),
                        tolerance = 0))) {
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
