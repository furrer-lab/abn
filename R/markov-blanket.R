#' Compute the Markov blanket
#'
#' This function computes the Markov blanket of a set of nodes given a DAG (Directed Acyclic Graph).
#'
#' @param dag a matrix or a formula statement (see details for format) defining the network structure, a directed acyclic graph (DAG).
#' @param node a character vector of the nodes for which the Markov Blanket should be returned.
#' @param data.dists a named list giving the distribution for each node in the network, see details.
#' @param data.df a data frame containing the data for the nodes in the network. Only needed if \code{dag} is a formula statement.
#'
#' @details
#' This function returns the Markov Blanket of a set of nodes given a DAG.
#'
#' The \code{dag} can be provided as a matrix where the rows and columns are the nodes names.
#' The matrix should be binary, where 1 indicates an edge from the column node (parent) to the row node (child).
#' The diagonal of the matrix should be 0 and the matrix should be acyclic.
#' The nodes names should be the same as the names of the distributions in \code{data.dists}.
#'
#' Alternatively, the \code{dag} can be provided using a formula statement (similar to glm).
#' This requires the \code{data.dists} and \code{data.df} arguments to be provided.
#' A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}.
#' The formula statement have to start with \code{~}.
#' In this example, \code{node1} has two parents (\code{parent1} and \code{parent2}).
#' \code{node2} and \code{node3} are children of the same parent (\code{parent3}).
#' The parents names have to exactly match those given in \code{name}.
#' \code{:} is the separator between either children or parents,
#' \code{|} separates children (left side) and parents (right side),
#' \code{+} separates terms, \code{.} replaces all the variables in \code{name}.
#'
#' @return character vector of node names from the Markov blanket.
#' @export
#'
#' @examples
#' ## Defining distribution and dag
#' dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian",
#'              e="binomial", f="binomial")
#' dag <- matrix(c(0,1,1,0,1,0,
#'                 0,0,1,1,0,1,
#'                 0,0,0,0,0,0,
#'                 0,0,0,0,0,0,
#'                 0,0,0,0,0,1,
#'                 0,0,0,0,0,0), nrow = 6L, ncol = 6L, byrow = TRUE)
#' colnames(dag) <- rownames(dag) <- names(dist)
#'
#' mb(dag, node = "b", data.dists = dist)
#' mb(dag, node = c("b","e"), data.dists = dist)
#' @keywords utilities
mb <- function(dag, node, data.dists=NULL, data.df=NULL) {
  # check if the dag is a valid dag
  dag <- check.valid.dag(dag = dag, data.df = data.df, is.ban.matrix = FALSE, group.var = NULL)

  # check if variables in data.dists are in the dag
  if (!is.null(data.dists)) {
    if (!all(names(data.dists)%in%colnames(dag))) {
      stop("Incorrect data.dists specification to compute the Markov Blanket")
    } else {
      # check order of variables in data.dists matches the order of variables in the dag
      if (!all(names(data.dists) == colnames(dag))) {
        warning("Incorrect order of variables in data.dists to compute the Markov Blanket. The order of variables in data.dists will be fixed to match the order of variables in the DAG.")

        # fix the order of variables in data.dists to match the order of variables in the dag
        data.dists <- data.dists[colnames(dag)]
      }
    }
  } else {
    stop("data.dists is required to compute the Markov Blanket.")
  }

  # check if the node is in the dag
  if (!all(node%in%colnames(dag))) {
    stop("Incorrect node specification 'node' to compute the Markov Blanket")
  }

  ## as of older version.
  # row parent column children mb.node.final <- list()
  mb.node.tmp <- list()
  for (n.element in node) {

    ## Parent + Children
    mb.children <- list()
    mb.parent <- list()
    for (i in 1:length(dag[1, ])) {
      if (dag[i, n.element] != 0) {
        mb.children[i] <- names(dag[, n.element])[i]
      }
      if (dag[n.element, i] != 0) {
        mb.parent[i] <- names(dag[n.element, ])[i]
      }
    }
    # delete NULL element
    mb.children <- unlist(mb.children[!sapply(mb.children, is.null)])
    mb.parent <- unlist(mb.parent[!sapply(mb.parent, is.null)])

    ## Parent of children
    mb.parent.children <- list()
    for (node.children in mb.children) {
      for (i in 1:length(dag[1, ])) {
        if (dag[node.children, i] != 0) {
          mb.parent.children[i] <- names(dag[node.children, ])[i]
        }
      }
    }
    # delete NULL element
    mb.parent.children <- unlist(mb.parent.children[!sapply(mb.parent.children, is.null)])

    # add all list
    mb.node <- unlist(list(mb.children, mb.parent, mb.parent.children))

    # unique element
    mb.node <- unique(mb.node)

    # delete index node
    mb.node.wo <- NULL
    if (length(mb.node) != 0) {
      for (i in 1:length(mb.node)) {
        if (mb.node[c(i)] == n.element) {
          mb.node.wo <- mb.node[-c(i)]
        }
      }
    }
    if (is.null(mb.node.wo)) {
      mb.node.wo <- mb.node
    }

    ## store out of loop
    mb.node.tmp <- unlist(list(mb.node.tmp, mb.node.wo))
    # unique element
    mb.node.tmp <- unique(mb.node.tmp)
    # delete NULL element
    mb.node.tmp <- unlist(mb.node.tmp[!sapply(mb.node.tmp, is.null)])

  }  #EOF loop through node

  return(mb.node.tmp)
}


# EOF
