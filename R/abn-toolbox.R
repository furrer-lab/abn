#' Logit of proportions
#'
#' See also the C implementation \code{?abn::logit_cpp()}.
#'
#' @param x numeric with values between \code{[0,1]}.
#'
#' @return numeric vector of same length as \code{x}.
#' @export
#' @returns numeric vector of same length as \code{x}.
#' @keywords utilities
logit <- function(x) {
  return(log(x/(1 - x)))
}

#' expit of proportions
#'
#' See also the C implementation \code{?abn::expit_cpp()}.
#'
#' @param x numeric with values between \code{[0,1]}.
#'
#' @return numeric vector of same length as \code{x}.
#' @export
#' @keywords utilities
expit <- function(x) {
  return(exp(x)/(1 + exp(x)))
}

#' Odds Ratio from a matrix
#'
#' Compute the odds ratio from a contingency table or a matrix.
#'
#' @param x a 2x2 table or matrix.
#'
#' @return A real value.
#' @export
#'
#' @keywords utilities
or <- function(x) {

    x <- as.matrix(x)

    if (dim(x)[1] != 2 || dim(x)[2] != 2) {
        stop("The contingency table should be 2-dimensionnal with a 2x2 formulation.")
    }
    if (any(x < 0)){
      stop("There are negative values in the contingency table.")
    }

    if (x[1, 2] == 0) {
        x[1, 2] <- 0.5
    }
    if (x[2, 1] == 0) {
        x[2, 1] <- 0.5
    }

    out <- (x[1, 1] * x[2, 2])/(x[1, 2] * x[2, 1])
    return(out)
}

#' Probability to odds
#'
#' @param x numeric vector of probabilities with values between \code{[0,1]}.
#'
#' @return numeric vector of same length as \code{x}.
#' @export
#' @keywords utilities
odds <- function(x) {
  if(x < 0 || x > 1){
    stop("x must be a probability in range of 0 to 1.")
  } else {
    return((x/(1 - x)))
  }
}

#' Compare two DAGs or EGs
#'
#' @description
#' Function that returns multiple graph metrics to compare two DAGs
#' or essential graphs, known as confusion matrix or error matrix.
#'
#' @usage compareDag(ref, test, node.names = NULL, checkDAG = TRUE)
#'
#' @param ref a matrix or a formula statement (see details for format) defining
#' the reference network structure, a directed acyclic graph (DAG).
#' Note that row names must be set or given in \code{node.names}
#' if the DAG is given via a formula statement.
#' @param test a matrix or a formula statement (see details for format) defining
#' the test network structure, a directed acyclic graph (DAG).
#'  Note that row names must be set or given in \code{node.names}
#'  if the DAG is given via a formula statement.
#' @param node.names a vector of names if the DAGs are given via formula, see details.
#' @param checkDAG should the DAGs be tested for DAGs (default).
#'
#' @details
#' This R function returns standard Directed Acyclic Graph comparison metrics.
#' In statistical classification, those metrics are known as a
#' confusion matrix or error matrix.
#'
#' Those metrics allows visualization of the difference between different DAGs.
#' In the case where comparing TRUTH to learned structure or two learned structures,
#' those metrics allow the user to estimate the performance of the learning algorithm.
#' In order to compute the metrics, a contingency table is computed of a
#' pondered difference of the adjacency matrices od the two graphs.
#'
#' The \code{ref} or \code{test} can be provided using a formula statement
#' (similar to GLM input).
#' A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}.
#' The formula statement have to start with \code{~}.
#' In this example, node1 has two parents (parent1 and parent2).
#' node2 and node3 have the same parent3.
#' The parents names have to exactly match those given in \code{node.names}.
#' \code{:} is the separtor between either children or parents,
#' \code{|} separates children (left side) and parents (right side),
#' \code{+} separates terms, \code{.} replaces all the variables in \code{node.names}.
#'
#' To test for essential graphs (or graphs) in general, the test for DAG
#' need to be switched off \code{checkDAG=FALSE}.
#' The function \code{compareEG()} is a wrapper to \code{compareDag(, checkDAG=FALSE)}.
#'
#' @returns
#' \describe{
#' \item{\code{TP}}{True Positive}
#' \item{\code{TN}}{True Negative}
#' \item{\code{FP}}{False Positive}
#' \item{\code{FN}}{False Negative}
#' \item{\code{CP}}{Condition Positive (ref)}
#' \item{\code{CN}}{Condition Negative (ref)}
#' \item{\code{PCP}}{Predicted Condition Positive (test)}
#' \item{\code{PCN}}{Predicted Condition Negative (test)}
#' \item{\code{True Positive Rate}}{\deqn{=\frac{\sum TP}{\sum CP}}{=\frac{\sum TP}{\sum CP}}}
#' \item{\code{False Positive Rate}}{\deqn{=\frac{\sum FP}{\sum CN}}{=\frac{\sum FP}{\sum CN}}}
#' \item{\code{Accuracy}}{\deqn{=\frac{\sum TP + \sum TN}{Total population}}{=\frac{\sum TP + \sum TN}{Total population}}}
#' \item{\code{G-measure}}{\deqn{\sqrt {{\frac {TP}{TP+FP}}\cdot {\frac {TP}{TP+FN}}}}{\sqrt {{\frac {TP}{TP+FP}}\cdot {\frac {TP}{TP+FN}}}}}
#' \item{\code{F1-Score}}{\deqn{\frac{2 \sum TP}{2 \sum TP + \sum FN + \sum FP}}{\frac{2 \sum TP}{2 \sum TP + \sum FN + \sum FP}}}
#' \item{\code{Positive Predictive Value}}{\deqn{\frac{\sum TP}{\sum PCP}}{\frac{\sum TP}{\sum PCP}}}
#' \item{\code{False Ommision Rate}}{\deqn{\frac{\sum FN}{\sum PCN}}{\frac{\sum FN}{\sum PCN}}}
#' \item{\code{Hamming-Distance}}{Number of changes needed to match the matrices.}
#' }
#'
#' @export
#' @references Sammut, Claude, and Geoffrey I. Webb. (2017). Encyclopedia of machine learning and data mining. Springer.
#' @keywords utilities
#' @examples
#' test.m <- matrix(data = c(0,1,0,
#'                           0,0,0,
#'                           1,0,0), nrow = 3, ncol = 3)
#' ref.m <- matrix(data = c(0,0,0,
#'                          1,0,0,
#'                          1,0,0), nrow = 3, ncol = 3)
#'
#' colnames(test.m) <- rownames(test.m) <- colnames(ref.m) <- colnames(ref.m) <- c("a", "b", "c")
#'
#' unlist(compareDag(ref = ref.m, test = test.m))
compareDag <- function(ref, test, node.names = NULL, checkDAG = TRUE) {

    ## check ref dag
    if (checkDAG) {
        ref <- validate_abnDag(  ref, data.df=node.names, returnDAG=TRUE)
        test <- validate_abnDag( test, data.df=node.names, returnDAG=TRUE)
    }

    if (any(dim(ref) != dim(test))) {
        stop("The reference and test DAG have not the same size")
    }

    n <- dim(ref)[1]

    ## unit matrix
    ref[ref != 0] <- 1
    test[test != 0] <- 1

    diff.matrix <- ref - (0.5 * test)

    diff.matrix.tab <- table(factor(diff.matrix, levels = c(-0.5, 0, 0.5, 1)))

    if(sum(ref == 1)==0 | sum(test == 1)==0){
        warning("If the test or reference DAG is empty, some of the estimates are not defined.")
    }

    ## output
    out <- list()

    out[["TPR"]] <- (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]))/(sum(ref == 1))
    out[["FPR"]] <- (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == -0.5]))/(sum(ref == 0))
    out[["Accuracy"]] <- (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]) + as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0]))/(dim(ref)[1]^2)
    out[["FDR"]] <- as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 1])/(sum(test == 1))
    out[["G-measure"]] <- sqrt(as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5])/(sum(test == 1)) * (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]))/(sum(ref == 1)))
    out[["F1-score"]] <- (2/((1/as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5])/(sum(test == 1))) + (1/(as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]))/(sum(ref == 1)))))
    out[["PPV"]] <- as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5])/(sum(test == 1))
    out[["FOR"]] <- as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 1])/(sum(test == 1))

    # transforming "reverse arc" into -0.5
    for( i in 1:n){
        for( j in i:n){
            if( diff.matrix[i,j]!=0 & diff.matrix[j,i]!=0){
                diff.matrix[i,j] <- -(diff.matrix[i,j] + diff.matrix[j,i])
                diff.matrix[j,i] <- 0
            }
        }
    }
    diff.matrix.tab <- table(factor(diff.matrix, levels = c(-0.5, 0, 0.5, 1)))

    out[["Hamming-distance"]] <- sum(as.numeric(diff.matrix.tab[names(diff.matrix.tab) %in% c(-0.5, 1)]))

    return(out)
}

#' @inherit compareDag
#' @usage compareEG(ref, test)
compareEG <- function(ref, test){
  compareDag(ref, test, checkDAG=FALSE)
}

#' Compute standard information for a DAG.
#'
#' This function returns standard metrics for DAG description. A list that
#' contains the number of nodes, the number of arcs,
#' the average Markov blanket size, the neighborhood average set size,
#' the parent average set size and children average set size.
#'
#' @usage infoDag(object, node.names = NULL)
#'
#' @param object an object of class \code{abnLearned}, \code{abnFit}.
#'  Alternatively, a matrix or a formula statement defining the network structure,
#'  a directed acyclic graph (DAG).
#'  Note that row names must be set up or given in \code{node.names}.
#' @param node.names a vector of names if the DAG is given via formula, see details.
#'
#' @details
#' This function returns a named list with the following entries:
#' the number of nodes, the number of arcs, the average Markov blanket size,
#' the neighborhood average set size, the parent average set size, and the
#' children's average set size.
#'
#' The \code{dag} can be provided using a formula statement (similar to glm).
#' A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}.
#' The formula statement have to start with \code{~}.
#' In this example, node1 has two parents (parent1 and parent2).
#' node2 and node3 have the same parent3.
#' The parents names have to exactly match those given in \code{node.names}.
#' \code{:} is the separator between either children or parents,
#' \code{|} separates children (left side) and parents (right side),
#' \code{+} separates terms, \code{.} replaces all the variables in \code{node.names}.
#'
#' @return A named list that contains following entries:
#' the number of nodes, the number of arcs,
#' the average Markov blanket size, the neighborhood average set size,
#' the parent average set size and children average set size.
#'
#' @references West, D. B. (2001). Introduction to graph theory. Vol. 2. Upper Saddle River: Prentice Hall.
#'
#' @export
#'
#' @examples
#' ## Creating a dag:
#' dag <- matrix(c(0,0,0,0, 1,0,0,0, 1,1,0,1, 0,1,0,0), nrow = 4, ncol = 4)
#' dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian")
#' colnames(dag) <- rownames(dag) <- names(dist)
#'
#' infoDag(dag)
#' plot(createAbnDag(dag = dag, data.dists = dist))
#' @keywords utilities
infoDag <- function(object, node.names = NULL) {

    if (inherits(object, "abnLearned")) {
      dag <- object$dag
    } else   if (inherits(object, "abnFit")) {
      dag <- object$abnDag$dag
    } else  if (is.matrix(object)) {
      dag <- abs(object)
      dag[dag != 0] <- 1
      #       diag(dag) <- 0  # RF: we want to test if proper dag!!
      dag <- check.valid.dag(dag = dag, is.ban.matrix = FALSE, group.var = NULL)
      ## naming
      if (is.null(colnames(dag))) {
          colnames(dag) <- rownames(dag) <- node.names
      }
    } else if (grepl("~", as.character(dag)[1], fixed = TRUE)) {
        dag <- formula_abn(f = dag, name = node.names)
        dag <- check.valid.dag(dag = dag, is.ban.matrix = FALSE, group.var = NULL)
    } else {
       stop("'object' must either be a matrix or a formula expression, or of class 'abnFit', 'abnLearned'.")
    }

    if (is.null(node.names)) {
        node.names <- colnames(dag)
    }

    ## ======================== test for conformity over! ========================

    out <- list()
    ## number of nodes
    n.nodes <- dim(dag)[1]

    ## number of arcs
    n.arcs <- sum(dag)

    ## =========================== average markov blanket size ===========================

    mb.size <- vector(mode = "numeric", length = length(node.names))

    for (i in 1:length(node.names)) {

        node <- node.names[i]
        # row children column parent

        ## Parent + Children
        mb.children <- list()
        mb.parent <- list()
        for (j in 1:length(dag[1, ])) {
            if (dag[node, j] != 0) {
                mb.children[j] <- names(dag[node, ])[j]
            }
            if (dag[j, node] != 0) {
                mb.parent[j] <- names(dag[, node])[j]
            }
        }
        # delete NULL element
        mb.children <- unlist(mb.children[!sapply(mb.children, is.null)])
        mb.parent <- unlist(mb.parent[!sapply(mb.parent, is.null)])

        ## Parent of children
        mb.parent.children <- list()
        for (node.children in mb.children) {
            for (k in 1:length(dag[1, ])) {
                if (dag[k, node.children] != 0) {
                  mb.parent.children[k] <- names(dag[, node.children])[k]
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
            for (l in 1:length(mb.node)) {
                if (mb.node[c(l)] == node) {
                  mb.node.wo <- mb.node[-c(l)]
                }
            }
        }
        if (is.null(mb.node.wo)) {
            mb.node.wo <- mb.node
        }

        mb.size[i] <- length(mb.node.wo)
    }

    mb.average <- mean(mb.size)

    ## average Neighborhood

    nh.size <- vector(mode = "numeric", length = length(node.names))
    parent.size <- vector(mode = "numeric", length = length(node.names))
    children.size <- vector(mode = "numeric", length = length(node.names))

    for (i in 1:length(node.names)) {
        nh.size[i] <- sum(dag[i, ]) + sum(dag[, i])
        parent.size[i] <- sum(dag[i, ])
        children.size[i] <- sum(dag[, i])
    }
    nh.average <- mean(nh.size)
    parent.average <- mean(parent.size)
    children.average <- mean(children.size)

    ## output
    out[["n.nodes"]] <- n.nodes
    out[["n.arcs"]] <- n.arcs
    out[["mb.average"]] <- mb.average
    out[["nh.average"]] <- nh.average
    out[["parent.average"]] <- parent.average
    out[["children.average"]] <- children.average

    return(out)
}

#' Simulate a DAG with with arbitrary arcs density
#'
#' Provided with node names, returns an \code{abnDAG}.
#' Arc density refers to the chance of a node being connected to the node before it.
#'
#'
#' @param node.name a vector of character giving the names of the nodes. It gives the size of the simulated DAG.
#' @param data.dists named list corresponding to the \code{node.name} specifying the distribution for each node. If not provided arbitrary distributions are assigned to the nodes.
#' @param edge.density number in \code{[0,1]} specifying the edge probability in the dag.
#' @param verbose print more information on the run.
#'
#' @return object of class \code{abnDag} consisting of a named matrix, a named list giving the distribution for each node and an empty element for the data.
#'
#' @details
#' This function generates DAGs by sampling triangular matrices and reorder columns and rows randomly.
#' The network density (\code{edge.density}) is used column-wise as binomial sampling probability.
#' Then the matrix is named using the user-provided names.
#' @importFrom stats rbinom
#' @export
#' @examples
#' simdag <- simulateDag(node.name = c("a", "b", "c", "d"),
#'                       edge.density = 0.5,
#'                       data.dists = list(a = "gaussian",
#'                                         b = "binomial",
#'                                         c = "poisson",
#'                                         d = "multinomial"))
#'
#' ## Example using Ozon entries:
#' dist <- list(Ozone="gaussian",   Solar.R="gaussian",  Wind="gaussian",
#'              Temp="gaussian",    Month="gaussian",    Day="gaussian")
#' out <- simulateDag(node.name = names(dist), data.dists = dist, edge.density = 0.8)
#' plot(out)
#' @keywords utilities
simulateDag <-
  function(node.name,
           data.dists = NULL,
           edge.density = 0.5,
           verbose = FALSE) {
    ## test
    if (is.null(node.name)) {
      stop("'node.name' is missing. Hence, I don't know the number of nodes in the DAG.")
    } else {
      nvars <- length(node.name)
    }

    if (edge.density > 1 | edge.density < 0) {
      stop("'edge.density' should be a real number in [0,1]")
    }

    if (is.null(data.dists)) {
      if(verbose){message("No 'data.dists' provided. Using arbitrary distributions.")}

      data.dists <-
        sample(
          x = c("gaussian", "binomial", "poisson", "multinomial"),
          size = nvars,
          replace = TRUE
        )
      names(data.dists) <- node.name
    } else {
      data.dists <-
        validate_dists(data.dists = data.dists, returnDists = TRUE)
    }

    ## simulate dag
    dag <- matrix(data = 0,
                  nrow = nvars,
                  ncol = nvars)

    for (j in 2:(nvars)) {
      dag[c(j:nvars), (j - 1)] <-
        rbinom(n = (nvars - j + 1),
               size = 1,
               prob = edge.density)
    }
    order.m <- sample(x = nvars,
                      size = nvars,
                      replace = FALSE)

    ## changing order
    dag <- dag[order.m, order.m]  #order.m

    ## naming
    colnames(dag) <- rownames(dag) <- node.name

    ## structure
    out <- createAbnDag(dag = dag, data.dists = data.dists)

    return(out)
  }

#' Computes skewness of a distribution
#'
#' @param x a numeric vector
#'
#' @return integer
#' @export
#' @keywords utilities
skewness <- function(x) {
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
}

#' Construct the essential graph
#'
#' Constructs different versions of the essential graph from a given DAG.
#' External function that computes essential graph of a dag Minimal PDAG:
#' The only directed edges are those who participate in v-structure Completed PDAG:
#' very directed edge corresponds to a compelled edge, and every undirected
#' edge corresponds to a reversible edge
#'
#' @usage essentialGraph(dag, node.names = NULL, PDAG = "minimal")
#'
#' @param dag a matrix or a formula statement (see \sQuote{Details} for format)
#' defining the network structure, a directed acyclic graph (DAG).
#' @param node.names a vector of names if the DAG is given via formula, see \sQuote{Details}.
#' @param PDAG a character value that can be: minimal or complete, see \sQuote{Details}.
#'
#' @details
#' This function returns an essential graph from a DAG,
#' aka acyclic partially directed graph (PDAG).
#' This can be useful if the learning procedure is defined up to a Markov class
#' of equivalence.
#' A minimal PDAG is defined as only directed edges are those who participate
#' in v-structure. Whereas the completed PDAG: every directed edge corresponds
#' to a compelled edge, and every undirected edge corresponds to a reversible edge.
#'
#' The \code{dag} can be provided using a formula statement (similar to glm).
#' A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}.
#' The formula statement have to start with \code{~}.
#' In this example, node1 has two parents (parent1 and parent2).
#' node2 and node3 have the same parent3.
#' The parents names have to exactly match those given in \code{node.names}.
#' \code{:} is the separator between either children or parents,
#' \code{|} separates children (left side) and parents (right side),
#' \code{+} separates terms, \code{.} replaces all the variables in \code{node.names}.
#'
#' @return A matrix giving the PDAG.
#'
#' @references
#' West, D. B. (2001). Introduction to Graph Theory. Vol. 2. Upper Saddle River: Prentice Hall.
#' Chickering, D. M. (2013) A Transformational Characterization of Equivalent Bayesian Network Structures, arXiv:1302.4938.
#'
#' @export
#' @keywords utilities
#' @examples
#' dag <- matrix(c(0,0,0, 1,0,0, 1,1,0), nrow = 3, ncol = 3)
#' dist <- list(a="gaussian", b="gaussian", c="gaussian")
#' colnames(dag) <- rownames(dag) <- names(dist)
#'
#' essentialGraph(dag)
essentialGraph <- function(dag, node.names = NULL, PDAG = "minimal") {

    ## dag transformation
    if (is.matrix(dag)) {
        # check.valid.dag(dag=dag,data.df=data.df,is.ban.matrix=FALSE,group.var=group.var) unit matrix
        dag[dag != 0] <- 1
        node.names <- colnames(dag)

    } else {
        if (grepl("~", as.character(dag)[1], fixed = TRUE)) {
            dag <- formula_abn(f = dag, name = node.names)
        } else {
            stop("Dag specification must either be a matrix or a formula expression")
        }
    }

    ## compute essential graph moral graph

    dim.dag <- dim(dag)[1]
    moral <- matrix(data = 0, nrow = dim.dag, ncol = dim.dag)

    for (i in 1:dim.dag) {
        for (j in 1:dim.dag) {
            if (dag[i, j] == 1) {
                moral[i, j] <- 1
                moral[j, i] <- 1
            }
        }
    }

    colnames(moral) <- rownames(moral) <- node.names

    ## essential arcs
    if (PDAG == "minimal") {
        for (i in 1:dim.dag) {
            if (sum(dag[i, ]) >= 2) {
                for (j in 1:dim.dag) {
                  if (dag[i, j] == 1) {
                    moral[j, i] <- 0
                  }
                }
            }
        }


        colnames(moral) <- rownames(moral) <- node.names
        return(moral)
    }

    if (PDAG == "completed") {
        for (i in 1:dim.dag) {
            if (sum(dag[i, ]) >= 2) {
                for (j in 1:dim.dag) {
                  if (dag[i, j] == 1) {
                    moral[j, i] <- 0
                  }
                  if (dag[j, i] == 1) {
                    moral[i, j] <- 0
                  }
                }
            }
        }

        colnames(moral) <- rownames(moral) <- node.names
        return(moral)
    }
}

#' Compute the score's contribution in a network of each observation.
#'
#' This function computes the score's contribution of each observation to the total network score.
#'
#' @usage scoreContribution(object = NULL,
#'                          dag = NULL, data.df = NULL, data.dists = NULL,
#'                          verbose = FALSE)
#'
#' @param object an object of class '\code{abnLearned}'
#' produced by \code{\link{mostProbable}},
#' \code{\link{searchHeuristic}} or \code{\link{searchHillClimber}}.
#' @param dag a matrix or a formula statement (see details) defining the network structure,
#' a directed acyclic graph (DAG), see details for format.
#' Note that colnames and rownames must be set.
#' @param data.df a data frame containing the data used for learning the network,
#' binary variables must be declared as factors and
#' no missing values all allowed in any variable.
#' @param data.dists a named list giving the distribution for each node
#' in the network, see details.
#' @param verbose if \code{TRUE} then provides some additional output.
#'
#' @details
#' This function computes the score contribution of each observation
#' to the total network score.
#' This function is available only in the `mle` settings.
#' To do so one uses the \code{\link{glm}} and \code{\link{predict}} functions.
#' This function is an attempt to perform diagnostic for an ABN analysis.
#'
#' @return A named list that contains the scores contributions:
#' maximum likelihood, aic, bic, mdl and diagonal values of the hat matrix.
#' @export
#' @keywords utilities
#' @importFrom stats glm predict.glm dbinom dnorm sigma dpois hatvalues
#' @examples
#' \dontrun{
#' ## Use a subset of a built-in simulated data set
#' mydat <- ex1.dag.data[,c("b1","g1","p1")]
#'
#' ## setup distribution list for each node
#' mydists <- list(b1="binomial", g1="gaussian", p1="poisson")
#'
#' ## now build cache
#' mycache <- buildScoreCache(data.df = mydat, data.dists = mydists, max.parents = 2, method = "mle")
#'
#' ## Find the globally best DAG
#' mp.dag <- mostProbable(score.cache=mycache, score="bic", verbose = FALSE)
#'
#' out <- scoreContribution(object = mp.dag)
#'
#' ## Observations contribution per network node
#' boxplot(out$bic)
#' }
scoreContribution <- function(object = NULL,
                              dag = NULL,
                              data.df = NULL,
                              data.dists = NULL,
                              verbose = FALSE){

    ## method abnCache
    if (!is.null(object)){
        if (inherits(object, "abnLearned")) {
            dag <- object$dag
            data.df <- object$score.cache$data.df
            data.dists <- object$score.cache$data.dists
        }}

    ## transform factor into 0/1

    node.ordering <- names(data.dists)
    nobs <- dim(data.df)[1]

    ll <- matrix(data = 0,nrow = dim(data.df)[1],ncol = dim(data.df)[2])
    nb.param <- matrix(data = 0,nrow = dim(data.df)[1],ncol = dim(data.df)[2])
    nb.parents <- matrix(data = 0,nrow = dim(data.df)[1],ncol = dim(data.df)[2])
    colnames(ll) <- colnames(nb.param) <- colnames(nb.parents) <- node.ordering
    nb.parents <- rowSums(dag)/nobs
    for (node in node.ordering) {

        if(as.character(data.dists[node])=="binomial"){
            Y <- as.factor(data.matrix(data.df[, node]))
        } else {
            Y <- data.matrix(data.df[, node])
        }
        X <- data.matrix(cbind(rep(1,nobs),data.df[, as.logical(dag[node,])]))
        if(as.character(data.dists[node])=="gaussian"){
            nb.param[,node] <- (dim(X)[2] + 1)/nobs
        }else{
            nb.param[,node] <- dim(X)[2]/nobs
        }

        x <- glm(formula = Y ~ -1 + X,family = as.character(data.dists[node]))
        yhat <- predict.glm(x,newdata = x$data, type = "response")

        switch (as.character(data.dists[node]),
                "binomial" = {
                    ll[,node] <- dbinom(as.numeric(Y)-1, size = 1L, prob = yhat, log = TRUE)
                },
                "gaussian" = {
                    ll[,node] <- dnorm(Y, mean = yhat, sd = sigma(x),log = TRUE)
                },
                "poisson" = {
                    ll[,node] <- dpois(Y, lambda = yhat, log = TRUE)
                }
        )
        hv <- hatvalues(x)

    }

    aic <- - 2*ll + 2*nb.param
    bic <- - 2*ll + nb.param*(log(nobs))
    mdl <- bic + (1/nobs + nb.parents) * log(nobs)

    out <- list("mlik" = ll, "aic" = aic, "bic" = bic, "mdl" = mdl, "hatvalues" = hv)

    return(out)

}
