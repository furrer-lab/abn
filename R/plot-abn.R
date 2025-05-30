#' Plot an ABN graphic
#'
#' Plot an ABN DAG using formula statement or a matrix in using Rgraphviz through the graphAM class.
#'
#' @usage plotAbn(dag, data.dists=NULL, markov.blanket.node=NULL, fitted.values=NULL,
#'                digits=2, edge.strength=NULL, edge.strength.lwd=5, edge.direction="pc",
#'                edge.color="black", edge.linetype="solid", edge.arrowsize=0.6,
#'                edge.fontsize=node.fontsize, node.fontsize=12,
#'                node.fillcolor=c("lightblue", "brown3", "chartreuse3"),
#'                node.fillcolor.list=NULL,
#'                node.shape=c("circle", "box", "ellipse", "diamond"),
#'                plot=TRUE,
#'                data.df=NULL, ...)
#'
#' @param dag a matrix or a formula statement (see details for format) defining
#' the network structure, a Directed Acyclic Graph (DAG).
#' Note that rownames must be set or given in \code{data.dists}.
#' @param data.dists a named list giving the distribution for each node in the network, see details.
#' @param markov.blanket.node name of variables to display its Markov blanket.
#' @param fitted.values modes or coefficents outputted from \code{\link{fitAbn}}.
#' @param digits number of digits to display the \code{fitted.values}.
#' @param edge.strength a named matrix containing evaluations of edge strength
#' which will change the arcs width (could be Mutual information, p-values,
#' number of bootstrap retrieve samples or the outcome of the \code{\link{linkStrength}}).
#' @param edge.strength.lwd maximum line width for \code{edge.strength}.
#' @param edge.direction character giving the direction in which arcs should
#' be plotted, \code{pc} (parent to child) or \code{cp} (child to parent) or \code{undirected}.
#' @param edge.color the colour of the edge.
#' @param edge.linetype the linetype of the edge. Defaults to \code{"solid"}.
#' Valid values are the same as for the R's base graphic parameter \code{lty}.
#' @param edge.arrowsize the thickness of the arrows. Not relevant if \code{arc.strength} is provided.
#' @param edge.fontsize the font size of the arcs fitted values.
#' @param node.fontsize the font size of the nodes names.
#' @param node.fillcolor the colour of the node. Second and third element is
#' used for the Markov blanket and node of the Markov blanket.
#' @param node.fillcolor.list the list of node that should be coloured.
#' @param node.shape the shape of the nodes according the Gaussian, binomial, Poisson and multinomial distributions.
#' @param plot logical variable, if set to \code{TRUE} then the graph is plotted.
#' @param data.df NULL or a data frame containing the data for the nodes in the network. Only needed if \code{dag} is a formula statement. If dag is an object of class \code{abnFit}, then \code{data.df} is used from there.
#' @param ... arguments passed to the plotting function.
#'
#' @details
#' By default binomial nodes are squares, multinomial nodes are empty, Gaussian nodes are circles and poison nodes are ellipses.
#'
#' The \code{dag} can be provided using a formula statement (similar to glm). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}.
#'
#' The construction is based on the \pkg{graph} package. Properties of the graph can be changed after the construction, see \sQuote{Examples}.
#'
#' @return A rendered graph, if \code{plot=TRUE}. The \code{graphAM} object is returned invisibly.
#'
#' @seealso \code{\link[graph]{graphAM-class}}, \code{\link[graph]{edgeRenderInfo}}
#' @keywords internal
#' @export
#' @importFrom Rgraphviz renderGraph layoutGraph
#' @importFrom graph edgeRenderInfo edgeRenderInfo<-
#' @importFrom methods new
#' @examples
#' # Define distribution list
#' dist <- list(a = "gaussian",
#'              b = "gaussian",
#'              c = "gaussian",
#'              d = "gaussian",
#'              e = "binomial",
#'              f = "binomial")
#'
#' # Define a matrix formulation
#' edge_strength <- matrix(c(0, 0.5, 0.5, 0.7, 0.1, 0,
#'                           0, 0, 0.3, 0.1, 0, 0.8,
#'                           0, 0, 0, 0.35, 0.66, 0,
#'                           0, 0, 0, 0, 0.9, 0,
#'                           0, 0, 0, 0, 0, 0.8,
#'                           0, 0, 0, 0, 0, 0),
#'                         nrow = 6L,
#'                         ncol = 6L,
#'                         byrow = TRUE)
#'
#' ## Naming of the matrix
#' colnames(edge_strength) <- rownames(edge_strength) <- names(dist)
#'
#' ## Random Data
#' df <- data.frame(a = rnorm(100),
#'                 b = rnorm(100),
#'                 c = rnorm(100),
#'                 d = rnorm(100),
#'                 e = rbinom(100, 1, 0.5),
#'                 f = rbinom(100, 1, 0.5))
#'
#' ## Plot form a matrix
#' plotAbn(dag = edge_strength,
#'         data.dists = dist)
#'
#' ## Edge strength
#' plotAbn(dag = ~ a | b:c:d:e + b | c:d:f + c | d:e + d | e + e | f,
#'         data.dists = dist,
#'         edge.strength = edge_strength,
#'         data.df = df)
#'
#' ## Plot from a formula for a different DAG!
#' plotAbn(dag = ~ a | b:c:e + b | c:d:f + e | f,
#'         data.dists = dist,
#'         data.df = df)
#'
#' ## Markov blanket
#' plotAbn(dag = ~ a | b:c:e + b | c:d:f + e | f,
#'         data.dists = dist,
#'         markov.blanket.node = "e",
#'         data.df = df)
#'
#' ## Change col for all edges
#' tmp <- plotAbn(dag = ~ a | b:c:e + b | c:d:f + e | f,
#'                data.dists = dist,
#'                plot = FALSE,
#'                data.df = df)
#' graph::edgeRenderInfo(tmp) <- list(col = "blue")
#' Rgraphviz::renderGraph(tmp)
#'
#' ## Change lty for individual ones. Named vector is necessary
#' tmp <- plotAbn(dag = ~ a | b:c:e + b | c:d:f + e | f,
#'                data.dists = dist,
#'                plot = FALSE,
#'                data.df = df)
#' edgelty <- rep(c("solid", "dotted"), c(6, 1))
#' names(edgelty) <- names(graph::edgeRenderInfo(tmp, "col"))
#' graph::edgeRenderInfo(tmp) <- list(lty = edgelty)
#' Rgraphviz::renderGraph(tmp)
plotAbn <- function(dag,
                    data.dists = NULL,
                    markov.blanket.node = NULL,
                    fitted.values = NULL,
                    digits = 2,
                    edge.strength = NULL,
                    edge.strength.lwd = 5,
                    edge.direction = "pc",
                    edge.color = "black",
                    edge.linetype = "solid",
                    edge.arrowsize = 0.6,
                    edge.fontsize = node.fontsize,
                    node.fontsize = 12,
                    node.fillcolor = c("lightblue", "brown3", "chartreuse3"),
                    node.fillcolor.list = NULL,
                    node.shape = c("circle", "box", "ellipse", "diamond"),
                    plot = TRUE,
                    data.df = NULL,
                    ...) {
  # Actually, the plot argument is wrong! i do not need the adjacency structure only. I need all but the plotting. i.e., all but the rendering of the graph.


  # The following is not relevant. The nodes are calculated via mb. They are not colored.
  #    if(!is.null(markov.blanket.node) & ("multinomial" %in% (data.dists))) warning("Multinomial nodes are excluded from Markov blanket computation.")

  ## for compatibility purpose
  if (inherits(x = dag, what = "abnLearned")) {
    data.dists <- dag$score.cache$data.dists
    dag <- dag$dag
  }
  if (inherits(x = dag, what = "abnFit")) {
    data.dists <- dag$abnDag$data.dists
    dag <- dag$abnDag$dag
    data.df <- dag$abnDag$data.df
  }
  if (is.null(data.dists)) stop("'data.dist' need to be provided.")
  name <- names(data.dists)


  ## dag transformation
  if (!is.null(dag)) {
    if (is.matrix(dag)) {
      ## run a series of checks on the DAG passed
      dag <- abs(dag)
      ## consistency checks
      diag(dag) <- 0
      dag[dag > 0] <- 1

      ## naming
      if (is.null(rownames(dag))) {
        colnames(dag) <- name
        rownames(dag) <- name
      }
      dag <- check.valid.dag(dag = dag, is.ban.matrix = FALSE, group.var = NULL)
    } else {
      if (grepl("~", as.character(dag)[1], fixed = TRUE)) {
        dag <- formula_abn(f = dag, name = name)
        ## run a series of checks on the DAG passed
        dag <- check.valid.dag(dag = dag, data.df = data.df, is.ban.matrix = FALSE, group.var = NULL)
      }
    }
  } else {
    stop("'dag' specification must either be a matrix or a formula expression")
  }

  # contains Rgraphviz
  if (edge.direction  == "undirected") {
    dag <- dag + t(dag)
    dag[dag !=  0] <- 1     # this should not be necessary!
  }


  ## create an object graph
  am.graph <- new(Class = "graphAM", adjMat = dag,
                  edgemode = ifelse(edge.direction == "undirected", "undirected", "directed"))

  ## ========= SHAPE =========
  ## Shape: plot differential depending on the distribution
  node.shape <- rep(node.shape, 4)
  shape <- rep(node.shape[1], length(data.dists))
  shape[data.dists == "binomial"] <- node.shape[2]
  shape[data.dists == "poisson"] <- node.shape[3]
  shape[data.dists == "multinomial"] <- node.shape[4]
  names(shape) <- names(data.dists)

  ## ================= NODE FILLED COLOR =================
  ## fill with default value, change if MB or fillcolor.list is requested
  fillcolor <- rep(node.fillcolor[1], length(data.dists))
  names(fillcolor) <- names(data.dists)

  ## =============== MARKOV BLANKET ===============
  ## Markov Blanket: plot the MB of a given node
  if (!is.null(markov.blanket.node)) {
    markov.blanket <- mb(dag, node =markov.blanket.node, data.dists = data.dists, data.df = data.df)
    fillcolor[names(data.dists) %in%  markov.blanket]  <- node.fillcolor[3]
    fillcolor[names(data.dists) %in%  markov.blanket.node]  <- node.fillcolor[2]

  } else if (!is.null(node.fillcolor.list)) {
    fillcolor[names(data.dists) %in%  node.fillcolor.list] <- node.fillcolor[2]
  }

  names.edges <- names(Rgraphviz::buildEdgeList(am.graph))
  ## =============== Fitted values ===============
  ## Plot the fitted values in abn as edges label
  if (!is.null(fitted.values)) {
    space <- "   "
    edge.label <- c()
    for (i in 1:length(fitted.values)) {
      if ((length(fitted.values[[i]]) > 1) & !(data.dists[names(fitted.values)[i]] %in% c("gaussian", "multinomial"))) {
        # If data distribution is binomial or poisson:
        # get parent names
        parent.names <- names(which(dag[i, ] == 1))
        # iterate over all parents
        p <- 1
        while (p <= length(parent.names)) {
          # If no parents, skip this.
          # if parent is multinomial
          if (data.dists[parent.names[p]] == "multinomial") {
            # get the number of categories of the parent
            ncat <- sum(stringi::stri_detect(str = colnames(fitted.values[[parent.names[p]]]), fixed = "intercept")) + 1
            # iterate over all fitted.values[[i]] corresponding to the parent categories
            # extract the values of the respective categories
            fitted.values_cat <- fitted.values[[i]][stringi::stri_detect(str = colnames(fitted.values[[i]]), regex = paste0("^", parent.names[p]))]
            # extract the labels of the respective categories
            fitted.values_cat_labels <- colnames(fitted.values[[i]])[stringi::stri_detect(str = colnames(fitted.values[[i]]), regex = paste0("^", parent.names[p]))]
            # paste the labels and the values together
            parent_catlabels <- c()
            for (k in 1:length(fitted.values_cat)) {
              parent_catlabels <- c(parent_catlabels, paste(fitted.values_cat_labels[k], signif(fitted.values_cat[k], digits = digits), sep = ": "))
            }
            parent_catlabels <- paste(stringi::stri_flatten(str = parent_catlabels, collapse = "\n"))
            # add it to the edge.label
            edge.label <- c(edge.label, paste(parent_catlabels))
          } else if (data.dists[parent.names[p]] != "multinomial") {
            # if parent is not multinomial
            # Workaround because the output of "bayes" and "mle" differ a little
            if (!is.null(names(fitted.values[[i]]))) {
              # "bayes" output has names because it's an array
              edge.label <- c(edge.label, paste(space, signif(fitted.values[[i]][stringi::stri_detect(str = names(fitted.values[[i]]), regex = paste0(parent.names[p], "$"))], digits = digits)))
            } else {
              # "mle" output has colnames because it's a matrix
              edge.label <- c(edge.label, paste(space, signif(fitted.values[[i]][stringi::stri_detect(str = colnames(fitted.values[[i]]), regex = paste0(parent.names[p], "$"))], digits = digits)))
            }
          }
          p <- p+1
        }
      } else if ((length(fitted.values[[i]]) > 1) & (data.dists[names(fitted.values)[i]] == "gaussian")) {
        # get parent names
        parent.names <- names(which(dag[i, ] == 1))
        # iterate over all parents
        p <- 1
        while (p <= length(parent.names)) {
          # If no parents, skip this.
          # if parent is multinomial
          if (data.dists[parent.names[p]] == "multinomial") {
            # get the number of categories of the parent
            ncat <- sum(stringi::stri_detect(str = colnames(fitted.values[[parent.names[p]]]), fixed = "intercept")) + 1
            # iterate over all fitted.values[[i]] corresponding to the parent categories
            # extract the values of the respective categories
            fitted.values_cat <- fitted.values[[i]][stringi::stri_detect(str = colnames(fitted.values[[i]]), regex = paste0("^", parent.names[p]))]
            # extract the labels of the respective categories
            fitted.values_cat_labels <- colnames(fitted.values[[i]])[stringi::stri_detect(str = colnames(fitted.values[[i]]), regex = paste0("^", parent.names[p]))]
            # paste the labels and the values together
            parent_catlabels <- c()
            for (k in 1:length(fitted.values_cat)) {
              parent_catlabels <- c(parent_catlabels, paste(fitted.values_cat_labels[k], signif(fitted.values_cat[k], digits = digits), sep = ": "))
            }
            parent_catlabels <- paste(stringi::stri_flatten(str = parent_catlabels, collapse = "\n"))
            # add it to the edge.label
            edge.label <- c(edge.label, paste(parent_catlabels))
          } else if (data.dists[parent.names[p]] != "multinomial") {
            # if parent is not multinomial do as always
            # Workaround because the output of "bayes" and "mle" differ a little
            if (!is.null(names(fitted.values[[i]]))) {
              # "bayes" output has names because it's an array
              edge.label <- c(edge.label, paste(space, signif(fitted.values[[i]][stringi::stri_detect(str = names(fitted.values[[i]]), regex = paste0(parent.names[p], "$"))], digits = digits)))
            } else {
              # "mle" output has colnames because it's a matrix
              edge.label <- c(edge.label, paste(space, signif(fitted.values[[i]][stringi::stri_detect(str = colnames(fitted.values[[i]]), regex = paste0(parent.names[p], "$"))], digits = digits)))
            }
          }
          p <- p + 1
        }
      } else if ((length(fitted.values[[i]]) > 1) & (data.dists[names(fitted.values)[i]] == "multinomial")) {
        # number of categories of node i
        ncat <- sum(stringi::stri_detect(str = colnames(fitted.values[[i]]), fixed = "intercept")) + 1
        for (j in 1:(ncat - 1) - 1) {
          # iterate over all parents of node i
          catlabels <- c()
          for (k in seq(ncat + j, length(fitted.values[[i]]), ncat - 1)){
            # iterate over all categories of node i with respect to the parent j
            catlabels <- c(catlabels, paste(colnames(fitted.values[[i]])[k], signif(fitted.values[[i]][k], digits = digits), sep = ": "))
          }
          edge.label <- c(edge.label, paste(stringi::stri_flatten(str = catlabels, collapse = "\n")))
        }
      }
    }
  } else {
    edge.label <- rep(" ", length(names.edges))
  }
  names(edge.label) <- names.edges
  ## =================== Arc Strength ===================
  ## Arc strength: plot the AS of the dag arcs
  #    if (is.matrix(edge.strength) & (edge.direction != "undirected")) {
  if (is.matrix(edge.strength)) {
    if (any(edge.strength < 0)) stop("'edge.strength' should be positive")
    if (any(edge.strength[dag == 0] > 0)) stop("'edge.strength' does not match dag")
    min.as <- min(edge.strength[edge.strength > 0])
    max.as <- max(edge.strength[edge.strength > 0])

    edge.strength.norm <- (edge.strength - min.as) / (max.as - min.as)
    edge.strength.norm[edge.strength.norm < 0] <- 0
    edge.lwd <- list()
    for (i in 1:length(dag[1, ])) {
      for (j in 1:length(dag[1, ])) {
        if (dag[i, j] == 1) {
          edge.lwd <- cbind(edge.lwd, round(edge.strength.lwd * edge.strength.norm[i, j]) + 1)
        }
      }
    }
  } else {
    edge.lwd <- rep(1, length(names.edges))
  }
  class(edge.lwd) <- "character"
  names(edge.lwd) <- names.edges

  ## ====== Plot ======
  attrs <- list(graph = list(rankdir = "BT"),
                node = list(fontsize = node.fontsize, fixedsize = FALSE),
                edge = list(arrowsize = edge.arrowsize, color = edge.color, lty = edge.linetype, fontsize = edge.fontsize))
  nodeAttrs <- list(fillcolor = fillcolor, shape = shape)
  edgeAttrs <- list(label = edge.label, lwd = edge.lwd)
  #     print(edgeAttrs)
  #    if (all(shape %in% c("circle","box","ellipse")))  {
  am.graph <- layoutGraph(am.graph, attrs = attrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)

  if (edge.direction == "pc")  {     # specify appropriate direction!
    edgeRenderInfo(am.graph) <- list(arrowtail = "open")
    edgeRenderInfo(am.graph) <- list(arrowhead = "none")
    #      edgeRenderInfo(am.graph) <- list(direction=NULL)# MESSES up!!! not needed.
  }
  edgeRenderInfo(am.graph) <- list(lwd = edge.lwd)

  #  if (plot) renderGraph(am.graph, attrs=attrs, nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs)
  if (plot) renderGraph(am.graph, ...)

  #   } else {
  #        am.graph <- layoutGraph(am.graph, attrs=attrs, nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs, ...)
  # the following does not work in R
  #        edgeRenderInfo(am.graph)[["direction"]] <- "back"
  # hence
  #        warning("edge.direction='pc' is not working with diamond shapes.")
  #        edgeRenderInfo(am.graph) <- list(lwd=edge.lwd)
  #        if (plot) renderGraph(am.graph,attrs=attrs, nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs)
  #    }

  invisible(am.graph)
}  #EOF
