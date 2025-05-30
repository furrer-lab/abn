% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-abn.R
\name{plotAbn}
\alias{plotAbn}
\title{Plot an ABN graphic}
\usage{
plotAbn(dag, data.dists=NULL, markov.blanket.node=NULL, fitted.values=NULL,
               digits=2, edge.strength=NULL, edge.strength.lwd=5, edge.direction="pc",
               edge.color="black", edge.linetype="solid", edge.arrowsize=0.6,
               edge.fontsize=node.fontsize, node.fontsize=12,
               node.fillcolor=c("lightblue", "brown3", "chartreuse3"),
               node.fillcolor.list=NULL,
               node.shape=c("circle", "box", "ellipse", "diamond"),
               plot=TRUE,
               data.df=NULL, ...)
}
\arguments{
\item{dag}{a matrix or a formula statement (see details for format) defining
the network structure, a Directed Acyclic Graph (DAG).
Note that rownames must be set or given in \code{data.dists}.}

\item{data.dists}{a named list giving the distribution for each node in the network, see details.}

\item{markov.blanket.node}{name of variables to display its Markov blanket.}

\item{fitted.values}{modes or coefficents outputted from \code{\link{fitAbn}}.}

\item{digits}{number of digits to display the \code{fitted.values}.}

\item{edge.strength}{a named matrix containing evaluations of edge strength
which will change the arcs width (could be Mutual information, p-values,
number of bootstrap retrieve samples or the outcome of the \code{\link{linkStrength}}).}

\item{edge.strength.lwd}{maximum line width for \code{edge.strength}.}

\item{edge.direction}{character giving the direction in which arcs should
be plotted, \code{pc} (parent to child) or \code{cp} (child to parent) or \code{undirected}.}

\item{edge.color}{the colour of the edge.}

\item{edge.linetype}{the linetype of the edge. Defaults to \code{"solid"}.
Valid values are the same as for the R's base graphic parameter \code{lty}.}

\item{edge.arrowsize}{the thickness of the arrows. Not relevant if \code{arc.strength} is provided.}

\item{edge.fontsize}{the font size of the arcs fitted values.}

\item{node.fontsize}{the font size of the nodes names.}

\item{node.fillcolor}{the colour of the node. Second and third element is
used for the Markov blanket and node of the Markov blanket.}

\item{node.fillcolor.list}{the list of node that should be coloured.}

\item{node.shape}{the shape of the nodes according the Gaussian, binomial, Poisson and multinomial distributions.}

\item{plot}{logical variable, if set to \code{TRUE} then the graph is plotted.}

\item{data.df}{NULL or a data frame containing the data for the nodes in the network. Only needed if \code{dag} is a formula statement. If dag is an object of class \code{abnFit}, then \code{data.df} is used from there.}

\item{...}{arguments passed to the plotting function.}
}
\value{
A rendered graph, if \code{plot=TRUE}. The \code{graphAM} object is returned invisibly.
}
\description{
Plot an ABN DAG using formula statement or a matrix in using Rgraphviz through the graphAM class.
}
\details{
By default binomial nodes are squares, multinomial nodes are empty, Gaussian nodes are circles and poison nodes are ellipses.

The \code{dag} can be provided using a formula statement (similar to glm). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}.

The construction is based on the \pkg{graph} package. Properties of the graph can be changed after the construction, see \sQuote{Examples}.
}
\examples{
# Define distribution list
dist <- list(a = "gaussian",
             b = "gaussian",
             c = "gaussian",
             d = "gaussian",
             e = "binomial",
             f = "binomial")

# Define a matrix formulation
edge_strength <- matrix(c(0, 0.5, 0.5, 0.7, 0.1, 0,
                          0, 0, 0.3, 0.1, 0, 0.8,
                          0, 0, 0, 0.35, 0.66, 0,
                          0, 0, 0, 0, 0.9, 0,
                          0, 0, 0, 0, 0, 0.8,
                          0, 0, 0, 0, 0, 0),
                        nrow = 6L,
                        ncol = 6L,
                        byrow = TRUE)

## Naming of the matrix
colnames(edge_strength) <- rownames(edge_strength) <- names(dist)

## Random Data
df <- data.frame(a = rnorm(100),
                b = rnorm(100),
                c = rnorm(100),
                d = rnorm(100),
                e = rbinom(100, 1, 0.5),
                f = rbinom(100, 1, 0.5))

## Plot form a matrix
plotAbn(dag = edge_strength,
        data.dists = dist)

## Edge strength
plotAbn(dag = ~ a | b:c:d:e + b | c:d:f + c | d:e + d | e + e | f,
        data.dists = dist,
        edge.strength = edge_strength,
        data.df = df)

## Plot from a formula for a different DAG!
plotAbn(dag = ~ a | b:c:e + b | c:d:f + e | f,
        data.dists = dist,
        data.df = df)

## Markov blanket
plotAbn(dag = ~ a | b:c:e + b | c:d:f + e | f,
        data.dists = dist,
        markov.blanket.node = "e",
        data.df = df)

## Change col for all edges
tmp <- plotAbn(dag = ~ a | b:c:e + b | c:d:f + e | f,
               data.dists = dist,
               plot = FALSE,
               data.df = df)
graph::edgeRenderInfo(tmp) <- list(col = "blue")
Rgraphviz::renderGraph(tmp)

## Change lty for individual ones. Named vector is necessary
tmp <- plotAbn(dag = ~ a | b:c:e + b | c:d:f + e | f,
               data.dists = dist,
               plot = FALSE,
               data.df = df)
edgelty <- rep(c("solid", "dotted"), c(6, 1))
names(edgelty) <- names(graph::edgeRenderInfo(tmp, "col"))
graph::edgeRenderInfo(tmp) <- list(lty = edgelty)
Rgraphviz::renderGraph(tmp)
}
\seealso{
\code{\link[graph]{graphAM-class}}, \code{\link[graph]{edgeRenderInfo}}
}
\keyword{internal}
