% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abn-toolbox.R
\name{simulateDag}
\alias{simulateDag}
\title{Simulate a DAG with with arbitrary arcs density}
\usage{
simulateDag(node.name, data.dists = NULL, edge.density = 0.5, verbose = FALSE)
}
\arguments{
\item{node.name}{a vector of character giving the names of the nodes. It gives the size of the simulated DAG.}

\item{data.dists}{named list corresponding to the \code{node.name} specifying the distribution for each node. If not provided arbitrary distributions are assigned to the nodes.}

\item{edge.density}{number in \code{[0,1]} specifying the edge probability in the dag.}

\item{verbose}{print more information on the run.}
}
\value{
object of class \code{abnDag} consisting of a named matrix, a named list giving the distribution for each node and an empty element for the data.
}
\description{
Provided with node names, returns an \code{abnDAG}.
Arc density refers to the chance of a node being connected to the node before it.
}
\details{
This function generates DAGs by sampling triangular matrices and reorder columns and rows randomly.
The network density (\code{edge.density}) is used column-wise as binomial sampling probability.
Then the matrix is named using the user-provided names.
}
\examples{
simdag <- simulateDag(node.name = c("a", "b", "c", "d"),
                      edge.density = 0.5,
                      data.dists = list(a = "gaussian",
                                        b = "binomial",
                                        c = "poisson",
                                        d = "multinomial"))

## Example using Ozon entries:
dist <- list(Ozone="gaussian",   Solar.R="gaussian",  Wind="gaussian",
             Temp="gaussian",    Month="gaussian",    Day="gaussian")
out <- simulateDag(node.name = names(dist), data.dists = dist, edge.density = 0.8)
plot(out)
}
\keyword{utilities}
