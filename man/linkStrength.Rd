% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/link-strength.R
\name{linkStrength}
\alias{linkStrength}
\title{Returns the strengths of the edge connections in a Bayesian Network learned from observational data}
\usage{
linkStrength(dag,
                    data.df = NULL,
                    data.dists = NULL,
                    method = c("mi.raw",
                               "mi.raw.pc",
                               "mi.corr",
                               "ls",
                               "ls.pc",
                               "stat.dist"),
                    discretization.method = "doane")
}
\arguments{
\item{dag}{a matrix or a formula statement (see details for format) defining
the network structure, a directed acyclic graph (DAG).
Note that rownames must be set or given in \code{data.dist} if the DAG is
given via a formula statement.}

\item{data.df}{a data frame containing the data used for learning each node,
binary variables must be declared as factors.}

\item{data.dists}{a named list giving the distribution for each node in the
network, see \sQuote{Details}.}

\item{method}{the method to be used. See \sQuote{Details}.}

\item{discretization.method}{a character vector giving the discretization
method to use. See \code{\link{discretization}}.}
}
\value{
The function returns a named matrix with the requested metric.
}
\description{
A flexible implementation of multiple proxy for strength measures useful for
visualizing the edge connections in a Bayesian Network learned from observational data.
}
\details{
This function returns multiple proxies for estimating the connection strength
of the edges of a possibly discretized Bayesian network's data set.
The returned connection strength measures are: the Raw Mutual Information
(\code{mi.raw}), the Percentage Mutual information (\code{mi.raw.pc}),
the Raw Mutual Information computed via correlation (\code{mi.corr}),
the link strength (\code{ls}), the percentage link strength (\code{ls.pc})
and the statistical distance (\code{stat.dist}).

The general concept of entropy is defined for probability distributions.
The probability is estimated from data using frequency tables.
Then the estimates are plug-in in the definition of the entropy to return
the so-called empirical entropy. A standard known problem of empirical entropy
is that the estimations are biased due to the sampling noise.
This is also known that the bias will decrease as the sample size increases.
The mutual information estimation is computed from the observed frequencies
through a plug-in estimator based on entropy.
For the case of an arc going from the node X to the node Y and the remaining
set of parent of Y is denoted as Z.

The mutual information is defined as I(X, Y) = H(X) + H(Y) - H(X, Y),
where H() is the entropy.

The Percentage Mutual information is defined as PI(X,Y) = I(X,Y)/H(Y|Z).

The Mutual Information computed via correlation is defined as
MI(X,Y) = -0.5 log(1-cor(X,Y)^2).

The link strength is defined as LS(X->Y) = H(Y|Z)-H(Y|X,Z).

The percentage link strength is defined as PLS(X->Y) = LS(X->Y) / H(Y|Z).

The statistical distance is defined as SD(X,Y) = 1- MI(X,Y) / max(H(X),H(Y)).
}
\examples{
# Gaussian
N <- 1000
mydists <- list(a="gaussian",
                b="gaussian",
                c="gaussian")
a <- rnorm(n = N, mean = 0, sd = 1)
b <- 1 + 2*rnorm(n = N, mean = 5, sd = 1)
c <- 2 + 1*a + 2*b + rnorm(n = N, mean = 2, sd = 1)
mydf <- data.frame("a" = a,
                   "b" = b,
                   "c" = c)
mycache.mle <- buildScoreCache(data.df = mydf,
                               data.dists = mydists,
                               method = "mle",
                               max.parents = 2)
mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
linkstr <- linkStrength(dag = mydag.mp$dag,
                        data.df = mydf,
                        data.dists = mydists,
                        method = "ls",
                        discretization.method = "sturges")
}
\references{
Boerlage, B. (1992).  Link strength in Bayesian networks. Diss. University of British Columbia.
Ebert-Uphoff, Imme. "Tutorial on how to measure link strengths in discrete Bayesian networks." (2009).
}
\keyword{utilities}
