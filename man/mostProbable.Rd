% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mostprobable.R
\name{mostProbable}
\alias{mostProbable}
\title{Find most probable DAG structure}
\usage{
mostProbable(score.cache, score="bic", prior.choice=1, verbose=TRUE, ...)
}
\arguments{
\item{score.cache}{object of class \code{abnCache} typically outputted by from \code{buildScoreCache()}.}

\item{score}{which score should be used to score the network. Possible choices are \code{aic, bic, mdl, mlik}.}

\item{prior.choice}{an integer, 1 or 2, where 1 is a uniform structural prior and 2 uses a weighted prior, see details}

\item{verbose}{if TRUE then provides some additional output.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
An object of class \code{abnMostprobable}, which is a list containing: a matrix giving the DAG definition of the most probable posterior structure, the cache of pre-computed scores and the score used for selection.
}
\description{
Find most probable DAG structure using exact order based approach due to Koivisto and Sood, 2004.
}
\details{
The procedure runs the exact order based structure discovery approach of Koivisto and Sood (2004) to find the most probable posterior network (DAG). The local.score is the node cache, as created using \code{\link{buildScoreCache}} (or manually provided the same format is used). Note that the scope of this search is given by the options used in local.score, for example, by restricting the number of parents or the ban or retain constraints given there.

This routine can take a long time to complete and is highly sensitive to the number of nodes in the network. It is recommended to use this on a reduced data set to get an idea as to the computational practicality of this approach.  In particular, memory usage can quickly increase to beyond what may be available. For additive models, problems comprising up to 20 nodes are feasible on most machines. Memory requirements can increase considerably after this, but then so does the run time making this less practical. It is recommended that some form of over-modeling adjustment is performed on this resulting DAG (unless dealing with vast numbers of observations), for example, using parametric bootstrapping, which is straightforward to implement in MCMC engines such as JAGS or WinBUGS. See the case studies at \url{https://r-bayesian-networks.org/}
or the files provided in the package directories \code{inst/bootstrapping_example} and \code{inst/old_vignette}
for details.

The parameter \code{prior.choice} determines the prior used within each node for a given choice of parent combination. In Koivisto and Sood (2004) p.554, a form of prior is used, which assumes that the prior probability for parent combinations comprising of the same number of parents are all equal. Specifically, that the prior probability for parent set G with cardinality |G| is proportional to \code{1/[n-1 choose |G|]} where there are n total nodes. Note that this favors parent combinations with either very low or very high cardinality, which may not be appropriate. This prior is used when \code{prior.choice=2}. When \code{prior.choice=1} an uninformative prior is used where parent combinations of all cardinalities are equally likely. The latter is equivalent to the structural prior used in the heuristic searches e.g., \code{searchHillclimber} or \code{searchHeuristic}.

Note that the network score (log marginal likelihood) of the most probable DAG is not returned as it can easily be computed using \code{\link{fitAbn}}, see examples below.
}
\examples{
\dontrun{
##############################
## Example 1
##############################
## This data comes with 'abn' see ?ex1.dag.data
mydat <- ex1.dag.data[1:5000, c(1:7, 10)]

## Setup distribution list for each node:
mydists <- list(b1 = "binomial",
                p1 = "poisson",
                g1 = "gaussian",
                b2 = "binomial",
                p2 = "poisson",
                b3 = "binomial",
                g2 = "gaussian",
                g3 = "gaussian")

## Parent limits, for speed purposes quite specific here:
max_par <- list("b1" = 0,
                "p1" = 0,
                "g1" = 1,
                "b2" = 1,
                "p2" = 2,
                "b3" = 3,
                "g2" = 3,
                "g3" = 2)
## Now build cache (no constraints in ban nor retain)
mycache <- buildScoreCache(data.df = mydat,
                           data.dists = mydists,
                           max.parents = max_par)

## Find the globally best DAG:
mp_dag <- mostProbable(score.cache = mycache)
myres <- fitAbn(object = mp_dag,
                create.graph = TRUE)
plot(myres) # plot the best model

## Fit the known true DAG (up to variables 'b4' and 'b5'):
true_dag <- matrix(data = 0, ncol = 8, nrow = 8)
colnames(true_dag) <- rownames(true_dag) <- names(mydists)

true_dag["p2", c("b1", "p1")] <- 1
true_dag["b3", c("b1", "g1", "b2")] <- 1
true_dag["g2", c("p1", "g1", "b2")] <- 1
true_dag["g3", c("g1", "b2")] <- 1

fitAbn(dag = true_dag,
       data.df = mydat,
       data.dists = mydists)$mlik

#################################################################
## Example 2 - models with random effects
#################################################################
## This data comes with abn see ?ex3.dag.data
mydat <- ex3.dag.data[, c(1:4, 14)]
mydists <- list(b1 = "binomial",
                b2 = "binomial",
                b3 = "binomial",
                b4 = "binomial")

## This takes a few seconds and requires INLA:
mycache_mixed <- buildScoreCache(data.df = mydat,
                                 data.dists = mydists,
                                 group.var = "group",
                                 max.parents = 2)

## Find the most probable DAG:
mp_dag <- mostProbable(score.cache = mycache_mixed)
## and get goodness of fit:
fitAbn(object = mp_dag,
       group.var = "group")$mlik
}
}
\references{
Koivisto, M. V. (2004). Exact Structure Discovery in Bayesian Networks, Journal of Machine Learning Research, vol 5, 549-573.
}
\keyword{models}
