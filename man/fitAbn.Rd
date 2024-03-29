% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitabn.R, R/fitabn_bayes.R, R/fitabn_mle.R
\name{fitAbn}
\alias{fitAbn}
\alias{fitAbn.bayes}
\alias{fitAbn.mle}
\alias{regressionLoop}
\title{Fit an additive Bayesian network model}
\usage{
fitAbn(object = NULL,
       dag = NULL,
       data.df = NULL,
       data.dists = NULL,
       method = NULL,
       group.var = NULL,
       adj.vars = NULL,
       cor.vars = NULL,
       centre = TRUE,
       compute.fixed = FALSE,
       control = NULL,
       verbose = FALSE,
       debugging = FALSE,
       ...)

fitAbn.bayes(
  dag = NULL,
  data.df = NULL,
  data.dists = NULL,
  group.var = NULL,
  cor.vars = NULL,
  centre = TRUE,
  compute.fixed = FALSE,
  control = fit.control(method = "bayes"),
  mylist = NULL,
  grouped.vars = NULL,
  group.ids = NULL,
  force.method = NULL,
  verbose = FALSE,
  debugging = FALSE
)

fitAbn.mle(
  dag = NULL,
  data.df = NULL,
  data.dists = NULL,
  group.var = NULL,
  grouped.vars = NULL,
  group.ids = NULL,
  adj.vars = NULL,
  cor.vars = NULL,
  centre = TRUE,
  control = fit.control(method = "mle"),
  verbose = FALSE,
  debugging = FALSE
)

regressionLoop(
  i = NULL,
  dag = NULL,
  data.df = NULL,
  data.df.multi = NULL,
  data.dists = NULL,
  group.var = NULL,
  grouped.vars = NULL,
  group.ids = NULL,
  control = NULL,
  nvars = NULL,
  nobs = NULL,
  dag.multi = NULL,
  verbose = NULL
)
}
\arguments{
\item{object}{an object of class \code{abnLearned} produced by \code{\link{mostProbable}}, \code{\link{searchHeuristic}} or \code{\link{searchHillClimber}}.}

\item{dag}{a matrix or a formula statement (see details) defining the network structure, a directed acyclic graph (DAG), see details for format. Note that column names and row names must be set up.}

\item{data.df}{a data frame containing the data used for learning the network, binary variables must be declared as factors, and no missing values all allowed in any variable.}

\item{data.dists}{a named list giving the distribution for each node in the network, see details.}

\item{method}{if \code{NULL}, takes method of \code{object}, otherwise \code{"bayes"} or \code{"mle"} for the method to be used, see details.}

\item{group.var}{only applicable for mixed models and gives the column name in \code{data.df} of the grouping variable (which must be a factor denoting group membership).}

\item{adj.vars}{a character vector giving the column names in \code{data.df} for which the network score has to be adjusted for, see details.}

\item{cor.vars}{a character vector giving the column names in data.df for which a mixed model should be used (\code{method = 'bayes'} only).}

\item{centre}{should the observations in each Gaussian node first be standardised to mean zero and standard deviation one.}

\item{compute.fixed}{a logical flag, set to \code{TRUE} for computation of marginal posterior distributions, see details.}

\item{control}{a list of control parameters. See \code{\link{fit.control}} for the names of the settable control values and their effect.}

\item{verbose}{if \code{TRUE} then provides some additional output, in particular the code used to call INLA, if applicable.}

\item{debugging}{if \code{TRUE} and \code{method = 'mle'} this enables to step into the for-loop.}

\item{\dots}{additional arguments passed for optimization.}

\item{mylist}{result returned from \code{\link{check.valid.data}}.}

\item{grouped.vars}{result returned from \code{\link{check.valid.groups}}. Column indexes of all variables which are affected from grouping effect.}

\item{group.ids}{result returned from \code{\link{check.valid.groups}}. Vector of group allocation for each observation (row) in 'data.df'.}

\item{force.method}{"notset", "INLA" or "C". This is specified in \code{\link{buildScoreCache}(control=list(max.mode.error=...))}.}

\item{i}{number of child-node (mostly corresponds to child node index e.g. in dag).}

\item{data.df.multi}{extended data.df for one-hot-encoded multinomial variables.}

\item{nvars}{number of variables in data.dists.}

\item{nobs}{number of observations in data.df.}

\item{dag.multi}{extended dag for one-hot-encoded multinomial variables.}
}
\value{
An object of class \code{abnFit}. A named list. One entry for each of the variables in \code{data.df} (excluding the grouping variable, if present) which contains an estimate of the log marginal likelihood for that individual node. An entry "mlik" which is the total log marginal likelihood for the full ABN model. A vector of \code{error.codes} - non-zero if a numerical error or warning occurred, and a vector error.code.desc giving a text description of the error. A list \code{modes}, which contains all the mode estimates for each parameter at each node. A vector called Hessian accuracy, which is the estimated local error in the log marginal likelihood for each node.  If \code{compute.fixed=TRUE} then a list entry called \code{marginals} which contains a named entry for every parameter in the ABN and each entry in this list is a two-column matrix where the first column is the value of the marginal parameter, say x, and the second column is the respective density value, pdf(x). Also, a list called \code{marginal.quantiles} is produced, giving the quantiles for each marginal posterior distribution.

list
}
\description{
Fits an additive Bayesian network to observed data and is equivalent to Bayesian or information-theoretic multi-dimensional regression modeling.
Two numerical options are available in the Bayesian settings, standard Laplace approximation or else an integrated nested Laplace approximation provided via a call to the R INLA library (see \href{https://www.r-inla.org/}{r-inla.org} - this is not hosted on CRAN).
}
\details{
\subsection{If \code{method="Bayes"}:}{

The procedure \code{fitAbn} fits an additive Bayesian network model to data where each node (variable - a column in data.df) can be either: presence/absence (Bernoulli); continuous (Gaussian); or an unbounded count (Poisson). Multinomial distributions are only supported with \code{method = "mle"} (see below).
The model comprises of a set of conditionally independent generalized linear regressions with or without random effects.
Internal code is used by default for numerical estimation in nodes without random effects, and INLA is the default for nodes with random effects.
This default behavior can be overridden using \code{control=list(max.mode.error=...)}. The default is \code{max.mode.error=10}, which means that the modes estimated from INLA output must be within 10\\% of those estimated using internal code.
Otherwise, the internal code is used rather than INLA.
To force the use of INLA on all nodes, use \code{max.mode.error=100}, which then ignores this check, to force the use of internal code then use \code{max.mode.error=0}.
For the numerical reliability and perform of \pkg{abn} see \url{https://r-bayesian-networks.org/}.
Generally speaking, INLA can be swift and accurate, but in several cases, it can perform very poorly and so some care is required (which is why there is an internal check on the modes).
Binary variables must be declared as factors with two levels, and the argument \code{data.dists} must be a list with named arguments, one for each of the variables in \code{data.df} (except a grouping variable - if present!), where each entry is either "poisson","binomial", "multinomial" or "gaussian", see examples below.
The "poisson" and "binomial" distributions use log and logit link functions, respectively.
Note that "binomial" here actually means only binary, one Bernoulli trial per row in \code{data.df}.

If the data are grouped into correlated blocks - wherein a standard regression context a mixed model might be used - then a network comprising of one or more nodes where a generalized linear mixed model is used (but limited to only a single random effect).
This is achieved by specifying parameters \code{group.var} and \code{cor.vars}.
Where the former defines the group membership variable, which should be a factor indicating which observations belong to the same grouping.
The parameter \code{cor.vars} is a character vector that contains the names of the nodes for which a mixed model should be used. This is not yet implemented with \code{method = 'mle'}.
For example, in some problems, it may be appropriate for all variables (except \code{group.var}) in data.df to be parametrized as a mixed model while in others it may only be a single variable for which grouping adjustment is required (as the remainder of variables are covariates measured at group level).

In the network structure definition, \code{dag}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.
The \code{dag} can be provided using a formula statement (similar to GLM).
A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}.
The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3.
The parents names must match those given in \code{data.df}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{data.df}.

If \code{compute.fixed=TRUE} then the marginal posterior distributions for all parameters are computed.
Note the current algorithm used to determine the evaluation grid is rather crude and may need to be manually refined using \code{variate.vec} (one parameter at a time) for publication-quality density estimates.
Note that a manual grid can only be used with internal code and not INLA (which uses its own grid).
The end points are defined as where the value of the marginal density drops below a given threshold \code{pdf.min}.
When estimating the log marginal likelihood in models with random effects (using internal code rather than INLA), an attempt is made to minimize the error by comparing the estimates given between a 3pt and 5pt rule when estimating the Hessian in the Laplace approximation.
The modes used in each case are identical. The first derivatives are computed using gsl's adaptive finite difference function, and this is embedding inside the standard 3pt and 5pt rules for the second derivatives.
In all cases, a central difference approximation is tried first with a forward difference being a fall back (as the precision parameters are strictly positive).
The error is minimized through choosing an optimal step size using gsl's Nelder-Mead optimization, and if this fails, (e.g., is larger than \code{max.hessian.error}) then the Brent-Dekker root bracketing method is used as a fallback.
If the error cannot be reduced to below \code{max.hessian.error}, then the step size, which gave the lowest error during the searches (across potentially many different initial bracket choices), is used for the final Hessian evaluations in the Laplace approximation.
}

\subsection{If \code{method="mle"}:}{

The procedure \code{fitAbn} with the argument \code{method= "mle"} fits an additive Bayesian network model to data where each node (variable - a column in data.df) can be either: presence/absence (Bernoulli); continuous (Gaussian); an unbounded count (Poisson); or a discrete variable (Multinomial).
The model comprises of a set of conditionally independent generalized linear regressions with or without adjustment.
Binary and discrete variables must be declared as factors and the argument \code{data.dists} must be a list with named arguments, one for each of the variables in \code{data.df}, where each entry is either "poisson","binomial", "multinomial" or "gaussian", see examples below.
The "poisson" and "binomial" distributions use log and logit link functions, respectively.
Note that "binomial" here actually means only binary, one Bernoulli trial per row in data.df.

If the data are grouped into correlated blocks - wherein a standard regression context a mixed-effect model might be used - then a network comprising of one or more nodes where a generalized linear mixed model is used (but limited to only a single random intercept).
This is achieved by specifying parameter \code{group.var} (\code{cor.vars} as with \code{method = "bayes"} is not yet implemented with \code{method = "mle"}).
The parameter \code{group.var} defines the group membership variable, which should be a factor indicating which observations belong to the same grouping.
This corresponds to \code{"1|group.var"} in the formula notation of e.g. \pkg{lme4}.

In the context of \code{fitAbn} adjustment means that irrespective to the adjacency matrix the adjustment variable set (\code{adj.vars}) will be add as covariate to every node defined by \code{cor.vars}.

In the network structure definition, \code{dag}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.
The \code{dag} can be provided using a formula statement (similar to GLM). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}.
The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3.
The parents names have to exactly match those given in \code{data.df}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{data.df}.

The Information-theoretic based network scores used in \code{fitAbn} with argument \code{method="mle"} are the maximum likelihood (mlik, called marginal likelihood in this context as it is computed node wise), the Akaike Information Criteria (aic), the Bayesian Information Criteria (bic) and the Minimum distance Length (mdl). The classical definitions of those metrics are given in Kratzer and Furrer (2018).

The numerical routine is based on an iterative scheme to estimate the regression coefficients. The Iterative Reweighed Least Square (IRLS) programmed using Rcpp/RcppArmadrillo. One hard coded feature of \code{fitAbn} with argument \code{method="mle"} is a conditional use of a bias reduced binomial regression when a classical Generalized Linear Model (GLM) fails to estimate the maximum likelihood of the given model accurately. Additionally, a QR decomposition is performed to check for rank deficiency. If the model is rank deficient and the BR GLM fails to estimate it, then predictors are sequentially removed. This feature aims at better estimating network scores when data sparsity is present.

A special care should be taken when interpreting or even displaying p-values computed with \code{fitAbn}. Indeed, the full model is already selected using goodness of fit metrics based on the (same) full dataset.

The \code{control} argument is a list with separate arguments for the Bayesian and MLE implementation. See  \code{\link{fit.control}} for details.
}
}
\section{Functions}{
\itemize{
\item \code{fitAbn.bayes()}: Internal function called by \code{fitAbn}.

\item \code{fitAbn.mle()}: Internal function called by \code{fitAbn}.

\item \code{regressionLoop()}: Internal function called by \code{fitAbn.mle()}.

}}
\examples{
\dontrun{
## Built-in dataset with a subset of cols
mydat <- ex0.dag.data[, c("b1", "b2", "b3", "g1", "b4", "p2", "p4")]

## setup distribution list for each node
mydists <- list(b1 = "binomial",
                b2 = "binomial",
                b3 = "binomial",
                g1 = "gaussian",
                b4 = "binomial",
                p2 = "poisson",
                p4 = "poisson")

## Null model - all independent variables
mydag_empty <- matrix(0, nrow = 7, ncol = 7)
colnames(mydag_empty) <- rownames(mydag_empty) <- names(mydat)

## Now fit the model to calculate its goodness-of-fit
myres <- fitAbn(dag = mydag_empty,
                data.df = mydat,
                data.dists = mydists)

## Log-marginal likelihood goodness-of-fit for complete DAG
print(myres$mlik)

## fitAbn accepts also the formula statement
myres <- fitAbn(dag = ~ b1 | b2 + b2 | p4:g1 + g1 | p2 + b3 | g1 + b4 | b1 + p4 | g1,
                data.df = mydat,
                data.dists = mydists)
print(myres$mlik) # a much weaker fit than full independence DAG

# Plot the DAG via Rgraphviz
plot(myres)

## Or equivalently using the formula statement, with plotting
## Now repeat but include some dependencies first
mydag <- mydag_empty
mydag["b1", "b2"] <- 1 # b1<-b2 and so on
mydag["b2", "p4"] <- mydag["b2", "g1"] <- mydag["g1", "p2"] <- 1
mydag["b3", "g1"] <- mydag["b4", "b1"] <- mydag["p4", "g1"] <- 1
myres_alt <- fitAbn(dag = mydag,
                    data.df = mydat,
                    data.dists = mydists)
plot(myres_alt)

## -----------------------------------------------------------------------------
## This function contains an MLE implementation accessible through a method
## parameter use built-in simulated data set
## -----------------------------------------------------------------------------
myres_mle <- fitAbn(dag = ~ b1 | b2 + b2 | p4 + g1 + g1 | p2 + b3 | g1 + b4 | b1 + p4 | g1,
                    data.df = mydat,
                    data.dists = mydists,
                    method = "mle")

## Print the output for mle first then for Bayes:
print(myres_mle)
plot(myres_mle)

print(myres)
plot(myres)

## This is a basic plot of some posterior densities. The algorithm used for
## selecting density points is quite straightforward, but it might result
## in a sparse distribution. Therefore, we also recompute the density over
## an evenly spaced grid of 50 points between the two endpoints that had
## a minimum PDF at f = min.pdf.
## Setting max.mode.error = 0 forces the use of the internal C code.
myres_c <- fitAbn(dag = mydag,
                  data.df = mydat,
                  data.dists = mydists,
                  compute.fixed = TRUE,
                  control = list(max.mode.error = 0))

print(names(myres_c$marginals)) # gives all the different parameter names

## Repeat but use INLA for the numerics using max.mode.error = 100
## as using internal code is the default here rather than INLA
myres_inla <- fitAbn(dag = mydag,
                     data.df = mydat,
                     data.dists = mydists,
                     compute.fixed = TRUE,
                     control = list(max.mode.error = 100))

## Plot posterior densities
default_par <- par(no.readonly = TRUE) # save default par settings
par(mfrow = c(2, 2), mai = c(.7, .7, .2, .1))
plot(myres_c$marginals$b1[["b1 | (Intercept)"]], type = "l", xlab = "b1 | (Intercept)")
lines(myres_inla$marginals$b1[["b1 | (Intercept)"]], col = "blue")
plot(myres_c$marginals$b2[["b2 | p4"]], type = "l", xlab = "b2 | p4")
lines(myres_inla$marginals$b2[["b2 | p4"]], col = "blue")
plot(myres_c$marginals$g1[["g1 | precision"]], type = "l", xlab = "g1 | precision")
lines(myres_inla$marginals$g1[["g1 | precision"]], col = "blue")
plot(myres_c$marginals$b4[["b4 | b1"]], type = "l", xlab = "b4 | b1")
lines(myres_inla$marginals$b4[["b4 | b1"]], col = "blue")
par(default_par) # reset par settings

## An elementary mixed model example using built-in data specify DAG,
## only two variables using a subset of variables from ex3.dag.data
## both variables are assumed to need (separate) adjustment for the
## group variable, i.e., a binomial GLMM at each node

mydists <- list(b1 = "binomial",
                b2 = "binomial")

## Compute marginal likelihood - use internal code via max.mode.error=0
## as using INLA is the default here.
## Model where b1 <- b2
myres_c <- fitAbn(dag = ~b1 | b2,
                  data.df = ex3.dag.data[, c(1, 2, 14)],
                  data.dists = mydists,
                  group.var = "group",
                  cor.vars = c("b1", "b2"),
                  control = list(max.mode.error = 0))
print(myres_c) # show all the output

## compare mode for node b1 with glmer(), lme4::glmer is automatically attached.

## Now for marginals - INLA is strongly preferable for estimating marginals for
## nodes with random effects as it is far faster, but may not be reliable
## see https://r-bayesian-networks.org/

## INLA's estimates of the marginals, using high n.grid = 500
## as this makes the plots smoother - see below.
myres_inla <- fitAbn(dag = ~b1 | b2,
                   data.df = ex3.dag.data[, c(1, 2, 14)],
                  data.dists = mydists,
                  group.var = "group",
                  cor.vars = c("b1", "b2"),
                  compute.fixed = TRUE,
                  n.grid = 500,
                  control = list(max.mode.error = 100,
                                 max.hessian.error = 10E-02))

## this is NOT recommended - marginal density estimation using fitAbn in
## mixed models is really just for diagnostic purposes, better to use
## fitAbn.inla() here; but here goes... be patient
myres_c <- fitAbn(dag = ~b1 | b2,
                  data.df = ex3.dag.data[, c(1, 2, 14)],
                  data.dists = mydists,
                  group.var = "group",
                  cor.vars = c("b1", "b2"),
                  compute.fixed = TRUE,
                  control = list(max.mode.error = 0,
                                 max.hessian.error = 10E-02))

## compare marginals between internal and INLA.
default_par <- par(no.readonly = TRUE) # save default par settings
par(mfrow = c(2, 3))
# 5 parameters - two intercepts, one slope, two group level precisions
plot(myres_inla$marginals$b1[[1]], type = "l", col = "blue")
lines(myres_c$marginals$b1[[1]], col = "brown", lwd = 2)
plot(myres_inla$marginals$b1[[2]], type = "l", col = "blue")
lines(myres_c$marginals$b1[[2]], col = "brown", lwd = 2)
# the precision of group-level random effects
plot(myres_inla$marginals$b1[[3]], type = "l", col = "blue", xlim = c(0, 2))
lines(myres_c$marginals$b1[[3]], col = "brown", lwd = 2)
plot(myres_inla$marginals$b2[[1]], type = "l", col = "blue")
lines(myres_c$marginals$b2[[1]], col = "brown", lwd = 2)
plot(myres_inla$marginals$b2[[1]], type = "l", col = "blue")
lines(myres_c$marginals$b2[[1]], col = "brown", lwd = 2)
# the precision of group-level random effects
plot(myres_inla$marginals$b2[[2]], type = "l", col = "blue", xlim = c(0, 2))
lines(myres_c$marginals$b2[[2]], col = "brown", lwd = 2)
par(default_par) # reset par settings

### these are very similar although not exactly identical

## use internal code but only to compute a single parameter over a specified
## grid.
## This can be necessary if the simple auto grid finding functions does
## a poor job.
myres_c <- fitAbn(dag = ~b1 | b2,
                  data.df = ex3.dag.data[, c(1, 2, 14)],
                  data.dists = mydists,
                  group.var = "group",
                  cor.vars = c("b1", "b2"),
                  centre = FALSE,
                  compute.fixed = TRUE,
                  control = list(marginal.node = 1,
                                 marginal.param = 3, # precision term in node 1
                                 variate.vec = seq(0.05, 1.5, len = 25),
                                 max.hessian.error = 10E-02))

default_par <- par(no.readonly = TRUE) # save default par settings
par(mfrow = c(1, 2))
plot(myres_c$marginals$b1[[1]], type = "l", col = "blue") # still fairly sparse
# An easy way is to use spline to fill in the density without recomputing other
# points provided the original grid is not too sparse.
plot(spline(myres_c$marginals$b1[[1]], n = 100), type = "b", col = "brown")
par(default_par) # reset par settings
}

}
\references{
Kratzer, G., Lewis, F.I., Comin, A., Pittavino, M. and Furrer, R. (2019). "Additive Bayesian Network Modelling with the R Package abn". arXiv preprint arXiv:1911.09006.

Kratzer, G., and Furrer, R., 2018. Information-Theoretic Scoring Rules to Learn Additive Bayesian Network Applied to Epidemiology. Preprint; Arxiv: stat.ML/1808.01126.

Lewis, F. I., and McCormick, B. J. J. (2012). Revealing the complexity of health determinants in resource poor settings. \emph{American Journal Of Epidemiology}. DOI:10.1093/aje/KWS183.

Further information about \pkg{abn} can be found at: \href{https://r-bayesian-networks.org/}{r-bayesian-networks.org}.
}
\seealso{
\code{\link{buildScoreCache}}

Other fitAbn: 
\code{\link{fit.control}()}

Other Bayes: 
\code{\link{buildScoreCache}()},
\code{\link{calc.node.inla.glm}()},
\code{\link{calc.node.inla.glmm}()},
\code{\link{getmarginals}()}
}
\author{
Fraser Iain Lewis and Gilles Kratzer.
}
\concept{Bayes}
\concept{abn}
\concept{fitAbn}
\keyword{internal}
\keyword{models}
