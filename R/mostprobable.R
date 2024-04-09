###############################################################################
##
##    Copyright (C) 2002  Fraser Lewis
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################
#' Find most probable DAG structure
#'
#' Find most probable DAG structure using exact order based approach due to Koivisto and Sood, 2004.
#'
#' @usage mostProbable(score.cache, score="bic", prior.choice=1, verbose=TRUE, ...)
#'
#' @param score.cache object of class \code{abnCache} typically outputted by from \code{buildScoreCache()}.
#' @param score which score should be used to score the network. Possible choices are \code{aic, bic, mdl, mlik}.
#' @param prior.choice an integer, 1 or 2, where 1 is a uniform structural prior and 2 uses a weighted prior, see details
#' @param verbose if TRUE then provides some additional output.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The procedure runs the exact order based structure discovery approach of Koivisto and Sood (2004) to find the most probable posterior network (DAG). The local.score is the node cache, as created using \code{\link{buildScoreCache}} (or manually provided the same format is used). Note that the scope of this search is given by the options used in local.score, for example, by restricting the number of parents or the ban or retain constraints given there.
#'
#' This routine can take a long time to complete and is highly sensitive to the number of nodes in the network. It is recommended to use this on a reduced data set to get an idea as to the computational practicality of this approach.  In particular, memory usage can quickly increase to beyond what may be available. For additive models, problems comprising up to 20 nodes are feasible on most machines. Memory requirements can increase considerably after this, but then so does the run time making this less practical. It is recommended that some form of over-modeling adjustment is performed on this resulting DAG (unless dealing with vast numbers of observations), for example, using parametric bootstrapping, which is straightforward to implement in MCMC engines such as JAGS or WinBUGS. See the case studies at \url{https://r-bayesian-networks.org/}
#' or the files provided in the package directories \code{inst/bootstrapping_example} and \code{inst/old_vignette}
#' for details.
#'
#' The parameter \code{prior.choice} determines the prior used within each node for a given choice of parent combination. In Koivisto and Sood (2004) p.554, a form of prior is used, which assumes that the prior probability for parent combinations comprising of the same number of parents are all equal. Specifically, that the prior probability for parent set G with cardinality |G| is proportional to \code{1/[n-1 choose |G|]} where there are n total nodes. Note that this favors parent combinations with either very low or very high cardinality, which may not be appropriate. This prior is used when \code{prior.choice=2}. When \code{prior.choice=1} an uninformative prior is used where parent combinations of all cardinalities are equally likely. The latter is equivalent to the structural prior used in the heuristic searches e.g., \code{searchHillclimber} or \code{searchHeuristic}.
#'
#' Note that the network score (log marginal likelihood) of the most probable DAG is not returned as it can easily be computed using \code{\link{fitAbn}}, see examples below.
#'
#' @return An object of class \code{abnMostprobable}, which is a list containing: a matrix giving the DAG definition of the most probable posterior structure, the cache of pre-computed scores and the score used for selection.
#' @export
#' @references Koivisto, M. V. (2004). Exact Structure Discovery in Bayesian Networks, Journal of Machine Learning Research, vol 5, 549-573.
#'
#' @examples
#' \dontrun{
#' ##############################
#' ## Example 1
#' ##############################
#' ## This data comes with 'abn' see ?ex1.dag.data
#' mydat <- ex1.dag.data[1:5000, c(1:7, 10)]
#'
#' ## Setup distribution list for each node:
#' mydists <- list(b1 = "binomial",
#'                 p1 = "poisson",
#'                 g1 = "gaussian",
#'                 b2 = "binomial",
#'                 p2 = "poisson",
#'                 b3 = "binomial",
#'                 g2 = "gaussian",
#'                 g3 = "gaussian")
#'
#' ## Parent limits, for speed purposes quite specific here:
#' max_par <- list("b1" = 0,
#'                 "p1" = 0,
#'                 "g1" = 1,
#'                 "b2" = 1,
#'                 "p2" = 2,
#'                 "b3" = 3,
#'                 "g2" = 3,
#'                 "g3" = 2)
#' ## Now build cache (no constraints in ban nor retain)
#' mycache <- buildScoreCache(data.df = mydat,
#'                            data.dists = mydists,
#'                            max.parents = max_par)
#'
#' ## Find the globally best DAG:
#' mp_dag <- mostProbable(score.cache = mycache)
#' myres <- fitAbn(object = mp_dag,
#'                 create.graph = TRUE)
#' plot(myres) # plot the best model
#'
#' ## Fit the known true DAG (up to variables 'b4' and 'b5'):
#' true_dag <- matrix(data = 0, ncol = 8, nrow = 8)
#' colnames(true_dag) <- rownames(true_dag) <- names(mydists)
#'
#' true_dag["p2", c("b1", "p1")] <- 1
#' true_dag["b3", c("b1", "g1", "b2")] <- 1
#' true_dag["g2", c("p1", "g1", "b2")] <- 1
#' true_dag["g3", c("g1", "b2")] <- 1
#'
#' fitAbn(dag = true_dag,
#'        data.df = mydat,
#'        data.dists = mydists)$mlik
#'
#' #################################################################
#' ## Example 2 - models with random effects
#' #################################################################
#' ## This data comes with abn see ?ex3.dag.data
#' mydat <- ex3.dag.data[, c(1:4, 14)]
#' mydists <- list(b1 = "binomial",
#'                 b2 = "binomial",
#'                 b3 = "binomial",
#'                 b4 = "binomial")
#'
#' ## This takes a few seconds and requires INLA:
#' mycache_mixed <- buildScoreCache(data.df = mydat,
#'                                  data.dists = mydists,
#'                                  group.var = "group",
#'                                  max.parents = 2)
#'
#' ## Find the most probable DAG:
#' mp_dag <- mostProbable(score.cache = mycache_mixed)
#' ## and get goodness of fit:
#' fitAbn(object = mp_dag,
#'        group.var = "group")$mlik
#' }
#' @keywords models
mostProbable <- function(score.cache, score="bic", prior.choice=1,
                         verbose=TRUE, ...) {

  if (!inherits(score.cache,"abnCache")) {
    stop("score.cache should be an object of class 'abnCache' ")
  }
  score <- c("mlik","aic","bic", "mdl")[pmatch(tolower(score), c("mlik","aic","bic","mdl"))][1]
  if (is.na(score)) stop("wrong specification of 'score'.")

  ### FIXME: error here? it seems that the following is only required for method='mle'
  if (score.cache$method=='mle') {
    if(score=="aic"){score.cache$mlik <- (-score.cache$aic)}
    if(score=="bic"){score.cache$mlik <- (-score.cache$bic)}
    if(score=="mdl"){score.cache$mlik <- (-score.cache$mdl)}
  }

  data.df <- score.cache$data.df[,names(score.cache$data.dists)]; ## n.b. this might be adjusted from original data.df ! when adjusting for random effect

  loc.numnodes <- as.integer(dim(score.cache$node.defn)[2]);
  loc.maxparents <- max(apply(score.cache$node.defn,1,sum)); ## maximum number of parents in any node
  score.cache$children <- as.integer(score.cache$children-1); ## since C indexes from 0

  ## check for missing values - check both NA and NaN - should be just the latter but it may be possible
  ## I guess for these to switch back and forth between R and C
  score.cache$mlik <- ifelse(is.nan(score.cache$mlik),-.Machine$double.xmax,score.cache$mlik);## if node calc gave a NaN
  score.cache$mlik <- ifelse(is.na(score.cache$mlik),-.Machine$double.xmax,score.cache$mlik);## if node calc gave a NA

  if(is.null(data.df)){stop("Must provide data.df - data used in call to mostprobable()");}

  ## need the number of combinations per node
  loc.num.subsets <- as.integer(table(score.cache$children));
  if(length(loc.num.subsets)!=dim(data.df)[2]){stop("At least one node has no valid parent combinations given constraints applied!");}
  ##now get indexes where end node starts and stops
  loc.end <- cumsum(c(table(score.cache$children)));
  loc.start <- c(1,loc.end[-length(loc.end)]+1);
  loc.end <- loc.end-1;#C from 0
  loc.end <- as.integer(loc.end);
  loc.start <- loc.start-1;#C from 0
  loc.start <- as.integer(loc.start);

  if(prior.choice != 1 && prior.choice != 2){stop("prior choice must be 1 or 2!\n");}

  res.prob <- .Call("mostprobable_C",score.cache,loc.numnodes,loc.start,loc.end, as.integer(prior.choice), verbose
                    ,PACKAGE="abn" ## uncomment to load as package not shlib
  )
  loc.res <- matrix(data=res.prob[[1]], ncol=loc.numnodes, byrow=TRUE)
  colnames(loc.res) <- rownames(loc.res) <- names(data.df)

  junk <- gc(FALSE) ## some garbage collection

  out <- list(dag=(loc.res), score.cache=score.cache, score=score)
  class(out) <- c("abnMostprobable","abnLearned")
  return(out)

}


