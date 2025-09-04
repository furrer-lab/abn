#' Control the iterations in \code{\link{buildScoreCache}}
#'
#' Allow the user to set restrictions in the \code{\link{buildScoreCache}} for both the Bayesian and the MLE approach.
#' Control function similar to \code{\link{fit.control}}.
#'
#' @param method a character that takes one of two values: "bayes" or "mle". Overrides \code{method} argument from \code{\link{buildScoreCache}}.
#' @param max.mode.error if the estimated modes from INLA differ by a factor of \code{max.mode.error} or more from those computed internally, then results from INLA are replaced by those computed internally. To force INLA always to be used, then \code{max.mode.error=100}, to force INLA never to be used \code{max.mod.error=0}. See also \code{\link{fitAbn}}.
#' @param mean the prior mean for all the Gaussian additive terms for each node. INLA argument \code{control.fixed=list(mean.intercept=...)} and \code{control.fixed=list(mean=...)}.
#' @param prec the prior precision (\eqn{\tau = \frac{1}{\sigma^2}}) for all the Gaussian additive term for each node. INLA argument \code{control.fixed=list(prec.intercept=...)} and \code{control.fixed=list(prec=...)}.
#' @param loggam.shape the shape parameter in the Gamma distribution prior for the precision in a Gaussian node. INLA argument \code{control.family=list(hyper = list(prec = list(prior="loggamma",param=c(loggam.shape, loggam.inv.scale))))}.
#' @param loggam.inv.scale the inverse scale parameter in the Gamma distribution prior for the precision in a Gaussian node. INLA argument \code{control.family=list(hyper = list(prec = list(prior="loggamma",param=c(loggam.shape, loggam.inv.scale))))}.
#' @param max.iters total number of iterations allowed when estimating the modes in Laplace approximation. passed to \code{.Call("fit_single_node", ...)}.
#' @param epsabs absolute error when estimating the modes in Laplace approximation for models with no random effects. Passed to \code{.Call("fit_single_node", ...)}.
#' @param error.verbose logical, additional output in the case of errors occurring in the optimization. Passed to \code{.Call("fit_single_node", ...)}.
#' @param trace Non-negative integer. If positive, tracing information on the progress of the "L-BFGS-B" optimization is produced. Higher values may produce more tracing information. (There are six levels of tracing.  To understand exactly what these do see the source code.). Passed to \code{.Call("fit_single_node", ...)}.
#' @param epsabs.inner absolute error in the maximization step in the (nested) Laplace approximation for each random effect term. Passed to \code{.Call("fit_single_node", ...)}.
#' @param max.iters.inner total number of iterations in the maximization step in the nested Laplace approximation. Passed to \code{.Call("fit_single_node", ...)}.
#' @param finite.step.size suggested step length used in finite difference estimation of the derivatives for the (outer) Laplace approximation when estimating modes. Passed to \code{.Call("fit_single_node", ...)}.
#' @param hessian.params a numeric vector giving parameters for the adaptive algorithm, which determines the optimal stepsize in the finite-difference estimation of the hessian. First entry is the initial guess, second entry absolute error. Passed to \code{.Call("fit_single_node", ...)}.
#' @param max.iters.hessian integer, maximum number of iterations to use when determining an optimal finite difference approximation (Nelder-Mead). Passed to \code{.Call("fit_single_node", ...)}.
#' @param max.hessian.error if the estimated log marginal likelihood when using an adaptive 5pt finite-difference rule for the Hessian differs by more than \code{max.hessian.error} from when using an adaptive 3pt rule then continue to minimize the local error by switching to the Brent-Dekker root bracketing method. Passed to \code{.Call("fit_single_node", ...)}.
#' @param factor.brent if using Brent-Dekker root bracketing method then define the outer most interval end points as the best estimate of \eqn{h} (stepsize) from the Nelder-Mead as \eqn{h/factor.brent,h*factor.brent)}. Passed to \code{.Call("fit_single_node", ...)}.
#' @param maxiters.hessian.brent maximum number of iterations allowed in the Brent-Dekker method. Passed to \code{.Call("fit_single_node", ...)}.
#' @param num.intervals.brent the number of initial different bracket segments to try in the Brent-Dekker method. Passed to \code{.Call("fit_single_node", ...)}.
#' @param n.grid recompute density on an equally spaced grid with \code{n.grid} points.
#' @param max.irls total number of iterations for estimating network scores using an Iterative Reweighed Least Square algorithm. Is this DEPRECATED?
#' @param tol real number giving the minimal tolerance expected to terminate the Iterative Reweighed Least Square algorithm to estimate network score. Passed to \code{irls_binomial_cpp_fast_br} and \code{irls_poisson_cpp_fast}.
#' @param tolPwrss numeric scalar passed to \code{\link[lme4]{glmerControl}} - the tolerance for declaring convergence in the penalized iteratively weighted residual sum-of-squares step. Similar to \code{tol}.
#' @param check.rankX character passed to \code{\link[lme4]{lmerControl}} and \code{\link[lme4]{glmerControl}} - specifying if \code{rankMatrix(X)} should be compared with \code{ncol(X)} and if columns from the design matrix should possibly be dropped to ensure that it has full rank. Defaults to \code{message+drop.cols}.
#' @param check.scaleX character passed to \code{\link[lme4]{lmerControl}} and \code{\link[lme4]{glmerControl}} - check for problematic scaling of columns of fixed-effect model matrix, e.g. parameters measured on very different scales. Defaults to \code{message+rescale}.
#' @param check.conv.grad character passed to \code{\link[lme4]{lmerControl}} and \code{\link[lme4]{glmerControl}} - checking the gradient of the deviance function for convergence. Defaults to \code{message} but can be one of "ignore" - skip the test; "warning" - warn if test fails; "message" - print a message if test fails; "stop" - throw an error if test fails.
#' @param check.conv.singular character passed to \code{\link[lme4]{lmerControl}} and \code{\link[lme4]{glmerControl}} - checking for a singular fit, i.e. one where some parameters are on the boundary of the feasible space (for example, random effects variances equal to 0 or correlations between random effects equal to +/- 1.0). Defaults to \code{message} but can be one of "ignore" - skip the test; "warning" - warn if test fails; "message" - print a message if test fails; "stop" - throw an error if test fails.
#' @param check.conv.hess character passed to \code{\link[lme4]{lmerControl}} and \code{\link[lme4]{glmerControl}} - checking the Hessian of the deviance function for convergence. Defaults to \code{message} but can be one of "ignore" - skip the test; "warning" - warn if test fails; "message" - print a message if test fails; "stop" - throw an error if test fails.
#' @param xtol_abs Defaults to 1e-6 stop on small change of parameter value. Only for \code{method='mle', group.var=...}. Default convergence tolerance for fitted \code{(g)lmer} models is reduced to the value provided here if default values did not fit. This value here is passed to the \code{optCtrl} argument of \code{(g)lmer} (see help of \code{\link[lme4:convergence]{lme4::convergence()}}).
#' @param ftol_abs  Defaults to 1e-6 stop on small change in deviance. Similar to \code{xtol_abs}.
#' @param trace.mblogit logical indicating if output should be produced for each iteration. Directly passed to \code{trace} argument in \code{\link[mclogit]{mclogit.control}}. Is independent of \code{verbose}.
#' @param catcov.mblogit Defaults to "free" meaning that there are no restrictions on the covariances of random effects between the logit equations. Set to "diagonal" if random effects pertinent to different categories are uncorrelated or "single" if random effect variances pertinent to all categories are identical.
#' @param epsilon Defaults to 1e-8. Positive convergence tolerance \eqn{\epsilon} that is directly passed to the \code{control} argument of \code{mclogit::mblogit} as \code{mclogit.control}. Only for \code{method='mle', group.var=...}.
#' @param ncores The number of cores to parallelize to, see \sQuote{Details}. If >0, the number of CPU cores to be used. -1 for all available -1 core. Only for \code{method="mle"}.
#' @param cluster.type The type of cluster to be used, see \code{?parallel::makeCluster}. \code{abn} then defaults to \code{"PSOCK"} on Windows and \code{"FORK"} on Unix-like systems. With "FORK" the child process are started with \code{rscript_args = "--no-environ"} to avoid loading the whole workspace into each child.
#' @param only_glmmTMB_poisson logical, if TRUE only use \code{glmmTMB} to fit Poisson nodes with random effects. This is useful if \code{glmer} fails due to convergence issues. Default is FALSE.
#' @param seed a non-negative integer which sets the seed in \code{set.seed(seed)}.
#'
#' @details
#' Parallelization over all children is possible via the function \code{foreach} of the package \pkg{doParallel}.  \code{ncores=0} or \code{ncores=1} use single threaded \code{foreach}. \code{ncores=-1} uses all available cores but one.
#'
#' @return Named list according the provided arguments.
#' @seealso \code{\link{fit.control}}.
#' @family buildScoreCache
#' @export build.control
#' @examples
#' ctrlmle <- abn::build.control(method = "mle",
#'                         ncores = 0,
#'                         cluster.type = "PSOCK",
#'                         max.irls = 100,
#'                         tol = 10^-11,
#'                         tolPwrss = 1e-7,
#'                         check.rankX = "message+drop.cols",
#'                         check.scaleX = "message+rescale",
#'                         check.conv.grad = "message",
#'                         check.conv.singular = "message",
#'                         check.conv.hess = "message",
#'                         xtol_abs = 1e-6,
#'                         ftol_abs = 1e-6,
#'                         trace.mblogit = FALSE,
#'                         catcov.mblogit = "free",
#'                         epsilon = 1e-6,
#'                         only_glmmTMB_poisson=FALSE,
#'                         seed = 9062019L)
#' ctrlbayes <- abn::build.control(method = "bayes",
#'                            max.mode.error = 10,
#'                            mean = 0, prec = 0.001,
#'                            loggam.shape = 1,
#'                            loggam.inv.scale = 5e-05,
#'                            max.iters = 100,
#'                            epsabs = 1e-07,
#'                            error.verbose = FALSE,
#'                            epsabs.inner = 1e-06,
#'                            max.iters.inner = 100,
#'                            finite.step.size = 1e-07,
#'                            hessian.params = c(1e-04, 0.01),
#'                            max.iters.hessian = 10,
#'                            max.hessian.error = 0.5,
#'                            factor.brent = 100,
#'                            maxiters.hessian.brent = 100,
#'                            num.intervals.brent = 100,
#'                            tol = 10^-8,
#'                            seed = 9062019L)
#'
build.control <-
  function (method = "bayes",
            max.mode.error = 10,
            mean = 0,
            prec = 0.001,
            loggam.shape = 1,
            loggam.inv.scale = 5e-05,
            max.iters = 100,
            epsabs = 1e-07,
            error.verbose = FALSE,
            trace = 0L,
            epsabs.inner = 1e-06,
            max.iters.inner = 100,
            finite.step.size = 1e-07,
            hessian.params = c(1e-04, 0.01),
            max.iters.hessian = 10,
            max.hessian.error = 0.5,
            factor.brent = 100,
            maxiters.hessian.brent = 100,
            num.intervals.brent = 100,
            n.grid = 250,
            ncores = 1,
            cluster.type = "FORK",
            max.irls = 100,
            tol = 1e-8,
            tolPwrss = 1e-7,
            check.rankX = "message+drop.cols",
            check.scaleX = "message+rescale",
            check.conv.grad = "message",
            check.conv.singular = "message",
            check.conv.hess = "message",
            xtol_abs = 1e-6,
            ftol_abs = 1e-6,
            trace.mblogit = FALSE,
            catcov.mblogit = "free",
            epsilon = 1e-6,
            only_glmmTMB_poisson=FALSE,
            seed = 9062019L) {
    if (method == "bayes") {
      ctrl <-
        list(
          max.mode.error = max.mode.error,
          mean = mean,
          prec = prec,
          loggam.shape = loggam.shape,
          loggam.inv.scale = loggam.inv.scale,
          max.iters = max.iters,
          epsabs = epsabs,
          error.verbose = error.verbose,
          trace = trace,
          epsabs.inner = epsabs.inner,
          max.iters.inner = max.iters.inner,
          finite.step.size = finite.step.size,
          hessian.params = hessian.params,
          max.iters.hessian = max.iters.hessian,
          max.hessian.error = max.hessian.error,
          factor.brent = factor.brent,
          maxiters.hessian.brent = maxiters.hessian.brent,
          num.intervals.brent = num.intervals.brent,
          n.grid = n.grid,
          ncores = ncores,
          cluster.type = cluster.type,
          seed = seed
        )
    } else if (method == "mle") {
      ctrl <-
        list(
          ncores = ncores,
          cluster.type = cluster.type,
          max.iters = max.iters,
          max.irls = max.irls,
          tol = tol,
          tolPwrss = tolPwrss,
          check.rankX = check.rankX,
          check.scaleX = check.scaleX,
          check.conv.grad = check.conv.grad,
          check.conv.singular = check.conv.singular,
          check.conv.hess = check.conv.hess,
          xtol_abs=xtol_abs,
          ftol_abs=ftol_abs,
          trace.mblogit = trace.mblogit,
          catcov.mblogit = catcov.mblogit,
          epsilon = epsilon,
          only_glmmTMB_poisson=only_glmmTMB_poisson,
          seed = seed
        )
    } else {
      stop("`method` unknown.")
    }
    return(ctrl)
  }

## precompute a cache of scores to data
#' Build a cache of goodness of fit metrics for each node in a DAG, possibly subject to user-defined restrictions
#' @keywords buildScoreCache.mle buildScoreCache.bayes calc.node.inla.glm calc.node.inla.glmm fitAbn.bayes fitAbn.mle
#'
#' @description Iterates over all valid parent combinations - subject to ban, retain, and \code{max.parent} limits - for each node, or a subset of nodes, and computes a cache of scores (AIC, BIC, log marginal likelihood).
#' This cache can then be used in different DAG structural search algorithms.
#'
#' @usage
#' buildScoreCache(data.df = NULL,
#' data.dists = NULL,
#' method = "bayes",
#' group.var = NULL,
#' adj.vars = NULL,
#' cor.vars = NULL,
#' dag.banned = NULL,
#' dag.retained = NULL,
#' max.parents = NULL,
#' which.nodes = NULL,
#' defn.res = NULL,
#' centre = TRUE,
#' dry.run = FALSE,
#' control = NULL,
#' verbose = FALSE,
#' debugging = FALSE,
#' ...)
#'
#' @param data.df a data frame containing the data used for learning each node. Binary variables must be declared as factors.
#' @param data.dists a named list giving the distribution for each node in the network, see \sQuote{Details}.
#' @param method should a "Bayes" or "mle" approach be used, see \sQuote{Details}.
#' @param group.var variable name for nodes to be fitted as variable intercept as in a mixed-effects model ("Bayes" and "mle") and gives the column name in \code{data.df} of the grouping variable which must be a factor denoting group membership.
#' @param adj.vars a character vector giving the column names in \code{data.df} for which the network score has to be adjusted for, see \sQuote{Details}.
#' @param cor.vars a character vector giving the column names in \code{data.df} for which a mixed model should be used to adjust for within group correlation or pure adjustment ("bayes" only).
#' @param dag.banned a matrix or a formula statement (see \sQuote{Details} for format) defining which arcs are not permitted - banned - see \sQuote{Details} for format. Note that colnames and rownames must be set, otherwise same row/column names as data.df will be assumed. If set as NULL an empty matrix is assumed.
#' @param dag.retained a matrix or a formula statement (see \sQuote{Details} for format) defining which arcs are must be retained in any model search, see \sQuote{Details} for format. Note that colnames and rownames must be set, otherwise same row/column names as data.df will be assumed. If set as NULL an empty matrix is assumed.
#' @param max.parents a constant or named list giving the maximum number of parents allowed, the list version allows this to vary per node (only for \code{method="bayes"}. A constant can be a single integer, a numeric vector of the length of variables with the same integer for all variable (e.g. \code{c(2,2)}) or a named list with all values being the same (e.g. \code{list("A"=2, "B"=2)}).
#' @param which.nodes a vector giving the column indices of the variables to be included, if ignored all variables are included. This is used to subset \code{data.df}.
#' @param defn.res an optional user-supplied list of child and parent combinations, see \sQuote{Details}.
#' @param centre should the observations in each Gaussian node first be standardized to mean zero and standard deviation one, defaults to TRUE.
#' @param dry.run if TRUE then a list of the child nodes and parent combinations are returned but without estimation of node scores (log marginal likelihoods).
#' @param control a list of control parameters. See \code{\link{build.control}} for the names of the settable control values and their effect.
#' @param verbose if TRUE then provides some additional output.
#' @param debugging if \code{TRUE} and \code{method = 'mle'} this enables to step into the for-loop.
#' @param ... additional arguments passed for optimization.
#'
#' @details The function computes a cache of scores based on possible restrictions (maximum complexity, retained and banned arcs).
#' This function is very similar to \code{\link{fitAbn}} - see that help page for details of the type of models used and in particular \code{data.dists} specification - but rather than fit a single complete DAG \code{buildScoreCache} iterates over all different parent combinations for each node, creating a cache of scores.
#' This cache of score could be used to select the optimal network in other function such as \code{\link{searchHeuristic}} or \code{\link{mostProbable}}.
#' \sQuote{dag.banned} and \sQuote{dag.retained} specify which arcs are forced to be absent or present in the DAG, respectively.
#' If provided as matrix, rows represent child nodes and columns their parents for elements with a value $=1$.
#'
#' Two very different approaches are implemented: a Bayesian and frequentist approaches. They can be selected using the \code{method} argument.
#' ## If \code{method="bayes"}:
#' This function is used to calculate all individual node scores (log marginal likelihoods).
#' Internal code is used by default for numerical estimation in nodes without random effects, and INLA is the default for nodes with random effects.
#' This default behavior can be overridden using \code{control=list(max.mode.error=...)}. The default is \code{max.mode.error=10}, which means that the modes estimated from INLA output must be within 10\% of those estimated using internal code.
#' Otherwise, the internal code is used rather than INLA.
#' To force the use of INLA on all nodes, use max.mode.error=100, which then ignores this check, to force the use of internal code then use \code{max.mode.error=0}. For more detials, see \link{fitAbn}.
#' The variable \code{which.nodes} is to allow the computation to be separated by node, for example, over different CPUs using say \code{R CMD BATCH}.
#' This may useful and indeed likely essential with larger problems or those with random effects.
#' Note that in this case, the results must then be combined back into a list of identical formats to that produced by an individual call to \code{buildScoreCache},
#' comprising of all nodes (in the same order as the columns in \code{data.df}) before sending it to any search routines. Using \code{dry.run} can be useful here.
#' ## If \code{method="mle"}:
#' This function is used to calculate all individual information-theoretic node scores. The possible information-theoretic based network scores computed in \code{buildScoreCache} are the maximum likelihood (mlik, called marginal likelihood in this context as it is computed node wise),
#' the Akaike Information Criteria (aic), the Bayesian Information Criteria (bic) and the Minimum distance Length (mdl). The classical definitions of those metrics are given in Kratzer and Furrer (2018). This function computes a cache that can be fed into a model search algorithm.
#' The numerical routines used here are identical to those in \code{\link{fitAbn}} and see that help page for further details and also the quality assurance section on the \href{https://r-bayesian-networks.org/}{r-bayesian-networks.org} of the \pkg{abn} website for more details.
#'
#' @return A named list of class \code{abnCache}.
#' \describe{
#' \item{\code{children}}{a vector of the child node indexes (from 1) corresponding to the columns in data.df (ignoring any grouping variable)}
#' \item{\code{node.defn}}{a matrix giving the parent combination}
#' \item{\code{mlik}}{log marginal likelihood value for each node combination. If the model cannot be fitted then NA is returned. }
#' \item{\code{error.code}}{if non-zero then either the root finding algorithm (glm nodes) or the maximisation algorithm (glmm nodes) terminated in an unusual way suggesting a possible unreliable result, or else the finite difference hessian estimation produced and error or warning (glmm nodes). NULL if \code{method="mle"}.}
#' \item{\code{error.code.desc}}{a textual description of the \code{error.code}. NULL if \code{method="mle"}}
#' \item{\code{hessian.accuracy}}{An estimate of the error in the final mlik value for each parent combination - this is the absolute difference between two different adaptive finite difference rules where each computes the mlik value. NULL if \code{method="mle"}}
#' \item{\code{data.df}}{a version of the original data (for internal use only in other functions such as \code{\link{mostProbable}}).}
#' \item{\code{data.dists}}{the named list of nodes distributions (for internal use only in other functions such as \code{\link{mostProbable}}).}
#' \item{\code{max.parents}}{the maximum number of parents (for internal use only in other functions such as \code{\link{mostProbable}}).}
#' \item{\code{dag.retained}}{the matrix encoding the retained arcs (for internal use only in other functions such as \code{\link{searchHeuristic}}).}
#' \item{\code{dag.banned}}{the matrix encoding the banned arcs (for internal use only in other functions such as \code{\link{searchHeuristic}}).}
#' \item{\code{aic}}{aic value for each node combination. If the model cannot be fitted then NaN is returned. NULL if \code{method="bayes"}.}
#' \item{\code{bic}}{bic value for each node combination. If the model cannot be fitted then NaN is returned. NULL if \code{method="bayes"}.}
#' \item{\code{mdl}}{mdl value for each node combination. If the model cannot be fitted then NaN is returned. NULL if \code{method="bayes"}.}
#' }
#'
#' @references
#' Kratzer, Gilles, Fraser Lewis, Arianna Comin, Marta Pittavino, and Reinhard Furrer. “Additive Bayesian Network Modeling with the R Package Abn.” Journal of Statistical Software 105 (January 28, 2023): 1–41. https://doi.org/10.18637/jss.v105.i08.
#'
#' Kratzer, G., Lewis, F.I., Comin, A., Pittavino, M., and Furrer, R. (2019). "Additive Bayesian Network Modelling with the R Package abn". arXiv:1911.09006.
#'
#' Kratzer, G., and Furrer, R., (2018). "Information-Theoretic Scoring Rules to Learn Additive Bayesian Network Applied to Epidemiology". arXiv:1808.01126.
#'
#' Lewis, F. I., and McCormick, B. J. J. (2012). "Revealing the complexity of health determinants in resource poor settings". \emph{American Journal Of Epidemiology}. doi:10.1093/aje/KWS183).
#'
#' Further information about \pkg{abn} can be found at: \href{https://r-bayesian-networks.org/}{r-bayesian-networks.org}.
#' @seealso \code{\link{fitAbn}}
#' @family buildScoreCache
#' @export
#'
#' @examples
#' ## Simple example
#' # Generate data
#' N <- 1e6
#' mydists <- list(a="gaussian",
#'                 b="gaussian",
#'                 c="gaussian")
#' a <- rnorm(n = N, mean = 0, sd = 1)
#' b <- 1 + 2*rnorm(n = N, mean = 5, sd = 1)
#' c <- 2 + 1*a + 2*b + rnorm(n = N, mean = 2, sd = 1)
#' mydf <- data.frame("a" = scale(a),
#'                    "b" = scale(b),
#'                    "c" = scale(c))
#'
#' # ABN with MLE
#' mycache.mle <- buildScoreCache(data.df = mydf,
#'                                data.dists = mydists,
#'                                method = "mle",
#'                             max.parents = 2)
#' dag.mle <- mostProbable(score.cache = mycache.mle,
#'                         max.parents = 2)
#' myfit.mle <- fitAbn(object = dag.mle,
#'                     method = "mle",
#'                     max.parents = 2)
#' plot(myfit.mle)
#'
#' \dontrun{
#' # ABN with Bayes
#' if(requireNamespace("INLA", quietly = TRUE)){
#'   # Run only if INLA is available
#'   mycache.bayes <- buildScoreCache(data.df = mydf,
#'                                    data.dists = mydists,
#'                                    method = "bayes",
#'                                    max.parents = 2)
#'   dag.bayes <- mostProbable(score.cache = mycache.bayes,
#'                             max.parents = 2)
#'   myfit.bayes <- fitAbn(object = dag.bayes,
#'                         method = "bayes",
#'                         max.parents = 2)
#'   plot(myfit.bayes)
#' }
#' # Compare MLE and Bayes with lm
#' mymod.lm <- lm(c ~ a + b, data = mydf)
#' summary(mymod.lm)
#'
#' ##################################################################################################
#' ## Example 1 - "mle" vs. "bayes" and the later with using the internal C routine compared to INLA
#' ##################################################################################################
#'
#' # Subset of the build-in dataset, see  ?ex0.dag.data
#' mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")] ## take a subset of cols
#'
#' # setup distribution list for each node
#' mydists <- list(b1="binomial", b2="binomial", g1="gaussian",
#'                 g2="gaussian", b3="binomial", g3="gaussian")
#'
#' # Structural constraints
#' ## ban arc from b2 to b1
#' ## always retain arc from g2 to g1
#' ## parent limits - can be specified for each node separately
#' max.par <- list("b1"=2, "b2"=2, "g1"=2, "g2"=2, "b3"=2, "g3"=2)
#'
#' # now build the cache of pre-computed scores accordingly to the structural constraints
#' res.c <- buildScoreCache(data.df=mydat,
#'                          data.dists=mydists,
#'                          dag.banned= ~b1|b2,
#'                          dag.retained= ~g1|g2,
#'                          max.parents=max.par,
#'                          method="bayes")
#'
#'
#' # repeat but using R-INLA. The mlik's should be virtually identical.
#' # Force using of INLA build.control(max.mode.error=100)
#' if(requireNamespace("INLA", quietly = TRUE)){
#'   res.inla <- buildScoreCache(data.df=mydat,
#'                               data.dists=mydists,
#'                               dag.banned= ~b1|b2, # ban arc from b2 to b1
#'                               dag.retained= ~g1|g2, # always retain arc from g2 to g1
#'                               max.parents=max.par,
#'                               method="bayes",
#'                               control=build.control(max.mode.error=100))
#'
#'   ## comparison - very similar
#'   difference <- res.c$mlik - res.inla$mlik
#'   summary(difference)
#' }
#'
#' # Comparison Bayes with MLE (unconstrained):
#' res.mle <- buildScoreCache(data.df=mydat, data.dists=mydists,
#'                            max.parents=3, method="mle")
#' res.abn <- buildScoreCache(data.df=mydat, data.dists=mydists,
#'                            max.parents=3, method="bayes")
#' # of course different, but same order:
#' plot(-res.mle$bic, res.abn$mlik)
#'
#' #################################################################
#' ## Example 2 - mle with several cores
#' #################################################################
#'
#' ## Many variables, few observations
#' mydat <- ex0.dag.data
#' mydists <- as.list(rep(c("binomial", "gaussian", "poisson"), each=10))
#' names(mydists) <- names(mydat)
#'
#' system.time({
#'   res.mle1 <- buildScoreCache(data.df=mydat,
#'                               data.dists=mydists,
#'                               max.parents=2,
#'                               method="mle",
#'                               control = build.control(method = "mle",
#'                                                       ncores=1))})
#' system.time({
#'   res.mle2 <- buildScoreCache(data.df=mydat,
#'                               data.dists=mydists,
#'                               max.parents=2,
#'                               method="mle",
#'                               control = build.control(method = "mle",
#'                                                       ncores=2))})
#'
#' #################################################################
#' ## Example 3 - grouped data - random effects example e.g. glmm
#' #################################################################
#'
#' ## this data comes with abn see ?ex3.dag.data
#' mydat <- ex3.dag.data[,c("b1","b2","b3","b4","b5","b6","b7",
#'                          "b8","b9","b10","b11","b12","b13", "group")]
#'
#' mydists <- list(b1="binomial", b2="binomial", b3="binomial",
#'                 b4="binomial", b5="binomial", b6="binomial", b7="binomial",
#'                 b8="binomial", b9="binomial", b10="binomial",b11="binomial",
#'                 b12="binomial", b13="binomial" )
#' max.par <- 2
#'
#' ## in this example INLA is used as default since these are glmm nodes
#' ## when running this at node-parent combination 71 the default accuracy check on the
#' ## INLA modes is exceeded (default is a max. of 10 percent difference from
#' ## modes estimated using internal code) and a message is given that internal code
#' ## will be used in place of INLA's results.
#'
#' mycache.bayes <- buildScoreCache(data.df=mydat,
#'                                  data.dists=mydists,
#'                                  group.var="group",
#'                                  method = "bayes",
#'                                  max.parents=max.par)
#' dag.bayes <- mostProbable(score.cache=mycache.bayes)
#' plot(dag.bayes)
#'
#' mycache.mle <- buildScoreCache(data.df=mydat,
#'                                data.dists=mydists,
#'                                group.var="group",
#'                                method = "mle",
#'                                max.parents=max.par)
#' dag.mle <- mostProbable(score.cache=mycache.mle)
#' plot(dag.mle)
#' }
#'
#' @keywords models
#' @concept abn
buildScoreCache <- function(data.df = NULL,
                            data.dists = NULL,
                            method = "bayes",
                            group.var = NULL,
                            adj.vars = NULL,
                            cor.vars = NULL,
                            dag.banned = NULL,
                            dag.retained = NULL,
                            max.parents = NULL,
                            which.nodes = NULL,
                            defn.res = NULL,
                            centre = TRUE,
                            dry.run = FALSE,
                            control = NULL,
                            verbose = FALSE,
                            debugging = FALSE,
                            ...) {

  ## start tests
  # Check verbose
  if(!any(verbose %in% c(TRUE, FALSE))){
    stop(paste("'verbose' is not provided but should be TRUE/FALSE."))
  }

  # Check data
  if(!is.null(data.df)){
    mylist <- check.valid.data(data.df = data.df, data.dists = data.dists, group.var = group.var)
  } else if (is.na(data.df)){
    stop("'data.df' is NA but must be provided.")
  } else {
    stop(paste("'data.df' is not provided."))
  }

  # Check dists
  if(!is.null(data.dists)){
    data.dists <- validate_dists(data.dists = data.dists, returnDists=TRUE)
  } else if (is.na(data.dists)){
    stop("'data.dists' is NA but must be provided.")
  } else {
    stop(paste("'data.dists' is not provided."))
  }

  # Check method
  if(!is.null(method)){
    method <- tolower(method)
    if (is.na(method)){
      stop("'method' is NA but must be provided.")
    } else if(!(method %in% c("bayes", "mle"))){
      stop("`method` is unknown.")
    }
  } else {
    stop(paste("'method' is not provided."))
  }

  # Check group.var
  val_groups <- check.valid.groups(group.var = group.var, data.df = data.df, cor.vars = cor.vars, verbose = verbose)

  # Check cor.vars
  val_corvars <- check.valid.groups(group.var = group.var, data.df = data.df, cor.vars = cor.vars, verbose = verbose)

  if (!is.null(group.var) && !is.null(cor.vars) && !is.null(which.nodes)){
    if (group.var %in% colnames(data.df)[which.nodes]){
      stop("`group.var` should not be among `which.nodes`.")
    } else {
      data.df <- val_groups[["data.df"]]
    }
  } else {
    #TODO: think about this situation...
  }

  # Check adj.vars
  if((!is.null(adj.vars) & !is.null(cor.vars)) & !(is.null(cor.vars[adj.vars]))){stop("cor.vars contains adj.vars, please remove them")}

  # Check dag.banned
  dag.banned <- check.valid.dag(dag = dag.banned, data.df = data.df, is.ban.matrix = TRUE, group.var = group.var)

  # Check dag.retained
  dag.retained <- check.valid.dag(dag = dag.retained, data.df = data.df, is.ban.matrix = FALSE, group.var = group.var) # is.ban.matrix = FALSE to check for cycles

  ###
  # Check max.parents
  # which.nodes: Limits the number of variables. Selects the variables for consideration. Integer vector of column indeces.
  # max.parents: Limits the number of allowed parent variables. single, Integer, list or NULL. Output: integer vector of max.parents corresponding to the column index.
  # defn.res: User specified child-parent combinations. Should not disagree with which.nodes and the max.number of parents.

  # make here max.parents a vector first? This would disentangle the handling of the different input types of max.parents?
  max.parents.orig <- max.parents # store original max.parents for later
  if(is.list(max.parents)){
    # make numeric vector out of list.
    max.parents <- unlist(max.parents, use.names = FALSE)

    # check if all elements were unlisted
    if(length(max.parents) != length(max.parents.orig)){
      stop("Something is wrong in the max.parents list. Any unallowed NA or NULL?")
    }
  }

  # max.parent is NULL or numeric (vector) from here on.
  if (!is.null(defn.res)) {
    if (!is.null(which.nodes)) {
      if (!is.null(max.parents)) {
        # check max.parents on trimmed data set and trimmed distribution
        max.parents <- check.valid.parents(data.df = data.df[which.nodes], max.parents = max.parents, group.var = group.var)
      } else {
        # no max.parents provided
        max.p <- length(data.dists[which.nodes]) # not (n-1) here, to raise warning downstream.

        # check max.p
        max.parents <- check.valid.parents(data.df = data.df[which.nodes], max.parents = max.p, group.var = group.var)
        if(verbose){message("`max.parents` was provided as NULL. Therefore, I assume max.parents equal to the number of possible nodes.")}
      }
      # Check consistency of which.nodes and max.parents
      if(length(which.nodes) < max(unique(max.parents))){
        stop("`which.nodes` selected less nodes than max.parents assumes.")
      }
    } else {
      # no which.nodes provided
      if (!is.null(max.parents)) {
        # check max.parents
        max.parents <- check.valid.parents(data.df = data.df, max.parents = max.parents, group.var = group.var) # this works if max.parents is a single numeric, a list or a numeric vector.
      } else {
        # no max.parents provided
        max.p <- length(data.dists) # not (n-1) here, to raise warning downstream.

        # check max.p
        max.parents <- check.valid.parents(data.df = data.df, max.parents = max.p, group.var = group.var)
        if(verbose){message("`max.parents` was provided as NULL. Therefore, I assume max.parents equal to the number of possible nodes.")}
      }
    }
    # check if defn.res agrees with which.nodes
    which.nodes.defnres <- unique(defn.res$children)
    if (which.nodes.defnres != which.nodes){
      stop("`defn.res` and `which.nodes` provided but can only use either or.")
    }
    # check if defn.res agrees with max.parents
    max.parents.defnres <- max(apply(defn.res[["node.defn"]],1,sum))
    if (max.parents.defnres != max.parents){
      stop("`defn.res` and `max.parents` provided but can only use one of each.")
    }
  } else {
    # no defn.res provided
    if (!is.null(which.nodes)) {
      if (!is.null(max.parents)) {
        # check max.parents on trimmed data set and trimmed distribution
        max.parents <- check.valid.parents(data.df = data.df[which.nodes], max.parents = max.parents, group.var = group.var)
      } else {
        # no max.parents provided
        max.p <- length(data.dists[which.nodes]) # not (n-1) here, to raise warning downstream.

        # check max.p
        max.parents <- check.valid.parents(data.df = data.df[which.nodes], max.parents = max.p, group.var = group.var)
        if(verbose){message("`max.parents` was provided as NULL. Therefore, I assume max.parents equal to the number of possible nodes.")}
      }
      # Check consistency of which.nodes and max.parents
      if(length(which.nodes) < max(unique(max.parents))){
        stop("`which.nodes` selected less nodes than max.parents assumes.")
      }
    } else {
      # no which.nodes provided
      if (!is.null(max.parents)) {
        # check max.parents
        max.parents <- check.valid.parents(data.df = data.df, max.parents = max.parents, group.var = group.var) # this works if max.parents is a single numeric, a list or a numeric vector.
      } else {
        # no max.parents provided
        max.p <- length(data.dists) # not (n-1) here, to raise warning downstream.

        # check max.p
        max.parents <- check.valid.parents(data.df = data.df, max.parents = max.p, group.var = group.var)
        if(verbose){message("`max.parents` was provided as NULL. Therefore, I assume max.parents equal to the number of possible nodes.")}
      }
    }
  } # EOF checking max.parents, which.nodes, defn.res
  # This outputs max.parents as numeric (vector) with a maximal value for max.parents equal to the number of nodes in data.df or data.df[which.nodes] respectively (max.parents <= n).

  # If max.parents = max number of nodes, correct for (n-1)
  if(!is.null(which.nodes)){
    nvars <- length(data.dists[which.nodes])
  } else if(is.null(which.nodes)){
    nvars <- length(data.dists)
  } else {
    stop("`which.nodes` is invalid.")
  }

  if (length(max.parents) == 1){
    # max.parents is single integer
    if(max.parents>=nvars){
      max.parents <- nvars-1
      warning(paste("`max.parents` >= no. of variables. I set it to (no. of variables - 1)=", max.parents))
    }
  } else {
    # max.parents is numeric vector
    counter <- 0
    for (i in 1:length(max.parents)){
      if (max.parents[i] >= nvars) {
        max.parents[i] <- max.parents[i]-1
        counter <- counter+1
      }
    }
    if(counter >0){warning("Some values of `max.parents` >= no. of variables. I set them to (no. of variables - 1): ", max.parents)}
  }

  # Check which.nodes and defn.res
  if (is.null(defn.res)){
    if (is.null(which.nodes)){
      # make which.nodes = all.nodes
      which.nodes <- check.which.valid.nodes(data.df = data.df, which.nodes = which.nodes, group.var = group.var)
    } else {
      # check the provided which.nodes
      which.nodes <- check.which.valid.nodes(data.df = data.df, which.nodes = which.nodes, group.var = group.var)
    }
  } else {
    if (is.null(which.nodes)){
      # If user supplied children and parent combinations, make which.nodes based on them.
      which.nodes <- unique(defn.res$children)
    } else {
      stop("`defn.res` and `which.nodes` provided but can only use either or.")
    }
  }

  # Check centre
  if(!any(centre %in% c(TRUE, FALSE))){
    if(is.null(centre)){
      centre <- TRUE
      warning("`centre` is not provided. I set it to TRUE.")
    } else {
      stop(paste("`centre` is not provided but should be TRUE/FALSE."))
    }
  } else {
    # centre is provided as TRUE/FALSE
    # Check if there are nodes to centre or if there are no gaussian nodes
    if(!any(data.dists == "gaussian")) {
      centre <- FALSE
      if(verbose){message("No gaussian nodes to centre. I set `centre` to FALSE.")}
    } else {
      if(verbose){message("I set `centre` to: ", centre)}
      centre <- centre
    }
  }

  # Check dry.run
  if(!any(dry.run %in% c(TRUE, FALSE))){
    stop(paste("'dry.run' is not provided but should be TRUE/FALSE."))
  }

  # Check control args
  # if any arg matches a possible control arg from build.control(), warn and use it.
  build.control.args <- names(formals(build.control))[-which(names(formals(build.control))=="method")] # allow method to be provided as it has a different meaning here.
  provided.args <- names(match.call()[-1]) # remove the first element of match.call() which is empty.
  if(any(provided.args %in% build.control.args)){
    warning(paste("Some arguments match possible control arguments from `build.control()`.
                  I will use the provided arguments. Please use `control=build.control(...)` instead in the future."))
    ambiguous.args <- provided.args[which(provided.args %in% build.control.args)]
    for (i in 1:length(ambiguous.args)){
      control[[ambiguous.args[i]]] <- match.call()[-1][[ambiguous.args[i]]]
    }
  }
  ctrl <- check.valid.buildControls(control = control, method = method, verbose = verbose)

  if("max.mode.error" %in% names(ctrl)){
    if(ctrl[["max.mode.error"]]==0) {
      force.method <- "C"
    } else if(ctrl[["max.mode.error"]]==100) {
      force.method <- "INLA"
    } else {
      force.method <- "notset"
    }
  } else {
    # no max.mode.error among control list
    force.method <- "notset"
  }

  ## Check consistency of dag.retain and dag.banned
  # check retain does not ask for more arcs to be retained than allowed in max.parents
  # number of parents per node to retain
  max.retain <- apply(dag.retained,1,sum)
  if(length(which( (max.retain>max(unique(max.parents))) == TRUE))>0){stop("'dag.retained' is inconsistent with max.parents!")}

  # check that arcs than are banned are also not retained
  if(length(which(which(as.integer(dag.banned)==1)%in%which(as.integer(dag.retained)==1)==TRUE))>0){stop("'dag.banned' and 'dag.retained' are inconsistent!")}

  # returns ammended data.df and suitable variables
  list.group.var <- check.valid.groups(group.var = group.var, data.df = data.df, cor.vars = cor.vars)
  ## int vect of variables to be treated as grouped indexed from 1
  grouped.vars <- list.group.var$grouped.vars
  ## int vector of group membership ids
  group.ids <- list.group.var$group.ids
  ## this has removed the grouping variable from data.df
  data.df <- list.group.var$data.df

  ## return a list with entries bin, gaus, pois, ntrials and exposure
  # mylist <- check.valid.data(data.df = data.df, data.dists = data.dists, group.var = group.var)

  ## Actual call of the respective function
  if (method == "bayes") {
    out <-
      buildScoreCache.bayes(
        data.df = data.df,
        data.dists = data.dists,
        group.var = group.var,
        cor.vars = cor.vars,
        dag.banned = dag.banned,
        dag.retained = dag.retained,
        max.parents = max.parents,
        which.nodes = which.nodes,
        defn.res = defn.res,
        dry.run = dry.run,
        verbose = verbose,
        centre = centre,
        force.method = force.method,
        mylist = mylist,
        grouped.vars = grouped.vars,
        group.ids = group.ids,
        control = ctrl,
        debugging = debugging
      )
  } else if (method == "mle") {
    out <-
      buildScoreCache.mle(
        data.df = data.df,
        data.dists = data.dists,
        adj.vars = adj.vars,
        cor.vars = cor.vars,
        dag.banned = dag.banned,
        dag.retained = dag.retained,
        max.parents = max.parents,
        which.nodes = which.nodes,
        defn.res = defn.res,
        dry.run = dry.run,
        verbose = verbose,
        debugging = debugging,
        centre = centre,
        force.method = force.method,
        group.var = group.var,
        grouped.vars = grouped.vars,
        group.ids = group.ids,
        control = ctrl
      )
  } else {
    stop("`method` unknown.")
  }
  class(out) <- c("abnCache")
  return(out)

}  #EOF




