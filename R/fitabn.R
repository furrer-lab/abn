#' Control the iterations in \code{\link{fitAbn}}
#'
#' Allow the user to set restrictions in the \code{\link{fitAbn}} for both the Bayesian and the MLE approach.
#' Control function similar to \code{\link{build.control}}.
#'
#' @inheritParams build.control
#' @param min.pdf the value of the posterior density function below which we stop the estimation only used when computing marginals, see details.
#' @param std.area logical, should the area under the estimated posterior density be standardized to exactly one, useful for error checking.
#' @param marginal.quantiles vector giving quantiles at which to compute the posterior marginal distribution at.
#' @param max.grid.iter gives number of grid points to estimate posterior density at when not explicitly specifying a grid used to avoid excessively long computation.
#' @param marginal.node used in conjunction with \code{marginal.param} to allow bespoke estimate of a marginal density over a specific grid. value from 1 to the number of nodes.
#' @param marginal.param used in conjunction with \code{marginal.node}. value of 1 is for intercept, see modes entry in results for the appropriate number.
#' @param variate.vec a vector containing the places to evaluate the posterior marginal density, must be supplied if \code{marginal.node} is not null.
#'
#' @inherit build.control details
#' @inherit build.control return
#' @export fit.control
#' @examples
#' ctrlmle <- abn::fit.control(method = "mle",
#'                        max.irls = 100,
#'                        tol = 10^-11,
#'                        tolPwrss = 1e-7,
#'                        xtol_abs = 1e-6,
#'                        ftol_abs = 1e-6,
#'                        epsilon = 1e-6,
#'                        ncores = 2,
#'                        seed = 9062019L)
#' ctrlbayes <- abn::fit.control(method = "bayes",
#'                          mean = 0,
#'                          prec = 0.001,
#'                          loggam.shape = 1,
#'                          loggam.inv.scale = 5e-05,
#'                          max.mode.error = 10,
#'                          max.iters = 100,
#'                          epsabs = 1e-07,
#'                          error.verbose = FALSE,
#'                          epsabs.inner = 1e-06,
#'                          max.iters.inner = 100,
#'                          finite.step.size = 1e-07,
#'                          hessian.params = c(1e-04, 0.01),
#'                          max.iters.hessian = 10,
#'                          max.hessian.error = 1e-04,
#'                          factor.brent = 100,
#'                          maxiters.hessian.brent = 10,
#'                          num.intervals.brent = 100,
#'                          min.pdf = 0.001,
#'                          n.grid = 100,
#'                          std.area = TRUE,
#'                          marginal.quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
#'                          max.grid.iter = 1000,
#'                          marginal.node = NULL,
#'                          marginal.param = NULL,
#'                          variate.vec = NULL,
#'                          ncores = 1,
#'                          seed = 9062019L)
#'
#' @seealso \code{\link{build.control}}.
#' @family fitAbn
fit.control <-
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
            max.hessian.error = 1e-04,
            factor.brent = 100,
            maxiters.hessian.brent = 10,
            num.intervals.brent = 100,
            min.pdf = 0.001,
            n.grid = 250,
            std.area = TRUE,
            marginal.quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
            max.grid.iter = 1000,
            marginal.node = NULL,
            marginal.param = NULL,
            variate.vec = NULL,
            ncores = 1,
            max.irls = 100,
            tol = 1e-11,
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
          min.pdf = min.pdf,
          n.grid = n.grid,
          std.area = std.area,
          marginal.quantiles = marginal.quantiles,
          max.grid.iter = max.grid.iter,
          marginal.node = marginal.node,
          marginal.param = marginal.param,
          variate.vec = variate.vec,
          ncores = ncores,
          seed = seed
        )
    } else if (method == "mle") {
      ctrl <-
        list(
          ncores = ncores,
          max.iters = max.iters,
          max.irls = max.irls,
          tol = tol,
          tolPwrss = tolPwrss,
          check.rankX = check.rankX,
          check.scaleX = check.scaleX,
          check.conv.grad = check.conv.grad,
          check.conv.singular = check.conv.singular,
          check.conv.hess = check.conv.hess,
          xtol_abs = xtol_abs,
          ftol_abs = ftol_abs,
          trace.mblogit = trace.mblogit,
          catcov.mblogit = catcov.mblogit,
          epsilon = epsilon,
          seed = seed
        )
    } else {
      stop("`method` unknown.")
    }
    return(ctrl)
  }


#' Fit an additive Bayesian network model
#'
#' @description Fits an additive Bayesian network to observed data and is equivalent to Bayesian or information-theoretic multi-dimensional regression modeling.
#' Two numerical options are available in the Bayesian settings, standard Laplace approximation or else an integrated nested Laplace approximation provided via a call to the R INLA library (see \href{https://www.r-inla.org/}{r-inla.org} - this is not hosted on CRAN).
#' @usage
#' fitAbn(object = NULL,
#'        dag = NULL,
#'        data.df = NULL,
#'        data.dists = NULL,
#'        method = NULL,
#'        group.var = NULL,
#'        adj.vars = NULL,
#'        cor.vars = NULL,
#'        centre = TRUE,
#'        compute.fixed = FALSE,
#'        control = NULL,
#'        verbose = FALSE,
#'        debugging = FALSE,
#'        ...)
#'
#' @param object an object of class \code{abnLearned} produced by \code{\link{mostProbable}}, \code{\link{searchHeuristic}} or \code{\link{searchHillClimber}}.
#' @param dag a matrix or a formula statement (see details) defining the network structure, a directed acyclic graph (DAG), see details for format. Note that column names and row names must be set up.
#' @param data.df a data frame containing the data used for learning the network, binary variables must be declared as factors, and no missing values all allowed in any variable.
#' @param data.dists a named list giving the distribution for each node in the network, see details.
#' @param method if \code{NULL}, takes method of \code{object}, otherwise \code{"bayes"} or \code{"mle"} for the method to be used, see details.
#' @param group.var only applicable for mixed models and gives the column name in \code{data.df} of the grouping variable (which must be a factor denoting group membership).
#' @param adj.vars a character vector giving the column names in \code{data.df} for which the network score has to be adjusted for, see details.
#' @param cor.vars a character vector giving the column names in data.df for which a mixed model should be used (\code{method = 'bayes'} only).
#' @param centre should the observations in each Gaussian node first be standardised to mean zero and standard deviation one.
#' @param compute.fixed a logical flag, set to \code{TRUE} for computation of marginal posterior distributions, see details.
#' @param control a list of control parameters. See \code{\link{fit.control}} for the names of the settable control values and their effect.
#' @param verbose if \code{TRUE} then provides some additional output, in particular the code used to call INLA, if applicable.
#' @param debugging if \code{TRUE} and \code{method = 'mle'} this enables to step into the for-loop.
#' @param \dots additional arguments passed for optimization.
#'
#' @details
#' ## If \code{method="Bayes"}:
#' The procedure \code{fitAbn} fits an additive Bayesian network model to data where each node (variable - a column in data.df) can be either: presence/absence (Bernoulli); continuous (Gaussian); or an unbounded count (Poisson). Multinomial distributions are only supported with \code{method = "mle"} (see below).
#' The model comprises of a set of conditionally independent generalized linear regressions with or without random effects.
#' Internal code is used by default for numerical estimation in nodes without random effects, and INLA is the default for nodes with random effects.
#' This default behavior can be overridden using \code{control=list(max.mode.error=...)}. The default is \code{max.mode.error=10}, which means that the modes estimated from INLA output must be within 10\% of those estimated using internal code.
#' Otherwise, the internal code is used rather than INLA.
#' To force the use of INLA on all nodes, use \code{max.mode.error=100}, which then ignores this check, to force the use of internal code then use \code{max.mode.error=0}.
#' For the numerical reliability and perform of \pkg{abn} see \url{http://r-bayesian-networks.org}.
#' Generally speaking, INLA can be swift and accurate, but in several cases, it can perform very poorly and so some care is required (which is why there is an internal check on the modes).
#' Binary variables must be declared as factors with two levels, and the argument \code{data.dists} must be a list with named arguments, one for each of the variables in \code{data.df} (except a grouping variable - if present!), where each entry is either "poisson","binomial", "multinomial" or "gaussian", see examples below.
#' The "poisson" and "binomial" distributions use log and logit link functions, respectively.
#' Note that "binomial" here actually means only binary, one Bernoulli trial per row in \code{data.df}.
#'
#' If the data are grouped into correlated blocks - wherein a standard regression context a mixed model might be used - then a network comprising of one or more nodes where a generalized linear mixed model is used (but limited to only a single random effect).
#' This is achieved by specifying parameters \code{group.var} and \code{cor.vars}.
#' Where the former defines the group membership variable, which should be a factor indicating which observations belong to the same grouping.
#' The parameter \code{cor.vars} is a character vector that contains the names of the nodes for which a mixed model should be used. This is not yet implemented with \code{method = 'mle'}.
#' For example, in some problems, it may be appropriate for all variables (except \code{group.var}) in data.df to be parametrized as a mixed model while in others it may only be a single variable for which grouping adjustment is required (as the remainder of variables are covariates measured at group level).
#'
#' In the network structure definition, \code{dag}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.
#' The \code{dag} can be provided using a formula statement (similar to GLM).
#' A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}.
#' The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3.
#' The parents names must match those given in \code{data.df}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{data.df}.
#'
#' If \code{compute.fixed=TRUE} then the marginal posterior distributions for all parameters are computed.
#' Note the current algorithm used to determine the evaluation grid is rather crude and may need to be manually refined using \code{variate.vec} (one parameter at a time) for publication-quality density estimates.
#' Note that a manual grid can only be used with internal code and not INLA (which uses its own grid).
#' The end points are defined as where the value of the marginal density drops below a given threshold \code{pdf.min}.
#' When estimating the log marginal likelihood in models with random effects (using internal code rather than INLA), an attempt is made to minimize the error by comparing the estimates given between a 3pt and 5pt rule when estimating the Hessian in the Laplace approximation.
#' The modes used in each case are identical. The first derivatives are computed using gsl's adaptive finite difference function, and this is embedding inside the standard 3pt and 5pt rules for the second derivatives.
#' In all cases, a central difference approximation is tried first with a forward difference being a fall back (as the precision parameters are strictly positive).
#' The error is minimized through choosing an optimal step size using gsl's Nelder-Mead optimization, and if this fails, (e.g., is larger than \code{max.hessian.error}) then the Brent-Dekker root bracketing method is used as a fallback.
#' If the error cannot be reduced to below \code{max.hessian.error}, then the step size, which gave the lowest error during the searches (across potentially many different initial bracket choices), is used for the final Hessian evaluations in the Laplace approximation.
#'
#' ## If \code{method="mle"}:
#' The procedure \code{fitAbn} with the argument \code{method= "mle"} fits an additive Bayesian network model to data where each node (variable - a column in data.df) can be either: presence/absence (Bernoulli); continuous (Gaussian); an unbounded count (Poisson); or a discrete variable (Multinomial).
#' The model comprises of a set of conditionally independent generalized linear regressions with or without adjustment.
#' Binary and discrete variables must be declared as factors and the argument \code{data.dists} must be a list with named arguments, one for each of the variables in \code{data.df}, where each entry is either "poisson","binomial", "multinomial" or "gaussian", see examples below.
#' The "poisson" and "binomial" distributions use log and logit link functions, respectively.
#' Note that "binomial" here actually means only binary, one Bernoulli trial per row in data.df.
#'
#' If the data are grouped into correlated blocks - wherein a standard regression context a mixed-effect model might be used - then a network comprising of one or more nodes where a generalized linear mixed model is used (but limited to only a single random intercept).
#' This is achieved by specifying parameter \code{group.var} (\code{cor.vars} as with \code{method = "bayes"} is not yet implemented with \code{method = "mle"}).
#' The parameter \code{group.var} defines the group membership variable, which should be a factor indicating which observations belong to the same grouping.
#' This corresponds to \code{"1|group.var"} in the formula notation of e.g. \pkg{lme4}.
#'
#' In the context of \code{fitAbn} adjustment means that irrespective to the adjacency matrix the adjustment variable set (\code{adj.vars}) will be add as covariate to every node defined by \code{cor.vars}.
#'
#' In the network structure definition, \code{dag}, each row represents a node in the network, and the columns in each row define the parents for that particular node, see the example below for the specific format.
#' The \code{dag} can be provided using a formula statement (similar to GLM). A typical formula is \code{ ~ node1|parent1:parent2 + node2:node3|parent3}.
#' The formula statement have to start with \code{~}. In this example, node1 has two parents (parent1 and parent2). node2 and node3 have the same parent3.
#' The parents names have to exactly match those given in \code{data.df}. \code{:} is the separator between either children or parents, \code{|} separates children (left side) and parents (right side), \code{+} separates terms, \code{.} replaces all the variables in \code{data.df}.
#'
#' The Information-theoretic based network scores used in \code{fitAbn} with argument \code{method="mle"} are the maximum likelihood (mlik, called marginal likelihood in this context as it is computed node wise), the Akaike Information Criteria (aic), the Bayesian Information Criteria (bic) and the Minimum distance Length (mdl). The classical definitions of those metrics are given in Kratzer and Furrer (2018).
#'
#' The numerical routine is based on an iterative scheme to estimate the regression coefficients. The Iterative Reweighed Least Square (IRLS) programmed using Rcpp/RcppArmadrillo. One hard coded feature of \code{fitAbn} with argument \code{method="mle"} is a conditional use of a bias reduced binomial regression when a classical Generalized Linear Model (GLM) fails to estimate the maximum likelihood of the given model accurately. Additionally, a QR decomposition is performed to check for rank deficiency. If the model is rank deficient and the BR GLM fails to estimate it, then predictors are sequentially removed. This feature aims at better estimating network scores when data sparsity is present.
#'
#' A special care should be taken when interpreting or even displaying p-values computed with \code{fitAbn}. Indeed, the full model is already selected using goodness of fit metrics based on the (same) full dataset.
#'
#' The \code{control} argument is a list with separate arguments for the Bayesian and MLE implementation. See  \code{\link{fit.control}} for details.
#'
#' @return An object of class \code{abnFit}. A named list. One entry for each of the variables in \code{data.df} (excluding the grouping variable, if present) which contains an estimate of the log marginal likelihood for that individual node. An entry "mlik" which is the total log marginal likelihood for the full ABN model. A vector of \code{error.codes} - non-zero if a numerical error or warning occurred, and a vector error.code.desc giving a text description of the error. A list \code{modes}, which contains all the mode estimates for each parameter at each node. A vector called Hessian accuracy, which is the estimated local error in the log marginal likelihood for each node.  If \code{compute.fixed=TRUE} then a list entry called \code{marginals} which contains a named entry for every parameter in the ABN and each entry in this list is a two-column matrix where the first column is the value of the marginal parameter, say x, and the second column is the respective density value, pdf(x). Also, a list called \code{marginal.quantiles} is produced, giving the quantiles for each marginal posterior distribution.
#'
#' @references
#' Kratzer, G., Lewis, F.I., Comin, A., Pittavino, M. and Furrer, R. (2019). "Additive Bayesian Network Modelling with the R Package abn". arXiv preprint arXiv:1911.09006.
#'
#' Kratzer, G., and Furrer, R., 2018. Information-Theoretic Scoring Rules to Learn Additive Bayesian Network Applied to Epidemiology. Preprint; Arxiv: stat.ML/1808.01126.
#'
#' Lewis, F. I., and McCormick, B. J. J. (2012). Revealing the complexity of health determinants in resource poor settings. \emph{American Journal Of Epidemiology}. DOI:10.1093/aje/KWS183.
#'
#' Further information about \pkg{abn} can be found at: \href{http://r-bayesian-networks.org/}{r-bayesian-networks.org}.
#'
#' @author Fraser Iain Lewis and Gilles Kratzer.
#'
#' @seealso \code{\link{buildScoreCache}}
#' @family fitAbn
#' @export
#' @examples
#' ## Built-in dataset with a subset of cols
#' mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")]
#'
#' ## setup distribution list for each node
#' mydists <- list(b1="binomial", b2="binomial", b3="binomial", g1="gaussian",
#'                b4="binomial", p2="poisson", p4="poisson")
#'
#' ## Null model - all independent variables
#' mydag.empty <- matrix(0, nrow=7, ncol=7)
#' colnames(mydag.empty) <- rownames(mydag.empty) <- names(mydat)
#'
#' ## Now fit the model to calculate its goodness-of-fit
#' myres <- fitAbn(dag=mydag.empty, data.df=mydat, data.dists=mydists)
#'
#' ## Log-marginal likelihood goodness-of-fit for complete DAG
#' print(myres$mlik)
#'
#' ## fit  using the formula statement
#' #  including the creation of the graph of the DAG via Rgraphviz
#' myres <- fitAbn(dag=~b1|b2+b2|p4:g1+g1|p2+b3|g1+b4|b1+p4|g1,
#'                 data.df=mydat, data.dists=mydists)
#' print(myres$mlik) ## a much weaker fit than full independence DAG
#'
#' plotAbn(dag=myres$abnDag$dag, data.dists=mydists, fitted.values=myres$modes)
#'
#' ## Or equivalentelly using the formula statement, with plotting
#' ## Now repeat but include some dependencies first
#' mydag <- mydag.empty
#' mydag["b1","b2"] <- 1 # b1<-b2 and so on
#' mydag["b2","p4"] <- mydag["b2","g1"] <- mydag["g1","p2"] <- 1
#' mydag["b3","g1"] <- mydag["b4","b1"] <- mydag["p4","g1"] <- 1
#' myresAlt <- fitAbn(dag=mydag, data.df=mydat, data.dists=mydists)
#' plot(myresAlt)
#'
#' ## -----------------------------------------------------------------------------------
#' ## This function contains an MLE implementation accessible through a method parameter
#' ## use built-in simulated data set
#' ## -----------------------------------------------------------------------------------
#'
#' myres.mle <- fitAbn(dag=~b1|b2+b2|p4+g1+g1|p2+b3|g1+b4|b1+p4|g1,
#'                     data.df=mydat, data.dists=mydists, method="mle")
#'
#' ## Print the output for mle first then for Bayes:
#' print(myres.mle)
#' print(myres)
#'
#'
#' \dontrun{
#'   ## A simple plot of some posterior densities the algorithm which chooses
#'   ## density points is very simple any may be rather sparse so also recompute
#'   ## the density over an equally spaced grid of 50 points between the two
#'   ## end points which had at f=min.pdf
#'   ## max.mode.error=0 foces to use the internal c code
#'   myres.c <- fitAbn(dag=mydag, data.df=mydat, data.dists=mydists,
#'                     compute.fixed=TRUE,
#'                     control=list(max.mode.error=0))
#'
#'   print(names(myres.c$marginals)) ## gives all the different parameter names
#'
#'   ## Repeat but use INLA for the numerics using max.mode.error=100
#'   ## as using internal code is the default here rather than INLA
#'   myres.inla <- fitAbn(dag=mydag, data.df=mydat, data.dists=mydists,
#'                        compute.fixed=TRUE,
#'                        control=list(max.mode.error=100))
#'
#'   ## Plot posterior densities
#'   par(mfrow=c(2,2), mai=c(.7,.7,.2,.1))
#'   plot(myres.c$marginals$b1[["b1|(Intercept)"]], type="l", xlab="b1|(Intercept)")
#'   lines(myres.inla$marginals$b1[["b1|(Intercept)"]], col="blue")
#'   plot(myres.c$marginals$b2[["b2|p4"]], type="l", xlab="b2|p4")
#'   lines(myres.inla$marginals$b2[["b2|p4"]], col="blue")
#'   plot(myres.c$marginals$g1[["g1|precision"]], type="l", xlab="g1|precision")
#'     lines(myres.inla$marginals$g1[["g1|precision"]], col="blue")
#'     plot(myres.c$marginals$b4[["b4|b1"]], type="l", xlab="b4|b1")
#'     lines(myres.inla$marginals$b4[["b4|b1"]], col="blue")
#'
#'  ## An elementary mixed model example using built-in data specify DAG,
#'  ## only two variables using a subset of variables from ex3.dag.data
#'  ## both variables are assumed to need (separate) adjustment for the
#'  ## group variable, i.e., a binomial GLMM at each node
#'
#'
#'  mydists <- list(b1="binomial",  b2="binomial")
#'
#'  ## Compute marginal likelihood - use internal code via max.mode.error=0
#'  ## as using INLA is the default here.
#'  ## Model where b1 <- b2
#'  myres.c <- fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2,14)], data.dists=mydists,
#'                    group.var="group", cor.vars=c("b1","b2"),
#'                    control=list(max.mode.error=0))
#'  print(myres.c) ## show all the output
#'
#'  ## compare mode for node b1 with glmer(), lme4::glmer is automatically attached.
#'
#'  ## Now for marginals - INLA is strongly preferable for estimating marginals for nodes
#'  ## with random effects as it is far faster, but may not be reliable
#'  ## see http://r-bayesian-networks.org
#'
#'  ## INLA's estimates of the marginals, using high n.grid=500
#'  ## as this makes the plots smoother - see below.
#'  ## myres.inla <- fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2,14)],
#'  #                      data.dists=mydists,
#'  #                      group.var="group", cor.vars=c("b1", "b2"),
#'  #                      compute.fixed=TRUE, n.grid=500,
#'  #                      control=list(max.mode.error=100, max.hessian.error=10E-02))
#'
#'  ## this is NOT recommended - marginal density estimation using fitAbn in mixed models
#'  ## is really just for diagnostic purposes, better to use fitAbn.inla() here
#'  ## but here goes...be patient
#'  # myres.c <- fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2,14)], data.dists=mydists,
#'  #                   group.var="group", cor.vars=c("b1", "b2"), compute.fixed=TRUE,
#'  #                   control=list(max.mode.error=0, max.hessian.error=10E-02))
#'
#'  ## compare marginals between internal and INLA.
#'  # par(mfrow=c(2,3))
#'  ## 5 parameters - two intercepts, one slope, two group level precisions
#'  # plot(myres.inla$marginals$b1[[1]], type="l", col="blue")
#'  # lines(myres.c$marginals$b1[[1]], col="brown", lwd=2)
#'  # plot(myres.inla$marginals$b1[[2]], type="l", col="blue")
#'  # lines(myres.c$marginals$b1[[2]], col="brown", lwd=2)
#'  ## the precision of group-level random effects
#'  # plot(myres.inla$marginals$b1[[3]],type="l", col="blue", xlim=c(0,2))
#'  # lines(myres.c$marginals$b1[[3]],col="brown",lwd=2)
#'  # plot(myres.inla$marginals$b2[[1]],type="l", col="blue")
#'  # lines(myres.c$marginals$b2[[1]],col="brown",lwd=2)
#'  # plot(myres.inla$marginals$b2[[1]], type="l", col="blue")
#'  # lines(myres.c$marginals$b2[[1]], col="brown", lwd=2)
#'  ## the precision of group-level random effects
#'  # plot(myres.inla$marginals$b2[[2]], type="l", col="blue", xlim=c(0,2))
#'  # lines(myres.c$marginals$b2[[2]], col="brown", lwd=2)
#'
#'  ### these are very similar although not exactly identical
#'
#'  ## use internal code but only to compute a single parameter over a specified grid
#'  ## This can be necessary if the simple auto grid finding functions does a poor job
#'
#'  #myres.c <- fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2,14)], data.dists=mydists,
#'  #                  group.var="group", cor.vars=c("b1", "b2"),
#'  #                  centre=FALSE, compute.fixed=TRUE,
#'  #                  control=list(marginal.node=1, marginal.param=3,## precision term in node 1
#'  #                  variate.vec=seq(0.05, 1.5, len=25), max.hessian.error=10E-02))
#'
#'  #par(mfrow=c(1,2))
#'  #plot(myres.c$marginals[[1]], type="l", col="blue")## still fairly sparse
#'  ## An easy way is to use spline to fill in the density without recomputing other
#'  ## points provided the original grid is not too sparse.
#'  #plot(spline(myres.c$marginals[[1]], n=100), type="b", col="brown")
#' }
#'
#' @keywords models
#' @concept abn
fitAbn <- function(object = NULL,
                   dag = NULL,
                   data.df = NULL,
                   data.dists = NULL,
                   method = NULL,
                   group.var = NULL,
                   adj.vars = NULL,
                   cor.vars = NULL,
                   centre = TRUE,
                   compute.fixed = FALSE,
                   control=NULL,
                   verbose = FALSE,
                   debugging = FALSE,
                   ...) {

  # Check verbose
  if(!any(verbose %in% c(TRUE, FALSE))){
    stop(paste("'verbose' is not provided but should be TRUE/FALSE."))
  }

  # Check data
  if(!is.null(data.df) && !is.null(object)){
    # data and object (with data) provided together but don't know on which data to rely
    stop("`data.df` and `object` provided but can only accept one of them.")
  } else if (!is.null(data.df) && is.null(object)){
    # only data.df provided
    mylist <- check.valid.data(data.df = data.df, data.dists = data.dists, group.var = group.var)
  } else if (is.null(data.df) && !is.null(object)) {
    # only object provided
    NA # This will be catched downstream when other arguments have passed tests.
  } else if (is.null(data.df) || is.na(data.df)) {
    stop("'data.df' is missing but must be provided in this situation.")
  } else {
    # no data.df and no object argument provided
    stop("`data.df` and `object` are missing but one of them is required.")
  }

  # Check dists
  if(!is.null(data.dists) && !is.null(object)){
    # data and object (with data) provided together but don't know on which data to rely
    stop("`data.dists` and `object` provided but can only accept one of them.")
  } else if (!is.null(data.dists) && is.null(object)){
    # only data.dists provided
    data.dists <- validate_dists(data.dists = data.dists, returnDists=TRUE)
  } else if (is.null(data.dists) && !is.null(object)) {
    # only object provided
    NA # This will be catched downstream when other arguments have passed tests.
  } else if (is.na(data.dists)) {
    stop("'data.dists' is NA but must be provided.")
  } else {
    # no data.dists and no object argument provided
    stop("`data.dists` and `object` are missing but one of them is required.")
  }

  # check object$group.var matches group.var argument here
  if(!is.null(group.var) && !is.null(object)){
    if (inherits(x=object, what="abnLearned") || inherits(x=object, what="abnCache")){
      # object is provided in correct form.
      if (group.var == object$score.cache$group.var) {
        if (verbose){message("Ok. 'group.var' in score.cache object equals 'group.var' argument provided to fitAbn().")}
      } else if (group.var != object$score.cache$group.var) {
        stop(paste("Fitting (", group.var, ") and learned (", object$score.cache$group.var, ") 'group.var' argument differ."))
      } else {
        stop("Checking coherence of 'group.var' argument failed with unknown error. I should never end up here.")
      }
    } else {
      warning("Unrecognised score cache 'object'. May lead to errors later on.")
    }
  } else {
    # more extensive tests are performed further downstream.
    if (verbose){message("Ok. No 'group.var' in score.cache object or no 'group.var' argument provided to fitAbn(). I assume there is no grouping.")}
  }


  # Resolve historic issue with dag and object arguments
  if (!is.null(dag) && !is.null(object)){
    # dag and object argument provided together
    stop("`dag` and `object` provided but can only accept one of them.")
  } else if (!is.null(dag) && is.null(object)){
    # only dag provided
    if(inherits(x = dag, what = "abnLearned")){
      # if dag is of class "abnLearned" make it the object. This is the old way.
      object <- dag
      message("Best practice with abn > 2.0 requires to pass 'dag' as 'object' parameter.")

      dag <- object$dag
      object <- object$score.cache
      data.df <- object$data.df
      data.dists <- object$data.dists
      group.var <- object$group.var
      cor.vars <- object$cor.vars
      # adj.vars <- object$adj.vars
      fitmethod <- object$method
      mylist <- object$mylist
    } else {
      # Only formula DAGs are accepted otherwise
      validdag <- check.valid.dag(dag = dag, data.df = data.df, group.var = group.var)
      fitmethod <- NULL
      if (!is.matrix(validdag)) {
        stop("Unknown `dag` argument")
      } else {
        dag <- validdag
      }
    }
  } else if (is.null(dag) && !is.null(object)){
    # only object provided
    if (inherits(x=object, what="abnLearned")) {
      # if object is of class "abnLearned" extract its stuff
      dag <- object$dag
      object <- object$score.cache
      data.df <- object$data.df
      data.dists <- object$data.dists
      group.var <- object$group.var
      cor.vars <- object$cor.vars
      # adj.vars <- object$adj.vars
      fitmethod <- object$method
      mylist <- object$mylist
    } else {
      stop("Unknown type of `object`. Must be of class 'abnLearned'.")
    }
  } else {
    # no dag and no object argument provided
    stop("`dag` and `object` are missing but one of them is required.")
  }

  # Check method
  if(!is.null(method)){
    method <- tolower(method)
    if (is.na(method)){
      stop("'method' is NA but must be provided.")
    } else if(!(method %in% c("bayes", "mle"))){
      stop("`method` argument is unknown. Must be NULL='bayes' or 'mle'.")
    }
  } else {
    # if no method provided, assume "bayes"
    method <- "bayes"
  }

  if(!is.null(fitmethod)){
    if(!(fitmethod %in% c("bayes", "mle"))){
      stop("Fitting `method` is unknown.")
    } else if (fitmethod != method){
      warning(paste("Fitting and learned methods differ. Continuing with fitting method: ", method))
    } else {
      # if method and fitmethod are ok, overwrite method argument with method from object$method
      method <- fitmethod
    }
  } else {
    if(!is.null(method)){
      if(verbose){message("using method from fitAbn() argument.")}
    } else {
      stop(paste("'method' is not provided."))
    }
  }

  # Check group.var
  if (!is.null(group.var)){
    # we have grouping
    if (inherits(x=object, what="abnLearned") || inherits(x=object, what="abnCache")) {
      # grouping was checked before. Extract only the variables that were not already above.
      grouped.vars <- object$grouped.vars
      group.ids <- object$group.ids
    } else {
      # check grouping
      val_groups <- check.valid.groups(group.var = group.var, data.df = data.df, cor.vars = cor.vars, verbose = verbose)
      val_corvars <- check.valid.groups(group.var = group.var, data.df = data.df, cor.vars = cor.vars, verbose = verbose)

      grouped.vars <- val_groups$grouped.vars ## int vect of variables to be treated as grouped indexed from 1
      group.ids <- val_groups$group.ids ## int vector of group membership ids
      data.df <- val_groups$data.df ## this has removed the grouping variable from data.df
      dag <- check.valid.dag(dag = dag[val_groups$grouped.vars, val_groups$grouped.vars], data.df = val_groups$data.df) ## remove the grouping variable from the DAG.
    }
  } else if (is.null(group.var)) {
    # no group.var provided to fitAbn(). Check if group.var is provided through object.
    if (inherits(x=object, what="abnLearned") || inherits(x=object, what="abnCache")) {
      # grouping was checked before. Extract only the variables that were not already above.
      grouped.vars <- object$grouped.vars
      group.ids <- object$group.ids
    } else {
      # we have really no grouping
      grouped.vars <- 0L
      group.ids <- 0L
    }
  } else {
    stop("I do not recognize the 'group.var' argument provided to fitAbn().")
  }

  # Check adj.vars
  if((!is.null(adj.vars) & !is.null(cor.vars)) & !(is.null(cor.vars[adj.vars]))){stop("cor.vars contains adj.vars, please remove them")}

  # Check centre
  if(!any(centre %in% c(TRUE, FALSE))){
    stop("'centre' should be either TRUE or FALSE but it was provided as: ", centre)
    centre <- centre
  } else {
    if(verbose){message("Assessing validity of argument 'centre'. ... OK.")}
  }

  # Check control args
  ctrl <- check.valid.fitControls(control = control, method = method, verbose = verbose)

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

  ## Actual call of the respective function
  if (method == "bayes") {
    out <- fitAbn.bayes(dag,
                        data.df = data.df,
                        data.dists = data.dists,
                        group.var = group.var,
                        cor.vars = cor.vars,
                        compute.fixed = compute.fixed,
                        control = ctrl,
                        mylist = mylist,
                        grouped.vars = grouped.vars,
                        group.ids = group.ids,
                        force.method = force.method)
  } else if (method == "mle") {
    out <- fitAbn.mle(dag,
                      data.df = data.df,
                      data.dists = data.dists,
                      group.var = group.var,
                      grouped.vars = grouped.vars,
                      group.ids = group.ids,
                      adj.vars = adj.vars,
                      cor.vars = cor.vars,
                      verbose = verbose,
                      centre = centre,
                      control = ctrl,
                      debugging = debugging)
  } else {
    stop("`method` unknown.")
  }
  class(out) <- c("abnFit")
  return(out)
}  #EOF
