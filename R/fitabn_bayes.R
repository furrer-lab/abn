#' @describeIn fitAbn Internal function called by \code{fitAbn}.
#' @param mylist result returned from \code{\link{check.valid.data}}.
#' @param grouped.vars result returned from \code{\link{check.valid.groups}}. Column indexes of all variables which are affected from grouping effect.
#' @param group.ids result returned from \code{\link{check.valid.groups}}. Vector of group allocation for each observation (row) in 'data.df'.
#' @param force.method "notset", "INLA" or "C". This is specified in \code{\link{buildScoreCache}(control=list(max.mode.error=...))}.
#' @family Bayes
#' @importFrom stats spline sd
fitAbn.bayes <- function(dag=NULL,
                         data.df=NULL,
                         data.dists=NULL,
                         group.var=NULL,
                         cor.vars=NULL,
                         centre=TRUE,
                         compute.fixed=FALSE,
                         control = fit.control(method = "bayes"),
                         mylist = NULL,
                         grouped.vars = NULL,
                         group.ids = NULL,
                         force.method = NULL,
                         verbose=FALSE,
                         debugging=FALSE){
  # Multinomial nodes not yet implemented for method "bayes"
  if (any(unlist(data.dists, use.names = F) == "multinomial")) {
    stop("Multinomial nodes are not yet implemented for method 'bayes'. Try with method='mle'.") # Specifically, in file fit_single_node.c, there is no case for multinomial distribution.
  }

  # Set seed for reproducibility
  set.seed(control[["seed"]])

  ## coerce binary factors to become 0/1 integers - the 0 is based on the first entry in levels()
  if (!is.null(mylist$bin)) {
    ## have at least one binary variable
    for (i in mylist$bin) {
      data.df[, i] <- as.numeric(data.df[, i]) - 1
    }
  }

  ## standardize gaussian variables to zero mean and sd=1
  if (centre &&
      !is.null(mylist$gaus)) {
    ## have at least one gaussian variable
    for (i in mylist$gaus) {
      data.df[, i] <- (data.df[, i] - mean(data.df[, i])) / sd(data.df[, i])
    }
  }

  ## coerce all data points to doubles
  ## coerce ALL cols to double equivalents
  for (i in 1:dim(data.df)[2]) {
    data.df[, i] <-
      as.double(data.df[, i])
  }

  ## get distributions in terms of a numeric code
  var.types <- get.var.types(data.dists)

  ########################################################################################
  ## All checking over
  ## get to here we have suitable data/variables and not do the actual fitting
  ########################################################################################
  ### loop through each node and find out what kind of model is to be fitted and then pass to appropriate
  ### separate call for each individual node
  ## setup some storage for entire DAG
  res.list <- list()
  error.code <- NULL
  hessian.accuracy <- NULL
  mymodes <- NULL
  mymargs <- list()
  INLA.marginals <- NULL

  ###########################################################
  ## Iterate over each node in the DAG separately
  ###########################################################
  if (verbose) cat("########################################################\n")
  if (verbose) cat("###### Fitting DAG to data\n")
  if (debugging) {
    out <- list() # Create empty list to store results from for loop
    for (i in 1:dim(dag)[1]){
      out[[i]] <- forLoopContentFitBayes(child = i,
                                         dag = dag,
                                         data.df = data.df,
                                         var.types = var.types,
                                         grouped.vars = grouped.vars,
                                         group.ids = group.ids,
                                         control = control,
                                         INLA.marginals = INLA.marginals,
                                         verbose = verbose,
                                         force.method = force.method,
                                         data.dists = data.dists,
                                         mymodes = mymodes,
                                         error.code = error.code,
                                         hessian.accuracy = hessian.accuracy,
                                         mymargs = mymargs)
    }
  } else {
    ncores <- control[["ncores"]]
    if (ncores > 1) {
      if (verbose) {
        path <- path.expand(paste0(getwd(), "/fitabn_bayes.out"))
        message("Writing cluster output to ", path)
        if (file.exists(path)) {
          file.remove(path)
          message(paste("File exists and will be overwritten:", path))
        }

        cl <- makeCluster(ncores,
                          type = control[["cluster.type"]],
                          rscript_args = "--no-environ", # only available for "FORK"
                          outfile = path)
      } else {
        cl <- makeCluster(ncores,
                          rscript_args = "--no-environ", # only available for "FORK"
                          type = control[["cluster.type"]])
      }
      registerDoParallel(cl)
      out <- foreach(i = 1:dim(dag)[1],
                     .inorder = TRUE,
                     .packages = c("INLA"),
                     .export = "forLoopContentFitBayes",
                     .verbose = verbose) %dopar% {
                       forLoopContentFitBayes(child = i,
                                              dag = dag,
                                              data.df = data.df,
                                              var.types = var.types,
                                              grouped.vars = grouped.vars,
                                              group.ids = group.ids,
                                              control = control,
                                              INLA.marginals = INLA.marginals,
                                              verbose = verbose,
                                              force.method = force.method,
                                              data.dists = data.dists,
                                              mymodes = mymodes,
                                              error.code = error.code,
                                              hessian.accuracy = hessian.accuracy,
                                              mymargs = mymargs)
                     }
      # clean up multi-threading
      stopCluster(cl)
    } else {
      out <- foreach(i = 1:dim(dag)[1],
                     .inorder = TRUE,
                     .packages = c("INLA"),
                     .export = "forLoopContentFitBayes",
                     .verbose = verbose) %do% {
                       forLoopContentFitBayes(child = i,
                                              dag = dag,
                                              data.df = data.df,
                                              var.types = var.types,
                                              grouped.vars = grouped.vars,
                                              group.ids = group.ids,
                                              control = control,
                                              INLA.marginals = INLA.marginals,
                                              verbose = verbose,
                                              force.method = force.method,
                                              data.dists = data.dists,
                                              mymodes = mymodes,
                                              error.code = error.code,
                                              hessian.accuracy = hessian.accuracy,
                                              mymargs = mymargs)
                     }
    }
  }

  #################################################
  ## Further tidy up of results across all nodes
  #################################################
  names(out) <- colnames(dag)

  res.list[["modes"]] <- lapply(out, function(x){return(x$mymodes)})
  res.list[["error.code"]] <- unlist(lapply(out, function(x){return(x$error.code)}))
  res.list[["hessian.accuracy"]] <- unlist(lapply(out, function(x){return(x$hessian.accuracy)}))
  res.list[["error.code.desc"]] <- as.character(res.list[["error.code"]])
  res.list[["error.code.desc"]] <- ifelse(res.list[["error.code.desc"]]==0,"success",res.list[["error.code"]])
  res.list[["error.code.desc"]] <- ifelse(res.list[["error.code.desc"]]==1,"warning: mode results may be unreliable (optimiser terminated unusually)",res.list[["error.code.desc"]])
  res.list[["error.code.desc"]] <- ifelse(res.list[["error.code.desc"]]==2,"error - logscore is NA - model could not be fitted",res.list[["error.code.desc"]])
  res.list[["error.code.desc"]] <- ifelse(res.list[["error.code.desc"]]==4,"warning: fin.diff hessian estimation terminated unusually ",res.list[["error.code.desc"]])

  res.list[["mliknode"]] <- unlist(lapply(out, function(x){return(x$child.mlik)}))
  res.list[["mlik"]] <- sum(res.list[["mliknode"]]) ## overall mlik
  res.list[["mse"]] <- getMSEfromModes(res.list[["modes"]], data.dists)
  res.list[["coef"]] <- modes2coefs(res.list[["modes"]])

  res.list[["used.INLA"]] <- unlist(lapply(out, function(x){return(x$INLA.marginals)})) ## vector - TRUE if INLA used false otherwise

  #####
  ##### Additional part only run if user wants marginal distributions
  #####
  if(compute.fixed){
    if (verbose) cat("Processing marginal distributions for non-INLA nodes...\n")
    ## might have some already computed from INLA
    if(length(which(res.list[["used.INLA"]]==FALSE))==0){
      ## ALL INLA so finished
      res.list[["marginals"]] <- lapply(out, function(x){return(x$mymargs)})
    } else {
      ## At least one C and this creates its own res.list[["marginals"]]
      ## now get rest using C
      max.parents <- max(apply(dag,1,sum)) ## over all nodes - different from above

      res.list <- getmarginals(res.list = res.list, ## rest of arguments as for call to C fitabn
                               data.df = data.df,
                               dag.m = dag,
                               var.types = var.types,
                               max.parents = max.parents,
                               mean = control[["mean"]],
                               prec = control[["prec"]],
                               loggam.shape = control[["loggam.shape"]],
                               loggam.inv.scale = control[["loggam.inv.scale"]],
                               max.iters = control[["max.iters"]],
                               epsabs = control[["epsabs"]],
                               verbose = verbose,
                               error.verbose = control[["error.verbose"]],
                               trace = as.integer(control[["trace"]]),
                               grouped.vars = as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C (changed from earlier fitabn)
                               group.ids = as.integer(group.ids),
                               epsabs.inner = control[["epsabs.inner"]],
                               max.iters.inner = control[["max.iters.inner"]],
                               finite.step.size = control[["finite.step.size"]],
                               hessian.params = control[["hessian.params"]],
                               max.iters.hessian = control[["max.iters.hessian"]],
                               min.pdf = control[["min.pdf"]],
                               marginal.node = control[["marginal.node"]],
                               marginal.param = control[["marginal.param"]],
                               variate.vec = control[["variate.vec"]],
                               n.grid = control[["n.grid"]],
                               INLA.marginals = res.list[["used.INLA"]],
                               iter.max = control[["max.grid.iter"]],
                               max.hessian.error = as.double(control[["max.hessian.error"]]),
                               factor.brent = as.double(control[["factor.brent"]]),
                               maxiters.hessian.brent = as.integer(control[["maxiters.hessian.brent"]]),
                               num.intervals.brent = as.double(control[["num.intervals.brent"]])
      )
    }

    ## at least one INLA node so we need to combine
    ## res.list[["inla.margs"]] with res.list[["marginals"]]
    if(length(which(out$INLA.marginals==TRUE))>0){
      names(mymargs) <- colnames(dag)[which(out$INLA.marginals==TRUE)]
      res.list[["inla.margs"]] <- mymargs
      ## also have res.list[["marginals"]] from getmarginalsC above

      masterlist <- list()
      for(i in rownames(dag)){
        attempt1 <- which(names(res.list[["inla.margs"]])==i)
        attempt2 <- which(names(res.list[["marginals"]])==i)
        if(length(attempt1)>0){
          masterlist[[i]] <- res.list[["inla.margs"]][[attempt1]]
        } else {masterlist[[i]] <- res.list[["marginals"]][[attempt2]] }
      }

      #res.list[["marginals.master"]] <- masterlist
      res.list[["marginals"]] <- masterlist
      res.list[["inla.margs"]] <- NULL

    }

    ## Three final optional operations
    ## 1. evaluate density across an equally spaced grid - this used spline interpolation
    ## 2. standardise the area to unity - this should alread be very close but will do no harm
    ## 3. compute quantiles - this is done after 1. and 2. (assuming they were turned on

    ## 1.
    if(!is.null(control[["n.grid"]])){## evaluated density across an equal grid - use spline interpolation
      res.list[["marginals"]] <- eval.across.grid(res.list[["marginals"]],control[["n.grid"]],control[["marginal.node"]])
    }
    ## 2.
    if(control[["std.area"]]){## want to standardize area under curve to unity - might be slightly adrift as is depending on accuracy of approx's
      if(is.null(control[["n.grid"]])){ stop("must provide n.grid if using std.area!") }
      res.list[["marginals"]] <- std.area.under.grid(res.list[["marginals"]],control[["marginal.node"]])
    }
    ## 3.
    if(!is.null(control[["marginal.quantiles"]])){
      res.list[["marginal.quantiles"]] <- get.quantiles(res.list[["marginals"]],control[["marginal.quantiles"]],control[["marginal.node"]])
    }

  } ## end of compute.fixed


  res.list[["method"]] <- "bayes"
  res.list[["abnDag"]] <- createAbnDag(dag, data.df = data.df, data.dists = data.dists)

  if (verbose) cat("########End of DAG fitting #############################\n")
  if (verbose) cat("########################################################\n")
  return(res.list)
}


#' Regress each node on its parents.#'
#' @param child integer of node to be regressed
#' @param var.types vector of numeric encoding of distribution types. See \code{get.var.types(data.dists)}
#' @param INLA.marginals vector of logicals indicating which nodes are to be fitted using INLA
#' @param mymodes Empty list of modes for each node
#' @param error.code Empty element of error codes for each node
#' @param hessian.accuracy Empty element of hessian accuracies for each node
#' @param mymargs Empty list of marginals for each node
#'
#' @return list of mlik, modes, marginals, error codes, hessian accuracies and a logical if INLA was used for each node.
#' @keywords internal
forLoopContentFitBayes <- function(child = NULL,
                                   dag = NULL,
                                   data.df = NULL,
                                   var.types = NULL,
                                   grouped.vars = NULL,
                                   group.ids = NULL,
                                   control = NULL,
                                   INLA.marginals = NULL,
                                   verbose = NULL,
                                   force.method = NULL,
                                   data.dists = NULL,
                                   mymodes = NULL,
                                   error.code = NULL,
                                   hessian.accuracy = NULL,
                                   mymargs = NULL){
  ###########################################################
  ###########################################################
  if (verbose) cat("###### Processing...Node ",rownames(dag)[child],"\n")
  orig.force.method <- NULL
  used.inla <- TRUE
  ####################################################
  ### First case is the node a GLM
  if( !(child%in%grouped.vars)){## only compute option here is C since fast and INLA slower and less reliable

    if(force.method=="notset" || force.method=="C"){
      if (verbose) cat("Using internal code (Laplace glm)\n")
      r <- try(res.c <- .Call("fit_single_node",
                              data.df,
                              as.integer(child),## childnode
                              as.integer(dag[child,]),## parent combination
                              as.integer(dim(dag)[1]),## number of nodes/variables
                              as.integer(var.types),## type of densities
                              as.integer(sum(dag[child,])),## max.parents
                              as.double(control[["mean"]]),as.double(1/sqrt(control[["prec"]])),as.double(control[["loggam.shape"]]),as.double(1/control[["loggam.inv.scale"]]),
                              as.integer(control[["max.iters"]]),as.double(control[["epsabs"]]),
                              as.integer(verbose),as.integer(control[["error.verbose"]]),as.integer(control[["trace"]]),
                              as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C
                              as.integer(group.ids),## group memberships - note indexed from 1
                              as.double(control[["epsabs.inner"]]),
                              as.integer(control[["max.iters.inner"]]),
                              as.double(control[["finite.step.size"]]),
                              as.double(control[["hessian.params"]]),
                              as.integer(control[["max.iters.hessian"]]),
                              as.integer(0),## Not applicable
                              as.double(control[["max.hessian.error"]]),## Not applicable
                              as.double(control[["factor.brent"]]),## Not applicable
                              as.integer(control[["maxiters.hessian.brent"]]),## Not applicable
                              as.double(control[["num.intervals.brent"]])## Not applicable
                              ,PACKAGE="abn" ))
      if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){
        cat("## !!! Laplace approximation failed at node ", rownames(dag)[child],
            "\n## The additive formulation at this node is perhaps over-parameterized?\n",
            "## Fitting the glm at this node using glm() may provide more information.",sep="")
        stop("")
      }
      used.inla <- FALSE ## flip
    } else {## use INLA for glm
      if(!requireNamespace("INLA", quietly = TRUE)){stop("library INLA is not available!\nR-INLA is available from https://www.r-inla.org/download-install.") }
      mean.intercept <- control[["mean"]] ## use same as for rest of linear terms
      prec.intercept <- control[["prec"]] ## use same as for rest of linear terms
      if (verbose) cat("Using INLA (glm)\n")
      # save(list = c("child", "dag", "data.df", "data.dists", "mean.intercept", "prec.intercept", "control", "verbose"), file = "./tests/testthat/testdata/calc.node.inla.glm_1.RData")
      res.inla <- calc.node.inla.glm(child, dag, data.df, data.dists,
                                     rep(1,dim(data.df)[1]),## ntrials
                                     rep(1,dim(data.df)[1]),## exposure
                                     TRUE, mean.intercept, prec.intercept, control[["mean"]], control[["prec"]],
                                     control[["loggam.shape"]],control[["loggam.inv.scale"]],
                                     verbose.loc = verbose,
                                     nthreads = control[["ncores"]])


      if(is.logical(res.inla)){if (verbose) cat("INLA failed....so reverting to internal code\n")
        r <- try(res.c <- .Call("fit_single_node",
                                data.df,
                                as.integer(child),## childnode
                                as.integer(dag[child,]),## parent combination
                                as.integer(dim(dag)[1]),## number of nodes/variables
                                as.integer(var.types),## type of densities
                                as.integer(sum(dag[child,])),## max.parents
                                as.double(control[["mean"]]),as.double(1/sqrt(control[["prec"]])),as.double(control[["loggam.shape"]]),as.double(1/control[["loggam.inv.scale"]]),
                                as.integer(control[["max.iters"]]),as.double(control[["epsabs"]]),
                                as.integer(verbose),as.integer(control[["error.verbose"]]),as.integer(control[["trace"]]),
                                as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C
                                as.integer(group.ids),## group memberships - note indexed from 1
                                as.double(control[["epsabs.inner"]]),
                                as.integer(control[["max.iters.inner"]]),
                                as.double(control[["finite.step.size"]]),
                                as.double(control[["hessian.params"]]),
                                as.integer(control[["max.iters.hessian"]]),
                                as.integer(0),## Not applicable
                                as.double(control[["max.hessian.error"]]),## Not applicable
                                as.double(control[["factor.brent"]]),## Not applicable
                                as.integer(control[["maxiters.hessian.brent"]]),## Not applicable
                                as.double(control[["num.intervals.brent"]])## Not applicable
                                ,PACKAGE="abn" ))
        if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){
          cat("## !!! Laplace approximation failed at node ", rownames(dag)[child],
              "\n## The additive formulation at this node is perhaps over-parameterized?\n",
              "## Fitting the glm at this node using glm() may provide more information.",sep="")
          stop("")
        }

        used.inla <- FALSE ## flip
      } ## INLA failed

    } ## use INLA
    ## End of GLM node
    ###########################################################

  } else if (child%in%grouped.vars) {
    ###########################################################
    ## Have a GLMM node
    ###########################################################
    ## have a glmm, so two options, INLA or C

    if(force.method=="notset" || force.method=="INLA"){##
      if(!requireNamespace("INLA", quietly = TRUE)){stop("library INLA is not available!\nR-INLA is available from https://www.r-inla.org/download-install.") }
      mean.intercept <- control[["mean"]] ## use same as for rest of linear terms
      prec.intercept <- control[["prec"]] ## use same as for rest of linear terms
      res.inla <- calc.node.inla.glmm(child, dag,
                                      data.frame(data.df,group=group.ids),
                                      data.dists,
                                      rep(1,dim(data.df)[1]),## ntrials
                                      rep(1,dim(data.df)[1]),## exposure
                                      TRUE,## always compute marginals - since only way to check results
                                      mean.intercept, prec.intercept, control[["mean"]], control[["prec"]],control[["loggam.shape"]],control[["loggam.inv.scale"]],
                                      verbose.loc = verbose,
                                      nthreads = control[["ncores"]])
      #return(res.inla)  ## to return RAW INLA object
      ## CHECK FOR INLA CRASH
      if(is.logical(res.inla)){
        if (verbose) cat("INLA failed....so reverting to internal code\n")
        orig.force.method <- force.method ## save original
        force.method="C"  ## Easiest way is just to force C for this node
      } else {
        res.inla.modes <- getModeVector(list.fixed=res.inla$marginals.fixed,list.hyper=res.inla$marginals.hyperpar)
      }

    }

    if (verbose) cat("fit a glmm at node ",rownames(dag)[child],"using C\n")
    if(force.method=="notset"){
      r <- try(res.c <- .Call("fit_single_node",
                              data.df,
                              as.integer(child),## childnode
                              as.integer(dag[child,]),## parent combination
                              as.integer(dim(dag)[1]),## number of nodes/variables
                              as.integer(var.types),## type of densities
                              as.integer(sum(dag[child,])),## max.parents
                              as.double(control[["mean"]]),as.double(1/sqrt(control[["prec"]])),as.double(control[["loggam.shape"]]),as.double(1/control[["loggam.inv.scale"]]),
                              as.integer(control[["max.iters"]]),as.double(control[["epsabs"]]),
                              as.integer(verbose),as.integer(control[["error.verbose"]]),as.integer(control[["trace"]]),
                              as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C
                              as.integer(group.ids),## group memberships - note indexed from 1
                              as.double(control[["epsabs.inner"]]),
                              as.integer(control[["max.iters.inner"]]),
                              as.double(control[["finite.step.size"]]),
                              as.double(control[["hessian.params"]]),
                              as.integer(control[["max.iters.hessian"]]),
                              as.integer(1),## turn on ModesONLY
                              as.double(control[["max.hessian.error"]]),
                              as.double(control[["factor.brent"]]),
                              as.integer(control[["maxiters.hessian.brent"]]),
                              as.double(control[["num.intervals.brent"]])
                              ,PACKAGE="abn")) #print(res)
      if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){
        cat("## !!! Laplace approximation failed at node ", rownames(dag)[child],
            "\n## The additive formulation at this node is perhaps over-parameterized?\n",
            "## Fitting the glmm at this node using glmer() in lme4 may provide more information",sep="")
        stop("")
      }

      res.c.modes <- res.c[[1]][-c(1:3)] ## remove mlik - this is first entry, and error code and hessian accuracy
      res.c.modes <- res.c.modes[which(res.c.modes!=.Machine$double.xmax)] ## this discards all "empty" parameters
      ## get difference in modes proportion relative to C
      diff.in.modes <- (res.inla.modes-res.c.modes)/res.c.modes
      error.modes <- max(abs(diff.in.modes))
    }

    if( force.method=="C" || (force.method=="notset" && error.modes>(control[["max.mode.error"]]/100))){ ## INLA might be unreliable to use C (slower)

      if(force.method=="notset"){
        if (verbose) cat("Using internal code (Laplace glmm)\n=>max. abs. difference (in %) with INLA is ")
        if (verbose) cat(formatC(100*error.modes,format="f",digits=1)," and exceeds tolerance\n")
      } else {
        if (verbose) cat("Using internal code (Laplace glmm)\n")
      }

      r <- try(res.c <- .Call("fit_single_node",
                              data.df,
                              as.integer(child),## childnode
                              as.integer(dag[child,]),## parent combination
                              as.integer(dim(dag)[1]),## number of nodes/variables
                              as.integer(var.types),## type of densities
                              as.integer(sum(dag[child,])),## max.parents
                              as.double(control[["mean"]]),as.double(1/sqrt(control[["prec"]])),as.double(control[["loggam.shape"]]),as.double(1/control[["loggam.inv.scale"]]),
                              as.integer(control[["max.iters"]]),as.double(control[["epsabs"]]),
                              as.integer(verbose),as.integer(control[["error.verbose"]]),as.integer(control[["trace"]]),
                              as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C
                              as.integer(group.ids),## group memberships - note indexed from 1
                              as.double(control[["epsabs.inner"]]),
                              as.integer(control[["max.iters.inner"]]),
                              as.double(control[["finite.step.size"]]),
                              as.double(control[["hessian.params"]]),
                              as.integer(control[["max.iters.hessian"]]),
                              as.integer(0),## turn off ModesONLY
                              as.double(control[["max.hessian.error"]]),
                              as.double(control[["factor.brent"]]),
                              as.integer(control[["maxiters.hessian.brent"]]),
                              as.double(control[["num.intervals.brent"]])
                              ,PACKAGE="abn")
      )
      if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){
        cat("## !!! Laplace approximation failed at node ", rownames(dag)[child],
            "\n## The additive formulation at this node is perhaps over-parameterized?\n",
            "## Fitting the glmm at this node using glmer() in lme4 may provide more information.",sep="")
        stop("")
      }

      used.inla <- FALSE ## flip

    } else {
      if (verbose) cat("Using INLA (glmm)\n")
    }## end of if inla bad

    ###########################################################
    ## End of GLMM node
    ###########################################################
  } else {          ## end of if GLMM
    stop("There is an issue in 'child' and/or 'grouped.vars'.")
  }

  ###########################################################
  ## End of all external computations
  ###########################################################
  ## computation for current node is all done so sort out the
  ## output into nicer form and give labels
  ###########################################################
  child.name <- colnames(dag)[child]

  if(used.inla==FALSE){
    ## organize output from C
    child.mlik <- res.c[[1]][1]
    mymodes <- res.c[[1]][-c(1:3)] ## remove mlik - this is first entry, and error code and hessian accuracy
    mymodes <- mymodes[which(mymodes!=.Machine$double.xmax)] ## this discards all "empty" parameters
    error.code <- res.c[[1]][2]
    hessian.accuracy <- res.c[[1]][3]
    INLA.marginals <- c(INLA.marginals,FALSE)
  } else {
    ## organize output from INLA
    child.mlik <- res.inla$mlik[2] ## [2] is for Gaussian rather than Integrated estimate
    mymodes <- getModeVector(list.fixed=res.inla$marginals.fixed,list.hyper=res.inla$marginals.hyperpar)
    mymargs <- getMargsINLA(list.fixed=res.inla$marginals.fixed,list.hyper=res.inla$marginals.hyperpar)
    error.code <- NA ## not available from INLA
    hessian.accuracy <- NA ## not available from INLA
    INLA.marginals <- c(INLA.marginals,TRUE)
  }

  nom <- colnames(dag)[which(dag[child,]==1)]
  if(var.types[child]=="1"){nom <- c("(Intercept)",nom,"group.precision") }## binomial : just some naming for use later
  if(var.types[child]=="2" && !(child%in%grouped.vars)){nom <- c("(Intercept)",nom,"precision","precision") } ## gaus and not grouped
  if(var.types[child]=="2" && (child%in%grouped.vars)){nom <- c("(Intercept)",nom,"group.precision","precision") }
  if(var.types[child]=="3"){nom <- c("(Intercept)",nom,"group.precision") } ## pois}
  mynom <- NULL

  for(j in 1:length(mymodes)){
    mynom <- c(mynom,paste(colnames(dag)[child],nom[j],sep="|"))
  }
  names(mymodes) <- mynom
  if(used.inla==TRUE){
    names(mymargs) <- mynom }

  if(!is.null(orig.force.method)){force.method <- orig.force.method } ## reset force.method after INLA crash
  ############################################################
  ## Finished with current node
  ############################################################
  return(
    list(
      child.mlik = child.mlik,
      mymodes = mymodes,
      mymargs = mymargs,
      error.code = error.code,
      hessian.accuracy = hessian.accuracy,
      INLA.marginals = INLA.marginals
    )
  )
}

#' function to extract the mode from INLA output
#'
#' @param list.fixed list of matrices of two cols x, y
#' @param list.hyper list of hyperparameters
#'
#' @return vector
#' @export
#' @keywords internal
getModeVector <- function(list.fixed,list.hyper){

  ## listfixed is a list of matrices of two cols x, y
  modes <- NULL
  for(i in 1:length(list.fixed)){
    mymarg <- list.fixed[[i]] ## matrix 2 cols
    mymarg <- spline(mymarg)
    modes <- c(modes,mymarg$x[which(mymarg$y==max(mymarg$y))]) ## find x value which corresponds to the max y value
  }

  ## now for hyperparam if there is one
  if(is.list(list.hyper)){## might be no hyperparam e.g. binomial or poisson glmm
    if(length(list.hyper)==1){
      mymarg <- list.hyper[[1]] ## matrix 2 cols
      mymarg <- spline(mymarg)
      modes <- c(modes,mymarg$x[which(mymarg$y==max(mymarg$y))])
    }## find x value which corresponds to the max y value

    if(length(list.hyper)==2){ ## gaussian glm
      mymarg <- list.hyper[[2]] ## matrix 2 cols - group level precision
      mymarg <- spline(mymarg)
      modes <- c(modes,mymarg$x[which(mymarg$y==max(mymarg$y))])
      mymarg <- list.hyper[[1]] ## matrix 2 cols - residual level precision
      mymarg <- spline(mymarg)
      modes <- c(modes,mymarg$x[which(mymarg$y==max(mymarg$y))])
      ### [[2]] then [[1]] to keep same order a C code
    }
  }

  return(modes)

}

#' function to extract marginals from INLA output
#'
#' @param list.fixed list of matrices of two cols x, y
#' @param list.hyper list of hyperparameters
#'
#' @return vector
#' @export
#' @keywords internal
getMargsINLA <- function(list.fixed,list.hyper){

  ## listfixed is a list of matrices of two cols x, y
  margs <- list()
  for(i in 1:length(list.fixed)){
    margs[[i]] <- list.fixed[[i]] ## matrix 2 cols
  }
  ## now for hyperparam if there is one
  if(is.list(list.hyper)){## might be no hyperparam e.g. binomial or poisson glmm
    if(length(list.hyper)==1){## poisson or binomial glmm
      margs[[i+1]] <- list.hyper[[1]] }## matrix 2 cols
    if(length(list.hyper)==2){## gaussian glmm two entries
      margs[[i+1]] <- list.hyper[[2]]  ## this is the group level precision
      margs[[i+2]] <- list.hyper[[1]] }## this is the residual precision
    ## note this is [[2]] then [[1]] to keep same order a C code

  }

  return(margs)

}

#' function to get marginal across an equal grid
#'
#' @param mylist list of matrices of two cols x, y
#' @param n.grid grid size
#' @param single NULL or TRUE if only a single node and parameter
#'
#' @return list
#' @export
#' @keywords internal
eval.across.grid <- function(mylist,n.grid,single){

  if(is.null(single)){
    for(i in 1:length(mylist)){## for each node
      q.inner.list <- mylist[[i]] ##copy
      for(j in 1:length(q.inner.list)){## for each parameter
        mymat <- q.inner.list[[j]] ## copy - this is a matrix
        q.mat <- matrix(data=rep(NA,2*n.grid),ncol=2) ## create new matrix
        colnames(q.mat) <- c("x","f(x)")
        interp <- spline(mymat,n=n.grid) ## interpolate over equal spaced grid of n points
        q.mat[,1] <- interp$x
        q.mat[,2] <- interp$y
        q.inner.list[[j]] <- q.mat ## overwrite
      }
      mylist[[i]] <- q.inner.list
    }
  } else {## have only a single node and parameter
    mymat <- mylist[[1]]
    q.mat <- matrix(data=rep(NA,2*n.grid),ncol=2) ## create new matrix
    colnames(q.mat) <- c("x","f(x)")
    interp <- spline(mymat,n=n.grid) ## interpolate over equal spaced grid of n points
    q.mat[,1] <- interp$x
    q.mat[,2] <- interp$y
    mylist[[1]] <- q.mat ## overwrite
  }
  return(mylist)
}

#' Standard Area Under the Marginal
#'
#' function to get std. are under marginal to exactly unity.
#' It should be very close to unity but in some cases due to numerical accuracy
#' differences (since each point is a separate estimate) this might be a little adrift
#' turn this option off to see how reliable the original estimation is
#'
#' @param mylist list of matrices of two cols x, y
#' @param single NULL or TRUE if only a single node and parameter
#'
#' @return list
#' @export
#' @keywords internal
std.area.under.grid <- function(mylist,single){

  if(is.null(single)){
    for(i in 1:length(mylist)){## for each node
      q.inner.list <- mylist[[i]] ##copy
      for(j in 1:length(q.inner.list)){## for each parameter
        mymat <- q.inner.list[[j]] ## copy - this is a matrix
        cur.area <- (mymat[2,1]-mymat[1,1])*sum(mymat[,2]) ## delta.x used - equal sized grid
        mymat[,2] <- mymat[,2]/cur.area ## std to ~ 1.0
        q.inner.list[[j]] <- mymat ## overwrite
      }
      mylist[[i]] <- q.inner.list
    }
  } else {
    mymat <- mylist[[1]] ## copy - this is a matrix
    cur.area <- (mymat[2,1]-mymat[1,1])*sum(mymat[,2]) ## delta.x used - equal sized grid
    mymat[,2] <- mymat[,2]/cur.area ## std to ~ 1.0
    mylist[[1]] <- mymat ## overwrite

  }

  return(mylist)
}

#' function to extract quantiles from INLA output
#'
#' function to get to extract quantiles
#'
#' @param mylist list of matrices of two cols x, y
#' @param quantiles vector with the desired quantiles
#' @param single NULL or TRUE if only a single node and parameter
#'
#' @return list
#' @export
#' @keywords internal
get.quantiles <- function(mylist,quantiles, single){

  if(is.null(single)){
    for(i in 1:length(mylist)){
      q.inner.list <- mylist[[i]] ##copy
      for(j in 1:length(q.inner.list)){
        mymat <- q.inner.list[[j]] ## copy - this is a matrix
        q.mat <- matrix(data=rep(NA,2*length(quantiles)),ncol=2) ## create new matrix
        colnames(q.mat) <- c("P(X<=x)","x")
        q.mat[,1] <- quantiles
        q.mat <- get.ind.quantiles(q.mat,mymat) ## the actual quantile arithmetic
        q.inner.list[[j]] <- q.mat ## overwrite
      }
      mylist[[i]] <- q.inner.list
    }
  } else {
    mymat <- mylist[[1]] ## copy - this is a matrix
    q.mat <- matrix(data=rep(NA,2*length(quantiles)),ncol=2) ## create new matrix
    colnames(q.mat) <- c("P(X<=x)","x")
    q.mat[,1] <- quantiles
    q.mat <- get.ind.quantiles(q.mat,mymat) ## the actual quantile arithmetic
    mylist[[1]] <- q.mat ## overwrite

  }
  return(mylist)
}

#' @describeIn get.quantiles helper function for get.quantiles
#' @param outmat matrix where the first col has the desired quantiles. We want to estimate this and out in into the second col
#' @param inmat is the actual x,f(x) matrix
#' @return matrix
get.ind.quantiles <- function(outmat,inmat){
  ##outmat is a matrix where the first col has the desired quantiles
  ## we want to estimate this and out in into the second col
  ## inmat is the actual x,f(x) matrix
  qs.vec <- outmat[,1]
  x <- inmat[,1]
  fx <- inmat[,2]
  cum.den <- cumsum(fx) ## cumulative of f(x) values
  sum.den <- sum(fx) ## total sum of f(x)
  cum.f <- cum.den/sum.den ## cumulative density function
  row <- 1
  for(qs in qs.vec){## for each quantile
    ## find row in inmat which has cumulative f(x)>q.s
    outmat[row,2] <- x[which(cum.f>qs)[1]] ## find first row to exceed quantile value q.s and get x value
    row <- row+1
  }

  class(outmat) <- c("abnFit")

  return(outmat)
}

#' Extract Standard Deviations from all Gaussian Nodes
#'
#' @param modes list of modes.
#' @param dists list of distributions.
#'
#' @return named numeric vector. Names correspond to node name. Value to standard deviations.
getMSEfromModes <- function(modes, dists){
  modes_gaus <- unlist(unname(modes[unname(which(dists == "gaussian"))]), use.names = TRUE)
  if (!is.null(modes_gaus)){
    # if there is at least one gaussian node, extract the stddev
    taus <- modes_gaus[stri_detect_fixed(str = names(modes_gaus), pattern = "precision")]
    taus_names <- as.character(stringi::stri_split_fixed(str = names(taus), pattern = "|precision", omit_empty = TRUE, simplify = TRUE))
    names(taus) <- taus_names
    mses <- 1/taus # mse stores stddevs
    return(mses)
  } else {
    # if there is no gaussian, return NULL
    return(NULL)
  }
}

#' Convert modes to fitAbn.mle$coefs structure
#'
#' @param modes list of modes.
#'
#' @return list of matrix arrays.
modes2coefs <- function(modes){
  newmodes <- modes
  for (child in names(newmodes)){
    childnames <- names(newmodes[[child]])
    childnames <- stringi::stri_replace_all_fixed(str = childnames, pattern = "|(Intercept)", replacement = "|intercept")
    for (childcoef in seq(1:length(childnames))) {
      # iterate through all coefficients from current node
      if (stringi::stri_detect_fixed(str = childnames[childcoef], pattern = "intercept", negate = TRUE)){
        # Rename all but child|intercept
        childnames[childcoef] <- stringi::stri_replace_all_fixed(str = childnames[childcoef], pattern = paste0(child, "|"), replacement = "")
      }
    }
    names(newmodes[[child]]) <- childnames

    # Remove Precision items from gaussians
    for (childcoefi in seq(1:length(childnames))) {
      if (stringi::stri_detect_fixed(str = childnames[childcoefi], pattern = "precision", negate = FALSE)){
        newmodes[[child]] <- newmodes[[child]][-childcoefi]
      }
    }

    # Convert to list of matrix arrays
    newmodes[[child]] <- as.array(t(as.matrix(newmodes[[child]])))
  }
  return(newmodes)
}
