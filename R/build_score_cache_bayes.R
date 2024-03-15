#' \code{buildScoreCache.mle} and \code{buildScoreCache.bayes} are internal functions called by \code{buildScoreCache}.
#'
#' @describeIn buildScoreCache Fit a given DAG to data with method="bayes".
#' @param force.method "notset", "INLA" or "C". This is specified in \code{\link{buildScoreCache}(control=list(max.mode.error=...))}.
#' @param mylist result returned from \code{\link{check.valid.data}}.
#' @param grouped.vars result returned from \code{\link{check.valid.groups}}.
#' @param group.ids result returned from \code{\link{check.valid.groups}}.
#' @importFrom stats sd
#' @family Bayes
buildScoreCache.bayes <-
  function(data.df = NULL,
           data.dists = NULL,
           group.var = NULL,
           cor.vars = NULL,
           dag.banned = NULL,
           dag.retained = NULL,
           max.parents = NULL,
           which.nodes = NULL,
           defn.res = NULL,
           dry.run = FALSE,
           centre = TRUE,
           force.method = NULL,
           mylist = NULL,
           grouped.vars = NULL,
           group.ids = NULL,
           control = build.control(method = "bayes"),
           verbose = FALSE,
           debugging = FALSE) {
    # Multinomial nodes not yet implemented for method "bayes"
    if (any(unlist(data.dists, use.names = F) == "multinomial")) {
      stop("Multinomial nodes are not yet implemented for method 'bayes'. Try with method='mle'.") # Specifically, in file fit_single_node.c, there is no case for multinomial distribution.
    }

    set.seed(control[["seed"]])

    data.df.original <- data.df
    n <- length(data.dists)

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
      data.df[, i] <- as.double(data.df[, i])
    }

    ## get distributions in terms of a numeric code
    var.types <- get.var.types(data.dists)

    ########################################################################################
    ## All checking over
    ## get to here we have suitable data/variables and now generate all the parent combinations
    ########################################################################################
    ## down to here we have all the data correct and now call C buildnodecache() to create all the node definitions.
    if (is.null(defn.res)) {
      ## pass to C the number (number_of_nodes,banned_arc_as_vector, retain_arcs_as_vector, max_parents_as_vector
      res <- .Call("buildcachematrix",
                   dim(dag.banned)[1],
                   as.integer(dag.banned),
                   as.integer(dag.retained),
                   as.integer(max.parents),
                   which.nodes,
                   PACKAGE = "abn"
      )
      # extract the results from the list
      children <- res[[1]] # this is a vector of node numbers for each node and its parent combinations
      node.defn <- matrix(data = res[[2]],
                          byrow = TRUE,
                          ncol = dim(data.df)[2],
                          dimnames = list(NULL, names(data.df))) # this is a matrix of 0/1 for each node and its parent combinations.
      rm(res)
    } else {
      ## some check since user has supplied defn.res
      if (!is.list(defn.res)) {
        stop("'defn.res' must be a list")
      }
      if (length(defn.res) != 2) {
        stop("'defn.res' must have two entries")
      }
      if (!(is.vector(defn.res[[1]]) &&
            is.matrix(defn.res[[2]]))) {
        stop("'defn.res' is wrong format")
      }
      if (!(max(defn.res[[2]]) == 1 &&
            min(defn.res[[2]]) == 0)) {
        stop("'defn.res' is wrong format - must only be 0,1 in node definitions")
      }
      # extract the results from the list
      children <- defn.res[["children"]]
      node.defn <- defn.res[["node.defn"]]
    }

    # dry run - don't do any computation just return the node definitions
    if(dry.run){
      if (verbose) cat("No computation - returning only the node combinations\n")
      return(list(children=children,node.defn=node.defn))
    }

    dag.m <- matrix(rep(NA,dim(data.df)[2]^2),ncol=dim(data.df)[2])
    colnames(dag.m) <- rownames(dag.m) <- names(data.df)
    ###########################################################
    ## Iterate over each node in the DAG separately
    ### loop through each node and find out what kind of model is to be fitted and then pass to appropriate
    ### separate call for each individual node-parentcombination
    ###########################################################
    out <- list()
    rows <- length(children) # number of child-parent combinations

    if (verbose) cat("Start estimation loop.")
    if (debugging){
      # store each child-parent combination in a separate row
      res <- matrix(NA,
                    nrow = rows,
                    ncol = 5,
                    dimnames = list(NULL, c("childParentCombNo", "mlik", "error.code", "hessian.accuracy", "used.INLA")))
      for (i in 1:rows) {
        res[i, ] <- forLoopContentBayes(row.no = i,
                                        children = children,
                                        node.defn = node.defn,
                                        dag.m = dag.m,
                                        force.method = force.method,
                                        data.df = data.df,
                                        data.dists = data.dists,
                                        var.types = var.types,
                                        control = control,
                                        grouped.vars = grouped.vars,
                                        group.ids = group.ids,
                                        verbose = verbose)
      }
    } else {
      # no debugging

      # Prepare multithreading
      ncores <- control[["ncores"]]

      if (ncores > 1){
        if (verbose){
          path <- path.expand(paste0(getwd(), "/build_score_cache_bayes.out"))
          message(paste("Writing output to", path))
          if(file.exists(path)){
            file.remove(path)
            message(paste("File exists and will be overwritten:", path))
          }

          cl <- makeCluster(ncores,
                            type = control[["cluster.type"]],
                            rscript_args = "--no-environ", # only available for "FORK"
                            outfile = path)

        } else {
          # no redirection of output
          cl <- makeCluster(ncores,
                            type = control[["cluster.type"]],
                            rscript_args = "--no-environ") # only available for "FORK"
        }

        registerDoParallel(cl)

        res <- foreach(i = 1:rows,
                       .combine = "rbind",
                       .packages = c("INLA"),
                       .export = "forLoopContentBayes",
                       .verbose = verbose) %dopar% {
          forLoopContentBayes(row.no = i,
                              children = children,
                              node.defn = node.defn,
                              dag.m = dag.m,
                              force.method = force.method,
                              data.df = data.df,
                              data.dists = data.dists,
                              var.types = var.types,
                              control = control,
                              grouped.vars = grouped.vars,
                              group.ids = group.ids,
                              verbose = verbose)
                       }
        # clean up multi-threading
        stopCluster(cl)
      } else {
        res <- foreach(i = 1:rows,
                       .combine = "rbind",
                       .packages = c("INLA"),
                       .export = "forLoopContentBayes",
                       .verbose = verbose) %do% {
                         forLoopContentBayes(row.no = i,
                                             children = children,
                                             node.defn = node.defn,
                                             dag.m = dag.m,
                                             force.method = force.method,
                                             data.df = data.df,
                                             data.dists = data.dists,
                                             var.types = var.types,
                                             control = control,
                                             grouped.vars = grouped.vars,
                                             group.ids = group.ids,
                                             verbose = verbose)
                       }
      }
    }
    if (verbose) cat("################# End of cache building ################### \n")
    # sort the results by child-parent combination number
    res <- res[order(res[, 1]), ] # important to match later on with the node.defn and children vectors (especially if processed in parallel)

    out$children <- children
    out$node.defn <- node.defn

    # store the results in a list
    # out$childParentCombNo <- res[, 1] # no need to store this. Just for reference.
    out$mlik <- as.numeric(res[, 2])
    out$error.code <- as.numeric(res[, 3])
    out$hessian.accuracy <- as.numeric(res[, 4])
    out$used.INLA <- as.logical(res[, 5])

    # add error code descriptions
    out$error.code.desc <- as.character(out$error.code)
    out$error.code.desc[out$error.code.desc == 0] <- "success"
    out$error.code.desc[out$error.code.desc == 1] <- "warning: mode results may be unreliable (optimiser terminated unusually)"
    out$error.code.desc[out$error.code.desc == 2] <- "error - logscore is NA. Model could not be fitted"
    # there is no error code 3!
    out$error.code.desc[out$error.code.desc == 4] <- "warning: fin.diff hessian estimation terminated unusually"

    # Finalise the results
    out$data.df <- data.df.original
    out$data.dists <- data.dists # used in searchHeuristic() and mostProbable()
    out$max.parents <- max.parents
    out$dag.retained <- dag.retained
    out$dag.banned <- dag.banned
    out$group.var <- group.var
    out$group.ids <- group.ids
    out$group.vars <- grouped.vars
    out$cor.vars <- cor.vars
    out$mylist <- mylist
    out$method <- "bayes"

    return(out)
  }

#' From each child-parent(s) combination, regress each child on its parents in buildScoreCache.bayes()
#' @describeIn buildScoreCache Internal function called by \code{buildScoreCache.bayes()}.
#' @param row.no The row number of the child-parent combination to be processed.
#' @param children vector of child node integers.
#' @param node.defn child-parent combination table.
#' @param dag.m Empty adjacency matrix.
#' @param var.types vector of numeric encoding of distribution types. See \code{get.var.types(data.dists)}
#'
#' @return Named vector of results from one child-parent combination subject to the \code{row.no}.
#' The names are:
#' \describe{
#' \item{childParentCombNo}{The row number of the child-parent combination in the \code{node.defn} table.
#' This must be the same as the row number in \code{node.defn}:
#' careful if \code{buildScoreCache.bayes()} is run in parallel!}
#' \item{mlik}{The marginal log-likelihood of the child-parent combination.}
#' \item{error.code}{The error code returned by \code{inla()}.}
#' \item{hessian.accuracy}{The accuracy of the Hessian matrix returned by \code{inla()}.}
#' \item{used.INLA}{A logical value indicating whether \code{inla()} was used to fit the model.}
#' }
#' @export
#' @keywords internal
forLoopContentBayes <- function(row.no = NULL, # i
                                children = NULL,
                                node.defn = NULL,
                                dag.m = NULL,
                                force.method = NULL,
                                data.df = NULL,
                                data.dists = NULL,
                                var.types = NULL,
                                control = NULL,
                                grouped.vars = NULL,
                                group.ids = NULL,
                                verbose = FALSE) {
  child <- children[row.no]

  FAILED <- FALSE
  ## to catch any crashes...
  ###########################################################
  if (verbose) cat("###### Processing...",row.no," of ", nrow(node.defn) ,"\n")
  dag.m[,] <- 0
  ## reset to empty
  dag.m[child,] <- node.defn[row.no,]

  ## set parent combination
  orig.force.method <- NULL
  used.inla <- TRUE

  ####################################################
  ### First case is the node a GLM
  ####################################################
  if( !(child%in%grouped.vars)){
    ## only compute option here is C since fast and INLA slower and less reliable
    if(force.method=="notset" || force.method=="C"){
      if (verbose) cat("Using internal code (Laplace glm)\n")
      r <- try(res.c <- .Call("fit_single_node",
                              data.df,
                              as.integer(child), ## childnode
                              as.integer(dag.m[child,]), ## parent combination
                              as.integer(dim(dag.m)[1]),  ## number of nodes/variables
                              as.integer(var.types), ## type of densities
                              as.integer(sum(dag.m[child,])),  ## max.parents
                              as.double(control[["mean"]]),as.double(1/sqrt(control[["prec"]])),as.double(control[["loggam.shape"]]),as.double(1/control[["loggam.inv.scale"]]),
                              as.integer(control[["max.iters"]]),as.double(control[["epsabs"]]),
                              as.integer(verbose),as.integer(control[["error.verbose"]]),as.integer(control[["trace"]]),
                              as.integer(grouped.vars-1), ## int.vector of variables which are mixed model nodes -1 for C
                              as.integer(group.ids),  ## group memberships - note indexed from 1
                              as.double(control[["epsabs.inner"]]),
                              as.integer(control[["max.iters.inner"]]),
                              as.double(control[["finite.step.size"]]),
                              as.double(control[["hessian.params"]]),
                              as.integer(control[["max.iters.hessian"]]),
                              as.integer(0),  ## modes only - false here as only applies to glmms
                              as.double(control[["max.hessian.error"]]),       ## Not applicable
                              as.double(control[["factor.brent"]]),            ## Not applicable
                              as.integer(control[["maxiters.hessian.brent"]]), ## Not applicable
                              as.double(control[["num.intervals.brent"]]),     ## Not applicable
                              PACKAGE="abn"))
      if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){
        if (verbose) cat("## !!! Laplace approximation failed\n")
        FAILED <- TRUE
      }
      used.inla <- FALSE

    } else {
      ## use INLA for glm
      if(!requireNamespace("INLA", quietly = TRUE)){stop("library INLA is not available!\nINLA is available from https://www.r-inla.org/download-install.")}
      mean.intercept <- control[["mean"]]
      ## use same as for rest of linear terms
      prec.intercept <- control[["prec"]]
      ## use same as for rest of linear terms
      if (verbose) cat("Using INLA (glm)\n")
      res.inla <- calc.node.inla.glm(child,
                                     dag.m,
                                     data.df,
                                     data.dists,
                                     rep(1,dim(data.df)[1]),
                                     ## ntrials
                                     rep(1,dim(data.df)[1]),
                                     ## exposure
                                     TRUE, mean.intercept, prec.intercept, control[["mean"]], control[["prec"]],control[["loggam.shape"]],control[["loggam.inv.scale"]],
                                     verbose.loc = verbose,
                                     nthreads = control[["ncores"]])

      if(is.logical(res.inla)){
        if (verbose) cat("INLA failed... so reverting to internal code.\n")
        r <- try(res.c <- .Call("fit_single_node",
                                data.df,
                                as.integer(child),    ## childnode
                                as.integer(dag.m[child,]), ## parent combination
                                as.integer(dim(dag.m)[1]), ## number of nodes/variables
                                as.integer(var.types),     ## type of densities
                                as.integer(sum(dag.m[child,])), ## max.parents
                                as.double(control[["mean"]]),as.double(1/sqrt(control[["prec"]])),as.double(control[["loggam.shape"]]),as.double(1/control[["loggam.inv.scale"]]),
                                as.integer(control[["max.iters"]]),as.double(control[["epsabs"]]),
                                as.integer(verbose),as.integer(control[["error.verbose"]]),as.integer(control[["trace"]]),
                                as.integer(grouped.vars-1),  ## int.vector of variables which are mixed model nodes -1 for C
                                as.integer(group.ids),                 ## group memberships - note indexed from 1
                                as.double(control[["epsabs.inner"]]),
                                as.integer(control[["max.iters.inner"]]),
                                as.double(control[["finite.step.size"]]),
                                as.double(control[["hessian.params"]]),
                                as.integer(control[["max.iters.hessian"]]),
                                as.integer(0),     ## modes only - false here as only applies to glmms
                                as.double(control[["max.hessian.error"]]),       ## Not applicable
                                as.double(control[["factor.brent"]]),            ## Not applicable
                                as.integer(control[["maxiters.hessian.brent"]]), ## Not applicable
                                as.double(control[["num.intervals.brent"]]),     ## Not applicable
                                PACKAGE="abn"))
        if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){
          if (verbose) cat("## !!! Laplace approximation failed\n")
          FAILED <- TRUE
        }
        used.inla <- FALSE
        ## flip
      }
      ## INLA failed

    }
    ## use INLA
    ## End of GLM node
    ###########################################################

  } else {
    ###########################################################
    ## Have a GLMM node
    ###########################################################
    ## have a glmm, so two options, INLA or C
    if(force.method=="notset" || force.method=="INLA"){##
      if(!requireNamespace("INLA", quietly = TRUE)){
        stop("library INLA is not available!\nINLA is available from https://www.r-inla.org/download-install.");
      }
      mean.intercept <- control[["mean"]]
      ## use same as for rest of linear terms
      prec.intercept <- control[["prec"]]
      ## use same as for rest of linear terms
      res.inla <- calc.node.inla.glmm(child,
                                      dag.m,
                                      data.frame(data.df,group=group.ids),
                                      data.dists,
                                      rep(1,dim(data.df)[1]),   ## ntrials
                                      rep(1,dim(data.df)[1]),   ## exposure
                                      TRUE,    ## always compute marginals - since only way to check results
                                      mean.intercept, prec.intercept, control[["mean"]], control[["prec"]],control[["loggam.shape"]],control[["loggam.inv.scale"]],
                                      verbose.loc = verbose,
                                      nthreads = control[["ncores"]])
      ## CHECK FOR INLA CRASH
      if(is.logical(res.inla)){
        if (verbose) cat("INLA failed...  so reverting to internal code\n");
        orig.force.method <- force.method;## save original
        force.method="C"; ## Easiest way is just to force C for this node
      } else {
        res.inla.modes <- getModeVector(list.fixed=res.inla$marginals.fixed,list.hyper=res.inla$marginals.hyperpar);
      }
    }
    if (verbose) cat("fit a glmm at node ",rownames(dag.m)[child],"using C\n");
    if(force.method=="notset"){
      r <- try(res.c <- .Call("fit_single_node",
                              data.df,
                              as.integer(child),        ## childnode
                              as.integer(dag.m[child,]),## parent combination
                              as.integer(dim(dag.m)[1]),## number of nodes/variables
                              as.integer(var.types),## type of densities
                              as.integer(sum(dag.m[child,])),## max.parents
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
                              as.integer(1), ## turn on ModesONLY
                              as.double(control[["max.hessian.error"]]),## Not applicable
                              as.double(control[["factor.brent"]]),     ## Not applicable
                              as.integer(control[["maxiters.hessian.brent"]]),## Not applicable
                              as.double(control[["num.intervals.brent"]])## Not applicable
                              ,PACKAGE="abn" ## uncomment to load as package not shlib
      ));

      if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){ if (verbose) cat(" Laplace approximation failed\n");
        FAILED <- TRUE;}
      if(!FAILED){
        res.c.modes <- res.c[[1]][-c(1:3)];## remove mlik - this is first entry, and error code and hessian accuracy
        res.c.modes <- res.c.modes[which(res.c.modes!=.Machine$double.xmax)];## this discards all "empty" parameters
        ## get difference in modes proportion relative to C
        diff.in.modes <- (res.inla.modes-res.c.modes)/res.c.modes;
        error.modes <- max(abs(diff.in.modes));
      }
    } ## end of notset

    if( !FAILED && force.method=="C" || (force.method=="notset" && error.modes>(control[["max.mode.error"]]/100))){ ## INLA might be unreliable so use C (slower)
      if(force.method=="notset"){
        if (verbose) cat("Using internal code (Laplace glmm)\n=>max. abs. difference (in %) with INLA is ");
        if (verbose) cat(formatC(100*error.modes,format="f",digits=1)," and exceeds tolerance\n");
      } else {
        if (verbose) cat("Using internal code (Laplace glmm)\n");}
      r <- try(res.c <- .Call("fit_single_node",
                              data.df,
                              as.integer(child), ## childnode
                              as.integer(dag.m[child,]),## parent combination
                              as.integer(dim(dag.m)[1]),## number of nodes/variables
                              as.integer(var.types),## type of densities
                              as.integer(sum(dag.m[child,])),## max.parents
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
                              as.integer(0), ## turn on ModesONLY
                              as.double(control[["max.hessian.error"]]),## Not applicable
                              as.double(control[["factor.brent"]]),## Not applicable
                              as.integer(control[["maxiters.hessian.brent"]]),## Not applicable
                              as.double(control[["num.intervals.brent"]])## Not applicable
                              ,PACKAGE="abn" ## uncomment to load as package not shlib
      ));
      if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){
        if (verbose) cat("## !!! Laplace approximation failed\n");
        FAILED <- TRUE;
      }
      used.inla <- FALSE;## flip
    } else {
      if (verbose) cat("Using INLA (glmm)\n");
    }## end of if inla bad
  } ## end of if GLMM
  ###########################################################
  ## End of GLMM node
  ###########################################################

  ###########################################################
  ## End of all external computations
  ###########################################################
  ## computation for current node is all done so sort out the
  ## output into nicer form and give labels
  ###########################################################
  if(!FAILED){
    if(used.inla==FALSE){## organize output from C
      mlik <- res.c[[1]][1]
      error.code <- res.c[[1]][2]
      hessian.accuracy <- res.c[[1]][3]
      used.INLA <- FALSE
    } else {
      ## organize output from INLA
      mlik <- res.inla$mlik[2]## [2] is for Gaussian rather than Integrated estimate
      error.code <- NA## not available from INLA
      hessian.accuracy <- NA## not available from INLA
      used.INLA <- TRUE
    }
  } else {## FAILED
    mlik <- NA
    error.code <- 2## model could not be fitted
    hessian.accuracy <- NA
    used.INLA <- FALSE
  }
  if(!is.null(orig.force.method)){
    force.method <- orig.force.method;
  } ## reset force.method after INLA crash
  ############################################################
  ## Finished with current node
  ############################################################
  return(c(childParentCombNo=row.no, mlik=mlik,error.code=error.code,hessian.accuracy=hessian.accuracy,used.INLA=used.INLA))
} ## end of nodes loop
