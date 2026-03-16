# Extracted from buildScoreCache.mle()
# Note that: group.var is assumed to be NULL
call_forLoopContent_orig <- function(adj.vars,
                              nvars,
                              data.df,
                              data.df.multi,
                              data.df.lvl,
                              max.parents,
                              data.dists,
                              mycache,
                              control,
                              verbose){

  # ========= END CACHE CREATION ===========

  ## EOF cache creation
  row.num <- NULL   # To avoid check comment: 'no visible binding for global variable
  out <- list()
  rows <- length(mycache[["children"]])

  ##-----------------------------
  ##start loop for the regression
  ##-----------------------------
  if(verbose){cat("Start estimation loop.")}
  {
    res <- foreach(row.num = 1:rows,
                   .combine='rbind',
                   .packages = c("stats", "lme4", "mclogit", "nnet"),
                   .export = 'forLoopContent',
                   .verbose = verbose) %do% {
                     forLoopContent(row.num = row.num,
                                    mycache = mycache,
                                    data.dists = data.dists,
                                    data.df.multi = data.df.multi,
                                    adj.vars = NULL, #adj.vars,
                                    data.df = data.df,
                                    data.df.lvl = data.df.lvl,
                                    group.var = NULL, #group.var,
                                    group.ids = NULL, #group.ids,
                                    control = control,
                                    n = nvars,
                                    verbose = FALSE) #verbose)
                   }

  }

  out[["children"]] <- mycache[["children"]]
  out[["node.defn"]] <- mycache$node.defn
  out[["mlik"]] <- as.numeric( res[,1] )
  out[["error.code"]] <- list()
  out[["hessian.accuracy"]] <- list()
  out[["used.INLA"]] <- list()
  out[["error.code.desc"]] <- list()
  out[["data.df"]] <- data.df.lvl
  out[["data.dists"]] <- data.dists
  out[["max.parents"]] <- max.parents
  out[["dag.retained"]] <- NULL #dag.retained
  out[["dag.banned"]] <- NULL #dag.banned
  out[["group.var"]] <- NULL #group.var
  out[["group.ids"]] <- NULL #group.ids
  out[["grouped.vars"]] <- NULL #grouped.vars
  out[["aic"]] <- as.numeric( res[,2] )
  out[["bic"]] <- as.numeric( res[,3] )
  out[["mdl"]] <- as.numeric( res[,4] )

  out[["method"]] <- "mle"

  return(out)
}


# Modified function to avoid the rbind-in-loop while using foreach
# Note that: group.var is assumed to be NULL
call_forLoopContent_docall <- function(adj.vars,
                                     nvars,
                                     data.df,
                                     data.df.multi,
                                     data.df.lvl,
                                     max.parents,
                                     data.dists,
                                     mycache,
                                     control,
                                     verbose){

  # ========= END CACHE CREATION ===========

  ## EOF cache creation
  row.num <- NULL   # To avoid check comment: 'no visible binding for global variable
  out <- list()
  rows <- length(mycache[["children"]])

  ##-----------------------------
  ##start loop for the regression
  ##-----------------------------
  if(verbose){cat("Start estimation loop.")}
  {
    rows <- foreach(row.num = 1:rows,
                   #.combine='rbind',
                   .packages = c("stats", "lme4", "mclogit", "nnet"),
                   .export = 'forLoopContent',
                   .verbose = verbose) %do% {
                     forLoopContent(row.num = row.num,
                                    mycache = mycache,
                                    data.dists = data.dists,
                                    data.df.multi = data.df.multi,
                                    adj.vars = NULL, #adj.vars,
                                    data.df = data.df,
                                    data.df.lvl = data.df.lvl,
                                    group.var = NULL, #group.var,
                                    group.ids = NULL, #group.ids,
                                    control = control,
                                    n = nvars,
                                    verbose = FALSE) #verbose)
                   }
    res <- do.call(rbind, rows)

  }

  out[["children"]] <- mycache[["children"]]
  out[["node.defn"]] <- mycache$node.defn
  out[["mlik"]] <- as.numeric( res[,1] )
  out[["error.code"]] <- list()
  out[["hessian.accuracy"]] <- list()
  out[["used.INLA"]] <- list()
  out[["error.code.desc"]] <- list()
  out[["data.df"]] <- data.df.lvl
  out[["data.dists"]] <- data.dists
  out[["max.parents"]] <- max.parents
  out[["dag.retained"]] <- NULL #dag.retained
  out[["dag.banned"]] <- NULL #dag.banned
  out[["group.var"]] <- NULL #group.var
  out[["group.ids"]] <- NULL #group.ids
  out[["grouped.vars"]] <- NULL #grouped.vars
  out[["aic"]] <- as.numeric( res[,2] )
  out[["bic"]] <- as.numeric( res[,3] )
  out[["mdl"]] <- as.numeric( res[,4] )

  out[["method"]] <- "mle"

  return(out)
}


# Modified function with pure matrix preallocation (no foreach)
# Note that: group.var is assumed to be NULL
call_forLoopContent_noforeach <- function(adj.vars,
                                         nvars,
                                         data.df,
                                         data.df.multi,
                                         data.df.lvl,
                                         max.parents,
                                         data.dists,
                                         mycache,
                                         control,
                                         verbose){

  # ========= END CACHE CREATION ===========

  ## EOF cache creation
  row.num <- NULL   # To avoid check comment: 'no visible binding for global variable
  out <- list()
  rows <- length(mycache[["children"]])

  ##-----------------------------
  ##start loop for the regression
  ##-----------------------------
  if(verbose){cat("Start estimation loop.")}
  {
    res <- matrix(NA_real_, rows, 4)
    for(row.num in 1:rows)
          res[row.num,] <- forLoopContent(row.num = row.num,
                                     mycache = mycache,
                                     data.dists = data.dists,
                                     data.df.multi = data.df.multi,
                                     adj.vars = NULL, #adj.vars,
                                     data.df = data.df,
                                     data.df.lvl = data.df.lvl,
                                     group.var = NULL, #group.var,
                                     group.ids = NULL, #group.ids,
                                     control = control,
                                     n = nvars,
                                     verbose = FALSE) #verbose)

  }

  out[["children"]] <- mycache[["children"]]
  out[["node.defn"]] <- mycache$node.defn
  out[["mlik"]] <- as.numeric( res[,1] )
  out[["error.code"]] <- list()
  out[["hessian.accuracy"]] <- list()
  out[["used.INLA"]] <- list()
  out[["error.code.desc"]] <- list()
  out[["data.df"]] <- data.df.lvl
  out[["data.dists"]] <- data.dists
  out[["max.parents"]] <- max.parents
  out[["dag.retained"]] <- NULL #dag.retained
  out[["dag.banned"]] <- NULL #dag.banned
  out[["group.var"]] <- NULL #group.var
  out[["group.ids"]] <- NULL #group.ids
  out[["grouped.vars"]] <- NULL #grouped.vars
  out[["aic"]] <- as.numeric( res[,2] )
  out[["bic"]] <- as.numeric( res[,3] )
  out[["mdl"]] <- as.numeric( res[,4] )

  out[["method"]] <- "mle"

  return(out)
}



# Modified function to avoid any in-loop concatenation:
# test the foreach without any .combine and do the combine into a matrix afterwords
# Note that: group.var is assumed to be NULL
call_forLoopContent_foreachrows <- function(adj.vars,
                                         nvars,
                                         data.df,
                                         data.df.multi,
                                         data.df.lvl,
                                         max.parents,
                                         data.dists,
                                         mycache,
                                         control,
                                         verbose){

  # ========= END CACHE CREATION ===========

  ## EOF cache creation
  row.num <- NULL   # To avoid check comment: 'no visible binding for global variable
  out <- list()
  rows <- length(mycache[["children"]])

  ##-----------------------------
  ##start loop for the regression
  ##-----------------------------
  if(verbose){cat("Start estimation loop.")}
  {
    lines <- foreach(row.num = 1:rows,
                   .packages = c("stats", "lme4", "mclogit", "nnet"),
                   .export = 'forLoopContent',
                   .verbose = verbose) %do% {
                    forLoopContent(row.num = row.num,
                                    mycache = mycache,
                                    data.dists = data.dists,
                                    data.df.multi = data.df.multi,
                                    adj.vars = NULL, #adj.vars,
                                    data.df = data.df,
                                    data.df.lvl = data.df.lvl,
                                    group.var = NULL, #group.var,
                                    group.ids = NULL, #group.ids,
                                    control = control,
                                    n = nvars,
                                    verbose = FALSE) #verbose)
                   }
    res <- matrix(unlist(lines), ncol = 4, byrow = TRUE)
  }

  out[["children"]] <- mycache[["children"]]
  out[["node.defn"]] <- mycache$node.defn
  out[["mlik"]] <- as.numeric( res[,1] )
  out[["error.code"]] <- list()
  out[["hessian.accuracy"]] <- list()
  out[["used.INLA"]] <- list()
  out[["error.code.desc"]] <- list()
  out[["data.df"]] <- data.df.lvl
  out[["data.dists"]] <- data.dists
  out[["max.parents"]] <- max.parents
  out[["dag.retained"]] <- NULL #dag.retained
  out[["dag.banned"]] <- NULL #dag.banned
  out[["group.var"]] <- NULL #group.var
  out[["group.ids"]] <- NULL #group.ids
  out[["grouped.vars"]] <- NULL #grouped.vars
  out[["aic"]] <- as.numeric( res[,2] )
  out[["bic"]] <- as.numeric( res[,3] )
  out[["mdl"]] <- as.numeric( res[,4] )

  out[["method"]] <- "mle"

  return(out)
}
