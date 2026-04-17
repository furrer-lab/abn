# Extracted from buildScoreCache.mle()
# Note that: group.var is assumed to be NULL, so only the corresponding code has been extracted
forLoopContent_copy <-
  function(row.num,
           mycache,
           data.dists,
           data.df.multi,
           adj.vars,
           data.df,
           data.df.lvl,
           group.var,
           group.ids,
           control,
           n,
           verbose) {
    # if (is.null(group.var))
    {
      # we have no grouping (do what was always done).
      if (verbose) message("we have no grouping (no mixed-effects model)")
      child <- mycache[["children"]][row.num] # child as integer
      distribution <- data.dists[child]
      Y <- data.matrix(data.df[, child])
      if (verbose) {
        # current node's name
        child.name <- colnames(mycache[["node.defn"]])[child]
        # current node's parents names
        parents.names <- names(which(mycache[["node.defn"]][row.num,] != 0))
        message(paste("regressing", child.name, "on", paste(parents.names, collapse = ", ")))
      }


      if (is.null(adj.vars)) {
        if ("multinomial" %in% data.dists[as.logical(mycache$node.defn[row.num, ])]) {
          #X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
          X <- data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ]), drop = FALSE]
        } else {
          X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
        }
      } else {
        if ("multinomial" %in% data.dists[as.logical(mycache$node.defn.adj[row.num, ])]) {
          #X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
          X <- data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ]), drop = FALSE]
        } else {
          X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
        }
      }
      ## Rank deficiency
      num.na <- 0

      R <- rank_cpp(X)
      r <- ncol(X)
      R_col <- R/r

      if (R_col != 1 & as.character(distribution) == "binomial") {
        Y1 <- if (is.factor(Y)) numeric(Y) else  Y

        while (rank_cpp(X)/ncol(X) != 1) {
          X <- X[, -1]
          num.na <- num.na + 1
          if (is.null(ncol(X)))
            X <- as.matrix(X)
        }

        tryCatch(fit <- irls_binomial_cpp_fast_br(A = X, b = Y1, maxit = control[["max.iters"]], tol = control[["tol"]]))
        # tryCatch(fit <- irls_binomial_cpp_fast_br(A = X, b = Y1, maxit = control[["max.iters"]], tol = control[["tol"]]),
        #  error = function(e) {
        #       while (rank_cpp(X)/ncol(X) != 1) {
        #         X <- X[, -1]
        #         num.na <- num.na + 1
        #         if (is.null(ncol(X)))
        #           X <- as.matrix(X)
        #       }
        #       fit <- irls_binomial_cpp_fast_br(A = X, b = Y1, maxit = control[["max.iters"]], tol = control[["tol"]])
        #   }, finally = fit)
      } else {
        switch(as.character(distribution),
               binomial = {
                 Y1 <- if (is.factor(Y)) numeric(Y) else Y
                 fit <- irls_binomial_cpp_fast_br(A = X, b = Y1, maxit = control[["max.iters"]], tol = control[["tol"]])
                 if (is.na(sum(fit[[1]]))) fit <- irls_binomial_cpp_fast_br(A = X, b = Y, maxit = control[["max.iters"]], tol = control[["tol"]])
               }, gaussian = {
                 suppressWarnings(
                   fit <- irls_gaussian_cpp_fast(A = X, b = Y, maxit = control[["max.iters"]], tol = control[["tol"]])
                  )
               }, poisson = {
                 suppressWarnings(
                   fit <- irls_poisson_cpp_fast(A = X, b = Y, maxit = control[["max.iters"]], tol = control[["tol"]])
                   )
               }, multinomial = {

                 Ymulti <- data.matrix(model.matrix(~-1 + data.df.lvl[, child]))

                 p <- ncol(Ymulti)
                 mask <- c(rep(FALSE, r + 1L), rep(c(FALSE, rep(TRUE, r)), p - 1L))

                 tmp <- nnet.default(x = X, y = Ymulti, mask = mask, size = 0,
                                     skip = TRUE, softmax = TRUE, rang = 0, trace = FALSE)

                 fit <- NULL
                 fit$loglik <- -tmp$value
                 edf <- ifelse(length(tmp$lev) == 2L, 1, length(tmp$lev) - 1) * R
                 fit$aic <- 2 * tmp$value + 2 * edf
                 fit$bic <- 2 * tmp$value + edf * log(dim(data.df)[1])
              })
      }

      # Prepare return values
      if (!is.null(fit)) {
        if (verbose) {
          message("Successfully fitted local model.")
        }
        c(fit$loglik,
          fit$aic,
          fit$bic,
          fit$bic + (1 + sum(mycache[["node.defn.multi"]][row.num, ]) - num.na) * log(n))
      } else if (is.null(fit)) {
        if (verbose) {
          message("Failed to fit local model. Returning very low scores.")
        }
        # no convergence, singularity, rank-deficiency, return very low scores
        c(rep(-Inf, 4))
      } else {
        stop("Unknown state of fit. I should never end up here.")
      }
    }
  } # End of forLoopContent()

# Modified version of the forLoopContent_copy function extracted from buildScoreCache.mle()
# X and Y are precomputed outside and passed as argument to the function
forLoopContent_precomputedXY <-
  function(row.num,
           X,
           Y,
           mycache,
           data.dists,
           data.df.multi,
           adj.vars,
           data.df,
           data.df.lvl,
           group.var,
           group.ids,
           control,
           n,
           verbose) {
    {
      # we have no grouping (do what was always done).
      if (verbose) message("we have no grouping (no mixed-effects model)")
      child <- mycache[["children"]][row.num] # child as integer
      distribution <- data.dists[child]
      #Y <- data.matrix(data.df[, child])
      #Y <- data.df[, child, drop = FALSE]
      #Y <- data.matrix(mycache$Y[[child]])
      if (verbose) {
        # current node's name
        child.name <- colnames(mycache[["node.defn"]])[child]
        # current node's parents names
        parents.names <- names(which(mycache[["node.defn"]][row.num,] != 0))
        message(paste("regressing", child.name, "on", paste(parents.names, collapse = ", ")))
      }


      #if (is.null(adj.vars)) {
      #  if ("multinomial" %in% data.dists[as.logical(mycache$node.defn[row.num, ])]) {
      #    #X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      #    X <- data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ]), drop = FALSE]
      #  } else {
      #X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      #X <- cbind(1, data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ]), drop = FALSE])
      #    p <- sum(as.logical(mycache[["node.defn.multi"]][row.num, ]))
      #    X <- matrix(NA_real_, nrow = length(data.df[, 1]), ncol = p + 1)
      #    X[, 1] <- 1
      #    X[, 2:(p + 1)] <- data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ]), drop = FALSE]
      #  }
      #} else {
      #  if ("multinomial" %in% data.dists[as.logical(mycache$node.defn.adj[row.num, ])]) {
      #X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      #    X <- data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ]), drop = FALSE]
      #  } else {
      #X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      #X <- cbind(1, data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ]), drop = FALSE])
      #    p <- sum(as.logical(mycache[["node.defn.multi"]][row.num, ]))
      #    X <- matrix(NA_real_, nrow = length(data.df[, 1]), ncol = p + 1)
      #    X[, 1] <- 1
      #    X[, 2:(p + 1)] <- data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ]), drop = FALSE]
      #  }
      #}
      ##input_binomial <- readRDS(file = "scripts/input_binomial_cpp.RData")
      ##X <- input_binomial$A
      ## Rank deficiency
      #X <- data.matrix(mycache$X[[row.num]])
      num.na <- 0

      R <- rank_cpp(X)
      r <- ncol(X)
      R_col <- R/r

      if (R_col != 1 & as.character(distribution) == "binomial") {
        Y1 <- if (is.factor(Y)) numeric(Y) else  Y

        while (rank_cpp(X)/ncol(X) != 1) {
          X <- X[, -1]
          num.na <- num.na + 1
          if (is.null(ncol(X)))
            X <- as.matrix(X)
        }

        #tryCatch(
        fit <- irls_binomial_cpp_fast_br(A = X, b = Y1, maxit = control[["max.iters"]], tol = control[["tol"]])
        #)
        # tryCatch(fit <- irls_binomial_cpp_fast_br(A = X, b = Y1, maxit = control[["max.iters"]], tol = control[["tol"]]),
        #  error = function(e) {
        #       while (rank_cpp(X)/ncol(X) != 1) {
        #         X <- X[, -1]
        #         num.na <- num.na + 1
        #         if (is.null(ncol(X)))
        #           X <- as.matrix(X)
        #       }
        #       fit <- irls_binomial_cpp_fast_br(A = X, b = Y1, maxit = control[["max.iters"]], tol = control[["tol"]])
        #   }, finally = fit)
      } else {

        switch(as.character(distribution),
               binomial = {
                 Y1 <- if (is.factor(Y)) numeric(Y) else Y
                 fit <- irls_binomial_cpp_fast_br(A = X, b = Y1, maxit = control[["max.iters"]], tol = control[["tol"]])
                 if (is.na(sum(fit[[1]]))) fit <- irls_binomial_cpp_fast_br(A = X, b = Y, maxit = control[["max.iters"]], tol = control[["tol"]])
               }, gaussian = {
                 suppressWarnings(
                   fit <- irls_gaussian_cpp_fast(A = X, b = Y, maxit = control[["max.iters"]], tol = control[["tol"]])
                 )
               }, poisson = {
                 suppressWarnings(
                   fit <- irls_poisson_cpp_fast(A = X, b = Y, maxit = control[["max.iters"]], tol = control[["tol"]])
                 )
               }, multinomial = {

                 #Ymulti <- data.matrix(model.matrix(~-1 + data.df.lvl[, child]))
                 levels <- levels(data.df.lvl[, child])
                 Ymulti <- data.df.multi[,paste0(names(data.df.lvl)[child], levels)]
                 #Ymulti <- model.matrix(~-1 + data.df.lvl[, child])

                 #p <- ncol(Ymulti)
                 p <- length(levels)
                 mask <- c(rep(FALSE, r + 1L), rep(c(FALSE, rep(TRUE, r)), p - 1L))

                 tmp <- nnet.default(x = X, y = Ymulti, mask = mask, size = 0,
                                     skip = TRUE, softmax = TRUE, rang = 0, trace = FALSE)

                 fit <- NULL
                 fit$loglik <- -tmp$value
                 edf <- ifelse(length(tmp$lev) == 2L, 1, length(tmp$lev) - 1) * R
                 fit$aic <- 2 * tmp$value + 2 * edf
                 fit$bic <- 2 * tmp$value + edf * log(dim(data.df)[1])
               })
      }

      # Prepare return values
      if (!is.null(fit)) {
        if (verbose) {
          message("Successfully fitted local model.")
        }
        c(fit$loglik,
          fit$aic,
          fit$bic,
          fit$bic + (1 + sum(mycache[["node.defn.multi"]][row.num, ]) - num.na) * log(n))
      } else if (is.null(fit)) {
        if (verbose) {
          message("Failed to fit local model. Returning very low scores.")
        }
        # no convergence, singularity, rank-deficiency, return very low scores
        c(rep(-Inf, 4))
      } else {
        stop("Unknown state of fit. I should never end up here.")
      }
    }
  } # End of forLoopContent()


# Function calling forLoopContent within the foreach loop
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
# Precompute X and Y
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
    for(row.num in 1:rows){
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
                                     verbose = FALSE)
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


# Help function to compute X
compute_X <- function(row, adj.vars, multinomial, data.df.multi) {
  if (is.null(adj.vars)) {
    if (multinomial) {
      #X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      X <- data.df.multi[, as.logical(row), drop = FALSE]
    } else {
      #X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      X <- cbind(1, data.df.multi[, as.logical(row), drop = FALSE])
    }
  } else {
    if (multinomial) {
      #X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      X <- data.df.multi[, as.logical(row), drop = FALSE]
    } else {
      #X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      X <- cbind(1, data.df.multi[, as.logical(row), drop = FALSE])
    }
  }
  X
}


# Modified function with pure matrix preallocation (no foreach)
# Precompute X and Y
# Note that: group.var is assumed to be NULL
call_forLoopContent_precomputedXY <- function(adj.vars,
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
    data.df.multi.num <- data.matrix(data.df.multi)
    data.df.num <- data.matrix(data.df)
    cache_X <- new.env(hash = TRUE, parent = emptyenv())
    Y_list <- vector("list", nvars)
    for (k in 1:nvars){
      Y_list[[k]] <- data.df.num[, k, drop = FALSE]
    }
    for(row.num in 1:rows){
      key <- bits_to_key(mycache$node.defn.multi[row.num, ])
      if (exists(key, envir = cache_X, inherits = FALSE)) {
        X <- get(key, envir = cache_X, inherits = FALSE)
      } else {
        if ("multinomial" %in% data.dists[as.logical(mycache$node.defn[row.num, ])])
          multinomial = TRUE
        else multinomial = FALSE
        obj <- compute_X(mycache$node.defn.multi[row.num, ], adj.vars, multinomial, data.df.multi.num)
        assign(key, obj, envir = cache_X)
        X <- obj
      }
      child <- mycache[["children"]][row.num]
      Y <- Y_list[[child]]
      res[row.num,] <- forLoopContent_precomputedXY(row.num = row.num,
                                                    X = X,
                                                    Y = Y,
                                                    mycache = mycache,
                                                    data.dists = data.dists,
                                                    data.df.multi = data.df.multi.num,
                                                    adj.vars = NULL, #adj.vars,
                                                    data.df = data.df.num,
                                                    data.df.lvl = data.df.lvl,
                                                    group.var = NULL, #group.var,
                                                    group.ids = NULL, #group.ids,
                                                    control = control,
                                                    n = nvars,
                                                    verbose = FALSE)
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



