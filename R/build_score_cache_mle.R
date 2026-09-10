#' From each child-parent(s) combination, regress each child on its parents in buildScoreCache.mle()
#' @describeIn buildScoreCache Internal function called by \code{buildScoreCache.mle()}.
#' @param row.num number of child-node (mostly corresponds to child node index e.g. in dag).
#' @param mycache prepared cache.
#' @param data.df.multi extended data.df for one-hot-encoded multinomial variables.
#' @param data.df.lvl copy of original \code{data.df}.
#' @param n corresponds to \code{nvars}, number of variables in data.dists.
#' @importFrom stats AIC BIC sd model.matrix as.formula lm glm logLik
#' @importFrom nnet nnet.default multinom
#' @importFrom mclogit mblogit
#' @returns A list that will be passed to \code{\link{buildScoreCache.mle}}.
#' @keywords internal
forLoopContent <-
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
    if (!is.null(group.var)) {
      if (verbose) {
        message("we have grouping (mixed-effects model) glmm")
      }
    ###
    # Build arguments for glmm.bic
    ###
    # current node's name
    child <- mycache[["children"]][row.num] # child as integer
    child.name <- colnames(mycache[["node.defn"]])[child]
    # current node's parents names
    parents.names <- names(which(mycache[["node.defn"]][row.num,] != 0))
    # current node's distribution
    child.dist <- data.dists[child] #       distribution <- data.dists[child]
    # data.df including group.var=group.ids
    data.df.grouping <- data.frame(data.df, grp = factor(group.ids))
    names(data.df.grouping) <- c(names(data.df), group.var)

    ###
    # Beginn of GLMM
    ###
    # Build formula statement
    if (length(parents.names) == 0){
      # no parent only random effect
      model <- as.formula(paste(child.name, "~ 1 + (1|", group.var, ")"))
    } else if (length(parents.names) > 0){
      model <- as.formula(paste(child.name, "~ (1|", group.var, ")+", paste(parents.names, collapse = "+")))
    } else {
      stop("Cannot build model formula. Unknown predictor names of parent nodes in 'parents.name'.")
    }

    # Main part: Depending on child's distribution, call the appropriate modelling function
    switch (as.character(child.dist),
            gaussian = {
              if (verbose) {message(paste("using lmer with model:", deparse1(model)))} else NA
                tryCatch({
                  fit <- lme4::lmer(model, data = data.df.grouping)
                }, error=function(e)NULL)

                if (is.null(fit)){
                  # relax tolerances for change in parameter values and objective function.
                  tryCatch({fit <- lme4::lmer(model, data = data.df.grouping,
                                              control = lme4::lmerControl(check.rankX = control[["check.rankX"]],
                                                                          check.scaleX = control[["check.scaleX"]],
                                                                          check.conv.grad = lme4::.makeCC(action = control[["check.conv.grad"]]),
                                                                          check.conv.singular = lme4::.makeCC(action = control[["check.conv.singular"]]),
                                                                          check.conv.hess = lme4::.makeCC(action = control[["check.conv.hess"]]),
                                                                          optCtrl = list(xtol_abs = control[["xtol_abs"]],
                                                                                         ftol = control[["ftol_abs"]])))
                  }, error=function(e)NULL)
                }

              # if fit is still NULL, do not modify model further as this would change the structure of the dag but return very low score (further down)
            },
            binomial = {
                if (verbose) {message(paste("using glmer with model:", deparse1(model)))} else NA
                tryCatch({
                  fit <- lme4::glmer(model, data = data.df.grouping, family = "binomial")
                }, error=function(e)NULL)

                if (is.null(fit)){
                  # relax tolerances for change in parameter values and objective function.
                  tryCatch({fit <- lme4::glmer(model, data = data.df.grouping, family = "binomial",
                                               control = lme4::glmerControl(tolPwrss = control[["tolPwrss"]],
                                                                            check.rankX = control[["check.rankX"]],
                                                                            check.scaleX = control[["check.scaleX"]],
                                                                            check.conv.grad = lme4::.makeCC(action = control[["check.conv.grad"]]),
                                                                            check.conv.singular = lme4::.makeCC(action = control[["check.conv.singular"]]),
                                                                            check.conv.hess = lme4::.makeCC(action = control[["check.conv.hess"]]),
                                                                            optCtrl = list(xtol_abs = control[["xtol_abs"]],
                                                                                           ftol = control[["ftol_abs"]])))
                  }, error=function(e)NULL)
                }

                # if fit is still NULL, do not modify model further as this would change the structure of the dag but return very low score (further down)
            },
            poisson = {
              if (verbose) {message(paste("using glmer with model:", deparse1(model)))} else NA

              fit <- NULL
              if(control[["only_glmmTMB_poisson"]] == FALSE) {
                tryCatch({
                  fit <- lme4::glmer(model, data = data.df.grouping, family = "poisson")
                }, error=function(e)NULL)

                if (is.null(fit)){
                  # relax tolerances for change in parameter values and objective function.
                  tryCatch({fit <- lme4::glmer(model, data = data.df.grouping, family = "poisson",
                                               control = lme4::glmerControl(tolPwrss = control[["tolPwrss"]],
                                                                            check.rankX = control[["check.rankX"]],
                                                                            check.scaleX = control[["check.scaleX"]],
                                                                            check.conv.grad = lme4::.makeCC(action = control[["check.conv.grad"]]),
                                                                            check.conv.singular = lme4::.makeCC(action = control[["check.conv.singular"]]),
                                                                            check.conv.hess = lme4::.makeCC(action = control[["check.conv.hess"]]),
                                                                            optCtrl = list(xtol_abs = control[["xtol_abs"]],
                                                                                           ftol = control[["ftol_abs"]])))
                  }, error=function(e)NULL)
                }
              } else if(control[["only_glmmTMB_poisson"]] == TRUE) {
                if (is.null(fit)){
                  # try glmmTMB as alternative
                  if (verbose) {message(paste("trying glmmTMB with model:", deparse1(model)))} else NA
                  tryCatch({
                    fit <- glmmTMB::glmmTMB(model, data = data.df.grouping, family = "poisson")
                  }, error=function(e)NULL)
                }
              } else {
                stop("Invalid 'only_glmmTMB_poisson' argument. Must be one of TRUE or FALSE.")
              }

              # if fit is still NULL, do not modify model further as this would change the structure of the dag but return very low score (further down)
            },
            multinomial = {
              if (length(parents.names) == 0){
                model_basic <- as.formula(paste(child.name, "~ 1"))
                model_random <- as.formula(paste("~ 1|", group.var, sep = ""))
              } else {
                model_basic <- as.formula(paste(child.name, "~ ", paste(parents.names, collapse = "+")))
                model_random <- as.formula(paste("~ 1|", group.var, sep = ""))
              }
              if(verbose) message(paste("using mblogit with fixed term:", deparse1(model_basic), "and random term:", deparse1(model_random))) else NA

              if (control[["catcov.mblogit"]] == "free"){
                tryCatch({
                  fit <- mclogit::mblogit(formula = model_basic, random = model_random, data = data.df.grouping, catCov = "free")
                }, error=function(e) NULL)
              } else if (control[["catcov.mblogit"]] == "diagonal"){
                tryCatch({
                  fit <- mclogit::mblogit(formula = model_basic, random = model_random, data = data.df.grouping, catCov = "diagonal")
                  # manipulate VarCov to bring in correct shape
                  fit_vcov <- matrix(0, nrow = length(fit$VarCov), ncol = length(fit$VarCov))
                  diag(fit_vcov) <- unlist(fit$VarCov, use.names = FALSE)
                  # col and rownames
                  cn <- c()
                  rn <- c()
                  for(i in 1:length(fit$VarCov)){
                    cn[i] <- colnames(fit$VarCov[[i]])
                    rn[i] <- rownames(fit$VarCov[[i]])
                  }
                  colnames(fit_vcov) <- cn
                  rownames(fit_vcov) <- rn
                  # replace
                  fit_vcov <- list(fit_vcov)
                  names(fit_vcov) <- group.var
                  fit$VarCov <- fit_vcov
                }, error=function(e) NULL)
              } else if (control[["catcov.mblogit"]] == "single"){
                stop("'catcov.mblogit' == 'single' is not yet implemented.")
                fit <- mclogit::mblogit(formula = model_basic, random = model_random, data = data.df.grouping, catCov = "single")
                # manipulate VarCov to bring in correct shape
              } else {
                stop("invalid 'catcov.mblogit' argument. Must be one of 'free', 'diagonal' or 'single'.")
              }

                if (is.null(fit)){
                  # relax tolerances for change in parameter values and objective function.
                  # tryCatch({fit <- mclogit::mblogit(formula = model_basic, random = model_random, data = data.df.grouping, catCov = control[["catcov.mblogit"]]
                                                    # control = mclogit::mclogit.control(epsilon = control[["epsilon"]],
                                                    #                                    trace = control[["trace.mblogit"]]))
                  # }, error=function(e)NULL)
                  if (control[["catcov.mblogit"]] == "free"){
                    tryCatch({
                      fit <- mclogit::mblogit(formula = model_basic, random = model_random, data = data.df.grouping, catCov = "free",
                                              control = mclogit::mclogit.control(epsilon = control[["epsilon"]],
                                                                                 trace = control[["trace.mblogit"]]))
                    }, error=function(e) NULL)
                  } else if (control[["catcov.mblogit"]] == "diagonal"){
                    tryCatch({
                      fit <- mclogit::mblogit(formula = model_basic, random = model_random, data = data.df.grouping, catCov = "diagonal",
                                              control = mclogit::mclogit.control(epsilon = control[["epsilon"]],
                                                                                 trace = control[["trace.mblogit"]]))
                      # manipulate VarCov to bring in correct shape
                      fit_vcov <- matrix(0, nrow = length(fit$VarCov), ncol = length(fit$VarCov))
                      diag(fit_vcov) <- unlist(fit$VarCov, use.names = FALSE)
                      # col and rownames
                      cn <- c()
                      rn <- c()
                      for(i in 1:length(fit$VarCov)){
                        cn[i] <- colnames(fit$VarCov[[i]])
                        rn[i] <- rownames(fit$VarCov[[i]])
                      }
                      colnames(fit_vcov) <- cn
                      rownames(fit_vcov) <- rn
                      # replace
                      fit_vcov <- list(fit_vcov)
                      names(fit_vcov) <- group.var
                      fit$VarCov <- fit_vcov
                    }, error=function(e) NULL)
                  } else if (control[["catcov.mblogit"]] == "single"){
                    stop("'catcov.mblogit' == 'single' is not yet implemented.")
                    fit <- mclogit::mblogit(formula = model_basic, random = model_random, data = data.df.grouping, catCov = "single",
                                            control = mclogit::mclogit.control(epsilon = control[["epsilon"]],
                                                                               trace = control[["trace.mblogit"]]))
                    # manipulate VarCov to bring in correct shape
                  } else {
                    stop("invalid 'catcov.mblogit' argument. Must be one of 'free', 'diagonal' or 'single'.")
                  }
                }
                # if fit is still NULL, do not modify model further as this would change the structure of the dag but return very low score (further down)
            }
    )
    ###
    # End of GLMM
    ###
    # collect values to return
    if(!is.null(fit)){
      if (verbose) {
        message("Successfully fitted local model.")
      }
      fit_loglik <- logLik(fit)
      fit_aic <- AIC(fit)
      fit_bic <- BIC(fit)
      # fit_mdl <- fit_bic + (1 + sum(mycache[["node.defn.multi"]][row.num, ]) - num.na) * log(ncol(data.df.grouping)) # num.na: number of "unused" variables determinded in removing columns if rank deficient
      fit_mdl <- fit_bic + (attributes(logLik(fit))$df-1) * log(ncol(data.df.grouping)) # logLik returns object of degrees of freedom which corresponds to the models parameters +1

      c(fit_loglik, fit_aic, fit_bic, fit_mdl)
    } else if(is.null(fit)){
      if (verbose) {
        message("Failed to fit local model. Returning very low scores.")
      }
      # no convergence, singularity, rank-deficiency, return very low scores
      c(rep(-Inf, 4))
    } else {
      stop("Unknown state of fit. I should never end up here.")
    }
  } else if (is.null(group.var)) {
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
        X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      } else {
        X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      }
    } else {
      if ("multinomial" %in% data.dists[as.logical(mycache$node.defn.adj[row.num, ])]) {
        X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
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
               suppressWarnings(fit <- irls_gaussian_cpp_fast(A = X, b = Y, maxit = control[["max.iters"]], tol = control[["tol"]]))
             }, poisson = {
               suppressWarnings(fit <- irls_poisson_cpp_fast(A = X, b = Y, maxit = control[["max.iters"]], tol = control[["tol"]]))
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
  } else {
    stop("Invalid `group.var`.")
  }
} # End of forLoopContent()

#' \code{buildScoreCache.mle} and \code{buildScoreCache.bayes} are internal functions called by \code{buildScoreCache}.
#'
#' @describeIn buildScoreCache Fit a given DAG to data with method="mle".
#' @param force.method "notset", "INLA" or "C". This is specified in \code{\link{buildScoreCache}(control=list(max.mode.error=...))}.
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom utils combn
#' @importFrom doParallel registerDoParallel
buildScoreCache.mle <-
  function(data.df = NULL,
           data.dists = NULL,
           max.parents = NULL,
           adj.vars = NULL,
           cor.vars = NULL,
           dag.banned = NULL,
           dag.retained = NULL,
           which.nodes = NULL,
           centre = TRUE,
           defn.res = NULL,
           dry.run = FALSE,
           verbose = FALSE,
           debugging = FALSE,
           force.method = NULL,
           group.var = NULL,
           grouped.vars = NULL,
           group.ids = NULL,
           control = build.control(method = "mle")) {

    set.seed(control[["seed"]])

    ## which.nodes
    if (!is.null(which.nodes)) {
        data.df <- data.df[, which.nodes]
        data.dists <- data.dists[which.nodes]
    }

    ## number of variable:
    nvars <- length(data.dists)
    nobs <- dim(data.df)[1]

    # formating factor
    data.df.lvl <- data.df

    ## standardize gaussian variables to zero mean and sd=1 have at least one gaussian variable
    if (centre && !is.null(data.dists == "gaussian")) {
        for (i in names(data.dists)[(data.dists == "gaussian")]) {
            data.df[, i] <- (data.df[, i] - mean(data.df[, i]))/sd(data.df[, i])
        }
    }

    for (i in 1:nvars) {
        if (data.dists[[i]] == "binomial" & !inherits(data.df[, i], "numeric")) {
            data.df[, i] <- as.numeric(factor(data.df[, i])) - 1
        }
        if (data.dists[[i]] == "multinomial") {
            data.df[, i] <- factor(data.df[, i])
        }
    }

    # adjustment: storing of df
    if (!is.null(adj.vars)) {
        data.df.adj <- data.df
        data.df <- data.df[, -adj.vars]
        nvars <- nvars - length(adj.vars)
    }

    ############################## Function to create the cache


    if (!is.null(defn.res)) {
        max.parents <- max(apply(defn.res[["node.defn"]], 1, sum))

    } else {
        ## Computing the cache
        fun.return <- function(x, n) {
          # x: vector or single integer
          # n: number of variables/columns in data.df
          # returns: vector of size n-1 with 1 for all x and zeros otherwise
            v <- rep(0, n - 1)
            v[x] <- 1
            return(v)
        }

        node.defn <- matrix(data = as.integer(0), nrow = 1L, ncol = nvars)
        children <- 1

        for (j in 1:nvars) {
            if (j != 1) {
                node.defn <- rbind(node.defn, matrix(data = as.integer(0),
                                                     nrow = 1L, ncol = nvars))
                children <- cbind(children, j)
            }
            # node.defn <- rbind(node.defn,matrix(data = 0,nrow = 1,ncol = n))

          if(is.list(max.parents)){
            stop("ISSUE: `max.parents` as list is not yet implemented further down here. Try with a single numeric value as max.parents instead.")
            if(!is.null(which.nodes)){
              stop("ISSUE: `max.parents` as list in combination with `which.nodes` is not yet implemented further down here. Try with single numeric as max.parents instead.")
            }
          } else if (is.numeric(max.parents) && length(max.parents)>1){
            if (length(unique(max.parents)) == 1){
              max.parents <- unique(max.parents)
            } else {
              stop("ISSUE: `max.parents` with node specific values that are not all the same, is not yet implemented further down here.")
            }
          }

          if(max.parents == nvars){
            max.parents <- max.parents-1
            warning(paste("`max.par` == no. of variables. I set it to (no. of variables - 1)=", max.parents)) #NOTE: This might cause differences to method="bayes"!
          }

            for (i in 1:(max.parents)) {
                tmp <- t(combn(x = (nvars - 1), m = i, FUN = fun.return, n = nvars, simplify = TRUE))
                tmp <- t(apply(X = tmp, MARGIN = 1, FUN = function(x) append(x = x, values = 0, after = j - 1)))

                node.defn <- rbind(node.defn, tmp)

                # children position
                children <- cbind(children, t(rep(j, length(tmp[, 1]))))
            }
        }

        # children <- rowSums(node.defn)
        colnames(node.defn) <- colnames(data.df)
        ## Coerce numeric matrix into integer matrix !!!
        node.defn <- apply(node.defn, c(1, 2), function(x) {
            (as.integer(x))
        })

        children <- as.integer(children)
        # node.defn_ <- node.defn

        ## DAG RETAIN/BANNED
        for (i in 1:nvars) {
            for (j in 1:nvars) {

                ## DAG RETAIN
                if (dag.retained[i, j] != 0) {
                  tmp.indices <- which(children == i & node.defn[, j] == 0)

                  if (length(tmp.indices) != 0) {
                    node.defn <- node.defn[-tmp.indices, ]
                    children <- children[-tmp.indices]
                  }
                }

                ## DAG BANNED
                if (dag.banned[i, j] != 0) {
                  tmp.indices <- which(children == i & node.defn[, j] == 1)

                  if (length(tmp.indices) != 0) {
                    node.defn <- node.defn[-tmp.indices, ]
                    children <- children[-tmp.indices]
                  }
                }

            }
        }

        mycache <- list(children = as.integer(children), node.defn = (node.defn))

        ###------------------------------###
        ### start limiting max.parent list###
        ###------------------------------###

        ###FIXME
        if (is.list(max.parents)) {
            for (z in 1:nvars) {
                tmp <- mycache[["node.defn"]][mycache[["children"]] == z, ]
                if (is.null(dim(tmp))) stop("Increase parents for node ",z," (due to retain)")

                if (any(diff(unlist(max.parents)) !=0))
                    stop("For method='mle', unique number of parents required")
                mycache[["node.defn"]][mycache[["children"]] == z, ] <- tmp[rowSums(tmp) <= unlist(max.parents[z]), ]
            }
        }


        ###----------------###
        ### start adjustment###
        ###----------------###

        if (!is.null(adj.vars)) {

            # mycache$node.defn.adj <- mycache$node.defn

            ## adding adjustment column set to zero mycache$node.defn.adj <- cbind(mycache$node.defn,matrix(data = 0,nrow = dim(mycache$node.defn)[1],ncol = length(adj.vars)))
            mycache$node.defn <- cbind(mycache$node.defn, matrix(data = 0, nrow = dim(mycache$node.defn)[1], ncol = length(adj.vars)))

            if (is.null(cor.vars)) {
                cor.vars <- colnames(data.df)
            }

            ## adjustment variables

            mycache$node.defn[mycache$children[match(cor.vars, colnames(data.df))], dim(data.df)[2]:dim(data.df.adj)] <- 1

            ## output
            colnames(mycache$node.defn) <- c(colnames(data.df), adj.vars)

            mycache$node.defn <- mycache$node.defn[, names(data.df.adj)]
            data.df <- data.df.adj
        }

        ##----------------------
        ## multinomial adaptation
        ##----------------------

        # unpacking the multinomial variables in the cache
        repetition.multi <- vector(length = nvars)

        for (i in 1:nvars) {
            if (data.dists[[i]] %in% c("binomial", "poisson", "gaussian")) {
                repetition.multi[i] <- 1
            } else {
                repetition.multi[i] <- nlevels(data.df.lvl[, i])
            }
        }

        if (!is.null(adj.vars)) {
            mycache$node.defn.multi <- mycache$node.defn.adj[, rep(1:nvars, repetition.multi)]
            data.df <- data.df.adj[, colnames(mycache$node.defn.adj)]
        } else {
            mycache$node.defn.multi <- mycache$node.defn[, rep(1:nvars, repetition.multi)]

        }

        # unpacking the multinomial variables in the data.df

        data.df.multi <- NULL

        for (i in 1:nvars) {
            if (data.dists[[i]] %in% c("binomial", "poisson", "gaussian")) {
                data.df.multi <- as.data.frame(cbind(data.df.multi, data.df[, i]))
                colnames(data.df.multi)[length(colnames(data.df.multi))] <- colnames(data.df)[i]
            } else {
                tmp <- model.matrix(~-1 + factor(data.df.lvl[, i]))
                colnames(tmp) <- paste0(names(data.df.lvl)[i], levels(factor(data.df.lvl[, i])))
                data.df.multi <- as.data.frame(cbind(data.df.multi, tmp))
            }
        }

    }
    if (dry.run) {
        return(mycache)
    }

    ## EOF cache creation
    row.num <- NULL   # To avoid check comment: 'no visible binding for global variable
    out <- list()
    rows <- length(mycache[["children"]])

    ##-----------------------------
    ##start loop for the regression
    ##-----------------------------
    if(verbose){cat("Start estimation loop.")}
    if(debugging){
      res <- matrix(nrow = rows, ncol = 4) # each score in one column
      for (i in 1:rows){
        # for each child ~ parents fit regression model
        res[i, ] <- forLoopContent(row.num = i,
                                   mycache = mycache,
                                   data.dists = data.dists,
                                   data.df.multi = data.df.multi,
                                   adj.vars = adj.vars,
                                   data.df = data.df,
                                   data.df.lvl = data.df.lvl,
                                   group.var = group.var,
                                   group.ids = group.ids,
                                   control = control,
                                   n = nvars,
                                   verbose = verbose)
      }
    } else {
      # no debugging

      # Prepare multithreading
      ncores <- control[["ncores"]]

      if (ncores > 1) {
        if (verbose){
          path <- path.expand(paste0(getwd(), "/build_score_cache_mle.out"))
          message(paste("Writing cluster output to: ", path))
          if(file.exists(path)){
            file.remove(path)
            message(paste("File exists and will be overwritten:", path))
          }

          cl <- makeCluster(ncores,
                            type = control[["cluster.type"]],
                            rscript_args = "--no-environ", # only available for "FORK"
                            outfile=path)
        } else {
          # no redirection
          cl <- makeCluster(ncores,
                            type = control[["cluster.type"]],
                            rscript_args = "--no-environ" # only available for "FORK"
                            )
        }

        registerDoParallel(cl)

        # NOTE for development: Make sure the recent package version is installed in the .libPath() available to the nodes.
        res <- foreach(row.num = 1:rows,
                       .combine='rbind',
                       .packages = c("stats", "lme4", "mclogit", "nnet"),
                       .export = 'forLoopContent',
                       .verbose = verbose) %dopar% {
                         forLoopContent(row.num = row.num,
                                        mycache = mycache,
                                        data.dists = data.dists,
                                        data.df.multi = data.df.multi,
                                        adj.vars = adj.vars,
                                        data.df = data.df,
                                        data.df.lvl = data.df.lvl,
                                        group.var = group.var,
                                        group.ids = group.ids,
                                        control = control,
                                        n = nvars,
                                        verbose = verbose)
                       }
        # clean up multithreading
        stopCluster(cl)

      } else {
        res <- foreach(row.num = 1:rows,
                       .combine='rbind',
                       .packages = c("stats", "lme4", "mclogit", "nnet"),
                       .export = 'forLoopContent',
                       .verbose = verbose) %do% {
                         forLoopContent(row.num = row.num,
                                        mycache = mycache,
                                        data.dists = data.dists,
                                        data.df.multi = data.df.multi,
                                        adj.vars = adj.vars,
                                        data.df = data.df,
                                        data.df.lvl = data.df.lvl,
                                        group.var = group.var,
                                        group.ids = group.ids,
                                        control = control,
                                        n = nvars,
                                        verbose = verbose)
                       }
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
    out[["dag.retained"]] <- dag.retained
    out[["dag.banned"]] <- dag.banned
    out[["group.var"]] <- group.var
    out[["group.ids"]] <- group.ids
    out[["grouped.vars"]] <- grouped.vars
    out[["aic"]] <- as.numeric( res[,2] )
    out[["bic"]] <- as.numeric( res[,3] )
    out[["mdl"]] <- as.numeric( res[,4] )

    out[["method"]] <- "mle"

    return(out)
}

