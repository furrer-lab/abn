#' @describeIn fitAbn Internal function called by \code{fitAbn}.
#' @param grouped.vars result returned from \code{\link{check.valid.groups}}. Column indexes of all variables which are affected from grouping effect.
#' @param group.ids result returned from \code{\link{check.valid.groups}}. Vector of group allocation for each observation (row) in 'data.df'.
#' @importFrom stats sd model.matrix
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
fitAbn.mle <- function(dag = NULL,
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
                       debugging = FALSE) {
  set.seed(control[["seed"]])

  ## number of variable:
  nvars <- length(data.dists)
  nobs <- dim(data.df)[1]

  ## standardize gaussian variables to zero mean and sd=1
  if(centre && !is.null(data.dists=="gaussian")){## have at least one gaussian variable
    for(i in names(data.dists)[(data.dists=="gaussian")]){data.df[,i] <- (data.df[,i]-mean(data.df[,i]))/sd(data.df[,i]);}
  }

  ##formatting
  for(i in 1:nvars){
    if(data.dists[[i]]=="binomial" & !inherits(data.df[, i], "numeric")){
      data.df[,i] <- as.numeric(factor(data.df[,i]))-1
    }
    if(data.dists[[i]]=="multinomial"){
      data.df[,i] <- factor(data.df[,i])
    }
  }

  ###----------------###
  ###start adjustment###
  ###----------------###
  if(!is.null(adj.vars)){
    if(is.null(cor.vars)){
      cor.vars <- colnames(data.df)
      cor.vars <- cor.vars[-adj.vars]
    }
    dag[cor.vars,adj.vars] <- 1
  }

  #unpacking the multinomial variables in the data.df
  data.df.multi <- NULL
  for(i in 1:nvars){
    if(data.dists[[i]] %in% c("binomial", "poisson", "gaussian")){
      data.df.multi <- as.data.frame(cbind(data.df.multi,data.df[,i]))
      colnames(data.df.multi)[length(colnames(data.df.multi))] <- colnames(data.df)[i]
    }else{
      tmp <- model.matrix(~-1+factor(data.df[,i]))
      colnames(tmp) <- paste0(names(data.df)[i], levels(factor(data.df[,i])))
      data.df.multi <- as.data.frame(cbind(data.df.multi, tmp))}
  }
  #extend dag to multinomial variables
  repetition.multi <- vector(length = nvars)
  for(i in 1:nvars){
    if(data.dists[[i]] %in% c("binomial", "poisson", "gaussian")){
      repetition.multi[i] <- 1
    }else{
      repetition.multi[i] <- nlevels(data.df[,i])
    }
  }
  dag.multi <- dag[,rep(1:nvars, repetition.multi)]

  ##-----------------------------
  ##start loop for the regression
  ##-----------------------------
  if(verbose){cat("Start estimation loop.")}
  if(debugging){
    res <- list()
    for (i in 1:length(dag[1,])){
      # for each child ~ parents fit regression model
      res[[i]] <- regressionLoop(i = i,
                               dag = dag,
                               data.df = data.df,
                               data.df.multi = data.df.multi,
                               data.dists = data.dists,
                               group.var = group.var,
                               grouped.vars = grouped.vars,
                               group.ids = group.ids,
                               control = control,
                               nvars = nvars,
                               nobs = nobs,
                               dag.multi = dag.multi,
                               verbose = verbose)
    }
  } else {
    # no debugging
    ncores <- control[["ncores"]]
    if (ncores > 1) {
      if (verbose){
        path <- path.expand(paste0(getwd(), "/fitAbn.out"))
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

      res <- foreach(i = 1:length(dag[1,]),
                     .packages = c("stats", "methods", "lme4", "mclogit", "nnet"),
                     .export = 'regressionLoop',
                     .verbose = verbose) %dopar% {
                       # for each child ~ parents fit regression model
                        regressionLoop(i = i,
                                       dag = dag,
                                       data.df = data.df,
                                       data.df.multi = data.df.multi,
                                       data.dists = data.dists,
                                       group.var = group.var,
                                       grouped.vars = grouped.vars,
                                       group.ids = group.ids,
                                       control = control,
                                       nvars = nvars,
                                       nobs = nobs,
                                       dag.multi = dag.multi,
                                       verbose = verbose)
                      }
    } else {
      res <- foreach(i = 1:length(dag[1,]),
                     .export = 'regressionLoop',
                     .verbose = verbose) %do% {
                       # for each child ~ parents fit regression model
                       regressionLoop(i = i,
                                   dag = dag,
                                   data.df = data.df,
                                   data.df.multi = data.df.multi,
                                   data.dists = data.dists,
                                   group.var = group.var,
                                   grouped.vars = grouped.vars,
                                   group.ids = group.ids,
                                   control = control,
                                   nvars = nvars,
                                   nobs = nobs,
                                   dag.multi = dag.multi,
                                   verbose = verbose)
                     }
    }
  }

  ## Prepare regressionLoop result for returning
  out <- list()
  for (childno in 1:nrow(dag)){
    out[["mliknode"]][childno] <- res[[childno]]$mliknode
    out[["mlik"]][childno] <- res[[childno]]$mlik
    out[["aicnode"]][childno] <- res[[childno]]$aicnode
    out[["aic"]][childno] <- res[[childno]]$aic
    out[["bicnode"]][childno] <- res[[childno]]$bicnode
    out[["bic"]][childno] <- res[[childno]]$bic
    out[["mdlnode"]][childno] <- res[[childno]]$mdlnode
    out[["mdl"]][childno] <- res[[childno]]$mdl
    out[["df"]][childno] <- tryCatch({
      res[[childno]]$df
    }, error=function(e){NA})
    out[["sse"]][childno] <- res[[childno]]$sse
    out[["mse"]][childno] <- tryCatch({
      res[[childno]]$mse
    }, error=function(e){NA})

    if(!is.null(group.var)){
      out[["mu"]][childno] <- tryCatch({
        list(res[[childno]]$mu)
      }, error=function(e){NA})
      out[["betas"]][childno] <- tryCatch({
        list(res[[childno]]$betas)
      }, error=function(e){NA})
      out[["sigma"]][childno] <- tryCatch({
        list(res[[childno]]$sigma)
      }, error=function(e){NA})
      out[["sigma_alpha"]][childno] <- tryCatch({
        list(res[[childno]]$sigma_alpha)
      }, error=function(e){NA})
    } else if (is.null(group.var)){
      out[["coef"]][[childno]] <- tryCatch({
        res[[childno]]$coef
      }, error=function(e){NA})
      out[["Stderror"]][[childno]] <-  tryCatch({
        res[[childno]]$var
      }, error=function(e){NA})
      out[["pvalue"]][[childno]] <- tryCatch({
        res[[childno]]$pvalue
      }, error=function(e){NA})
    } else {
      stop("Invalid `group.var`.")
    }
  }

  for (i in 1:length(out)){
    names(out[[i]]) <- names(data.df)
  }

  # Reduce node-level score to an overall model score
  out[["mlik"]] <- sum(out[["mliknode"]], na.rm = TRUE)
  out[["aic"]] <- sum(out[["aicnode"]], na.rm = TRUE)
  out[["bic"]] <- sum(out[["bicnode"]], na.rm = TRUE)

  out[["method"]] <- "mle"
  out[["abnDag"]] <- createAbnDag(dag, data.df = data.df, data.dists = data.dists)
  if(!is.null(group.var)){
    out[["group.var"]] <- group.var
    out[["group.ids"]] <- group.ids
    out[["grouped.vars"]] <- grouped.vars
  }
  return(out)
}

#' Regress each child on its parents in fitAbn()
#' @describeIn fitAbn Internal function called by \code{fitAbn.mle()}.
#' @param i number of child-node (mostly corresponds to child node index e.g. in dag).
#' @param nvars number of variables in data.dists.
#' @param nobs number of observations in data.df.
#' @param data.df.multi extended data.df for one-hot-encoded multinomial variables.
#' @param dag.multi extended dag for one-hot-encoded multinomial variables.
#' @importFrom stats coefficients coef residuals df.residual vcov AIC BIC sd model.matrix as.formula lm glm logLik
#' @importFrom methods is
#' @return list
#' @keywords internal
regressionLoop <- function(i = NULL, # number of child-node (mostly corresponds to child node index e.g. in dag)
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
                           verbose = NULL){

  child.dist <- data.dists[[i]] # current nodes/response variable's distribution
  parent.dists <- data.dists[as.logical(dag[i,])] # distributions of parent nodes from current child node/response
  child.name <- names(data.df)[i] # current nodes/response variable's name
  parents.names <- names(data.df)[as.logical(dag[i,])] # names of parent nodes from current child node/response
  parents.names.multi <- names(data.df.multi)[as.logical(dag.multi[i,])] # names of parent nodes from current child node/response
  USED_NNET <- FALSE # If TRUE, naming of output is done slightly different.
  USED_MBLOGIT <- FALSE # If TRUE, naming of output is done slightly different.

  if (!is.null(group.var)) {
    # we have grouping (mixed-effects model) GLMM

    ###
    # Build arguments for GLMM
    ###
    # data.df including group.var=group.ids
    data.df.grouping <- data.frame(data.df, grp = factor(group.ids))
    names(data.df.grouping) <- c(names(data.df), group.var)

    ###
    # Begin of GLMM
    ###
    # Build formula statement
    if (length(parents.names) == 0){
      model <- as.formula(paste(child.name, "~ 1 + (1|", group.var, ")"))
    } else if (length(parents.names) > 0){
      model <- as.formula(paste(child.name, "~ (1|", group.var, ")+", paste(parents.names, collapse = "+")))
    } else {
      stop("Cannot build model formula. Unknown predictor names of parent nodes in 'parents.name'.")
    }

    if(verbose) message(paste("\n", deparse1(model), "\n")) else NA

    # Main part GLMM: Depending on child's distribution, call the appropriate modelling function
    switch (as.character(child.dist),
            gaussian = {
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

                if (is.null(fit)) {
                  # if fit is still NULL, try other (all available) optimizer:
                  # fit same as above (not very elegant)
                  fit <- lme4::lmer(model, data = data.df.grouping,
                                    control = lme4::lmerControl(check.rankX = control[["check.rankX"]],
                                                                check.scaleX = control[["check.scaleX"]],
                                                                check.conv.grad = lme4::.makeCC(action = control[["check.conv.grad"]]),
                                                                check.conv.singular = lme4::.makeCC(action = control[["check.conv.singular"]]),
                                                                check.conv.hess = lme4::.makeCC(action = control[["check.conv.hess"]]),
                                                                optCtrl = list(xtol_abs = control[["xtol_abs"]],
                                                                               ftol = control[["ftol_abs"]])))
                  # refit with all available optimizers
                  fit_all <- lme4::allFit(fit)
                  # keep only results from optimizers that were reported as OK
                  fit_all_OK <- fit_all[sapply(fit_all, methods::is, "merMod")]
                  # extract messages from each optimizer
                  convergence_results <- lapply(fit_all_OK, function(x) x@optinfo$conv$lme4$messages)
                  # Get only results from optimizer without any message (meaning they did converge etc.)
                  converged_idx <- sapply(convergence_results, is.null)
                  # Keep only results from first optimizer (arbitrarily) that did converge without any message
                  if(sum(converged_idx) == 0){
                    if(verbose){message("No algorithms from allFit converged. You may still be able to use the results, but proceed with extreme caution.")}
                    fit <- NULL
                  } else {
                    fit <- fit_all[converged_idx][[1]]
                  }
                }
            },
            binomial = {
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

                if (is.null(fit)) {
                  # if fit is still NULL, try other (all available) optimizer:
                  # fit same as above (not very elegant)
                  fit <- lme4::glmer(model, data = data.df.grouping, family = "binomial",
                                     control = lme4::glmerControl(tolPwrss = control[["tolPwrss"]],
                                                                  check.rankX = control[["check.rankX"]],
                                                                  check.scaleX = control[["check.scaleX"]],
                                                                  check.conv.grad = lme4::.makeCC(action = control[["check.conv.grad"]]),
                                                                  check.conv.singular = lme4::.makeCC(action = control[["check.conv.singular"]]),
                                                                  check.conv.hess = lme4::.makeCC(action = control[["check.conv.hess"]]),
                                                                  optCtrl = list(xtol_abs = control[["xtol_abs"]],
                                                                                 ftol = control[["ftol_abs"]])))
                  # refit with all available optimizers
                  fit_all <- lme4::allFit(fit)
                  # keep only results from optimizers that were reported as OK
                  fit_all_OK <- fit_all[sapply(fit_all, methods::is, "merMod")]
                  # extract messages from each optimizer
                  convergence_results <- lapply(fit_all_OK, function(x) x@optinfo$conv$lme4$messages)
                  # Get only results from optimizer without any message (meaning they did converge etc.)
                  converged_idx <- sapply(convergence_results, is.null)
                  # Keep only results from first optimizer (arbitrarily) that did converge without any message
                  if(sum(converged_idx) == 0){
                    if(verbose){message("No algorithms from allFit converged. You may still be able to use the results, but proceed with extreme caution.")}
                    fit <- NULL
                  } else {
                    fit <- fit_all[converged_idx][[1]]
                  }
                }
            },
            poisson = {
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

                if (is.null(fit)) {
                  # if fit is still NULL, try other (all available) optimizer:
                  # fit same as above (not very elegant)
                  fit <- lme4::glmer(model, data = data.df.grouping, family = "poisson",
                                     control = lme4::glmerControl(tolPwrss = control[["tolPwrss"]],
                                                                  check.rankX = control[["check.rankX"]],
                                                                  check.scaleX = control[["check.scaleX"]],
                                                                  check.conv.grad = lme4::.makeCC(action = control[["check.conv.grad"]]),
                                                                  check.conv.singular = lme4::.makeCC(action = control[["check.conv.singular"]]),
                                                                  check.conv.hess = lme4::.makeCC(action = control[["check.conv.hess"]]),
                                                                  optCtrl = list(xtol_abs = control[["xtol_abs"]],
                                                                                 ftol = control[["ftol_abs"]])))
                  # refit with all available optimizers
                  fit_all <- lme4::allFit(fit)
                  # keep only results from optimizers that were reported as OK
                  fit_all_OK <- fit_all[sapply(fit_all, methods::is, "merMod")]
                  # extract messages from each optimizer
                  convergence_results <- lapply(fit_all_OK, function(x) x@optinfo$conv$lme4$messages)
                  # Get only results from optimizer without any message (meaning they did converge etc.)
                  converged_idx <- sapply(convergence_results, is.null)
                  # Keep only results from first optimizer (arbitrarily) that did converge without any message
                  if(sum(converged_idx) == 0){
                    if(verbose){message("No algorithms from allFit converged. You may still be able to use the results, but proceed with extreme caution.")}
                    fit <- NULL
                  } else {
                    fit <- fit_all[converged_idx][[1]]
                  }
                }
            },
            multinomial = {
              if (length(parents.names) == 0){
                model_basic <- as.formula(paste(child.name, "~ 1"))
              } else {
                model_basic <- as.formula(paste(child.name, "~ ", paste(parents.names, collapse = "+")))
              }
              model_random <- as.formula(paste("~ 1|", group.var, sep = ""))
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
                  tryCatch({fit <- mclogit::mblogit(formula = model_basic, random = model_random, data = data.df.grouping,
                                                    control = mclogit::mclogit.control(epsilon = control[["epsilon"]],
                                                                                       trace = control[["trace.mblogit"]]))
                  }, error=function(e)NULL)
                }
                if (is.null(fit)) {
                  # if fit is still NULL, try other (all available) optimizer:
                  # fit same as above (not very elegant)
                  fit <- lme4::lmer(model, data = data.df.grouping,
                                    control = lme4::lmerControl(optCtrl = list(xtol_abs=1e-6, ftol_abs=1e-6)))
                  # refit with all available optimizers
                  fit_all <- lme4::allFit(fit)
                  # keep only results from optimizers that were reported as OK
                  fit_all_OK <- fit_all[sapply(fit_all, methods::is, "merMod")]
                  # extract messages from each optimizer
                  convergence_results <- lapply(fit_all_OK, function(x) x@optinfo$conv$lme4$messages)
                  # Get only results from optimizer without any message (meaning they did converge etc.)
                  converged_idx <- sapply(convergence_results, is.null)
                  # Keep only results from first optimizer (arbitrarily) that did converge without any message
                  if(sum(converged_idx) == 0){
                    if(verbose){message("No algorithms from allFit converged. You may still be able to use the results, but proceed with extreme caution.")}
                    fit <- NULL
                  } else {
                    fit <- fit_all[converged_idx][[1]]
                  }
                }
              USED_MBLOGIT <- TRUE
            }
    )
    ###
    # End of GLMM
    ###
  } else if (is.null(group.var)) {
    # we have no grouping (do what was always done).

    # Separate outcome and predictor variables
    Y <- data.matrix(data.df[,i])

    if("multinomial" %in% data.dists[as.logical(dag[i,])]){
      # if predictor is multinomial, use one-hot-encoded categories
      X <- data.matrix(data.df.multi[,as.logical(dag.multi[i,])])
    }else{
      # predictor is not multinomial, keep predictors as is
      X <- data.matrix(cbind(rep(1,length(data.df[,1])),data.df[,as.logical(dag[i,])]))
    }

    num.na <- 0

    # Check for collinearity
    if(qr(X)$rank/ncol(X)!=1 & as.character(data.dists[i])=="binomial"){
      # if collinear, iteratively remove last variable from predictors and retry with internal IRLS.
      Y <- as.numeric(as.character(Y))
      fit <- tryCatch(irls_binomial_cpp_br(A = X, b = Y, maxit = control[["max.irls"]],tol = control[["tol"]]),error=function(e){
        while((qr(X)$rank/ncol(X))!=1){
          X <- X[,-1]
          num.na <- num.na+1
          if(is.null(ncol(X))) X <- as.matrix(X)
        }
        list(irls_binomial_cpp_br(A = X, b = Y, maxit = control[["max.irls"]],tol = control[["tol"]]),num.na)
      })
      num.na <- fit[[2]]
      fit <- fit[[1]]
    }else{
      # no collinearity detected.
      switch(as.character(data.dists[i]),
             "binomial"={
               Y <- as.numeric(as.character(Y))
               fit <- tryCatch(irls_binomial_cpp(A = X, b = Y, maxit = control[["max.irls"]],tol = control[["tol"]]),
                               error=function(e){
                                 irls_binomial_cpp_br(A = X, b = Y, maxit = control[["max.irls"]],tol = control[["tol"]])
                               }
               )
               if(is.na(sum(fit[[1]]))) {
                 fit <- irls_binomial_cpp_br(A = X, b = Y, maxit = control[["max.irls"]],tol = control[["tol"]])
               }
             },
             "gaussian"={
               fit <- irls_gaussian_cpp(A = X, b = Y, maxit = control[["max.irls"]],tol = control[["tol"]])
             },
             "poisson"={
               fit <- irls_poisson_cpp(A = X, b = Y, maxit = control[["max.irls"]],tol = control[["tol"]])
             },
             "multinomial"={
               tmp <- multinom(formula = Y~-1+X,Hess = FALSE,trace=FALSE)

               # Calculate scores and prepare output
               fit <- list()
               fit$coefficients <- as.matrix(as.vector(coef(tmp)))
               fit$names.coef <- row.names((coef(tmp)))
               fit$loglik <- - tmp$value
               edf <- ifelse(length(tmp$lev) == 2L, 1, length(tmp$lev) - 1) * qr(X)$rank
               fit$aic <- 2 * tmp$value + 2 * edf
               fit$bic <- 2 * tmp$value + edf * log(nobs)
               fit$sse <- sum(residuals(tmp)^2)
               fit$var.out <- as.matrix(as.vector(summary(tmp)$standard.errors))
               USED_NNET <- TRUE
             }
      )}
  } else {
    stop("Invalid `group.var`.")
  } # EOF regression loop main part

  ## Collect outputs from regLoop main part
  res <- list() # collect results for current child node
  if (!is.null(group.var)) {
    if(!is.null(fit)){
      res[["mliknode"]] <- logLik(fit)
      res[["mlik"]] <- Reduce("+", res[["mliknode"]])
      res[["aicnode"]] <- AIC(fit)
      res[["aic"]] <- Reduce("+", res[["aicnode"]])
      res[["bicnode"]] <- BIC(fit)
      res[["bic"]] <- Reduce("+", res[["bicnode"]])
      res[["mdlnode"]] <-  res[["bic"]] + (attributes(logLik(fit))$df-1) * log(ncol(data.df.grouping)) # logLik returns object of degrees of freedom which corresponds to the models parameters +1
      res[["mdl"]] <- Reduce("+", res[["mdlnode"]])
      res[["sse"]] <- sum(residuals(fit)^2) # residual sum of squares (sum_1^n(y-f(x_i)))^2 -> ||y-X*beta_hat||^2
      res[["mse"]] <- res[["sse"]] / stats::df.residual(fit) # the residual sum of squares divided by the residuals number of degrees of freedom.

      if (verbose){warning("Mixed-model degrees of freedom, coefficients, variance and p-values should be considered carefully.")}
      # Number of DFs is hotly debated (Bates et al. 2008, https://doi.org/10.1016/j.jml.2007.12.005) and situation specific.
      # MCMC or parametric bootstrapping approaches were proposed and implemented in the afex package but are computationally too heavy to work out here.
      # Therefore, no df and pvalues for GLMMs at the moment.
      res[["df"]] <- NA
      if(USED_MBLOGIT){
        # multinomial response
        res[["mu"]] <- fit$coefficients # = fit$coefmat[,1]; vector or matrix of size kx1 or kx(number of predictor variables) respectively
        res[["mu"]] <- fit$coefmat[,1]
        if(!is.na(dim(fit$coefmat)[2]) && dim(fit$coefmat)[2] >1){
          res[["betas"]] <- matrix(as.matrix(fit$coefmat)[, which(colnames(fit$coefmat) != "(Intercept)")],
                                   nrow = nrow(fit$coefmat), byrow = FALSE, # keep same shape as fit$coefmat
                                   dimnames = list(rownames(fit$coefmat), colnames(fit$coefmat)[-1]))
        } else {
          res[["betas"]] <- NA
        }
        res[["sigma"]] <- NA
        res[["sigma_alpha"]] <- as.matrix(fit$VarCov[[group.var]]) # matrix size kxk for k = (number of factor levels of response)
      } else {
        # not multinomial response
        res[["mu"]] <- lme4::fixef(fit)[[1]] # fixed-effect intercept (beta_0)
        res[["betas"]] <- lme4::fixef(fit)[-1] # fixed-effect slopes (beta_1, beta_2, ...)
        res[["sigma"]] <- as.data.frame(lme4::VarCorr(fit))$sdcor[-1] # random-effect residual
        res[["sigma_alpha"]] <- as.data.frame(lme4::VarCorr(fit))$sdcor[1] # random-effect intercept
      }
    } else if(is.null(fit)){
      # no convergence, return very low score
      stop("fit is NULL. I don't yet know how to handle this case.")
    } else {
      stop("Unknown state of fit. I should never end up here.")
    }
  } else if (is.null(group.var)) {
    res[["mliknode"]] <- fit$loglik
    res[["mlik"]] <- Reduce("+", res[["mliknode"]])
    res[["aicnode"]] <- fit$aic
    res[["aic"]] <- Reduce("+", res[["aicnode"]])
    res[["bicnode"]] <- fit$bic
    res[["bic"]] <- Reduce("+", res[["bicnode"]])
    res[["mdlnode"]] <- fit$bic + (1 + sum(dag.multi[i,] - num.na)) * log(nvars)
    res[["mdl"]] <- Reduce("+", res[["mdlnode"]])
    res[["df"]] <- (length(data.df[,1])-(sum(dag.multi[i,])+1))
    res[["sse"]] <- fit$sse  # residual sum of squares (sum_1^n(y-f(x_i)))^2 -> ||y-X*beta_hat||^2
    res[["mse"]] <- res[["sse"]]/res[["df"]] # the residual sum of squares divided by the residuals number of degrees of freedom.

    switch(as.character(data.dists[i]),
           "gaussian"={
             res[["coef"]] <- matrix(data = c(rep(NA,num.na),fit$coefficients),nrow = 1)

             tmp.catch <- tryCatch(expr = solve(t(X) %*% X),error = function(e) NaN)
             if(is.nan(tmp.catch[1])){
               res[["var"]] <- matrix(data = rep(NaN,length.out=ncol(X)),nrow = 1)
             }else{
               res[["var"]] <- matrix(data = sqrt(diag(res[["mse"]]*tmp.catch)),nrow = 1) # sqrt of diagonal of covariance matrix.
               res[["var"]] <- matrix(data = c(rep(NA,num.na),res[["var"]]),nrow = 1)
             }
           },
           "binomial"={
             res[["coef"]] <- matrix(data = c(rep(NA,num.na),fit$coefficients),nrow = 1)
             res[["var"]] <- tryCatch(as.matrix(sqrt(diag(solve(fit$varcov)))),
                                      error=function(e){
                                        as.matrix(sqrt(svd(fit$varcov)$d))
                                      }
             )
             res[["var"]] <- matrix(data = c(rep(NA,num.na),res[["var"]]),nrow = 1)
           },
           "poisson"={
             fit$coefficients <- cbind(t(rep(NA,num.na)),t(fit$coefficients))
             res[["coef"]] <- fit$coefficients

             res[["var"]] <- as.matrix(sqrt(diag(solve(fit$varcov))))
             res[["var"]] <- matrix(data = cbind(t(rep(NA,num.na)),t(res[["var"]])),nrow = 1)
           },
           "multinomial"={
             fit$coefficients <- rbind(as.matrix(rep(NA,num.na*length(fit$names.coef))),fit$coefficients)
             res[["coef"]] <- t(fit$coefficients)

             res[["var"]] <- matrix(data = rbind(as.matrix(rep(NA,num.na*length(fit$names.coef))),fit$var.out),nrow = 1)
           }
    )
  } else {
    stop("Invalid `group.var`.")
  }

  ## Naming of results
  if(!is.null(group.var)){
    # we have grouping and therefore a slightly different output of res
    if(child.dist!="multinomial"){
      # current response is not multinomial distributed and has therefore a slightly different output
      # DO NOT CHANGE NAMES HERE.
      # if("multinomial" %in% parent.dists){
      #   # we have one or more multinomial predictors
      #   colnames(res[["mu"]]) <- XX
      #   colnames(res[["betas"]]) <- XX
      #   colnames(res[["sigma"]]) <- XX
      #   colnames(res[["sigma_alpha"]]) <- XX
      # } else if (!("multinomial") %in% parent.dists){
      #   # we have no multinomial predictors
      #   colnames(res[["mu"]]) <- XX
      #   colnames(res[["betas"]]) <- XX
      #   colnames(res[["sigma"]]) <- XX
      #   colnames(res[["sigma_alpha"]]) <- XX
      # }
    } else if(child.dist=="multinomial"){
      # current response is multinomial distributed
      if("multinomial" %in% parent.dists){
        # we have one or more multinomial predictors
        # Naming convention: <childname>.<child_factorlevel_name>~<parentname>.<parent_factorlevel_name>
        names(res[["mu"]]) <- paste(child.name, names(res[["mu"]]), sep = ".")
        # colnames(res[["betas"]]) <- NA
        # colnames(res[["sigma"]]) <- NA
        names(res[["sigma_alpha"]]) <- paste(child.name, names(res[["sigma_alpha"]]), sep = ".")
      } else if (!("multinomial") %in% parent.dists){
        # we have no multinomial predictors
        # Naming convention: <childname>.<child_factorlevel_name>~<parentname>
        names(res[["mu"]]) <- paste(child.name, names(res[["mu"]]), sep = ".")
        # colnames(res[["betas"]]) <- NA
        # colnames(res[["sigma"]]) <- NA
        colnames(res[["sigma_alpha"]]) <- rownames(res[["sigma_alpha"]]) <- paste(child.name, colnames(res[["sigma_alpha"]]), sep = ".")
      }
    }
  } else if(is.null(group.var)){
    # no grouping and res has classic output format
    if(child.dist!="multinomial"){
      if("multinomial" %in% parent.dists){
        colnames(res[["coef"]]) <- c(parents.names.multi)
        colnames(res[["var"]]) <- c(parents.names.multi)
      }else{
        colnames(res[["coef"]]) <- c(paste(child.name,"|intercept",sep = ""),parents.names)
        colnames(res[["var"]]) <- c(paste(child.name,"|intercept",sep = ""),parents.names)
      }
    } else if(child.dist=="multinomial"){
      separator = ""
      if("multinomial" %in% parent.dists){
        colnames(res[["coef"]]) <- c(as.vector(outer(parents.names.multi, fit$names.coef, paste, sep=separator)))
        colnames(res[["var"]]) <- c(as.vector(outer(parents.names.multi, fit$names.coef, paste, sep=separator)))
      }else{
        if (USED_NNET){
          colnames(res[["coef"]]) <- c(paste(child.name,"|intercept.",fit$names.coef,sep = ""),as.vector(outer(parents.names, fit$names.coef, paste, sep=separator)))
          colnames(res[["var"]]) <- c(paste(child.name,"|intercept.",fit$names.coef,sep = ""),as.vector(outer(parents.names, fit$names.coef, paste, sep=separator)))
        } else if (USED_MBLOGIT){
          colnames(res[["coef"]]) <- c(paste(child.name,"|intercept.",rownames(fit$coefmat),sep = ""),as.vector(outer(parents.names, rownames(fit$coefmat), paste, sep=separator)))
          colnames(res[["var"]]) <-c(paste(child.name,"|intercept.",rownames(fit$coefmat),sep = ""),as.vector(outer(parents.names, rownames(fit$coefmat), paste, sep=separator)))
        } else {
          colnames(res[["coef"]]) <- tryCatch({
            c(paste(child.name,"|intercept.",names(stats::coefficients(fit)),sep = ""),as.vector(outer(parents.names, names(stats::coefficients(fit)), paste, sep=separator)))
          }, error=function(e){
            c(paste(child.name,"|intercept.",fit$names.coef,sep = ""),as.vector(outer(parents.names, fit$names.coef, paste, sep=separator)))
          })
          colnames(res[["var"]]) <- tryCatch({
            c(paste(child.name,"|intercept.",names(stats::coefficients(fit)),sep = ""),as.vector(outer(parents.names, names(stats::coefficients(fit)), paste, sep=separator)))
          }, error=function(e){
            c(paste(child.name,"|intercept.",fit$names.coef,sep = ""),as.vector(outer(parents.names, fit$names.coef, paste, sep=separator)))
          })
        }
      }
    }
  } else {
    stop("Invalid `group.var`.")
  }#naming

  return(res)
}
