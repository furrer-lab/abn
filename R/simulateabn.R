#' Bugs code for Gaussian response
#'
#' Bugs model for a normal distributed response variable
#' \eqn{X \sim \mathcal{N}(\mu,\,\sigma^{2})}.
#'
#' @param nodename character string of response variable name.
#' @param parentnames single character string (for one parent) or vector of characters (for multiple parent nodes) with parent node (predictor variables) names.
#' @param nodesintercept numeric intercept value for current response variable.
#' @param parentcoefs list named with parent node names of numeric vectors with coefficient values.
#' @param std integer with standard deviation of response variable that will be
#' converted to precision (see Details).
#'
#' @details
#' The variance of the normal distribution is \eqn{\frac{1}{\tau}}.
#'
#' @return Bugs model returned as stdout.
#'
#' @examples
#' gauss_bugs(nodename = "a",
#'            parentnames = c("b", "c"),
#'            nodesintercept = c(0.318077),
#'            parentcoefs = list("b"=c(b=0.3059395),
#'                               "c"=c(c=0.5555)),
#'            std = c(0.05773503))
#' @importFrom stringi stri_flatten
#' @seealso [makebugs] [simulateAbn]
#' @export
#' @keywords utilities internal
gauss_bugs <- function(nodename, nodesintercept, parentnames, parentcoefs, std){
  cat(nodename, " ~ dnorm(mu.", nodename, ", precision.", nodename, ") # Gaussian response", sep = "", fill = TRUE)
  # mu
  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    predictors <- c()
    for (i in 1:length(parentnames)){
      for (j in 1:length(parentcoefs[[i]])){
        predictors <- append(predictors, paste0(parentcoefs[[i]][j], "*", names(parentcoefs[[i]][j])))
      }
    }
    flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

    cat("mu.", nodename, " <- ", nodesintercept, " + ", flatpredictors, " # Linear regression", sep = "", fill = TRUE)
  } else {
    cat("mu.", nodename, " <- ", nodesintercept, " # Linear regression", sep = "", fill = TRUE)
  }

  # tau (precision)
  cat("precision.", nodename, " <- inverse(", std, ") # precision tau = 1/standard_dev", sep = "", fill = TRUE)
}

#' @describeIn gauss_bugs Bugs code for Gaussian response with varying intercept
#' @param nodesintercept overall mean of response. Parameter from fixed-effects intercept.
#' @param parentcoefs overall slope for each predictor (parent node) variable (fixed-effects).
#' @param sigma within-group variance. Parameter from random-effects residual.
#' @param sigma_alpha between-group variance. Parameter from random-effects intercept.
gauss_bugsGroup <- function(nodename, nodesintercept, parentnames, parentcoefs, sigma, sigma_alpha){
  ## Node as response
  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    # With parent nodes
    predictors <- c()
    for (i in 1:length(parentnames)){
      for (j in 1:length(parentcoefs[[i]])){
        predictors <- append(predictors, paste0(parentcoefs[[i]][j], "*", parentnames[i]))
      }
    }
    flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

    cat(nodename, " <- mu_", nodename, " + ", flatpredictors, " + alpha_", nodename, " + e_", nodename, sep = "", fill = TRUE)
  } else {
    # No parent nodes
    cat(nodename, " <- mu_", nodename, " + alpha_", nodename, " + e_", nodename, sep = "", fill = TRUE)
  }
  cat("mu_", nodename, " <- ", nodesintercept, sep = "", fill = TRUE)
  cat("e_", nodename, " ~ dnorm(mu_e_", nodename, ", tau_", nodename, ")", sep = "", fill = TRUE)
  cat("mu_e_", nodename, " <- 0", sep = "", fill = TRUE)
  cat("tau_", nodename, " <- inverse(sigma_", nodename, ")", sep = "", fill = TRUE)
  cat("sigma_", nodename, " <- ", sigma, sep = "", fill = TRUE)
  ## Grouping of node
  cat("alpha_", nodename, " ~ dnorm(mu_alpha_", nodename, ", tau_alpha_", nodename, ")", sep = "", fill = TRUE)
  cat("mu_alpha_", nodename, " <- 0", sep = "", fill = TRUE)
  cat("tau_alpha_", nodename, " <- inverse(sigma_alpha_", nodename, ")", sep = "", fill = TRUE)
  cat("sigma_alpha_", nodename, " <- ", sigma_alpha, sep = "", fill = TRUE)
}

#' Bugs code for Bernoulli response
#'
#' Bugs model for a Binomial response \eqn{X} in a single trial:
#' \eqn{X \sim \mathcal{B}(n=1, p) = \mathcal{Bernoulli}(p)}.
#'
#' @param nodename character string of response variable name.
#' @param parentnames single character string (for one parent) or vector of characters (for multiple parent nodes) with parent node (predictor variables) names.
#' @param nodesintercept numeric intercept value for current response variable.
#' @param parentcoefs list named with parent node names of numeric vectors with coefficient values.
#'
#' @return Bugs model returned as stdout.
#'
#' @examples
#' bern_bugs(nodename = "a",
#'           parentnames = c("b", "c"),
#'           nodesintercept = c(0.318077),
#'           parentcoefs = list("b"=c(b=0.3059395),
#'                              "c"=c(c=0.5555)))
#' @importFrom stringi stri_flatten
#' @seealso [makebugs] [simulateAbn]
#' @export
#' @keywords utilities internal
bern_bugs <- function(nodename, nodesintercept, parentnames, parentcoefs){
  cat(nodename, " ~ dbern(p.", nodename, ") # Bernoulli response", sep = "", fill = TRUE)
  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    predictors <- c()
    for (i in 1:length(parentnames)){
      for (j in 1:length(parentcoefs[[i]])){
        predictors <- append(predictors, paste0(parentcoefs[[i]][j], "*", names(parentcoefs[[i]][j])))
      }
    }
    flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

    cat("logit(p.", nodename, ") <- ", nodesintercept, " + ", flatpredictors, " # logistic regression", sep = "", fill = TRUE)
  } else {
    cat("logit(p.", nodename, ") <- ", nodesintercept, " # logistic regression", sep = "", fill = TRUE)
  }
}

#' @describeIn bern_bugs Bugs code for Bernoulli response with varying intercept
#' @param nodesintercept overall mean of response. Parameter from fixed-effects intercept.
#' @param parentcoefs overall slope for each predictor (parent node) variable (fixed-effects).
#' @param sigma_alpha between-group variance. Parameter from random-effects intercept.
bern_bugsGroup <- function(nodename, nodesintercept, parentnames, parentcoefs, sigma_alpha){
  ## Node as response
  cat(nodename, " ~ dbern(p_", nodename, ")", sep = "", fill = TRUE)
  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    # With parent nodes
    predictors <- c()
    for (i in 1:length(parentnames)){
      for (j in 1:length(parentcoefs[[i]])){
        predictors <- append(predictors, paste0(parentcoefs[[i]][j], "*", parentnames[i]))
      }
    }
    flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

    cat("logit(p_", nodename, ") <- mu_", nodename, " + ", flatpredictors, " + alpha_", nodename, sep = "", fill = TRUE)
  } else {
    # No parent nodes
    cat("logit(p_", nodename, ") <- mu_", nodename, " + alpha_", nodename, sep = "", fill = TRUE)
  }
  cat("mu_", nodename, " <- ", nodesintercept, sep = "", fill = TRUE)
  ## Grouping of node
  cat("alpha_", nodename, " ~ dnorm(mu_alpha_", nodename, ", tau_alpha_", nodename, ")", sep = "", fill = TRUE)
  cat("mu_alpha_", nodename, " <- 0", sep = "", fill = TRUE)
  cat("tau_alpha_", nodename, " <- inverse(sigma_alpha_", nodename, ")", sep = "", fill = TRUE)
  cat("sigma_alpha_", nodename, " <- ", sigma_alpha, sep = "", fill = TRUE)
}

#' Bugs code for Categorical response
#'
#' @param nodename character string of response variable name.
#' @param nodesCatIdx integer vector of length \eqn{|K-1|} and starting at \eqn{k+1} (see Examples).
#' @param parentnames single character string (for one parent) or vector of characters (for multiple parent nodes) with parent node (predictor variables) names.
#' @param nodesintercepts numeric vector with the intercept value for each category (see Details).
#' @param parentcoefs list named with parent node names of numeric vectors with coefficient values for each category (see Details).
#'
#' @details
#' The output of \code{\link[abn]{fitAbn}} with \code{method = "mle"} is based on
#' the output of logistic regression models fit with either \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}}, \code{\link[lme4]{glmer}}, \code{\link[nnet]{multinom}},
#' \code{\link[mclogit]{mblogit}} or internal \code{irls} methods.
#' They all use the first factor level as reference level.
#' Therefore, \code{nodesCatIdx} starts with index \eqn{2} and not \eqn{1}.
#' \code{nodesintercepts} and \code{parentcoefs} refer to the values of
#' \code{(Intercept)} and \code{Estimate} of the respective model output.
#' Predictor names build the keys in \code{parentcoef}.
#'
#' @return Bugs model returned as stdout.
#'
#' @examples
#' # A -> B
#' # Where B is a categorical variable with 4 levels.
#' categorical_bugs(nodename = "b",
#'                  nodesCatIdx = c(2, 3, 4),
#'                  parentnames = "a",
#'                  nodesintercepts = c(2.188650, 3.133928, 3.138531),
#'                  parentcoefs = list("a"=c(a=1.686432, a=3.134161, a=5.052104)))
#' @importFrom stringi stri_flatten
#' @seealso [makebugs] [simulateAbn]
#' @export
#' @keywords utilities internal
categorical_bugs <- function(nodename, nodesCatIdx, parentnames, nodesintercepts, parentcoefs){
  cat(nodename, " ~ dcat(p.", nodename, ") # Categorical response", sep = "", fill = TRUE)
  cat("p.", nodename, "[1] <- phi.", nodename, "[1]/sum(phi.", nodename, ") # soft-max", sep = "", fill = TRUE)
  cat("log(phi.", nodename, "[1]) <- 0 # Reference category", sep = "", fill = TRUE)

  for(i in 1:length(nodesCatIdx)){
    cat("p.", nodename, "[", nodesCatIdx[i], "] <- phi.", nodename, "[", nodesCatIdx[i], "]/sum(phi.", nodename, ") # soft-max", sep = "", fill = TRUE)
    if (!any(is.na(parentnames)) && length(parentnames) >0){
      predictors <- c()
      for (j in 1:length(parentnames)){
        predictors[j] <- paste0(parentcoefs[[j]][i], "*", parentnames[j])
      }
      flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

      cat("log(phi.", nodename, "[", nodesCatIdx[i], "]) <- ", nodesintercepts[i], " + ", flatpredictors, sep = "", fill = TRUE)
    } else {
      # no parent nodes
      cat("log(phi.", nodename, "[", nodesCatIdx[i], "]) <- ", nodesintercepts[i], sep = "", fill = TRUE)
    }
  }
}

#' @describeIn categorical_bugs Bugs code for Categorical response with varying intercept
#' @param nodesintercepts overall mean of response. Parameter from fixed-effects intercept.
#' @param parentcoefs overall slope for each predictor (parent node) variable (fixed-effects).
#' @param sigma within-group variance. Parameter from random-effects residual.
#' @param sigma_alpha between-group variance-covariance matrix. Parameters from random-effects intercept.
categorical_bugsGroup <- function(nodename, nodesCatIdx, nodesintercepts, parentnames, parentcoefs, sigma, sigma_alpha){
  cat(nodename, " ~ dcat(p_", nodename, ")", sep = "", fill = TRUE)
  cat("p_", nodename, "[1] <- phi_", nodename, "[1]/sum(phi_", nodename, ")", sep = "", fill = TRUE)
  cat("log(phi_", nodename, "[1]) <- 0 + alpha_", nodename, "[1]", sep = "", fill = TRUE)

  for(i in 1:length(nodesCatIdx)){
    cat("p_", nodename, "[", nodesCatIdx[i], "] <- phi_", nodename, "[", nodesCatIdx[i], "]/sum(phi_", nodename, ")", sep = "", fill = TRUE)
    if (!any(is.na(parentnames)) && length(parentnames) >0){
      predictors <- c()
      for (j in 1:length(parentnames)){
        predictors[j] <- paste0(parentcoefs[i,j], "*", parentnames[j])
      }
      flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

      cat("log(phi_", nodename, "[", nodesCatIdx[i], "]) <- ", nodesintercepts[i], " + ", flatpredictors, " + alpha_", nodename, "[", nodesCatIdx[i]-1, "]", sep = "", fill = TRUE)
    } else {
      # no parent nodes
      cat("log(phi_", nodename, "[", nodesCatIdx[i], "]) <- ", nodesintercepts[i], " + alpha_", nodename, "[", nodesCatIdx[i]-1, "]", sep = "", fill = TRUE)
    }
  }
  ## Grouping of node
  cat("alpha_", nodename, " ~ dmnorm.vcov(mu_alpha_", nodename, ", sigma_alpha_", nodename, ")", sep = "", fill = TRUE)
  for (k in 1:dim(sigma_alpha)[1]){
    cat("mu_alpha_", nodename, "[", k, "] <- 0", sep = "", fill = TRUE)
    for (j in 1:dim(sigma_alpha)[2]){
      cat("sigma_alpha_", nodename, "[", k, ", ", j, "] <- ", sigma_alpha[k,j], sep = "", fill = TRUE)
    }
  }
}

#' Bugs code for Poisson response
#'
#' Bugs model for count response variable
#' \eqn{X \sim \mathcal{Pois}(\lambda)}.
#'
#' @param nodename character string of response variable name.
#' @param parentnames single character string (for one parent) or vector of characters (for multiple parent nodes) with parent node (predictor variables) names.
#' @param nodesintercept numeric intercept value for current response variable.
#' @param parentcoefs list named with parent node names of numeric vectors with coefficient values.
#'
#' @return Bugs model returned as stdout.
#'
#' @examples
#' pois_bugs(nodename = "a",
#'           parentnames = c("b", "c"),
#'           nodesintercept = c(0.318077),
#'           parentcoefs = list("b"=c(b=0.3059395),
#'                              "c"=c(c=0.5555)))
#' @importFrom stringi stri_flatten
#' @seealso [makebugs] [simulateAbn]
#' @export
#' @keywords utilities internal
pois_bugs <- function(nodename, nodesintercept, parentnames, parentcoefs){
  cat(nodename, " ~ dpois(lambda.", nodename, ") # Poisson response", sep = "", fill = TRUE)
  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    predictors <- c()
    for (i in 1:length(parentnames)){
      for (j in 1:length(parentcoefs[[i]])){
        predictors <- append(predictors, paste0(parentcoefs[[i]][j], "*", names(parentcoefs[[i]][j])))
      }
    }
    flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

    cat("log(lambda.", nodename, ") <- ", nodesintercept, " + ", flatpredictors, " # logistic regression", sep = "", fill = TRUE)
  } else {
    cat("log(lambda.", nodename, ") <- ", nodesintercept, " # logistic regression", sep = "", fill = TRUE)
  }
}

#' @describeIn pois_bugs Bugs code for Poisson response with varying intercept
#' @param nodesintercept overall mean of response. Parameter from fixed-effects intercept.
#' @param parentcoefs overall slope for each predictor (parent node) variable (fixed-effects).
#' @param sigma_alpha between-group variance. Parameter from random-effects intercept.
pois_bugsGroup <- function(nodename, nodesintercept, parentnames, parentcoefs, sigma_alpha){
  ## Node as response
  cat(nodename, " ~ dpois(lambda_", nodename, ")", sep = "", fill = TRUE)

  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    # With parent nodes
    predictors <- c()
    for (i in 1:length(parentnames)){
      for (j in 1:length(parentcoefs[[i]])){
        predictors <- append(predictors, paste0(parentcoefs[[i]][j], "*", parentnames[i]))
      }
    }
    flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

    cat("log(lambda_", nodename, ") <- mu_", nodename, " + ", flatpredictors, " + alpha_", nodename, sep = "", fill = TRUE)
  } else {
    # No parent nodes
    cat("log(lambda_", nodename, ") <- mu_", nodename, " + alpha_", nodename, sep = "", fill = TRUE)
  }
  cat("mu_", nodename, " <- ", nodesintercept, sep = "", fill = TRUE)
  ## Grouping of node
  cat("alpha_", nodename, " ~ dnorm(mu_alpha_", nodename, ", tau_alpha_", nodename, ")", sep = "", fill = TRUE)
  cat("mu_alpha_", nodename, " <- 0", sep = "", fill = TRUE)
  cat("tau_alpha_", nodename, " <- inverse(sigma_alpha_", nodename, ")", sep = "", fill = TRUE)
  cat("sigma_alpha_", nodename, " <- ", sigma_alpha, sep = "", fill = TRUE)
}

#' Make BUGS model from fitted DAG
#'
#' @param dag named adjacency matrix representing the DAG. Names correspond to node names.
#' @param data.dists list of node distributions.
#' @param coefs a list named by the node names containing for each element a matrix with the nodes' coefficients.
#' @param stderrors a list named by the node names containing for each element a matrix with the nodes' standard errors
#'
#' @return Bugs model returned as stdout.
#'
#' @examples
#' ## Prepare data and arguments
#' mydists <- list(a="gaussian",
#'                 b="multinomial",
#'                 c="binomial",
#'                 d="poisson")
#' mydag <- matrix(0, 4, 4, byrow = TRUE,
#'                 dimnames = list(c("a", "b", "c", "d"),
#'                                 c("a", "b", "c", "d")))
#' mydag[2,1] <- mydag[3,2] <- mydag[4,3] <- 1
#' # plotAbn(mydag, data.dists = mydists)
#' mycoefs <- list("a"=matrix(-6.883383e-17, byrow = TRUE,
#'                            dimnames = list(NULL,
#'                                            "a|intercept")),
#'                 "b"=matrix(c(2.18865, 3.133928, 3.138531, 1.686432, 3.134161, 5.052104),
#'                            nrow= 1, byrow = TRUE,
#'                            dimnames = list(c(NULL),
#'                                       c("b|intercept.2", "b|intercept.3", "b|intercept.4",
#'                                       "a.2", "a.3", "a.4"))),
#'                 "c"=matrix(c(1.11, 2.22, 3.33, 4.44, 5.55),
#'                            nrow= 1, byrow = TRUE,
#'                            dimnames = list(c(NULL),
#'                                       c("c|intercept", "b1", "b2", "b3", "b4"))),
#'                 "d"=matrix(c(3.33, 4.44),
#'                            nrow= 1, byrow = TRUE,
#'                            dimnames = list(c(NULL),
#'                                       c("d|intercept", "c"))))
#' mymse <- c("a"=0,"b"=1,"c"=2,"d"=3)
#' ## Make BUGS model
#' makebugs(dag = mydag, data.dists = mydists, coefs = mycoefs, stderrors = mymse)
#'
#' @importFrom stringi stri_detect_fixed stri_split_fixed
#' @seealso [simulateAbn] [gauss_bugs] [bern_bugs] [categorical_bugs] [pois_bugs]
#' @export
#' @keywords utilities internal
makebugs <- function(dag, data.dists, coefs, stderrors){
  cat("model{", sep = "", fill = TRUE)
  # main loop
  for (n in 1:length(data.dists)){
    nodename <- names(data.dists[n])
    nodedist <- data.dists[[n]]
    nodesintercepts <- coefs[[nodename]][stringi::stri_detect_fixed(str = colnames(coefs[[nodename]]), pattern = "|intercept")]
    if (length(nodesintercepts) == 0){
      nodesintercepts <- 0
    }
    parentnames <- colnames(dag)[which(dag[nodename, ] == 1)]
    if (length(parentnames) < 1) {
      parentnames <- NA
    }

    parentcoefs <- list()
    for (p in 1:length(parentnames)){
      # iterate through all parent nodes
      pname <- parentnames[p]
      if(!is.na(pname) && pname != nodename){
        if (data.dists[[pname]] == "multinomial"){
          # if parent is multinomial we need a different coefficient extraction pipeline.
          # This is because of the different variable naming in fitAbn results and data.df.mult.
          parentsCatNames <- sapply(X = colnames(coefs[[nodename]]),
                                    FUN = function(x){
                                      lookup <- stringi::stri_detect_regex(str = x, pattern = paste0(pname))
                                      if(lookup){
                                        if(stringi::stri_length(x[lookup])>length(pname)){
                                          stringi::stri_split_fixed(str = x, pattern = paste0(pname), omit_empty = TRUE, simplify = TRUE, tokens_only = FALSE)
                                        } else if(stringi::stri_length(x[lookup])==length(pname)){
                                          pname
                                        }
                                      }
                                    })
          parentsCatNames <- stringi::stri_omit_empty(unlist(parentsCatNames))

          parentsCatIdx <- c()
          for (i in 1:length(parentsCatNames)) {
            parentsCatIdx[i] <- i
          }

          pcoefs <- coefs[[nodename]][stringi::stri_detect_regex(str = colnames(coefs[[nodename]]), pattern = paste0("^", pname))]
          pcoefs <- exp(pcoefs)/sum(exp(pcoefs)) # from log-odds to probability
          names(pcoefs) <- rep(pname, length(parentsCatNames))

        } else {
          # parent is not multinomial we have a slightly different extraction pipeline
          parentsCatNames <- sapply(X = colnames(coefs[[nodename]]),
                                    FUN = function(x){
                                      lookup <- stringi::stri_detect_regex(str = x, pattern = paste0("^", pname))
                                      if(lookup){
                                        if(stringi::stri_length(x[lookup])>length(pname)){
                                          stringi::stri_split_fixed(str = x, pattern = paste0("^", pname), omit_empty = TRUE, simplify = TRUE, tokens_only = FALSE)
                                        } else if(stringi::stri_length(x[lookup])==length(pname)){
                                          pname
                                        }
                                      }
                                    })
          parentsCatNames <- stringi::stri_omit_empty(unlist(parentsCatNames))

          pcoefs <- coefs[[nodename]][stringi::stri_detect_regex(str = colnames(coefs[[nodename]]), pattern = paste0("^", pname))]
          names(pcoefs) <- rep(pname, length(parentsCatNames))
        }
        parentcoefs[[pname]] <- pcoefs
      }
    }

    switch (data.dists[[n]],
            gaussian = {
              nodestderror <- stderrors[nodename]
              gauss_bugs(nodename = nodename,
                         nodesintercept = nodesintercepts,
                         parentnames = parentnames,
                         parentcoefs = parentcoefs,
                         std = nodestderror)
            },
            multinomial = {
              nodesCatNames <- sapply(X = colnames(coefs[[nodename]]),
                                      FUN = function(x){
                                        if(stringi::stri_detect_fixed(str = x, pattern = "intercept")){
                                          stringi::stri_split_fixed(str = x, pattern = "|intercept.")
                                        }
                                      }
              )
              nodesCatNames <- unlist(nodesCatNames, use.names = FALSE)[which(unlist(nodesCatNames) != nodename)]
              nodesCatIdx <- 0L
              for (i in 1:length(nodesCatNames)) {
                nodesCatIdx[i] <- i+1 # +1 because reference category is not among the nodesCatnames
              }
              categorical_bugs(nodename = nodename,
                               nodesCatIdx = nodesCatIdx,
                               parentnames = parentnames,
                               nodesintercepts = nodesintercepts,
                               parentcoefs = parentcoefs)
            },
            binomial = {
              bern_bugs(nodename = nodename,
                        nodesintercept = nodesintercepts,
                        parentnames = parentnames,
                        parentcoefs = parentcoefs)
            },
            poisson = {
              pois_bugs(nodename = nodename,
                        nodesintercept = nodesintercepts,
                        parentnames = parentnames,
                        parentcoefs = parentcoefs)
            }
    )
  }
  # print EOF
  cat("}", sep = "", fill = TRUE)
}

#' Make BUGS model from fitted DAG with grouping
#'
#' @inheritParams makebugs
#' @inheritParams fitAbn.mle
#' @param mu Standard deviation of fixed effects.
#' @param betas Coefficients/slopes of fixed effects .
#' @param sigma variance of random effects.
#' @param sigma_alpha variance-covariance matrix corresponding to covariances output from \code{\link[mclogit]{mblogit}}.
#'
#' @return Bugs model returned as stdout.
#'
#' @importFrom stringi stri_detect_fixed stri_split_fixed
#' @seealso [simulateAbn] [gauss_bugsGroup] [bern_bugsGroup] [categorical_bugsGroup] [pois_bugsGroup]
#' @export
#' @keywords utilities internal
makebugsGroup <- function(dag,
                          data.dists,
                          # coefs,
                          # var,
                          stderrors,
                          group.var,
                          mu,
                          betas,
                          sigma,
                          sigma_alpha){
  cat("model{", sep = "", fill = TRUE)
  # main loop
  for (n in 1:length(data.dists)){
    ###
    # Extract the arguments. This is slightly different from the non-group function due to the different structure in 'abnFit'.
    ###
    nodename <- names(data.dists[n])
    nodedist <- data.dists[[n]]
    nodesintercepts <- mu[[n]]
    parentnames <- colnames(dag)[which(dag[nodename, ] == 1)]
    parentcoefs <- betas[[nodename]]
    if (nodedist %in% c("binomial", "poisson", "multinomial")){
      sig <- NA
      sig_alpha <- as.matrix(sigma_alpha[[nodename]])
    } else {
      # gaussian
      sig <- sigma[[nodename]]
      sig_alpha <- sigma_alpha[[nodename]]
    }

    ###
    # Transform the arguments
    ###
    # If sigma or sigma_alpha contain 0 or negative values, JAGS can't invert them to the precision (matrix).
    # This is circumvented with assigning a very small value instead of zero.
    if (!any(is.na(sig)) && any(sig<=0)){
      sig[which(sig<=0)] <- 1e-50
    }
    if (!any(is.na(sig_alpha)) && any(sig_alpha<=0)){
      sig_alpha[which(sig_alpha<=0)] <- 1e-50
    }

    ###
    # Generate BUGS model
    ###
    switch (data.dists[[n]],
            gaussian = {
              nodestderror <- stderrors[nodename]
              gauss_bugsGroup(nodename = nodename, # response
                              nodesintercept = nodesintercepts, # fixed-effect intercept (fit@beta[1])
                              parentnames = parentnames, # predictor(s)
                              parentcoefs = parentcoefs, # predictor(s) coefficients: fixed-effects (fit@beta[>1])
                              sigma = sig, # random-effect residual (fit@theta[>1])
                              sigma_alpha = sig_alpha) # random-effect group intercept (fit@theta[1])
            },
            multinomial = {
              nodesCatNames <- sapply(X = names(mu[[nodename]]),
                                      FUN = function(x){
                                        if(stringi::stri_detect_fixed(str = x, pattern = "~")){
                                          tmp <- stringi::stri_split_fixed(str = x, pattern = "~", simplify = T)
                                          unique(stringi::stri_split_fixed(str = tmp[,1], pattern = paste0(nodename, "."), simplify = TRUE, omit_empty = TRUE))
                                        }
                                      }
              )
              nodesCatIdx <- 0L
              for (i in 1:length(nodesCatNames)) {
                nodesCatIdx[i] <- i+1 # +1 because reference category is not among the nodesCatnames
              }
              categorical_bugsGroup(nodename = nodename,
                                    nodesCatIdx = nodesCatIdx,
                                    parentnames = parentnames,
                                    nodesintercepts = nodesintercepts,
                                    parentcoefs = parentcoefs,
                                    sigma = sig,
                                    sigma_alpha = sig_alpha)
            },
            binomial = {
              bern_bugsGroup(nodename = nodename,
                             nodesintercept = nodesintercepts,
                             parentnames = parentnames,
                             parentcoefs = parentcoefs,
                             sigma_alpha = sig_alpha)
            },
            poisson = {
              pois_bugsGroup(nodename = nodename,
                             nodesintercept = nodesintercepts,
                             parentnames = parentnames,
                             parentcoefs = parentcoefs,
                             sigma_alpha = sig_alpha)
            }
    )
  }
  # print EOF
  cat("}", sep = "", fill = TRUE)
}

#' Simulate data from a fitted additive Bayesian network.
#'
#' @param object of type \code{abnFit}.
#' @param run.simulation call JAGS to simulate data (default is \code{TRUE}).
#' @param bugsfile A path to a valid file or \code{NULL} (default)  to delete the bugs file after simulation.
#' @param n.chains number of parallel chains for the model.
#' @param n.adapt number of iteration for adaptation. If \code{n.adapt} is set to zero, then no adaptation takes place.
#' @param n.thin thinning interval for monitors.
#' @param n.iter number of iteration to monitor.
#' @param seed by default set to 42.
#' @param verbose if TRUE prints additional output
#' @param debug if TRUE prints bug file content to stdout and does not run simulations.
#'
#' @return data.frame
#'
#' @examples
#' df <- FCV[, c(12:15)]
#' mydists <- list(Outdoor="binomial",
#'                 Sex="multinomial",
#'                 GroupSize="poisson",
#'                 Age="gaussian")
#'
#' ## buildScoreCache -> mostProbable() -> fitAbn()
#' suppressWarnings({
#'   mycache.mle <- buildScoreCache(data.df = df, data.dists = mydists, method = "mle",
#'                                  adj.vars = NULL, cor.vars = NULL,
#'                                  dag.banned = NULL, dag.retained = NULL,
#'                                  max.parents = 1,
#'                                  which.nodes = NULL, defn.res = NULL)
#' }) # ignore non-convergence warnings
#' mp.dag.mle <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
#' myres.mle <- fitAbn(object = mp.dag.mle, method = "mle")
#'
#' myres.sim <- simulateAbn(object = myres.mle,
#'                              run.simulation = TRUE,
#'                              bugsfile = NULL,
#'                              verbose = FALSE)
#' str(myres.sim)
#' prop.table(table(myres.sim$Outdoor))
#' prop.table(table(df$Outdoor))
#'
#' @importFrom rjags jags.model coda.samples
#' @importFrom utils capture.output
#' @seealso [makebugs]
#' @export
simulateAbn <- function(object = NULL,
                        run.simulation = TRUE,
                        bugsfile = NULL,
                        n.chains = 10L,
                        n.adapt = 1000L,
                        n.thin = 100L,
                        n.iter = 10000L,
                        seed = 42L,
                        verbose = FALSE,
                        debug = FALSE){
  ###
  # Validation of arguments
  ###
  if (!inherits(verbose, "logical")){
    stop("'verbose' must be either TRUE or FALSE.")
  }

  if (!inherits(seed, "integer") | seed<0){
    stop("'seed' must be a non-negative integer")
  }

  if (!inherits(n.chains, "integer") | n.chains<0){
    stop("'n.chains' must be a non-negative integer")
  }
  if (!inherits(n.adapt, "integer") | n.adapt<0){
    stop("'n.adapt' must be a non-negative integer")
  }
  if (!inherits(n.thin, "integer") | n.thin<0){
    stop("'n.thin' must be a non-negative integer")
  }
  if (!inherits(n.iter, "integer") | n.iter<0){
    stop("'n.iter' must be a non-negative integer")
  }

  if(!is.null(bugsfile)){
    if(file.exists(bugsfile) && verbose){
      warning("'bugs.file' already exists and will be overwritten.")
    } else {
      outfile <- file.path(bugsfile)
    }
  } else if (is.null(bugsfile)){
    # create temporary file that will be deleted afterwards
    outfile <- tempfile(pattern = "model.bugs")
  } else {
    stop("'bugs.file' must be NULL or an existing path to a file.")
  }

  if (!is.null(object)){
    # Check object
    if(!inherits(object, "abnFit")){
      stop("'object' is not of class 'abnFit'.")
    }

    # extract relevant stuff
    dag <- object$abnDag$dag
    data.dists <- object$abnDag$data.dists
    mse <- object$mse
    if(!is.null(object$group.var)){
      group.var <- object$group.var
      group.ids <- object$group.ids
      grouped.vars <- object$grouped.vars
      mu <- object$mu
      betas <- object$betas
      sigm <- object$sigma
      sigm_alpha <- object$sigma_alpha
    } else {
      coefs <- object$coef
      group.var <- NULL
    }

  } else {
    # relevant stuff is provided otherwise
    stop("Please provide an object of class 'abnFit'.")
  }

  if(!inherits(run.simulation, "logical")){
    stop("'run.simulation' must be either TRUE or FALSE.")
  }

  ###
  # Create BUGS Model
  ###
  # increase console width to avoid line breaks that cause JAGS errors.
  origwidth <- getOption("width")
  options(width=1000)

  if(!is.null(group.var)){
    # save(file = "tests/testthat/testdata/makebugsGauss_data.Rdata",
    #      list = c("dag", "data.dists", "mse", "group.var", "mu", "betas", "sigm", "sigm_alpha"))
    if (debug){
      makebugsGroup(dag = dag,
                    data.dists = data.dists,
                    stderrors = mse,
                    group.var = group.var,
                    mu = mu,
                    betas = betas,
                    sigma = sigm,
                    sigma_alpha = sigm_alpha)
      stop("JAGS is not run in debug mode")
    } else {
      capture.output(
        makebugsGroup(dag = dag,
                      data.dists = data.dists,
                      stderrors = mse,
                      mu = mu,
                      betas = betas,
                      sigma = sigm,
                      sigma_alpha = sigm_alpha),
        file = outfile,
        append = FALSE
      )
    }
  } else {
    if (debug){
      makebugs(dag = dag,
               data.dists = data.dists,
               coefs = coefs,
               stderrors = mse)
      stop("JAGS is not run in debug mode")
    } else {
      capture.output(
        makebugs(dag = dag,
                 data.dists = data.dists,
                 coefs = coefs,
                 stderrors = mse),
        file = outfile,
        append = FALSE
      )
    }
  }
  # reset console width to original value
  options(width = origwidth)

  if(verbose){
    system(paste0("cat ", outfile))
  }

  ###
  # Run JAGS
  ###
  if(run.simulation){
    jagsverbose <- options(jags.pb = ifelse(verbose, "text", "none"))
    on.exit(jagsverbose)

    # initialise JAGS
    init <- list(.RNG.name = "base::Mersenne-Twister",
                 .RNG.seed = seed)

    # call JAGS
    jgs <- jags.model(file = outfile,
                      inits = init,
                      n.chains = n.chains,
                      n.adapt = n.adapt,
                      quiet = !verbose)
    res <- coda.samples(model = jgs,
                        variable.names = as.character(names(data.dists)),
                        thin = n.thin,
                        n.iter = n.iter)
    res <- do.call(rbind.data.frame, res)

    # format output
    # make factor if binomial and multinomial
    for (i in 1:length(data.dists)) {
      if (data.dists[i] %in% c("binomial", "multinomial")) {
        res[,names(data.dists[i])] <- as.factor(res[,names(data.dists[i])])
      }
    }

    return(res)
  } else {
    if(verbose){
      message(paste0("No simulation is run but bugfile is written to: ", outfile))
    }
    return(NULL)
  }
}
