#' Bugs code for Gaussian response
#'
#' Bugs model for a normal distributed response variable
#' \eqn{X \sim \mathcal{N}(\mu,\,\sigma^{2})}.
#'
#' @param nodename character string of response variable name.
#' @param parentnames single character string (for one parent) or vector of characters (for multiple parent nodes) with parent node (predictor variables) names.
#' @param nodesintercept numeric intercept value for current response variable.
#' @param parentcoefs list named with parent node names of numeric vectors with coefficient values.
#' @param data.dists list named with node names of distribution types.
#' @param variance integer with variance of response variable that will be
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
#'            data.dists = list("a = "gaussian",
#'                              "b" = "gaussian",
#'                              "c" = "gaussian"),
#'            variance = c(0.05773503))
#' @importFrom stringi stri_flatten
#' @seealso [makebugs] [simulateAbn]
#' @export
#' @keywords utilities internal
gauss_bugs <- function(nodename, nodesintercept, parentnames, parentcoefs, data.dists, variance){
  cat(nodename, " ~ dnorm(mu.", nodename, ", precision.", nodename, ") # Gaussian response", sep = "", fill = TRUE)
  # mu
  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    predictors <- vapply(parentnames, function(pname) {
      build_predictor_terms(nodename, pname, parentcoefs[[pname]], data.dists)
    }, character(1))
    flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

    cat("mu.", nodename, " <- ", nodesintercept, " + ", flatpredictors, " # Linear regression", sep = "", fill = TRUE)
  } else {
    cat("mu.", nodename, " <- ", nodesintercept, " # Linear regression", sep = "", fill = TRUE)
  }

  # tau (precision)
  cat("precision.", nodename, " <- 1/", variance, " # precision tau = 1/variance",
      sep = "", fill = TRUE)
}

#' @describeIn gauss_bugs Bugs code for Gaussian response with varying intercept
#' @param nodesintercept overall mean of response. Parameter from fixed-effects intercept.
#' @param parentcoefs overall slope for each predictor (parent node) variable (fixed-effects).
#' @param data.dists list named with node names of distribution types
#' @param sigma within-group variance. Parameter from random-effects residual.
#' @param sigma_alpha between-group variance. Parameter from random-effects intercept.
gauss_bugsGroup <- function(nodename, nodesintercept, parentnames, parentcoefs, data.dists, sigma, sigma_alpha){
  ## Node as response
  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    # With parent nodes
    predictors <- vapply(parentnames, function(pname) {
      build_predictor_terms(nodename, pname, parentcoefs[[pname]], data.dists)
    }, character(1))
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
#' @param data.dists list named with node names of distribution types.
#' @return Bugs model returned as stdout.
#'
#' @examples
#' bern_bugs(nodename = "a",
#'           parentnames = c("b", "c"),
#'           nodesintercept = c(0.318077),
#'           parentcoefs = list("b"=c(b=0.3059395),
#'                              "c"=c(c=0.5555)),
#'           data.dists = c("a" = "binomial",
#'                         "b" = "binomial",
#'                         "c" = "binomial"))
#' @importFrom stringi stri_flatten
#' @seealso [makebugs] [simulateAbn]
#' @export
#' @keywords utilities internal
bern_bugs <- function(nodename, nodesintercept, parentnames, parentcoefs, data.dists){
  cat(nodename, " ~ dbern(p.", nodename, ") # Bernoulli response", sep = "", fill = TRUE)
  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    predictors <- vapply(parentnames, function(pname) {
      build_predictor_terms(nodename, pname, parentcoefs[[pname]], data.dists)
    }, character(1))
    flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

    cat("logit(p.", nodename, ") <- ", nodesintercept, " + ", flatpredictors, " # logistic regression", sep = "", fill = TRUE)
  } else {
    cat("logit(p.", nodename, ") <- ", nodesintercept, " # logistic regression", sep = "", fill = TRUE)
  }
}

#' @describeIn bern_bugs Bugs code for Bernoulli response with varying intercept
#' @param nodesintercept overall mean of response. Parameter from fixed-effects intercept.
#' @param parentcoefs overall slope for each predictor (parent node) variable (fixed-effects).
#' @param data.dists list named with node names of distribution types
#' @param sigma_alpha between-group variance. Parameter from random-effects intercept.
bern_bugsGroup <- function(nodename, nodesintercept, parentnames, parentcoefs, data.dists, sigma_alpha){
  ## Node as response
  cat(nodename, " ~ dbern(p_", nodename, ")", sep = "", fill = TRUE)
  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    # With parent nodes
    predictors <- vapply(parentnames, function(pname) {
      build_predictor_terms(nodename, pname, parentcoefs[[pname]], data.dists)
    }, character(1))
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
#' @param data.dists list named with node names of distribution types.
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
#'                  parentcoefs = list("a"=c(a=1.686432, a=3.134161, a=5.052104)),
#'                  data.dists = list("a" = "gaussian", 
#'                                     "b" = "multinomial"))
#' @importFrom stringi stri_flatten
#' @seealso [makebugs] [simulateAbn]
#' @export
#' @keywords utilities internal
categorical_bugs <- function(nodename, nodesCatIdx, parentnames, nodesintercepts, parentcoefs, data.dists){
  cat(nodename, " ~ dcat(p.", nodename, ") # Categorical response", sep = "", fill = TRUE)
  cat("p.", nodename, "[1] <- phi.", nodename, "[1]/sum(phi.", nodename, ") # soft-max", sep = "", fill = TRUE)
  cat("log(phi.", nodename, "[1]) <- 0 # Reference category", sep = "", fill = TRUE)

  for(i in 1:length(nodesCatIdx)){
    cat("p.", nodename, "[", nodesCatIdx[i], "] <- phi.", nodename, "[", nodesCatIdx[i], "]/sum(phi.", nodename, ") # soft-max", sep = "", fill = TRUE)
    if (!any(is.na(parentnames)) && length(parentnames) >0){
      predictors <- vapply(parentnames, function(pname) {
        build_predictor_terms(nodename, pname, parentcoefs[[pname]], data.dists,
                              child_cat_i = i)
      }, character(1))
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
#' @param data.dists list named with node names of distribution types
#' @param sigma within-group variance. Parameter from random-effects residual.
#' @param sigma_alpha between-group variance-covariance matrix. Parameters from random-effects intercept.
categorical_bugsGroup <- function(nodename, nodesCatIdx, nodesintercepts, parentnames, parentcoefs, data.dists, sigma, sigma_alpha){
  cat(nodename, " ~ dcat(p_", nodename, ")", sep = "", fill = TRUE)
  cat("p_", nodename, "[1] <- phi_", nodename, "[1]/sum(phi_", nodename, ")", sep = "", fill = TRUE)
  cat("log(phi_", nodename, "[1]) <- 0 + alpha_", nodename, "[1]", sep = "", fill = TRUE)

  for(i in 1:length(nodesCatIdx)){
    cat("p_", nodename, "[", nodesCatIdx[i], "] <- phi_", nodename, "[", nodesCatIdx[i], "]/sum(phi_", nodename, ")", sep = "", fill = TRUE)
    if (!any(is.na(parentnames)) && length(parentnames) >0){
      predictors <- vapply(parentnames, function(pname) {
        build_predictor_terms(nodename, pname, parentcoefs[[pname]], data.dists,
                              child_cat_i = i)
      }, character(1))
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
#' @param data.dists list named with node names of distribution types.
#' @return Bugs model returned as stdout.
#'
#' @examples
#' pois_bugs(nodename = "a",
#'           parentnames = c("b", "c"),
#'           nodesintercept = c(0.318077),
#'           parentcoefs = list("b"=c(b=0.3059395),
#'                              "c"=c(c=0.5555)),
#'           data.dists = list("a" = "poisson",
#'                              "b" = "poisson", 
#'                              "c" = "poisson"))
#' @importFrom stringi stri_flatten
#' @seealso [makebugs] [simulateAbn]
#' @export
#' @keywords utilities internal
pois_bugs <- function(nodename, nodesintercept, parentnames, parentcoefs, data.dists){
  cat(nodename, " ~ dpois(lambda.", nodename, ") # Poisson response", sep = "", fill = TRUE)
  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    predictors <- vapply(parentnames, function(pname) {
      build_predictor_terms(nodename, pname, parentcoefs[[pname]], data.dists)
    }, character(1))
    flatpredictors <- stringi::stri_flatten(predictors, collapse = " + ")

    cat("log(lambda.", nodename, ") <- ", nodesintercept, " + ", flatpredictors, " # logistic regression", sep = "", fill = TRUE)
  } else {
    cat("log(lambda.", nodename, ") <- ", nodesintercept, " # logistic regression", sep = "", fill = TRUE)
  }
}

#' @describeIn pois_bugs Bugs code for Poisson response with varying intercept
#' @param nodesintercept overall mean of response. Parameter from fixed-effects intercept.
#' @param parentcoefs overall slope for each predictor (parent node) variable (fixed-effects).
#' @param data.dists list named with node names of distribution types.
#' @param sigma_alpha between-group variance. Parameter from random-effects intercept.
pois_bugsGroup <- function(nodename, nodesintercept, parentnames, parentcoefs, data.dists, sigma_alpha){
  ## Node as response
  cat(nodename, " ~ dpois(lambda_", nodename, ")", sep = "", fill = TRUE)

  if (!any(is.na(parentnames)) && length(parentnames > 0)){
    # With parent nodes
    predictors <- vapply(parentnames, function(pname) {
      build_predictor_terms(nodename, pname, parentcoefs[[pname]], data.dists)
    }, character(1))
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

build_predictor_terms <- function(nodename, pname, pcoef, data.dists, child_cat_i = NULL) {
  if (data.dists[[pname]] == "multinomial") {
    if (data.dists[[nodename]] == "multinomial") {
      # multinomial child of multinomial parent:
      # pcoef[i, k] = coef for child cat i, parent level k (all levels, 1-indexed)
      n.parent.cats <- ncol(pcoef)
      paste0(pcoef[child_cat_i, ], "*equals(", pname, ",",
             seq(2, n.parent.cats+1), ")", collapse = " + ")
    } else {
      # non-multinomial child of multinomial parent:
      # flat vector, one coef per parent level
      n.parent.cats <- length(pcoef)
      paste0(pcoef, "*equals(", pname, ",",
             seq(2, n.parent.cats+1), ")", collapse = " + ")
    }
  } else {
    # continuous/binomial/poisson parent: single coefficient
    if (!is.null(child_cat_i)) {
      paste0(pcoef[child_cat_i], "*", pname)
    } else {
      paste0(pcoef, "*", pname)
    }
  }
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
    
    intercept_cols <- colnames(coefs[[nodename]])[
      stringi::stri_detect_fixed(
        str     = colnames(coefs[[nodename]]),
        pattern = "|intercept"
      )
    ]
    if (length(intercept_cols) > 0) {
      nodesintercepts <- coefs[[nodename]][, intercept_cols, drop = TRUE]
    } else {
      nodesintercepts <- 0
    }
    
    parentnames <- colnames(dag)[which(dag[nodename, ] == 1)]
    if (length(parentnames) < 1) {
      parentnames <- NA
    }

    parentcoefs <- list()
    if (any(!is.na(parentnames))){
      parentnames <- parentnames[parentnames != nodename]
      
      if (length(parentnames) > 0) {
        parentnames <- parentnames[order(nchar(parentnames), decreasing = TRUE)]
        
        for (pname in parentnames){
          # iterate through all parent nodes
          p_dist <- data.dists[[pname]]
          all_cols <- colnames(coefs[[nodename]])
          
          pattern <- paste0("^", pname)
          matched_cols <- all_cols[stringi::stri_detect_regex(all_cols, pattern)]
          if (length(matched_cols) > 0) {
            other_parents <- setdiff(parentnames, pname)
            longer_parents <- other_parents[startsWith(other_parents, pname) & nchar(other_parents) > nchar(pname)]
            
            for (lp in longer_parents) {
              matched_cols <- matched_cols[!startsWith(matched_cols, lp)]
            }
          }
          coef_subset <- coefs[[nodename]][, matched_cols, drop = FALSE]

          parentcoefs[[pname]] <- coef_subset
        }
      }
    }

    switch (data.dists[[n]],
            gaussian = {
              nodestderror <- stderrors[nodename]
              gauss_bugs(nodename = nodename,
                         nodesintercept = nodesintercepts,
                         parentnames = parentnames,
                         parentcoefs = parentcoefs,
                         variance = nodestderror,
                         data.dists = data.dists)
            },
            multinomial = {
              if (length(intercept_cols) > 0) {
                # root multinomial: child category indices from intercept col names
                # e.g. "AGE|intercept.2", "AGE|intercept.3" -> c(2L, 3L)
                raw_cats <- stringi::stri_replace_first_regex(
                  str         = intercept_cols,
                  pattern     = paste0("^", nodename, "\\|intercept\\."),
                  replacement = ""
                )
                
                # Check if they are numeric digits or text labels
                if (all(grepl("^[0-9]+$", raw_cats))) {
                  nodesCatIdx <- as.integer(raw_cats)
                } else {
                  # BUGS requires integer indices. Map text factor levels to sequential integers (2, 3, 4...)
                  nodesCatIdx <- seq_len(length(raw_cats)) + 1 
                  
                  # Save the original text labels as an attribute so they are never lost
                  attr(nodesCatIdx, "level_names") <- raw_cats
                }
              } else {
                # multinomial child of multinomial parent: no intercept columns.
                # Child categories are encoded as the trailing digit of coef names,
                # e.g. "AGE12","AGE13" -> child cats 2, 3.
                all_col_names <- colnames(coefs[[nodename]])
                extracted_suffixes <- sub("^.*\\.", "", all_col_names) # if dot-separated
                
                if (all(grepl("^[0-9]+$", extracted_suffixes))) {
                  nodesCatIdx <- sort(unique(as.integer(extracted_suffixes)))
                } else {
                  unique_sufs <- sort(unique(extracted_suffixes))
                  nodesCatIdx <- seq_len(length(unique_sufs)) + 1
                  attr(nodesCatIdx, "level_names") <- unique_sufs
                }
                nodesintercepts <- rep(0, length(nodesCatIdx))
              }
              n.child.cats <- length(nodesCatIdx)
              
              # Reshape parentcoefs into matrix only for multinomial parents.
              # Multinomial parent coef names: {pname}{parent_cat}{child_cat}
              #   -> n.child.cats blocks of n.parent.cats coefficients -> matrix.
              # Non-multinomial parent coef names: {pname}{parent_cat}
              #   -> flat vector used directly with equals() in BUGS.
              for (pname in names(parentcoefs)) {
                if (data.dists[[pname]] == "multinomial") {
                  flat <- parentcoefs[[pname]]
                  n.coefs.per.child.cat <- length(flat) / n.child.cats
                  parentcoefs[[pname]] <- matrix(
                    flat,
                    nrow  = n.child.cats,
                    ncol  = n.coefs.per.child.cat,
                    byrow = TRUE
                  )
                }
                # non-multinomial parents: leave as flat vector
              }
              categorical_bugs(nodename = nodename,
                               nodesCatIdx = nodesCatIdx,
                               parentnames = parentnames,
                               nodesintercepts = nodesintercepts,
                               parentcoefs = parentcoefs,
                               data.dists = data.dists)
            },
            binomial = {
              bern_bugs(nodename = nodename,
                        nodesintercept = nodesintercepts,
                        parentnames = parentnames,
                        parentcoefs = parentcoefs,
                        data.dists = data.dists)
            },
            poisson = {
              pois_bugs(nodename = nodename,
                        nodesintercept = nodesintercepts,
                        parentnames = parentnames,
                        parentcoefs = parentcoefs,
                        data.dists = data.dists)
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
                              data.dists = data.dists,
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
                                    data.dists = data.dists,
                                    sigma = sig,
                                    sigma_alpha = sig_alpha)
            },
            binomial = {
              bern_bugsGroup(nodename = nodename,
                             nodesintercept = nodesintercepts,
                             parentnames = parentnames,
                             parentcoefs = parentcoefs,
                             data.dists = data.dists,
                             sigma_alpha = sig_alpha)
            },
            poisson = {
              pois_bugsGroup(nodename = nodename,
                             nodesintercept = nodesintercepts,
                             parentnames = parentnames,
                             parentcoefs = parentcoefs,
                             data.dists = data.dists,
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
#' @param n number of samples to generate. If \code{NULL} (default), n=n.iter/(n.chains*n.thin). 
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
                        n = NULL,
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

  if (!is.null(n)){
    n.chains = 1L
    n.adapt = 1000L
    n.thin = 1L
    n.iter = n/n.chains*n.thin
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
    if (!is.list(object) || !all(c("abnDag", "mse", "coef") %in% names(object))) {
      stop("'object' must be a list containing at least 'abnDag', 'mse', and 'coef'.")
    }
    if (!is.list(object$abnDag) || !all(c("dag", "data.dists") %in% names(object$abnDag))) {
      stop("'object$abnDag' must be a list containing 'dag' and 'data.dists'.")
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
  # save current options and set new ones
  origoptions <- options()
  on.exit(options(origoptions)) # reset options to original values after function exits
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
    # get back factors
    for (nodename in names(data.dists)) {
      if (data.dists[[nodename]] == "multinomial" && nodename %in% names(res)) {
        
        # Find the intercept columns in your coefs object for this node
        intercept_cols <- grep(paste0("^", nodename, "\\|intercept\\."), colnames(coefs[[nodename]]), value = TRUE)
        
        if (length(intercept_cols) > 0) {
          # Extract the text suffixes (e.g. "fc", "m", "mc")
          raw_cats <- stringi::stri_replace_first_regex(
            str         = intercept_cols,
            pattern     = paste0("^", nodename, "\\|intercept\\."),
            replacement = ""
          )
          
          if (!all(grepl("^[0-9]+$", raw_cats))) {
            if (!is.null(object$abnDag$data.df) && nodename %in% names(object$abnDag$data.df)) {
              df_levels <- levels(object$abnDag$data.df[[nodename]])
              
              # Find which level in the dataframe is the reference (the one NOT in raw_cats)
              ref_level <- setdiff(df_levels, raw_cats)
              
              # Reconstruct the exact order: Reference category is index 1, followed by raw_cats
              all_levels <- c(ref_level, raw_cats)
            } else {
              # Fallback if no dataset is present
              all_levels <- c("reference", raw_cats)
            }
            
            # Convert sampled integers to factors using the extracted names
            res[[nodename]] <- factor(
              res[[nodename]],
              levels = seq_along(all_levels),
              labels = all_levels
            )
          }
        }
      } else if (data.dists[[nodename]] %in% "binomial" && nodename %in% names(res)) {
        if (!is.null(object$abnDag$data.df) && nodename %in% names(object$abnDag$data.df)) {
          if (length(bin_levels) == 2){
            bin_levels <- levels(object$abnDag$data.df[[nodename]])
            res[[nodename]] <- factor(
               res[[nodename]],
               levels = c(0, 1),
               labels = bin_levels
            )
          }
        }
      }
    }
    
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

#' Generate random Directed Acyclic Graphs (DAGs) with specified density, number of nodes, and max number of parents
#'
#' This function generates a list of random DAG adjacency matrices based on the
#' node distribution and desired number of nodes, density level, and max number of parents
#'
#' @param p_distrib List specifying the number of nodes per distribution
#'    and multinomial levels.
#' @param max_parents. Numeric. Maximal number of parents in the graph.
#' @param d Numeric. Density level (proportion of edges to keep in the graph, between 0 and 1).
#' @param nbr_graphs Integer. Number of random DAGs to generate.
#' @param MB List. Contains `node` the name of the node which Markov blanket size should be controlled and `MB_size` the desired size of the Markov blanket.
#'
#' @return A list of adjacency matrices, each representing a DAG. Each matrix has
#'   dimension equal to the total number of nodes.
#'
#' @examples
#' p_distrib <- list(distrib = c(3, 3, 2, 2), levels_multi = c(3, 4))
#' names(p_distrib$distrib) <- c("gaussian","poisson","binomial","multinomial")
#' density <- 0.3
#' n_graphs <- 5
#' max_parents <- 2
#' dags <- generate_RandomDAGs(p_distrib, max_parents, d = density, n_graphs)
#' length(dags)  # Should be 5
#' sapply(dags, dim)  # Dimensions of adjacency matrices
#'
#' @importFrom igraph graph_from_adjacency_matrix is_connected
#' @export
generate_RandomDAGs <- function(p_distrib, max_parents, d, nbr_graphs, MB = list(MB_size = NULL, node = NULL)){
  if (!is.list(p_distrib) || length(p_distrib) != 2 || length(p_distrib$levels_multi) != p_distrib$distrib["multinomial"] || length(p_distrib$distrib) <= 0 || !all(names(p_distrib$distrib) %in% c("gaussian","binomial","multinomial","poisson"))) {
    stop("'p_distrib' should be a list containing two elements: number of nodes per distribution (only Gaussian, binomial, multinomial and Poisson are supported) and multinomial levels for multinomial nodes.")
  }
  
  if (!is.numeric(d) || length(d) != 1 || d <= 0 || d >= 1) {
    stop("'d' must be a numeric value strictly between 0 and 1.")
  }
  
  if (!is.numeric(max_parents) || length(max_parents) != 1 || max_parents < 0 || max_parents > (sum(p_distrib$distrib)-1)) {
    stop(paste0("'max_parents' must be a numeric value between 0 and ",(sum(p_distrib$distrib)-1),"."))
  }
  
  if (!is.numeric(nbr_graphs) || nbr_graphs <= 0 || floor(nbr_graphs) != nbr_graphs) {
    stop("'nbr_graphs' must be a positive integer.")
  }
  
  p <- sum(p_distrib$distrib)
  node_names <- paste0(substr(rep(names(p_distrib$distrib), p_distrib$distrib), 1, 1), # First letter of each distribution name
                       sequence(p_distrib$distrib))     # Sequence of numbers for each distribution
  node_names_short <- substr(rep(names(p_distrib$distrib), p_distrib$distrib), 1, 1)
  
  adj_matrices <- list()
  count <- 1
  
  T <- matrix(0, ncol=p, nrow=p)
  T[lower.tri(T, diag = TRUE)] <- 1
  T <- T-diag(diag(T))
  
  target_edges <- round(p * (p - 1) / 2 * d)
  
  while (length(adj_matrices)<nbr_graphs && count < 1e5){
    T_tmp <- T
    
    for (j in 1:p) {
      potential_parents <- which(T_tmp[, j] == 1)
      
      if (length(potential_parents) > max_parents) {
        drop_rows <- sample(potential_parents, length(potential_parents) - max_parents)
        T_tmp[drop_rows, j] <- 0
      }
    }
    current_edges <- sum(T_tmp)
    
    if (current_edges > target_edges) {
      # --- CASE 1: Too many edges, trim excess ---
      active_edges <- which(T_tmp == 1, arr.ind = TRUE)
      num_to_drop <- current_edges - target_edges
      drop_indices <- sample(nrow(active_edges), num_to_drop)
      
      for (k in drop_indices) {
        T_tmp[active_edges[k, 1], active_edges[k, 2]] <- 0
      }
      
    } else if (current_edges < target_edges) {
      # --- CASE 2: Too few edges, add missing ones safely ---
      while (sum(T_tmp) < target_edges) {
        # Find all inactive slots in the lower triangle
        inactive_edges <- which(T_tmp == 0 & lower.tri(T_tmp), arr.ind = TRUE)
        if (nrow(inactive_edges) == 0) break # Safety check if grid is full
        
        # Filter to only keep slots whose target column hasn't hit max_parents yet
        col_sums <- colSums(T_tmp)
        valid_slots <- inactive_edges[col_sums[inactive_edges[, 2]] < max_parents, , drop = FALSE]
        
        if (nrow(valid_slots) == 0) break # No more slots can accept parents without breaking max_parents
        
        # Randomly pick one valid slot and turn it on
        add_idx <- sample(nrow(valid_slots), 1)
        T_tmp[valid_slots[add_idx, 1], valid_slots[add_idx, 2]] <- 1
      }
    }
    
    if (sum(T_tmp) != target_edges) {
      count <- count + 1
      next
    }
    
    #edges_to_remove <- sample(which(T==1),round(p*(p-1)/2 - p*(p-1)/2*d))
    #T_tmp[edges_to_remove] <- 0
    perm <- sample(1:p)
    perm <- node_names[perm]
    perm_matrix <- matrix(0, ncol=p, nrow=p)
    colnames(perm_matrix) <- node_names
    row_indices <- seq_along(perm)
    col_indices <- match(colnames(perm_matrix), perm)
    perm_matrix[cbind(row_indices, col_indices)] <- 1
    adj_matrix <- perm_matrix %*% T_tmp %*% t(perm_matrix)
    rownames(adj_matrix) = colnames(adj_matrix) <- node_names
    graph <- graph_from_adjacency_matrix(adj_matrix)
    
    is_valid <- (length(which(rowSums(adj_matrix) + colSums(adj_matrix) == 0)) == 0 && 
                   igraph::is_connected(graph, mode = "weak"))
    
    if (is_valid && !is.null(max_parents)) {
      if (max(colSums(adj_matrix)) > max_parents) {
        is_valid <- FALSE
      }
    }
    
    if (is_valid && !is.null(MB$MB_size)) {
      mb_nodes <- find_MB(graph, MB$node)
      if (length(mb_nodes) != MB$MB_size) {
        is_valid <- FALSE
      }
    }
    
    if (is_valid) {
      adj_matrices <- c(adj_matrices, list(adj_matrix))
    }
    count <- count +1
  }
  return(adj_matrices)
}
