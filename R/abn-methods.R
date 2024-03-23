#' Print objects of class \code{abnDag}
#'
#' @param x Object of class \code{abnDag}
#' @param digits number of digits of the adjacency matrix.
#' @param ... additional parameters. Not used at the moment.
#'
#' @return outputs adjacency matrix and statement of the class of \code{x}.
#' @concept DAG
#' @export
#' @examples
#' mydag <- createAbnDag(dag = ~a+b|a, data.df = data.frame("a"=1, "b"=1))
#' print(mydag)
print.abnDag <- function(x, digits = 3L, ...){
  print(x$dag, digits = digits)
  cat("Class 'abnDag'.\n")
  invisible(x)
}

#' Prints summary statistics from an object of class \code{abnDag}
#'
#' @inherit infoDag params
#' @param ... additional parameters. Not used at the moment.
#' @concept DAG
#' @export
#' @examples
#' mydag <- createAbnDag(dag = ~a+b|a, data.df = data.frame("a"=1, "b"=1))
#' summary(mydag)
#' @returns List with summary statistics of the DAG.
summary.abnDag <- function(object, ...) {
  su <- infoDag(object$dag)
  return(su)
}

#' Plots DAG from an object of class \code{abnDag}
#'
#' @param x Object of class \code{abnDag}
#' @param ... additional parameters. Not used at the moment.
#'
#' @return \code{Rgraphviz::plot}
#' @concept DAG
#' @export
#' @importClassesFrom graph graphAM
#' @importFrom grDevices dev.new dev.flush
#' @examples
#' mydag <- createAbnDag(dag = ~a+b|a,
#'                       data.df = data.frame("a"=1, "b"=1),
#'                       data.dists = list(a="binomial", b="gaussian"))
#' plot(mydag)
plot.abnDag <- function(x, ...){
  # Check if the graph is a "abnDag" object
  if (!inherits(x, "abnDag")) stop('Function type not implemented yet. Use which="abnDag"')

  # Plot the graph
  g <- plotAbn(dag = x$dag, data.dists = x$data.dists, ...)
  invisible(g)
}


##-------------------------------------------------------------------------
## abnCache
##-------------------------------------------------------------------------
#' Print objects of class \code{abnCache}
#'
#' @param x Object of class \code{abnCache}
#' @param digits number of digits of the results.
#' @param ... additional parameters. Not used at the moment.
#'
#' @return summary statement of the class of \code{abnCache}.
#' @concept DAG
#' @export
#' @examples
#' ## Subset of the build-in dataset, see  ?ex0.dag.data
#' mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")] ## take a subset of cols
#'
#' ## setup distribution list for each node
#' mydists <- list(b1="binomial", b2="binomial", g1="gaussian",
#'                 g2="gaussian", b3="binomial", g3="gaussian")
#'
#' # Structural constraints
#' # ban arc from b2 to b1
#' # always retain arc from g2 to g1
#'
#' ## parent limits
#' max.par <- list("b1"=2, "b2"=2, "g1"=2, "g2"=2, "b3"=2, "g3"=2)
#'
#' ## now build the cache of pre-computed scores accordingly to the structural constraints
#' if(requireNamespace("INLA", quietly = TRUE)){
#'   # Run only if INLA is available
#' res.c <- buildScoreCache(data.df=mydat, data.dists=mydists,
#'                          dag.banned= ~b1|b2, dag.retained= ~g1|g2, max.parents=max.par)
#' print(res.c)
#' }
print.abnCache <- function(x, digits = 3, ...){

  cat("Number of nodes in the network: ",max(x$children), ".\n\n", sep='')
  if(x$method=="bayes"){
    cat("Distribution of the marginal likelihood: \n")
    print(summary(x[["mlik"]]), digits=digits)
  }

  if(x$method=="mle"){
    cat(" Distribution of the aic: \n")
    print(summary(x[["aic"]]), digits=digits)

    cat("\n Distribution of the bic: \n")
    print(summary(x[["bic"]]), digits=digits)

    cat("\n Distribution of the mdl: \n")
    print(summary(x[["mdl"]]), digits=digits)
  }
  invisible(x)
}

##-------------------------------------------------------------------------
## abnHeuristic
##-------------------------------------------------------------------------
#' Print objects of class \code{abnHeuristic}
#' @param x Object of class \code{abnHeuristic}
#' @param digits number of digits of the results.
#' @param ... additional parameters. Not used at the moment.
#' @return prints the best score found and the distribution of the scores.
#' @export
print.abnHeuristic <- function(x, digits = 2L, ...){
  cat("Best DAG' score found with",x$algo,"algorithm with", x$num.searches,"different searches limited to" , x$max.steps,"steps:\n")
  print(max(unlist(x$scores)), digits=digits)

  cat("\n Score distribution: \n")
  print(summary(unlist(x[["scores"]])), digits=digits)

  invisible(x)
}

#' Plot objects of class \code{abnHeuristic}
#' @param x Object of class \code{abnHeuristic}
#' @param ... additional parameters. Not used at the moment.
#' @importFrom graphics par plot points title lines
#' @importFrom grDevices rgb
#' @return plot of the scores of the heuristic search.
#' @export
plot.abnHeuristic <- function(x, ...){
  # Keep old par() settings and restore them at the end
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  # Plot the scores
  df <- unlist(x$scores)

  par(mfrow=c(1,2))
  plot(NULL, lty=1, xlab="Index of heuristic search", ylab="BN score", ylim = range(df), xlim = c(1,length(df)))
  for(i in 1:length(df)){
    if(sum(i==order(df, decreasing = FALSE)[1:10])){
      points(x=i,y=df[i], type="p", pch=19, col=rgb(0,0,1, 0.8),lwd = 2)
    } else {
      points(x=i,y=df[i], type="p", pch=19, col=rgb(0,0,0, 0.3))
    }
  }
  points(x = which.max(df), y = df[which.max(df)], col="red", pch=19)
  title("Networks final score")


  L <- (x$detailed.score)

  test <- array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))

  plot(NULL,lty=1, xlab="Number of Steps",ylab="BN score", ylim = range(test), xlim = c(1,length(test[,,1])))
  for(i in 1:length(L)){
    if(sum(i==order(df,decreasing = FALSE)[1:10])){
      points(x=1:(length(test[,,1])),y=test[1,,i], type="l", lty=1, col=rgb(0,0,1, 0.8),lwd = 2)
    } else {
      points(x=1:(length(test[,,1])),y=test[1,,i], type="l", lty=1, col=rgb(0,0,0, 0.17))
    }
  }
  lines(x=1:(length(test[,,1])),y=test[1,,which.max(df)], type="l", col="red", lwd=3)
  title("Networks score trajectory")
  invisible(x)
}

##-------------------------------------------------------------------------
## abnHillClimber
##-------------------------------------------------------------------------
#' Print objects of class \code{abnHillClimber}
#' @param x Object of class \code{abnHillClimber}
#' @param digits number of digits of the results.
#' @param ... additional parameters. Not used at the moment.
#' @returns prints the consensus DAG and the class of the object.
#' @export
print.abnHillClimber <- function(x, digits = 3L, ...){
  print(x$consensus, digits = digits)
  cat("Consensus DAG from 'search.hillclimber'  (class 'abnHillClimber').\n")
  invisible(x)
}

#' Plot objects of class \code{abnHillClimber}
#' @param x Object of class \code{abnHillClimber}
#' @param ... additional parameters. Not used at the moment.
#' @importFrom grDevices dev.new dev.flush
#' @returns plot of the consensus DAG.
#' @export
plot.abnHillClimber <- function(x, ...){
  # Check if the object is of class abnHillClimber
  if(!inherits(x, "abnHillClimber")){
    stop("The object is not of class 'abnHillClimber'")
  }

  # Plot the consensus DAG
  g <- plotAbn(dag = x$dag, data.dists = x$score.cache$data.dists)

  # Return the plot
  invisible(g)
}


##-------------------------------------------------------------------------
## abnMostprobable
##-------------------------------------------------------------------------

#' Print objects of class \code{abnMostprobable}
#' @param x Object of class \code{abnMostprobable}
#' @param digits number of digits of the results.
#' @param ... additional parameters. Not used at the moment.
#' @returns prints the mostprobable consensus DAG.
#' @export
print.abnMostprobable <- function(x, digits = 3L, ...){

  print(x$dag, digits = digits)
  cat("Consensus DAG from 'mostprobable', can be use with 'fitAbn'.\n")
  invisible(x)
}

#' Print summary of objects of class \code{abnMostprobable}
#' @param object Object of class \code{abnMostprobable}
#' @param ... additional parameters. Not used at the moment.
#' @returns prints the mostprobable consensus DAG and the number of observations used to calculate it.
#' @export
summary.abnMostprobable <- function(object, ...){
  cat("Optimal DAG from 'mostProbable':\n")
  print(object$dag)
  cat( paste0("Calculated on ", dim(object$score.cache$data.df)[1], " observations.\n"))
  cat( paste0("(Cache length ", length(object$score.cache$mlik), '.)\n'))
  invisible( object)
}

#' Plot objects of class \code{abnMostprobable}
#' @param x Object of class \code{abnMostprobable}
#' @param ... additional parameters. Not used at the moment.
#' @importFrom grDevices dev.new dev.flush
#' @returns plot of the mostprobable consensus DAG.
#' @export
plot.abnMostprobable <- function(x, ...){
  # Check if the object is of class abnMostprobable
  if(!inherits(x, "abnMostprobable")){
    stop("The object is not of class 'abnMostprobable'")
  }

  # Plot the DAG
  g <- plotAbn(dag = x$dag, data.dists = x$score.cache$data.dists)

  # Return the plot
  invisible(g)
}

##-------------------------------------------------------------------------
## abnFit
##-------------------------------------------------------------------------
#' Print objects of class \code{abnFit}
#' @param x Object of class \code{abnFit}
#' @param digits number of digits of the results.
#' @param ... additional parameters. Not used at the moment.
#' @returns prints the parameters of the fitted model.
#' @export
print.abnFit <- function(x, digits = 3L, ...){

  if(x$method=="mle"){
    cat("The ABN model was fitted using a Maximum Likelihood Estimation (MLE) approach.\n", fill = TRUE)
    if(!is.null(x$group.var)){
      cat("The model is a Generalized Linear Mixed Model (GLMM) with the following grouping variable: \n", x$group.var, "\n", fill = TRUE, sep = "")
      cat("Fixed-effect parameters (mu):\n")
      print(x$mu, digits=digits)
      cat("Fixed-effect coefficients (betas):\n")
      print(x$betas, digits=digits)
      cat("Random-effects residuals (sigma):\n")
      print(x$sigma, digits=digits)
      cat("Random-effects intercepts (sigma_alpha):\n")
      print(x$sigma_alpha, digits=digits)
    } else {
      cat("Coefficients:\n")
      print(x$coef, digits=digits)
    }
  } else if(x$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. The estimated modes (the highest posterior density values of the parameters) are:\n\n")
    print(x$modes, digits=digits)
  } else {
    stop("The method used to fit the model is not recognized.")
  }
  cat("Number of nodes in the network: ", ifelse(x$method=="mle", length(x$coef), length(x$modes)), "\n")
  invisible(x)
}

#' Print summary of objects of class \code{abnFit}
#' @param object Object of class \code{abnFit}
#' @param digits number of digits of the results.
#' @param ... additional parameters. Not used at the moment.
#' @returns prints summary statistics of the fitted model.
#' @export
summary.abnFit <- function(object, digits = 3L, ...){

  if(object$method=="mle"){
    cat("The ABN model was fitted using an mle approach. The estimated coefficients are:\n")
    print(object$coef, digits=3)

    cat("Number of nodes in the network: ",length(object$modes), ".\n", sep='')

    cat("The AIC network score per node is: \n")
    print(unlist(object[["aicnode"]]), digits=digits)

    cat("\n The BIC network score per node is: \n")
    print(unlist(object[["bicnode"]]), digits=digits)

    cat("\n The MDL network score per node is: \n")
    print(unlist(object[["mdlnode"]]), digits=digits)
  }

  if(object$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. The estimated modes are:\n")
    print(object$modes, digits=digits)

   cat("Number of nodes in the network: ",length(object$modes), ".\n\n", sep='')

   cat("The network score per node is:\n")
   print(unlist(object[1:length(object$modes)]))
  }

  invisible(object)
}

#' Print coefficients of objects of class \code{abnFit}
#' @param object Object of class \code{abnFit}
#' @param digits number of digits of the results.
#' @param verbose print additional output.
#' @param ... additional parameters. Not used at the moment.
#' @returns prints the coefficients of the fitted model.
#' @export
coef.abnFit <- function(object, digits = 3L, verbose = TRUE, ...){
  if(object$method=="mle"){
    cat("The ABN model was fitted using an mle approach. The estimated coefficients are:\n")
    print(object$coef, digits=digits)
  }

  if(object$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. The estimated modes are:\n")
    print(object$modes, digits=digits)
  }

  invisible(object)

}

#' Print AIC of objects of class \code{abnFit}
#' @param object Object of class \code{abnFit}
#' @param digits number of digits of the results.
#' @param verbose print additional output.
#' @param ... additional parameters. Not used at the moment.
#' @returns prints the AIC of the fitted model.
#' @export
AIC.abnFit <- function(object, digits = 3L, verbose = TRUE, ...){

  if(object$method=="mle"){

    cat("The ABN model was fitted using an mle approach. The AIC network score per node is: \n")
    print(unlist(object[["aicnode"]]), digits=digits)

  }

  if(object$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. AIC does not make sense but the network score per node is is is:\n")
    print(unlist(object[1:length(object$modes)]))
  }

  invisible(object)
}

#' Print BIC of objects of class \code{abnFit}
#' @param object Object of class \code{abnFit}
#' @param digits number of digits of the results.
#' @param verbose print additional output.
#' @param ... additional parameters. Not used at the moment.
#' @returns prints the BIC of the fitted model.
#' @export
BIC.abnFit <- function(object, digits = 3L, verbose = TRUE, ...){

  if(object$method=="mle"){

    cat("The ABN model was fitted using an mle approach. The BIC network score per node is: \n")
    print(unlist(object[["bicnode"]]), digits=3)

  }

  if(object$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. BIC does not make sense but the network score per node is is is:\n")
    print(unlist(object[1:length(object$modes)]))
  }

  invisible(object)
}

#' Print logLik of objects of class \code{abnFit}
#' @param object Object of class \code{abnFit}
#' @param digits number of digits of the results.
#' @param verbose print additional output.
#' @param ... additional parameters. Not used at the moment.
#' @returns prints the logLik of the fitted model.
#' @export
logLik.abnFit <- function(object, digits = 3L, verbose = TRUE, ...){

  if(object$method=="mle"){

    cat("The ABN model was fitted using an mle approach. The loglikelihood network score per node is: \n")
    print(unlist(object[["mliknode"]]), digits=3)

  }

  if(object$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. Loglikelihood does not make sense but the network score per node is is is:\n")
    print(unlist(object[1:length(object$modes)]))
  }

  invisible(object)
}

#' Print family of objects of class \code{abnFit}
#' @param object Object of class \code{abnFit}
#' @param ... additional parameters. Not used at the moment.
#' @importFrom stats family
#' @returns prints the distributions for each variable of the fitted model.
#' @exportS3Method abn::family abnFit
family.abnFit <- function(object, ...){

  cat("All link functions are canonical: \n
      gaussian node = identy, binomial node = logit, Poisson node = log and multinomial node = logit.\n\n")

  print(unlist(object$abnDag$data.dists))

  invisible(object)
}

#' Print number of observations of objects of class \code{abnFit}
#' @param object Object of class \code{abnFit}
#' @param ... additional parameters. Not used at the moment.
#' @importFrom stats nobs
#' @returns prints the number of observations of the fitted model.
#' @exportS3Method abn::nobs abnFit
nobs.abnFit <- function(object, ...){
  nrow(object$abnDag$data.df)
}

#' Plot objects of class \code{abnFit}
#' @param x Object of class \code{abnFit}
#' @param ... additional parameters. Not used at the moment.
#' @importFrom methods hasArg
#' @importFrom stats fitted.values
#' @returns a plot of the fitted model.
#' @export
plot.abnFit <- function(x, ...){
  # Check if the object is of class abnFit
  if (!inherits(x, "abnFit")) stop("x must be an object of class abnFit")

  if (hasArg(fitted.values)) {
    g <- plotAbn(x$abnDag$dag, data.dists = x$abnDag$data.dists, ...)
  } else {
    if(x$method=="mle"){
      g <- plotAbn(x$abnDag$dag, data.dists = x$abnDag$data.dists, fitted.values = x$coef, ...)
    } else {
      g <- plotAbn(x$abnDag$dag, data.dists = x$abnDag$data.dists, fitted.values = x$modes, ...)
    }
  }
  invisible(g)
}
