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
summary.abnDag <- function(object, ...) {
  su <- infoDag(object$dag)
  return(su)
}

#' Plots DAG from an object of class \code{abnDag}
#'
#' @param x Object of class \code{abnDag}
#' @param new defaults to \code{TRUE} for opening a new plot.
#' @param ... additional parameters. Not used at the moment.
#'
#' @return \code{Rgraphviz::plot}
#' @concept DAG
#' @export
#' @importClassesFrom graph graphAM
#' @importFrom grDevices dev.new dev.flush
#' @examples
#' mydag <- createAbnDag(dag = ~a+b|a, data.df = data.frame("a"=1, "b"=1))
#' plot(mydag)
plot.abnDag <- function(x, new=TRUE, ...){
  invisible({
    if (new) dev.new()
    on.exit(dev.flush())

    # Rgraphviz:
    mygraph <- new("graphAM", adjMat = t(x$dag), edgemode = "directed")
    g <- Rgraphviz::plot(x = mygraph)
  })
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
#'
#' res.c <- buildScoreCache(data.df=mydat, data.dists=mydists,
#'                          dag.banned= ~b1|b2, dag.retained= ~g1|g2, max.parents=max.par)
#' print(res.c)
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
#' @export
plot.abnHeuristic <- function(x, ...){

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
#' @export
print.abnHillClimber <- function(x, digits = 3L, ...){
  print(x$consensus, digits = digits)
  cat("Consensus DAG from 'search.hillclimber'  (class 'abnHillClimber').\n")
  invisible(x)
}

#' Plot objects of class \code{abnHillClimber}
#' @param x Object of class \code{abnHillClimber}
#' @param new defaults to \code{TRUE} for opening a new plot.
#' @param ... additional parameters. Not used at the moment.
#' @importFrom grDevices dev.new dev.flush
#' @export
plot.abnHillClimber <- function(x, new=TRUE, ...){

  if (new) dev.new()
  on.exit(dev.flush())

  # Rgraphviz
  mygraph <- new("graphAM", adjMat = x$consensus, edgemode = "directed")
  g <- Rgraphviz::plot(x = mygraph)
  invisible(g)
}


##-------------------------------------------------------------------------
## abnMostprobable
##-------------------------------------------------------------------------

#' Print objects of class \code{abnMostprobable}
#' @param x Object of class \code{abnMostprobable}
#' @param digits number of digits of the results.
#' @param ... additional parameters. Not used at the moment.
#' @export
print.abnMostprobable <- function(x, digits = 3L, ...){

  print(x$dag, digits = digits)
  cat("Consensus DAG from 'mostprobable', can be use with 'fitabn'.\n")
  invisible(x)
}

#' Print summary of objects of class \code{abnMostprobable}
#' @param object Object of class \code{abnMostprobable}
#' @param ... additional parameters. Not used at the moment.
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
#' @param new defaults to \code{TRUE} for opening a new plot.
#' @param ... additional parameters. Not used at the moment.
#' @importFrom grDevices dev.new dev.flush
#' @export
plot.abnMostprobable <- function(x, new=TRUE, ...){

  if (new) dev.new()
  on.exit(dev.flush())

  # Rgraphviz:
  mygraph <- new("graphAM", adjMat = t(x$dag), edgemode = "directed")
  g <- Rgraphviz::plot(x = mygraph)
  invisible(g)
}

##-------------------------------------------------------------------------
## abnFit
##-------------------------------------------------------------------------
#' Print objects of class \code{abnFit}
#' @param x Object of class \code{abnFit}
#' @param digits number of digits of the results.
#' @param ... additional parameters. Not used at the moment.
#' @export
print.abnFit <- function(x, digits = 3L, ...){

  if(x$method=="mle"){
    cat("The ABN model was fitted using an mle approach.", fill = TRUE)
    if(!is.null(x$group.var)){
      cat("GLMM with the grouping variable ", x$group.var, ".", fill = TRUE, sep = "")
      cat("The estimated fixed-effect parameters are:\n\n")
      print(x$mu, digits=digits)
      print(x$betas, digits=digits)
      cat("The estimated random-effect parameters are:\n\n")
      print(x$sigma, digits=digits)
      print(x$sigma_alpha, digits=digits)
      cat("Number of nodes in the network: ",length(x$abnDag$data.dists), " plus ", length(x$group.var), " grouping variable.", fill = TRUE, sep = "")
    } else {
      # no grouping
      cat("The estimated coefficients are:\n\n")
      print(x$coef, digits=digits)
      cat("Number of nodes in the network: ",length(x$coef), ".\n", sep='')
    }
  }

  if(x$method=="bayes"){
    cat("The ABN model was fitted using a Bayesian approach. The estimated modes are:\n\n")
    print(x$modes, digits=digits)
    cat("Number of nodes in the network: ",length(x$modes), ".\n", sep='')
  }

  invisible(x)
}

#' Print summary of objects of class \code{abnFit}
#' @param object Object of class \code{abnFit}
#' @param digits number of digits of the results.
#' @param ... additional parameters. Not used at the moment.
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
#' @exportS3Method abn::nobs abnFit
nobs.abnFit <- function(object, ...){
  nrow(object$abnDag$data.df)
}

#' Plot objects of class \code{abnFit}
#' @param x Object of class \code{abnFit}
#' @param which defaults to "abnFit".
#' @param ... additional parameters. Not used at the moment.
#' @importFrom methods hasArg
#' @importFrom stats fitted.values
#' @export
plot.abnFit <- function(x, which ="abnFit", ...){

  if (which != "abnFit") stop('Function type not implemented yet. Use which="abnFit"')

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
