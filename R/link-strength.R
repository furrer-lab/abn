#' Returns the strengths of the edge connections in a Bayesian Network learned from observational data
#'
#' A flexible implementation of multiple proxy for strength measures useful for
#' visualizing the edge connections in a Bayesian Network learned from observational data.
#'
#' @usage linkStrength(dag,
#'                     data.df = NULL,
#'                     data.dists = NULL,
#'                     method = c("mi.raw",
#'                                "mi.raw.pc",
#'                                "mi.corr",
#'                                "ls",
#'                                "ls.pc",
#'                                "stat.dist"),
#'                     discretization.method = "doane")
#'
#' @param dag a matrix or a formula statement (see details for format) defining
#' the network structure, a directed acyclic graph (DAG).
#' Note that rownames must be set or given in \code{data.dist} if the DAG is
#' given via a formula statement.
#' @param data.df a data frame containing the data used for learning each node,
#' binary variables must be declared as factors.
#' @param data.dists a named list giving the distribution for each node in the
#' network, see \sQuote{Details}.
#' @param method the method to be used. See \sQuote{Details}.
#' @param discretization.method a character vector giving the discretization
#' method to use. See \code{\link{discretization}}.
#'
#' @details
#' This function returns multiple proxies for estimating the connection strength
#' of the edges of a possibly discretized Bayesian network's data set.
#' The returned connection strength measures are: the Raw Mutual Information
#' (\code{mi.raw}), the Percentage Mutual information (\code{mi.raw.pc}),
#' the Raw Mutual Information computed via correlation (\code{mi.corr}),
#' the link strength (\code{ls}), the percentage link strength (\code{ls.pc})
#' and the statistical distance (\code{stat.dist}).
#'
#' The general concept of entropy is defined for probability distributions.
#' The probability is estimated from data using frequency tables.
#' Then the estimates are plug-in in the definition of the entropy to return
#' the so-called empirical entropy. A standard known problem of empirical entropy
#'  is that the estimations are biased due to the sampling noise.
#'  This is also known that the bias will decrease as the sample size increases.
#'  The mutual information estimation is computed from the observed frequencies
#'  through a plug-in estimator based on entropy.
#'  For the case of an arc going from the node X to the node Y and the remaining
#'  set of parent of Y is denoted as Z.
#'
#'  The mutual information is defined as I(X, Y) = H(X) + H(Y) - H(X, Y),
#'  where H() is the entropy.
#'
#'  The Percentage Mutual information is defined as PI(X,Y) = I(X,Y)/H(Y|Z).
#'
#'  The Mutual Information computed via correlation is defined as
#'  MI(X,Y) = -0.5 log(1-cor(X,Y)^2).
#'
#'  The link strength is defined as LS(X->Y) = H(Y|Z)-H(Y|X,Z).
#'
#'  The percentage link strength is defined as PLS(X->Y) = LS(X->Y) / H(Y|Z).
#'
#'  The statistical distance is defined as SD(X,Y) = 1- MI(X,Y) / max(H(X),H(Y)).
#'
#' @return The function returns a named matrix with the requested metric.
#' @references
#' Boerlage, B. (1992).  Link strength in Bayesian networks. Diss. University of British Columbia.
#' Ebert-Uphoff, Imme. "Tutorial on how to measure link strengths in discrete Bayesian networks." (2009).
#'
#' @export
#' @keywords utilities
#' @importFrom stats cor
#' @examples
#' # Gaussian
#' N <- 1000
#' mydists <- list(a="gaussian",
#'                 b="gaussian",
#'                 c="gaussian")
#' a <- rnorm(n = N, mean = 0, sd = 1)
#' b <- 1 + 2*rnorm(n = N, mean = 5, sd = 1)
#' c <- 2 + 1*a + 2*b + rnorm(n = N, mean = 2, sd = 1)
#' mydf <- data.frame("a" = a,
#'                    "b" = b,
#'                    "c" = c)
#' mycache.mle <- buildScoreCache(data.df = mydf,
#'                                data.dists = mydists,
#'                                method = "mle",
#'                                max.parents = 2)
#' mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
#' linkstr <- linkStrength(dag = mydag.mp$dag,
#'                         data.df = mydf,
#'                         data.dists = mydists,
#'                         method = "ls",
#'                         discretization.method = "sturges")
linkStrength <- function(dag, data.df = NULL, data.dists = NULL, method = c("mi.raw", "mi.raw.pc", "mi.corr", "ls", "ls.pc", "stat.dist"), discretization.method = "doane") {

    group.var <- NULL    # hard coded, no choice here
#    if (is.matrix(dag)) {
#        name <- names(dag)
#    } else {
#        name <- names(data.dists)
#    }
    # Contains Rgraphviz stuff

    ## dag transformation
#    if (!is.null(dag)) {
        if (is.matrix(dag)) {
            ## run a series of checks on the DAG passed
 #           dag <- abs(dag)
 #            diag(dag) <- 0
            dag <- check.valid.dag(dag = dag, data.df = data.df, is.ban.matrix = FALSE, group.var = group.var)
            ## naming
 #            if (is.null(colnames(dag))) {
 #                colnames(dag) <- name
 #                rownames(dag) <- name
 #           }
        } else {
          if (grepl("~", as.character(dag)[1], fixed = T)) {
                dag <- formula_abn(f = dag, name = names(data.dists))
                ## run a series of checks on the DAG passed
                dag <- check.valid.dag(dag = dag, is.ban.matrix = FALSE, group.var = group.var)
#            }
#        }
         } else {
            stop("Dag specification must either be a matrix or a formula expression")
         }
        }


    ## rows and columns
    n.row <- length(dag[1, ])
    n.col <- length(dag[, 1])


    ## MUTUAL INFORMATION

    if (method == "mi.raw" | method == "mi.raw.pc")
        {

            mi.m <- matrix(data = 0, nrow = n.row, ncol = n.col)
            for (i in 1:n.col) {
                for (j in 1:n.row) {
                  if (dag[j, i] != 0) {

                    mi.m[j, i] <- miData(freqs.table = discretization(data.df = data.df[, c(i, j)], data.dists = data.dists[c(i, j)], discretization.method), method = method)

                  }
                }
            }

            return(mi.m)
        }  #EOF: mi.raw

    if (method == "mi.corr")
        {

            mi.m <- matrix(data = 0, nrow = n.row, ncol = n.col)
            for (i in 1:n.col) {
                for (j in 1:n.row) {
                  if (dag[j, i] != 0) {

                    mi.m[j, i] <- -0.5 * log(1 - cor(x = data.df[, j], y = data.df[, i])^2)
                  }
                }
            }

            return(mi.m)
        }  #EOF: mi.theo

    if (method == "ls")
        {
            ## Formula: I(X;Y|Z)=H(X,Z)+H(Y,Z)-H(X,Y,Z)-H(Z)=H(X|Z)-H(X|Y,Z) with Z all other parents of Y i -> j = p -> c
            ls.m <- matrix(data = 0, nrow = n.row, ncol = n.col)
            for (i in 1:n.col) {
                for (j in 1:n.row) {
                  if (dag[j, i] != 0) {
                    parent.list <- t(dag[j, ])
                    parent.list[parent.list != 0] <- 1

                    parent.list[j] <- 0
                    parent.list[i] <- 0

                    parent.list.x <- parent.list
                    parent.list.x[i] <- 1
                    parent.list.y <- parent.list
                    parent.list.y[j] <- 1

                    if (sum(parent.list) != 0) {
                      # ls.m <- rbind(parent.list)
                      ls.m[j, i] <- entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.x)], data.dists = data.dists[as.logical(parent.list.x)], discretization.method = discretization.method)) +
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method)) -
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list)], data.dists = data.dists[as.logical(parent.list)], discretization.method = discretization.method))
                    } else {
                      ls.m[j, i] <- entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.x)], data.dists = data.dists[as.logical(parent.list.x)], discretization.method = discretization.method)) +
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method))
                    }
                  }
                }
            }
            return(ls.m)
        }  #EOF: link strength


    if (method == "ls.pc")
        {
            ## Formula: I(X;Y|Z)/H(Y|Z)=(H(X,Z)+H(Y,Z)-H(X,Y,Z)-H(Z))/H(Y|Z)=(H(X|Z)-H(X|Y,Z))/H(Y|Z). H(Y|Z)=H(Y,Z)-H(Z)

            ls.pc.m <- matrix(data = 0, nrow = n.row, ncol = n.col)
            for (i in 1:n.col) {
                for (j in 1:n.row) {
                  if (dag[j, i] != 0) {
                    parent.list <- t(dag[j, ])
                    parent.list[parent.list != 0] <- 1

                    parent.list[j] <- 0
                    parent.list[i] <- 0

                    parent.list.x <- parent.list
                    parent.list.x[i] <- 1
                    parent.list.y <- parent.list
                    parent.list.y[j] <- 1
                    parent.list.xy <- parent.list.y
                    parent.list.xy[i] <- 1

                    if (sum(parent.list) != 0) {
                      # ls.m <- rbind(parent.list)
                      ls.pc.m[j, i] <- entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.x)], data.dists = data.dists[as.logical(parent.list.x)], discretization.method = discretization.method)) +
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method)) -
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.xy)], data.dists = data.dists[as.logical(parent.list.xy)], discretization.method = discretization.method)) -
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list)], data.dists = data.dists[as.logical(parent.list)], discretization.method = discretization.method))
                      tmp <- (entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method)) -
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list)], data.dists = data.dists[as.logical(parent.list)], discretization.method = discretization.method)))
                      if (tmp != 0) {
                        ls.pc.m[j, i] <- ls.pc.m[j, i]/tmp
                      }

                    } else {
                      ls.pc.m[j, i] <- entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.x)], data.dists = data.dists[as.logical(parent.list.x)], discretization.method = discretization.method)) +
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method)) -
                        entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.xy)], data.dists = data.dists[as.logical(parent.list.xy)], discretization.method = discretization.method))
                      tmp <- (entropyData(freqs.table = discretization(data.df = data.df[, as.logical(parent.list.y)], data.dists = data.dists[as.logical(parent.list.y)], discretization.method = discretization.method)))
                      if (tmp != 0) {
                        ls.pc.m[j, i] <- ls.pc.m[j, i]/tmp
                      }
                    }
                  }
                }
            }
            return(ls.pc.m)
        }  #EOF: link strength


    if (method == "stat.dist")
        {
            ## Formula: statistical distance = 1- MI(X,Y)/max(H(X), H(Y))

            stat.dist <- matrix(data = 0, nrow = n.row, ncol = n.col)
            for (i in 1:n.col) {
                for (j in 1:n.row) {
                  if (dag[j, i] != 0) {

                    stat.dist[j, i] <- miData(freqs.table = discretization(data.df = data.df[, c(i, j)], data.dists = data.dists[c(i, j)], discretization.method), method = "mi.raw")

                  }
                }
            }
            return(stat.dist)
        }  #EOF: statistical distance


    # stat.dist <- matrix(data = 0,nrow = n.row,ncol = n.col) for(i in 1:n.col){ for(j in 1:n.row){ if(dag[j,i]!=0){ stat.dist[j,i] <- mi.data(freqs.table = discretization(data.df =
    # data.df[,c(i,j)],data.dists = data.dists[c(i,j)],discretization.method), method = method) #stat.dist[j,i] <- stat.dist[j,i]/max(entropy.data(freqs.table = discretization(data.df =
    # data.df[,i],data.dists=data.dists[i],discretization.method = discretization.method)), entropy.data(freqs.table = discretization(data.df =
    # data.df[,j],data.dists=data.dists[j],discretization.method = discretization.method))) #stat.dist[j,i] <- 1-stat.dist[j,i] }}} return(mi.m) }#EOF: statistical distance


}  #EOF
