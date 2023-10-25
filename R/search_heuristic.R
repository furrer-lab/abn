#' A family of heuristic algorithms that aims at finding high scoring directed acyclic graphs
#'
#' A flexible implementation of multiple greedy search algorithms to find high scoring network (DAG)
#'
#' @usage searchHeuristic(score.cache, score = "mlik",
#'                        num.searches = 1, seed = 42L, start.dag = NULL,
#'                        max.steps = 100,
#'                        algo = "hc", tabu.memory = 10, temperature = 0.9,
#'                        verbose = FALSE, ...)
#'
#' @param score.cache output from \code{buildScoreCache()}.
#' @param score which score should be used to score the network. Possible choices are \code{aic, bic, mdl, mlik}.
#' @param num.searches a positive integer giving the number of different search to run, see details.
#' @param seed a non-negative integer which sets the seed.
#' @param start.dag a DAG given as a matrix, see details for format, which can be used to explicity provide a starting point for the structural search.
#' @param max.steps a constant giving the number of search steps per search, see details.
#' @param algo which heuristic algorithm should be used. Possible choices are: \code{hc, tabu, sa}.
#' @param tabu.memory a non-negative integer number to set the memory of the \code{tabu} search.
#' @param temperature a real number giving the update in temperature for the \code{sa} (simulated annealing) search algorithm.
#' @param verbose if TRUE then provides some additional output.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' This function is a flexible implementation of multiple greedy heuristic algorithms,
#' particularly well adapted to the \code{abn} framework.
#' It targets multi-random restarts heuristic algorithms.
#' The user can select the \code{num.searches} and the maximum number of steps
#' within by \code{max.steps}. The optimization algorithm within each search is
#' relatively rudimentary.
#'
#' The function \code{searchHeuristic} is different from the
#' \code{\link{searchHillClimber}} in the sense that this function is fully
#' written in R, whereas the \code{\link{searchHillClimber}} is written in C
#' and thus expected to be faster. The function \code{\link{searchHillClimber}}
#' at each hill-climbing step consider every information from the pre-computed
#' scores cache while the function \code{\link{searchHeuristic}} performs a local
#' stepwise optimization. This function chooses a structural move (or edge move)
#' and compute the score's change. On this point, it is closer to the MCMCMC
#' algorithm from Madigan and York (1995) and Giudici and Castelo (2003)
#' with a single edge move.
#'
#' If the user select \code{random}, then a random valid DAG is selected.
#' The routine used favourise low density structure.
#' The function implements three algorithm selected with the
#' parameter \code{algo}: \code{hc}, \code{tabu} or \code{sa}.
#'
#' If \code{algo=hc}:
#' The Hill-climber algorithm (\code{hc}) is a single move algorithm.
#' At each Hill-climbing step within a search an arc is attempted to be added.
#' The new score is computed and compared to the previous network's score.
#'
#' If \code{algo=tabu}:
#' The same algorithm is as with \code{hc} is used, but a list of banned moves
#' is computed. The parameter \code{tabu.memory} controls the length of the tabu
#' list. The idea is that the classical Hill-climber algorithm is inefficient
#' when it should cross low probability regions to unblock from a local maximum
#' and reaching a higher score peak. By forcing the algorithm to choose some not
#' already used moves, this will force the algorithm to escape the local maximum.
#'
#' If \code{algo=sa}:
#' This variant of the heuristic search algorithm is based on simulated annealing
#' described by Metropolis et al. (1953).
#' Some accepted moves could result in a decrease of the network score.
#' The acceptance rate can be monitored with the parameter \code{temperature}.
#'
#' @return
#' An object of class \code{abnHeuristic} (which extends the class \code{abnLearnd}) and contains list with entires:
#' \describe{
#' \item{dags}{a list of DAGs}
#' \item{scores}{a vector giving the network score for the locally optimal network for each search}
#' \item{detailed.score}{a vector giving the evolution of the network score for the all the random restarts}
#' \item{score}{a number giving the network score for the locally optimal network}
#' \item{score.cache}{the pre-computed cache of scores}
#' \item{num.searches}{a numeric giving the number of random restart}
#' \item{max.steps}{a numeric giving the maximal number of optimization steps within each search}
#' \item{algorithm}{a character for indicating the algorithm used}
#' }
#' @export
#' @importFrom utils tail
#' @references
#' Heckerman, D., Geiger, D. and Chickering, D. M. (1995). Learning Bayesian networks: The combination of knowledge and statistical data. \emph{Machine Learning}, 20, 197-243.
#' Madigan, D. and York, J. (1995) "Bayesian graphical models for discrete data". International Statistical Review, 63:215232.
#' Giudici, P. and Castelo, R. (2003). "Improving Markov chain Monte Carlo model search for data mining". Machine Learning, 50:127158.
#' Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., & Teller, E. (1953). "Equation of state calculations by fast computing machines". The journal of chemical physics, 21(6), 1087-1092.
#'
#' @examples
#' \dontrun{
#' ##############################################
#' ## example: use built-in simulated data set
#' ##############################################
#'
#' mydat <- ex1.dag.data ## this data comes with abn see ?ex1.dag.data
#'
#' ## setup distribution list for each node
#' mydists<-list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
#'               p2="poisson", b3="binomial", g2="gaussian", b4="binomial",
#'               b5="binomial", g3="gaussian")
#'
#' mycache <- buildScoreCache(data.df = mydat, data.dists = mydists, max.parents = 2)
#'
#' ## Now peform 10 greedy searches
#' heur.res <- searchHeuristic(score.cache = mycache, data.dists = mydists,
#'                             start.dag = "random", num.searches = 10,
#'                             max.steps = 50)
#'
#' ## Plot (one) dag
#' plotAbn(heur.res$dags[[1]], data.dists = mydists)
#' }
searchHeuristic <- function(score.cache, score="mlik",
                            num.searches=1, seed=42L, start.dag=NULL,
                            max.steps=100,
                            algo="hc", tabu.memory=10, temperature=0.9,
                            verbose=FALSE,  ...) {

    if (!inherits(score.cache,"abnCache")) {
        stop("score.cache should be an object of class 'abnCache' ")
    }
    score <- c("mlik","aic","bic",
               "mdl")[pmatch(tolower(score), c("mlik","aic","bic","mdl"))][1]
    if (is.na(score)) stop("wrong specification of 'score'.")



    data.dists <- score.cache$data.dists
    if(is.vector(score.cache$max.parents)){
        max.parents <- max(score.cache$max.parents)
    } else {
    max.parents <- score.cache$max.parents
    }
    dag.retained <- score.cache$dag.retained
    dag.banned <- score.cache$dag.banned

    ## function
    resample <- function(x, ...) x[sample.int(length(x), ...)]
    n.var <- length(data.dists)

    ## output
    out.dags <- NULL
    out.scores <- list()


    out.detailed <- NULL


    ## seeding
    set.seed(seed=seed)

    for (searchIndex in 1:num.searches) {
        if (verbose)
            cat("processing search...", searchIndex, "\n")

        temperature.update <- 1

        ## Initializing matrix
        dag.tmp <- matrix(data=0, nrow=n.var, ncol=n.var)

        ## start zero matrix
        if (!is.null(start.dag)) {
            colnames(dag.tmp) <- rownames(dag.tmp) <- sample(1:n.var)
        }

        ## start random matrix
        if (is.null(start.dag)) {
            start.dag <- "random"
        }
        if (start.dag == "random") {
            vec.tmp <- c(rep(1, max.parents), rep(0, 2 * n.var - max.parents))
            for (lines in 1:(n.var - 1)) {
                dag.tmp[1 + lines, 1:lines] <- sample(vec.tmp)[1:lines]
            }
            colnames(dag.tmp) <- rownames(dag.tmp) <- sample(1:n.var)
        }
        if (is.matrix(start.dag))
            dag.tmp <- start.dag

        ## score init dag

        score.init <- vector(mode="numeric", length=n.var)

        if (score %in% c("bic", "aic", "mdl")) {
            sc <- cbind(score.cache$node.defn[, as.numeric(colnames(dag.tmp))], -score.cache[[score]])
        } else {
            sc <- cbind(score.cache$node.defn[, as.numeric(colnames(dag.tmp))], score.cache[[score]])
        }

        for (lines in 1:n.var) {

            sc.tmp <- sc[score.cache$children == as.numeric(colnames(dag.tmp)[lines]), ]

            score.init[lines] <- min(sc.tmp[which(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[lines, ])))), n.var + 1])

        }

        # dag.init <- dag.tmp start hc
        steps <- 1
        score.tmp <- score.init
        dag <- dag.tmp

        ## hill climbing algorithm
        if (algo == "hc") {


                out <- NULL


            while (steps < max.steps) {

                dag.tmp <- dag

                # change dag
                y <- resample(2:n.var, 1)
                x <- resample(1:(y - 1), 1)

                if (dag.tmp[y, x] == 0) {
                  if (sum(dag.tmp[y, ]) == max.parents) {
                    x <- which.max(unname(dag.tmp[y, ]))
                    dag.tmp[y, x] <- 0
                  }
                  dag.tmp[y, x] <- 1

                } else {

                  dag.tmp[y, x] <- 0
                }

                # compute score
                sc.tmp <- sc[score.cache$children == as.numeric(colnames(dag.tmp)[y]), ]
                score.test <- min(sc.tmp[which(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[y, ])))), n.var + 1])
                if (score.tmp[y] < score.test) {
                  score.tmp[y] <- score.test
                  dag <- dag.tmp
                }
                steps <- steps + 1


                  out <- cbind(out, sum(score.tmp))


            }  #eow

        }

        ## tabu algorithm
        if (algo == "tabu") {


                out <- NULL


            memory <- matrix(data=0, nrow=tabu.memory, ncol=2)
            prob <- c(0, 0)

            while (steps < max.steps) {

                dag.tmp <- dag

                ## tabu step
                while (sum(apply(tail(memory, tabu.memory), 1, function(x) identical((x), as.numeric(prob)))) > 0) {
                  # change dag
                  y <- resample(2:n.var, 1)
                  x <- resample(1:(y - 1), 1)

                  prob <- c(x, y)

                }

                memory <- rbind(memory, prob)

                if (dag.tmp[y, x] == 0) {
                  if (sum(dag.tmp[y, ]) == max.parents) {
                    x <- which.max(unname(dag.tmp[y, ]))
                    dag.tmp[y, x] <- 0
                  }
                  dag.tmp[y, x] <- 1

                } else {

                  dag.tmp[y, x] <- 0
                }

                # compute score
                sc.tmp <- sc[score.cache$children == as.numeric(colnames(dag.tmp)[y]), ]
                score.test <- min(sc.tmp[which(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[y, ])))), n.var + 1])
                if (score.tmp[y] < score.test) {
                  score.tmp[y] <- score.test
                  dag <- dag.tmp
                }
                steps <- steps + 1


                  out <- cbind(out, sum(score.tmp))


            }  #eow
        }

        ## simulated annealing (threashold acceptance)
        if (algo == "sa") {


                out <- NULL


            temperature.update <- 1


            while (steps < max.steps) {

                dag.tmp <- dag

                # change dag
                y <- resample(2:n.var, 1)
                x <- resample(1:(y - 1), 1)

                if (dag.tmp[y, x] == 0) {
                  if (sum(dag.tmp[y, ]) == max.parents) {
                    x <- which.max(unname(dag.tmp[y, ]))
                    dag.tmp[y, x] <- 0
                  }
                  dag.tmp[y, x] <- 1

                } else {

                  dag.tmp[y, x] <- 0
                }

                # compute score
                sc.tmp <- sc[score.cache$children == as.numeric(colnames(dag.tmp)[y]), ]
                score.test <- min(sc.tmp[which(apply(sc.tmp, 1, function(x) identical(unname(x[1:n.var]), unname(dag.tmp[y, ])))), n.var + 1])

                if (score.tmp[y] < score.test) {
                  score.tmp[y] <- score.test
                  dag <- dag.tmp

                }
                ## threshold acceptance (Dueck and Scheuer 1990) Dueck, G. and Scheuer, T. 'Threshold Accepting: A General Purpose Optimization Algorithm Appearing Superior to Simulated Annealing.' J. Comp.
                ## Phys. 90, 161-175, 1990.
                if (((score.tmp[y] - score.test) > 0) && (abs(score.tmp[y] - score.test) < (temperature.update * abs(score.init)/n.var))) {

                  if (verbose) {
                    print("Accepting negative move")
                  }

                  score.tmp[y] <- score.test
                  dag <- dag.tmp
                }

                steps <- steps + 1
                temperature.update <- temperature.update * temperature


                  out <- cbind(out, sum(score.tmp))


            }  #eow

        }

        ## output of one search
        dag <- dag[as.character(1:n.var), as.character(1:n.var)]
        colnames(dag) <- rownames(dag) <- names(data.dists)

        out.dags[[searchIndex]] <- (dag)
        out.scores[[searchIndex]] <- sum(score.tmp)


        out.detailed[[searchIndex]] <- out


    }  #eos


    ## return

    if (score %in% c("bic", "aic", "mdl")) {

            out <- list(dags=out.dags,
                 scores=lapply(X=out.scores, FUN=function(x) {
                     -x
                 }), detailed.score=lapply(X=out.detailed, FUN=function(x) {
                     -x
                 }), score=score, score.cache=score.cache, num.searches=num.searches, max.steps=max.steps,algorithm=algo
            )

            class(out) <- c("abnHeuristic", "mle")

            return(out)
    } else {
        out <- list(dags=out.dags, scores=out.scores, detailed.score=out.detailed,
                    score=score, score.cache=score.cache, num.searches=num.searches,
                    max.steps=max.steps,algorithm=algo)
        class(out) <- c("abnHeuristic", "bayes")
        return(out)
    }


}  #eof
