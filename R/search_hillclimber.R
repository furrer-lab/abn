#' Find high scoring directed acyclic graphs using heuristic search.
#'
#' Find high scoring network (DAG) structures using a random re-starts greedy hill-climber heuristic search.
#'
#' @usage searchHillClimber(score.cache, score = "mlik", num.searches = 1, seed = 42,
#'                          start.dag = NULL, support.threshold = 0.5, timing.on = TRUE,
#'                          dag.retained = NULL, verbose = FALSE, ...)
#'
#' @param score.cache output from \code{buildScoreCache()}.
#' @param score character giving which network score should be used to
#' select the structure. Currently \code{'mlik'} only.
#' @param num.searches number of times to run the search.
#' @param seed non-negative integer which sets the seed in the GSL random number generator.
#' @param start.dag a DAG given as a matrix, see details for format, which can be used to provide a starting point for the structural search explicitly.
#' @param support.threshold the proportion of search results - each locally optimal DAG - in which each arc must appear to be a part of the consensus network.
#' @param timing.on extra output in terms of duration computation.
#' @param dag.retained a DAG given as a matrix, see details for format. This is necessary if the score.cache was created using an explicit retain matrix, and the same retain matrix should be used here. dag.retained is used by the algorithm which generates the initial random DAG to ensure that the necessary arcs are retained.
#' @param verbose extra output.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The procedure runs a greedy hill-climbing search similar, but not identical,
#' to the method presented initially in Heckerman et al. 1995.
#' (Machine Learning, 20, 197-243).
#' Each search begins with a randomly chosen DAG structure where this is
#' constructed in such a way as to attempt to choose a DAG uniformly from
#' the vast landscape of possible structures. The algorithm used is as follows:
#' given a node cache (from \code{\link{buildScoreCache}},
#' then within this set of all allowed local parent combinations,
#' a random combination is chosen for each node.
#' This is then combined into a full DAG, which is then checked for cycles,
#' where this check iterates over the nodes in a random order.
#' If all nodes have at least one child (i.e., at least one cycle is present),
#' then the first node examined has all its children removed, and the check for
#' cycles is then repeated. If this has removed the only cycle present,
#' then this DAG is used at the starting point for the search,
#' if not, a second node is chosen (randomly) and the process is then
#' repeated until a DAG is obtained.
#'
#' The actual hill-climbing algorithm itself differs slightly from that presented
#' in Heckerman et al. as a full cache of all possible local combinations are available.
#' At each hill-climbing step, everything in the node cache is considered.
#' In other words, all possible single swaps between members of cache currently present
#' in the DAG and those in the full cache. The single swap, which provides the greatest
#' increase in goodness of fit is chosen. A single swap here is the removal
#' or addition of any one node-parent combination present in the cache while
#' avoiding a cycle. This means that as well as all single arc changes
#' (addition or removal), multiple arc changes are also considered at each same step,
#' note however that arc reversal in this scheme takes two steps (as this requires
#' first removal of a parent arc from one node and then addition of a parent arc
#' to a different node). The original algorithm perturbed the current DAG by only
#' a single arc at each step but also included arc reversal.
#' The current implementation may not be any more efficient than the original
#' but is arguably more natural given a pre-computed cache of local scores.
#'
#' A start DAG may be provided in which case num.searches must
#' equal 1 - this option is really just to provide a local search
#' around a previously identified optimal DAG.
#'
#' This function is designed for two different purposes:
#' i) interactive visualization; and
#' ii) longer batch processing. The former provides easy visual "eyeballing" of
#' data in terms of its majority consensus network (or similar threshold),
#' which is a graphical structure which comprises of all arcs which feature in
#' a given proportion (\code{support.threshold}) of locally optimal DAGs already
#' identified during the run. The general hope is that this structure will
#' stabilize - become fixed - relatively quickly, at least for problems with
#' smaller numbers of nodes.
#'
#' @return A list with entries:
#' \describe{
#' \item{init.score}{a vector giving network score for initial network from which the search commenced}
#' \item{final.score}{a vector giving the network score for the locally optimal network}
#' \item{init.dag}{list of matrices, initial DAGs}
#' \item{final.dag}{list of matrices, locally optimal DAGs}
#' \item{consensus}{a matrix holding a binary graph, not necessary a DAG!}
#' \item{support.threshold}{percentage supported used to create consensus matrix}
#' }
#'
#' @references Lewis, F. I., and McCormick, B. J. J. (2012). Revealing the complexity of health determinants in resource poor settings. \emph{American Journal Of Epidemiology}. DOI:10.1093/aje/KWS183).
#' @export searchHillClimber
searchHillClimber <- function(score.cache, score="mlik",
                              num.searches=1, seed=42, start.dag=NULL,
                              support.threshold=0.5, timing.on=TRUE, dag.retained=NULL,
                              verbose=FALSE, ...) {

    if (!inherits(score.cache,"abnCache")) {
        stop("score.cache should be an object of class 'abnCache' ")
    }
    score <- c("mlik")[pmatch(tolower(score), c("mlik"))][1]
    if (is.na(score)) stop("not implemented 'score'.")


    data.df <- score.cache$data.df  ## n.b. this might be adjusted from original data.df depending on the group variable
    ## check that all nodes are included in the cache - easy to forget if parallelising
    #if (length(unique(score.cache$children)) != dim(data.df)[2]) {
    #    stop("all nodes must be included in the cache!")
    #}
    ## check valid arguments
    if (!is.numeric(num.searches) || length(num.searches) != 1) {
        stop("invalid num.searches")
    }
    if (seed < 0) {
        stop("seed must be a positive integer")
    }
    if (!is.logical(verbose)) {
        stop("invalid verbose - must be logical")
    }
    if (!is.logical(timing.on)) {
        stop("invalid timing.on")
    }
    # if(!is.logical(trace)){ stop('invalid trace value must be logical')}
    if (!is.numeric(support.threshold) || support.threshold < 0 || support.threshold > 1) {
        stop("invalid support.threshold value")
    }
    # if(trace==TRUE && !interactive()){stop('asked for graphical trace but not running interactively')} create a random network from which to commence searching. This is a random combination chosen
    # from all valid parent combination i.e. in the node cache and then checked for cyclicity

    if (!is.null(dag.retained)) {
        check.valid.dag(dag.retained, data.df=score.cache$data.df, is.ban.matrix=FALSE, group.var=NULL)
    } else {
        dag.retained <- check.valid.dag(dag.retained, data.df=score.cache$data.df, is.ban.matrix=FALSE, group.var=NULL)
    }  ##if null just create empty mat and return


    ## check that the cache has no invalid nodes - these should be removed
    myNA <- which(is.na(score.cache$mlik))
    if (length(myNA) > 0) {
        if (verbose)
            cat("### NOTE: the score.cache has missing values in mlik - assigning to -.Machine$double.xmax ###\n")
        score.cache$mlik[myNA] <- -.Machine$double.xmax  ## the most negative number possible, i.e. -infinity
    }

    use.start.dag <- FALSE
    ## user supplied start network must also have data passed
    if (!is.null(start.dag)) {
        if (is.null(data.df)) {
            stop("must supply data.df when using explicit start.dag")
        }
        if (num.searches != 1) {
            stop("num.searches must equal 1 when using explicit start dag - otherwise an identical search is repeated\n => the final DAG identified depends only on the initial search location\n")
        }
        check.valid.dag(dag=start.dag, data.df=data.df, is.ban.matrix=FALSE, NULL)
        use.start.dag <- TRUE
    }

    cache.defn <- as.integer(score.cache[["node.defn"]])  ## into one long vector filled by col
    children <- as.integer(score.cache[["children"]])
    nodescores <- as.double(score.cache[["mlik"]])
    numVars <- as.integer(dim(score.cache[["node.defn"]])[2])  ## number of variables in a DAG
    numRows <- as.integer(dim(score.cache[["node.defn"]])[1])  ## total number of different parent combinations

    numparents.per.node <- as.integer(table(children))  ## number of parent combinations per variable

    ## all inside C get a random DAG from which to start - N.B. last argument is verbose
    res <- .Call("searchhill", children, cache.defn, nodescores, numVars, numRows, numparents.per.node, as.integer(seed), verbose, as.integer(timing.on), as.integer(use.start.dag), as.integer(start.dag),
        as.integer(num.searches), as.integer(dag.retained), PACKAGE="abn"  ## uncomment to load as package not shlib
)

    ## results format is res[[1]] is a vector of network scores in order: init score1, final score1, init score2, final score2,....etc with res[[2]]=init network 1, res[[3]] final network 1,
    ## res[[4]] init network2, res[[5]] init network2,... etc
    for (i in 2:length(res)) {
        colnames(res[[i]]) <- colnames(score.cache[[2]])
        rownames(res[[i]]) <- colnames(score.cache[[2]])
    }
    ## now re-organise res into something easier to analyse - a list of three lists list 1 - vector of scores and list of initial matrices, list 2 -vector of score and list of final matrices
    scores <- res[[1]]
    init.indexes <- seq(1, 2 * num.searches, by=2)
    fin.indexes <- seq(2, 2 * num.searches, by=2)

    init.scores <- scores[init.indexes]  ## scores from initial networks
    fin.scores <- scores[fin.indexes]
    init.mat <- res[init.indexes + 1]  ##offset for score vector
    fin.mat <- res[fin.indexes + 1]  ##offset for score vector
    rm(res)

    con.dag <- fin.mat[[1]]
    if (num.searches > 1)
        {
            for (i in 2:num.searches) {
                con.dag <- con.dag + fin.mat[[i]]
            }
        }  ## get total arc freq
    ## now make binary according to threshold passed
    con.dag.binary <- ifelse(con.dag >= ceiling(num.searches * support.threshold), 1, 0)
    ######################################################### create graph object part
        out <- list(init.score=init.scores,
                    final.score=fin.scores,
                    init.dag=init.mat,
                    final.dag=fin.mat,
                    # consensus=(con.dag.binary),
                    dag=(con.dag.binary),
                    support.threshold=support.threshold,
                    score.cache=score.cache)
        class(out) <- c("abnHillClimber", "abnLearned")
        return(out)


}
