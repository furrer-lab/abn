# Refactoring cache creation with the idea to include it into the forLoopContent()
computeCache_maria <- function(nvars,
                               data.df,
                               max.parents){
    for (child in 1:nvars){
      res <- forLoopNode_maria(child = child,
                               max.parents = max.parents,
                               nvars = nvars)
      if (child == 1){
        node.defn_all <- res$node.defn
        children_all <- res$children
      } else {
        node.defn_all <- rbind(node.defn_all, res$node.defn)
        children_all <- cbind(children_all, res$children)
      }


    }
    colnames(node.defn_all) <- colnames(data.df)
    mode(node.defn_all) <- "integer"
    children_all <- as.integer(children_all)
    ## HERE CODE for DAG RETAIN/BANNED
    mycache_all <- list(children = (children_all), node.defn = (node.defn_all))
    mycache_all

  }


forLoopNode_maria <- function(child,
                              max.parents,
                              nvars){

    node.defn_temp <- matrix(data = as.integer(0), nrow = 1L, ncol = nvars)
    children_temp <- child
    j <- child
    {
      for (i in 1:(max.parents)) {
        tmp <- t(combn(x = (nvars - 1), m = i, FUN = fun.return, n = nvars, simplify = TRUE))
        tmp <- t(apply(X = tmp, MARGIN = 1, FUN = function(x) append(x = x, values = 0, after = j - 1)))

        node.defn_temp <- rbind(node.defn_temp, tmp)
        children_temp <- cbind(children_temp, t(rep(j, length(tmp[, 1]))))
      }
    }
    mycache_temp <- list(children = (children_temp), node.defn = (node.defn_temp))
    mycache_temp
  }

fun.return <- function(x, n) {
  # x: vector or single integer
  # n: number of variables/columns in data.df
  # returns: vector of size n-1 with 1 for all x and zeros otherwise
  v <- rep(0, n - 1)
  v[x] <- 1
  return(v)
}


# Extracted from buildScoreCache.mle()
computeCache_orig <- function(nvars,
                              data.df,
                              max.parents){

    node.defn <- matrix(data = as.integer(0), nrow = 1L, ncol = nvars)
    children <- 1
    for (j in 1:nvars) {
      if (j != 1) {
        node.defn <- rbind(node.defn, matrix(data = as.integer(0),
                                             nrow = 1L, ncol = nvars))
        children <- cbind(children, j)
      }
      for (i in 1:(max.parents)) {
        tmp <- t(combn(x = (nvars - 1), m = i, FUN = fun.return, n = nvars, simplify = TRUE))
        tmp <- t(apply(X = tmp, MARGIN = 1, FUN = function(x) append(x = x, values = 0, after = j - 1)))
        node.defn <- rbind(node.defn, tmp)
        # children position
        children <- cbind(children, t(rep(j, length(tmp[, 1]))))
      }
    }
    colnames(node.defn) <- colnames(data.df)
    mode(node.defn) <- "integer"
    children <- as.integer(children)
    mycache <- list(children = (children), node.defn = (node.defn))
    mycache
  }

