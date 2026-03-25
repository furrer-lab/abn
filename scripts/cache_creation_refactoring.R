# Refactoring cache creation with the idea to include it into the forLoopContent()
computeCache_inForLoop <- function(adj.vars,
                                   nvars,
                               data.df,
                               data.df.lvl,
                               max.parents,
                               data.dists){
    node.defn_all.list <- vector('list', nvars)
    children_all.list <- vector('list', nvars)
    for (child in 1:nvars){
      res <- forLoopNode_withCache(child = child,
                               max.parents = max.parents,
                               nvars = nvars)
      node.defn_all.list[[child]] <- res$node.defn
      children_all.list[[child]] <- res$children
    }
    node.defn_all <- do.call('rbind', node.defn_all.list)
    children_all <- do.call('cbind', children_all.list)
    colnames(node.defn_all) <- colnames(data.df)
    mode(node.defn_all) <- "integer"
    children_all <- as.integer(children_all)
    ## HERE CODE for DAG RETAIN/BANNED
    mycache_all <- list(children = (children_all), node.defn = (node.defn_all))

    # ASSUMING adj=NULL and ignore adjustment
    mycache_all$node.defn.adj <- mycache_all$node.defn

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
      mycache_all$node.defn.multi <- mycache_all$node.defn.adj[, rep(1:nvars, repetition.multi)]
      data.df <- data.df.adj[, colnames(mycache_all$node.defn.adj)]
    } else {
      mycache_all$node.defn.multi <- mycache_all$node.defn[, rep(1:nvars, repetition.multi)]

    }

    mycache_all

  }


forLoopNode_withCache <- function(child,
                              max.parents,
                              nvars){

    node.defn_temp.list <- vector("list", max.parents+1)
    children_temp.list <- vector("list", max.parents+1)
    node.defn_temp.list[[1]] <- matrix(data = as.integer(0), nrow = 1L, ncol = nvars)
    children_temp.list[[1]] <- child
    j <- child
    for (i in 1:(max.parents)) {
      tmp <- t(combn(x = (nvars - 1), m = i, FUN = fun.return, n = nvars, simplify = TRUE))
      tmp <- t(apply(X = tmp, MARGIN = 1, FUN = function(x) append(x = x, values = 0, after = j - 1)))

      node.defn_temp.list[[i+1]] <- tmp
      children_temp.list[[i+1]] <- t(rep(j, length(tmp[, 1])))
    }
    node.defn_temp <- do.call('rbind', node.defn_temp.list)
    children_temp <- do.call('cbind', children_temp.list)
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
# Note: adjustment is ignored (mycache$node.defn.adj <- mycache$node.defn)
computeCache_orig <- function(adj.vars,
                              nvars,
                              data.df,
                              data.df.lvl,
                              max.parents,
                              data.dists){
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
  colnames(node.defn) <- colnames(data.df)
  ## Coerce numeric matrix into integer matrix !!!
  #node.defn <- apply(node.defn, c(1, 2), function(x) {
  #  (as.integer(x))
  #})
  mode(node.defn) <- "integer"
  children <- as.integer(children)
  mycache <- list(children = (children), node.defn = (node.defn))

  # ASSUMING adj=NULL and ignore adjustment
  mycache$node.defn.adj <- mycache$node.defn

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

  return(mycache)
}

# In buildScoreCache.mle(): replace the growing objects within the loop with preallocation
# and the do.Call outside the loop
computeCache_doCall <- function(adj.vars,
                                nvars,
                              data.df,
                              data.df.lvl,
                              max.parents,
                              data.dists){
  node.defn_all.list <- vector('list', nvars)
  children_all.list <- vector('list', nvars)
  for (j in 1:nvars) {
    node.defn_temp.list <- vector('list', max.parents+1)
    children_temp.list <- vector('list', max.parents+1)
    node.defn_temp.list[[1]] <- matrix(data = as.integer(0), nrow = 1L, ncol = nvars)
    children_temp.list[[1]] <- j
    for (i in 1:(max.parents)) {
      tmp <- t(combn(x = (nvars - 1), m = i, FUN = fun.return, n = nvars, simplify = TRUE))
      tmp <- t(apply(X = tmp, MARGIN = 1, FUN = function(x) append(x = x, values = 0, after = j - 1)))
      node.defn_temp.list[[i+1]] <- tmp
      children_temp.list[[i+1]] <- t(rep(j, length(tmp[, 1])))
    }
    node.defn_all.list[[j]] <- do.call('rbind', node.defn_temp.list)
    children_all.list[[j]] <- do.call('cbind', children_temp.list)
  }
  node.defn <- do.call('rbind', node.defn_all.list)
  children <- do.call('cbind', children_all.list)
  colnames(node.defn) <- colnames(data.df)
  mode(node.defn) <- "integer"
  children <- as.integer(children)
  mycache <- list(children = (children), node.defn = (node.defn))

  # ASSUMING adj=NULL and ignore adjustment
  mycache$node.defn.adj <- mycache$node.defn

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

  return(mycache)
}

# In buildScoreCache.mle(): replace the growing node.defn objects within the loop with a matrix preallocation
# Optimize the tmp matrix calculation
computeCache_prealloc <- function(adj.vars,
                                nvars,
                                data.df,
                                data.df.lvl,
                                max.parents,
                                data.dists){
  tmp.list <- vector('list', max.parents)
  cache_nrow <- 1
  for (i in 1:(max.parents)) {
    tmp.list[[i]] <- t(combn(x = (nvars - 1), m = i, FUN = fun.return, n = nvars, simplify = TRUE))
    cache_nrow <- cache_nrow + nrow(tmp.list[[i]])
  }
  cache_nrow <- nvars * (cache_nrow)
  n_filled <- 0
  node.defn_all.list <- matrix(data = as.integer(0), nrow = cache_nrow, ncol = nvars)
  children.list <- vector('list', cache_nrow)
  for (j in 1:nvars) {
    n_filled <- n_filled + 1
    node.defn_all.list[n_filled, ] <- matrix(data = as.integer(0), nrow = 1L, ncol = nvars)
    children.list[[n_filled]] <- j
    for (i in 1:(max.parents)) {
      tmp <- t(apply(X = tmp.list[[i]], MARGIN = 1, FUN = function(x) append(x = x, values = 0, after = j - 1)))
      node.defn_all.list[(n_filled+1):(n_filled+nrow(tmp)),] <- tmp
      children.list[(n_filled+1):(n_filled+nrow(tmp))] <- t(rep(j, length(tmp[, 1])))
      n_filled <- n_filled + nrow(tmp)
    }
  }
  node.defn <- node.defn_all.list
  children <- unlist(children.list)
  colnames(node.defn) <- colnames(data.df)
  mode(node.defn) <- "integer"
  children <- as.integer(children)
  mycache <- list(children = (children), node.defn = (node.defn))

  # ASSUMING adj=NULL and ignore adjustment
  mycache$node.defn.adj <- mycache$node.defn

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

  return(mycache)
}
