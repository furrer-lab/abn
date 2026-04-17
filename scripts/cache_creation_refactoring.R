# Refactoring cache creation with the idea to include it into the forLoopContent()
computeCache_inForLoop <- function(adj.vars,
                                   nvars,
                                   data.df,
                                   data.df.lvl,
                                   max.parents,
                                   data.dists){
  tmp.list <- vector('list', max.parents)
  tmp_nrow <- 1
  for (i in 1:(max.parents)) {
    tmp.list[[i]] <- t(combn(x = (nvars - 1), m = i, FUN = fun.return, n = nvars, simplify = TRUE))
    tmp_nrow <- tmp_nrow + nrow(tmp.list[[i]])
  }
  cache_nrow <- nvars * (tmp_nrow)
  n_filled <- 0
  node.defn_all.list <- matrix(data = as.integer(0), nrow = cache_nrow, ncol = nvars)
  children.list <- vector('list', cache_nrow)
  for (child in 1:nvars){
    res <- forLoopNode_withCache_mod(child = child,
                                 tmp.list = tmp.list,
                                 tmp_nrow = tmp_nrow,
                                 max.parents = max.parents,
                                 nvars = nvars)
    node.defn_all.list[(n_filled+1):(n_filled+tmp_nrow),] <- res$node.defn
    children.list[(n_filled+1):(n_filled+tmp_nrow)] <- res$children
    n_filled <- n_filled + tmp_nrow
  }
  node.defn <- node.defn_all.list
  children <- unlist(children.list)
  colnames(node.defn) <- colnames(data.df)
  mode(node.defn) <- "integer"
  children <- as.integer(children)
  ## HERE CODE for DAG RETAIN/BANNED
  mycache_all <- list(children = (children), node.defn = (node.defn))

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
                                  tmp.list,
                                  tmp_nrow,
                                  max.parents,
                                  nvars){
  n_filled <- 0
  node.defn_temp.list <- matrix(data = as.integer(0), nrow = tmp_nrow, ncol = nvars)
  children_temp.list <- vector("list", tmp_nrow)
  n_filled <- n_filled + 1
  node.defn_temp.list[n_filled, ] <- matrix(data = as.integer(0), nrow = 1L, ncol = nvars)
  children_temp.list[[n_filled]] <- child
  for (i in 1:(max.parents)) {
    tmp <- t(apply(X = tmp.list[[i]], MARGIN = 1, FUN = function(x) append(x = x, values = 0, after = child - 1)))
    node.defn_temp.list[(n_filled+1):(n_filled+nrow(tmp)),] <- tmp
    children_temp.list[(n_filled+1):(n_filled+nrow(tmp))] <- t(rep(child, length(tmp[, 1])))
    n_filled <- n_filled + nrow(tmp)
  }
  mycache_temp <- list(children = (children_temp.list), node.defn = (node.defn_temp.list))
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


# In buildScoreCache.mle(): replace the growing node.defn objects within the loop with a matrix preallocation
# Optimize the tmp matrix calculation
# Include pre-calculation of the X and Y matrices to be passed to the forLoopContent function (only for performance comparison)
computeCache_precomputeX <- function(adj.vars,
                                  nvars,
                                  data.df,
                                  data.df.multi,
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

  cache_X <- new.env(hash = TRUE, parent = emptyenv())
  X_list <- vector("list", nrow(mycache$node.defn.multi))
  for (k in seq_len(nrow(mycache$node.defn.multi))) {
    key <- bits_to_key(mycache$node.defn.multi[k, ])
    if (exists(key, envir = cache_X, inherits = FALSE)) {
      X_list[[k]] <- get(key, envir = cache_X, inherits = FALSE)
    } else {
      if ("multinomial" %in% data.dists[as.logical(mycache$node.defn[k, ])])
        multinomial = TRUE
      else multinomial = FALSE
      obj <- compute_X(mycache$node.defn.multi[k, ], adj.vars, multinomial, data.df.multi)
      assign(key, obj, envir = cache_X)
      X_list[[k]] <- obj
    }
  }
  mycache$X <- X_list

  Y_list <- vector("list", nvars)
  for (k in 1:nvars){
   Y_list[[k]] <- data.df[, k, drop = FALSE]
  }
  mycache$Y <- Y_list

  return(mycache)
}

bits_to_key <- function(x) {
  as.character(sum(x * 2^seq.int(from = length(x) - 1, to = 0)))
}

compute_X <- function(row, adj.vars, multinomial, data.df.multi) {
  if (is.null(adj.vars)) {
    if (multinomial) {
      #X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      X <- data.df.multi[, as.logical(row), drop = FALSE]
    } else {
      #X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      X <- cbind(1, data.df.multi[, as.logical(row), drop = FALSE])
    }
  } else {
    if (multinomial) {
      #X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      X <- data.df.multi[, as.logical(row), drop = FALSE]
    } else {
      #X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
      X <- cbind(1, data.df.multi[, as.logical(row), drop = FALSE])
    }
  }
   X
}
