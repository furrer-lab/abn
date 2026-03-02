# Code extracted from buildScoreCache.mle() to get only the parameters
# needed to build the cache matrix

setupScoreCache.mle <-
  function(data.df = NULL,
           data.dists = NULL,
           max.parents = NULL,
           adj.vars = NULL,
           cor.vars = NULL,
           dag.banned = NULL,
           dag.retained = NULL,
           which.nodes = NULL,
           centre = TRUE,
           defn.res = NULL,
           dry.run = FALSE,
           verbose = FALSE,
           debugging = FALSE,
           force.method = NULL,
           group.var = NULL,
           grouped.vars = NULL,
           group.ids = NULL,
           control = build.control(method = "mle")) {

    set.seed(control[["seed"]])

    ## which.nodes
    if (!is.null(which.nodes)) {
        data.df <- data.df[, which.nodes]
        data.dists <- data.dists[which.nodes]
    }

    ## number of variable:
    nvars <- length(data.dists)
    nobs <- dim(data.df)[1]

    # formating factor
    data.df.lvl <- data.df

    ## standardize gaussian variables to zero mean and sd=1 have at least one gaussian variable
    if (centre && !is.null(data.dists == "gaussian")) {
        for (i in names(data.dists)[(data.dists == "gaussian")]) {
            data.df[, i] <- (data.df[, i] - mean(data.df[, i]))/sd(data.df[, i])
        }
    }

    for (i in 1:nvars) {
        if (data.dists[[i]] == "binomial" & !inherits(data.df[, i], "numeric")) {
            data.df[, i] <- as.numeric(factor(data.df[, i])) - 1
        }
        if (data.dists[[i]] == "multinomial") {
            data.df[, i] <- factor(data.df[, i])
        }
    }

    # unpacking the multinomial variables in the data.df
    data.df.multi <- NULL
    for (i in 1:nvars) {
      if (data.dists[[i]] %in% c("binomial", "poisson", "gaussian")) {
        data.df.multi <- as.data.frame(cbind(data.df.multi, data.df[, i]))
        colnames(data.df.multi)[length(colnames(data.df.multi))] <- colnames(data.df)[i]
      } else {
        tmp <- model.matrix(~-1 + factor(data.df.lvl[, i]))
        colnames(tmp) <- paste0(names(data.df.lvl)[i], levels(factor(data.df.lvl[, i])))
        data.df.multi <- as.data.frame(cbind(data.df.multi, tmp))
      }
    }

    return(list(nvars = nvars, nobs = nobs, data.df.lvl = data.df.lvl, data.df = data.df, data.df.multi = data.df.multi))

  }

# Modified version with preallocation of data.df.multi and colnames
setupScoreCache.mle_mod <-
  function(data.df = NULL,
           data.dists = NULL,
           max.parents = NULL,
           adj.vars = NULL,
           cor.vars = NULL,
           dag.banned = NULL,
           dag.retained = NULL,
           which.nodes = NULL,
           centre = TRUE,
           defn.res = NULL,
           dry.run = FALSE,
           verbose = FALSE,
           debugging = FALSE,
           force.method = NULL,
           group.var = NULL,
           grouped.vars = NULL,
           group.ids = NULL,
           control = build.control(method = "mle")) {

    set.seed(control[["seed"]])

    ## which.nodes
    if (!is.null(which.nodes)) {
      data.df <- data.df[, which.nodes]
      data.dists <- data.dists[which.nodes]
    }

    ## number of variable:
    nvars <- length(data.dists)
    nobs <- dim(data.df)[1]

    # formating factor
    data.df.lvl <- data.df

    ## standardize gaussian variables to zero mean and sd=1 have at least one gaussian variable
    if (centre && !is.null(data.dists == "gaussian")) {
      for (i in names(data.dists)[(data.dists == "gaussian")]) {
        data.df[, i] <- (data.df[, i] - mean(data.df[, i]))/sd(data.df[, i])
      }
    }

    for (i in 1:nvars) {
      if (data.dists[[i]] == "binomial" & !inherits(data.df[, i], "numeric")) {
        data.df[, i] <- as.numeric(factor(data.df[, i])) - 1
      }
      if (data.dists[[i]] == "multinomial") {
        data.df[, i] <- factor(data.df[, i])
      }
    }

    # unpacking the multinomial variables in the data.df

    #data.df.multi <- NULL
    #data.df.multi.list <- vector('list', nvars)
    #colnames.list <- vector('list', nvars)
    #for (i in 1:nvars) {
    #  if (data.dists[[i]] %in% c("binomial", "poisson", "gaussian")) {
    #    data.df.multi.list[[i]] <- data.df[, i]
    #    colnames.list[[i]] <- colnames(data.df)[i]
    #    #data.df.multi <- as.data.frame(cbind(data.df.multi, data.df[, i]))
    #    #colnames(data.df.multi)[length(colnames(data.df.multi))] <- colnames(data.df)[i]
    #  } else {
    #    tmp <- model.matrix(~-1 + factor(data.df.lvl[, i]))
    #    colnames.list[[i]] <- paste0(names(data.df.lvl)[i], levels(factor(data.df.lvl[, i])))
    #    data.df.multi.list[[i]] <- tmp
    #    #colnames(tmp) <- paste0(names(data.df.lvl)[i], levels(factor(data.df.lvl[, i])))
    #    #data.df.multi <- as.data.frame(cbind(data.df.multi, tmp))
    #  }
    #}
    #data.df.multi <- as.data.frame(do.call('cbind', data.df.multi.list))
    #colnames(data.df.multi) <- unlist(colnames.list)

    n_rows <- nrow(data.df)
    data.df.multi.list <- matrix(NA_real_, nrow = n_rows, ncol = nvars)
    n_filled <- 0
    chunk_size <- 10    # grow by this many columns
    colnames.list <- vector('list', nvars)
    for (i in 1:nvars) {
      if (data.dists[[i]] %in% c("binomial", "poisson", "gaussian")) {
        new_col <- data.df[, i]
        n_filled <- n_filled + 1
        if (n_filled > ncol(data.df.multi.list)) {
          new_ncol <- ncol(data.df.multi.list) + chunk_size
          new_mat <- matrix(NA_real_, nrow = n_rows, ncol = new_ncol)
          new_mat[, 1:n_filled-1] <- data.df.multi.list[, 1:n_filled-1]
          data.df.multi.list <- new_mat
        }
        data.df.multi.list[, n_filled] <- new_col
        colnames.list[[i]] <- colnames(data.df)[i]
      } else {
        tmp <- model.matrix(~-1 + factor(data.df.lvl[, i]))
        colnames.list[[i]] <- paste0(names(data.df.lvl)[i], levels(factor(data.df.lvl[, i])))
        new_col <- tmp
        cols_needed <- n_filled + ncol(tmp)
        if (cols_needed > ncol(data.df.multi.list)) {
          new_ncol <- ncol(data.df.multi.list) + chunk_size
          new_mat <- matrix(NA_real_, nrow = n_rows, ncol = new_ncol)
          new_mat[, 1:n_filled] <- data.df.multi.list[, 1:n_filled]
          data.df.multi.list <- new_mat
        }
        data.df.multi.list[, (n_filled + 1):cols_needed] <- new_col
        n_filled <- cols_needed
      }
    }
    data.df.multi <- as.data.frame(data.df.multi.list[, 1:n_filled, drop = FALSE])
    colnames(data.df.multi) <- unlist(colnames.list)

    return(list(nvars = nvars, nobs = nobs, data.df.lvl = data.df.lvl, data.df = data.df, data.df.multi = data.df.multi))

  }
