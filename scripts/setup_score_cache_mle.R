# Code extracted from buildScoreCache.mle() to get the parameters
# used by the main code to be eventually refactored

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

    # adjustment: storing of df
    if (!is.null(adj.vars)) {
        data.df.adj <- data.df
        data.df <- data.df[, -adj.vars]
        nvars <- nvars - length(adj.vars)
    }

    return(list(nvars = nvars, nobs = nobs, data.df.lvl = data.df.lvl, data.df = data.df))

  }
