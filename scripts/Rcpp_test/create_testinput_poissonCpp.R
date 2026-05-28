source("scripts/setup_score_cache_mle.R")
source("scripts/cache_creation_refactoring.R")
params <- readRDS(file = "tests/testthat/testdata/input_score_cache_mle.RData")
params$max.parents <- 4 #unique(params$max.parents)

mleparams <- setupScoreCache.mle(data.df = params$data.df,
                                 data.dists = params$data.dists,
                                 adj.vars = params$adj.vars,
                                 cor.vars = params$cor.vars,
                                 dag.banned = params$dag.banned,
                                 dag.retained = params$dag.retained,
                                 max.parents = params$max.parents,
                                 which.nodes = params$which.nodes,
                                 defn.res = params$defn.res,
                                 dry.run = params$dry.run,
                                 verbose = params$verbose,
                                 debugging = params$debugging,
                                 centre = params$centre,
                                 force.method = params$force.method,
                                 group.var = params$group.var,
                                 grouped.vars = params$grouped.vars,
                                 group.ids = params$group.ids,
                                 control = params$control)

cache_orig = computeCache_orig(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.lv, params$max.parents, params$data.dists)


rows <- length(cache_orig[["children"]])
row.num <- 26059
child <- cache_orig[["children"]][row.num]
distribution <- params$data.dists[child]
Y <- data.matrix(mleparams$data.df[, child])
Y1 <- if (is.factor(Y)) numeric(Y) else Y
# Note: adjustment is ignored (mycache$node.defn.adj <- mycache$node.defn)
if (is.null(params$adj.vars)) {
  if ("multinomial" %in% params$data.dists[as.logical(cache_orig[["node.defn"]][row.num, ])]) {
    X <- data.matrix(cbind(mleparams$data.df.multi[, as.logical(cache_orig[["node.defn.multi"]][row.num, ])]))
  } else {
    X <- data.matrix(cbind(rep(1, length(mleparams$data.df[, 1])), mleparams$data.df.multi[, as.logical(cache_orig[["node.defn.multi"]][row.num, ])]))
  }
} else {
  if ("multinomial" %in% params$data.dists[as.logical(mycache$node.defn.adj[row.num, ])]) {
    X <- data.matrix(cbind(data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
  } else {
    X <- data.matrix(cbind(rep(1, length(data.df[, 1])), data.df.multi[, as.logical(mycache[["node.defn.multi"]][row.num, ])]))
  }
}


saveRDS(list(A = X, b1 = Y1, b = Y, maxit = control[["max.iters"]], tol = control[["tol"]]), file = "scripts/input_poisson_cpp.RData")


# For testing the fit results uncomment lines below:
# fit <- irls_poisson_cpp_fast(A = X, b = Y, maxit = params$control[["max.iters"]], tol = params$control[["tol"]])
# input_poisson <- readRDS(file = "scripts/input_poisson_cpp.RData")
# fit2 <- irls_poisson_cpp_fast(input_poisson$A, input_poisson$b, input_poisson$maxit, input_poisson$tol)
# identical(fit, fit2)

