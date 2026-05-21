library(parallel)
library(doParallel)
load("data/n_1000_p_16_nodedists_Balanced_s_02_graph1.Rdata")
mydf <- n_1000_p_16_nodedists_Balanced_s_02_graph1$data
mydists <- n_1000_p_16_nodedists_Balanced_s_02_graph1$dists

source("scripts/setup_score_cache.R")
params <- readRDS(file = "tests/testthat/testdata/input_score_cache_mle.RData") # with control = NULL
params_c1 <- setupScoreCache(data.df = mydf,
                          data.dists = mydists,
                          method = "mle",
                          max.parents = 3,
                          group.var = NULL,
                          adj.vars = NULL,
                          cor.vars = NULL,
                          dag.banned = NULL,
                          dag.retained = NULL,
                          which.nodes = NULL,
                          defn.res = NULL,
                          centre = TRUE,
                          dry.run = FALSE,
                          control = build.control(method = "mle",  ncores=1),
                          verbose = FALSE,
                          debugging = FALSE)

params_c2 <- setupScoreCache(data.df = mydf,
                             data.dists = mydists,
                             method = "mle",
                             max.parents = 3,
                             group.var = NULL,
                             adj.vars = NULL,
                             cor.vars = NULL,
                             dag.banned = NULL,
                             dag.retained = NULL,
                             which.nodes = NULL,
                             defn.res = NULL,
                             centre = TRUE,
                             dry.run = FALSE,
                             control = build.control(method = "mle",  ncores=2),
                             verbose = FALSE,
                             debugging = FALSE)

identical(params, params_c1)

# Compare results using the original code buildScoreCache():
lb <- bench::mark(
bsc.orig = buildScoreCache(data.df = params_c1$data.df,
                           data.dists = params_c1$data.dists,
                           method = 'mle',
                           adj.vars = params_c1$adj.vars,
                           cor.vars = params_c1$cor.vars,
                           dag.banned = params_c1$dag.banned,
                           dag.retained = params_c1$dag.retained,
                           max.parents = params_c1$max.parents,
                           which.nodes = params_c1$which.nodes,
                           defn.res = params_c1$defn.res,
                           dry.run = params_c1$dry.run,
                           verbose = params_c1$verbose,
                           debugging = params_c1$debugging,
                           centre = params_c1$centre,
                           force.method = params_c1$force.method,
                           group.var = params_c1$group.var,
                           grouped.vars = params_c1$grouped.vars,
                           group.ids = params_c1$group.ids,
                           control = params_c1$control),
bsc.orig_c2 = buildScoreCache(data.df = params_c2$data.df,
                              data.dists = params_c2$data.dists,
                              method = 'mle',
                              adj.vars = params_c2$adj.vars,
                              cor.vars = params_c2$cor.vars,
                              dag.banned = params_c2$dag.banned,
                              dag.retained = params_c2$dag.retained,
                              max.parents = params_c2$max.parents,
                              which.nodes = params_c2$which.nodes,
                              defn.res = params_c2$defn.res,
                              dry.run = params_c2$dry.run,
                              verbose = params_c2$verbose,
                              debugging = params_c2$debugging,
                              centre = params_c2$centre,
                              force.method = params_c2$force.method,
                              group.var = params_c2$group.var,
                              grouped.vars = params_c2$grouped.vars,
                              group.ids = params_c2$group.ids,
                              control = params_c2$control)
)

# Compare results using the buildScoreCache.mle():
# - code extracted from the buildScoreCache() function
# - it only runs all checks and data cleaning before calling the buildScoreCache.mle()
lb <- bench::mark(
  bsc.mle.orig_c1 = buildScoreCache.mle(data.df = params_c1$data.df,
                                     data.dists = params_c1$data.dists,
                                     adj.vars = params_c1$adj.vars,
                                     cor.vars = params_c1$cor.vars,
                                     dag.banned = params_c1$dag.banned,
                                     dag.retained = params_c1$dag.retained,
                                     max.parents = params_c1$max.parents,
                                     which.nodes = params_c1$which.nodes,
                                     defn.res = params_c1$defn.res,
                                     dry.run = params_c1$dry.run,
                                     verbose = params_c1$verbose,
                                     debugging = params_c1$debugging,
                                     centre = params_c1$centre,
                                     force.method = params_c1$force.method,
                                     group.var = params_c1$group.var,
                                     grouped.vars = params_c1$grouped.vars,
                                     group.ids = params_c1$group.ids,
                                     control = params_c1$control),
  bsc.mle.orig_c2 = buildScoreCache.mle(data.df = params_c2$data.df,
                                        data.dists = params_c2$data.dists,
                                        adj.vars = params_c2$adj.vars,
                                        cor.vars = params_c2$cor.vars,
                                        dag.banned = params_c2$dag.banned,
                                        dag.retained = params_c2$dag.retained,
                                        max.parents = params_c2$max.parents,
                                        which.nodes = params_c2$which.nodes,
                                        defn.res = params_c2$defn.res,
                                        dry.run = params_c2$dry.run,
                                        verbose = params_c2$verbose,
                                        debugging = params_c2$debugging,
                                        centre = params_c2$centre,
                                        force.method = params_c2$force.method,
                                        group.var = params_c2$group.var,
                                        grouped.vars = params_c2$grouped.vars,
                                        group.ids = params_c2$group.ids,
                                        control = params_c2$control)
)


# Same results, only different output class
class(bsc.mle.orig_c1) <- c("abnCache")
identical(bsc.orig, bsc.mle.orig_c1)
