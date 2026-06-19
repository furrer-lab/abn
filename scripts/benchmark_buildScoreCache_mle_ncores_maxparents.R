# Performance comparison for the original buildScoreCache() function when multi cores are used.
# This scripts also include the identity check between the original code and the help function
# Run: Rscript scripts/benchmark_buildScoreCache_mle_ncores_maxparents.R
# Note: Failing case removed from the computation: child 14, poisson;  parents:     "b1" "g2" "g3" "p1"; function: `/src/irls_poisson_fast.cpp`

library(abn)
library(parallel)
library(doParallel)
load("data/n_1000_p_16_nodedists_Balanced_s_02_graph1.Rdata")
mydf <- n_1000_p_16_nodedists_Balanced_s_02_graph1$data[,-c(14)]
mydists <- n_1000_p_16_nodedists_Balanced_s_02_graph1$dists[-c(14)]

source("scripts/setup_score_cache.R")
# The code lines in the setupScoreCache function are extracted from buildScoreCache().
# setupScoreCache only runs all checks and data cleaning without calling the buildScoreCache.mle() function
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


get_params <- function(mydf, mydists, mp, ncores) {
	params <- setupScoreCache(data.df = mydf,
                             data.dists = mydists,
                             method = "mle",
                             max.parents = mp,
                             group.var = NULL,
                             adj.vars = NULL,
                             cor.vars = NULL,
                             dag.banned = NULL,
                             dag.retained = NULL,
                             which.nodes = NULL,
                             defn.res = NULL,
                             centre = TRUE,
                             dry.run = FALSE,
                             control = build.control(method = "mle",  ncores=ncores),
                             verbose = FALSE,
                             debugging = FALSE)
	params
}


# When benchmarking parallel executions (ncores > 1) run bench::mark with memory = False
# Memory profiling would be otherwise unreliable.
# Runtime benchmarking is possible, but not the comparison of allocation behaviour.
grid <- expand.grid(ncores = c(1, 2, 3, 4), maxp = c(3, 4, 5, 6))
results <- bench::press(
  ncores = grid$ncores,
  maxp = grid$maxp,
  {
    params <- get_params(mydf, mydists, maxp, ncores)
    gc()
    bench::mark(
		bsc.mle.new = newBuildScoreCache.mle(data.df = params$data.df,
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
                                 control = params$control),

		bsc.mle.orig = buildScoreCache.mle(data.df = params$data.df,
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
                                        control = params$control),
		iterations = 1,
                memory = FALSE)
  }
)

saveRDS(results, "bench_press_results.rds")
flat <- results[, sapply(results, function(x) !is.list(x))]
utils::write.csv(flat, "bench_press_results_flat.csv", row.names = FALSE)

