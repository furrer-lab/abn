params <- readRDS(file = "tests/testthat/testdata/input_score_cache_mle.RData")
params$max.parents <- unique(params$max.parents)

### Compare changes to the buildScoreCache.mle with the original code

# Tests
testthat::test_file("tests/testthat/test-build_score_cache_mle_newcache.R")

# Microbenchmarking
(lb <- bench::mark(
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
                                 control = params$control)
		   ))
# Results of the benchmarking:
# A tibble: 2 × 13
#   expression    min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory     time       gc
#   <bch:expr>  <bch> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>     <list>     <list>
# 1 bsc.mle.or… 22.2s  22.2s    0.0449    1.88GB    1.35      1    30      22.2s <named list> <Rprofmem> <bench_tm> <tibble>
# 2 bsc.mle.new 20.9s  20.9s    0.0478    1.84GB    0.717     1    15      20.9s <named list> <Rprofmem> <bench_tm> <tibble>
#

# Profiling
library(profvis)
params <- readRDS(file = "tests/testthat/testdata/input_score_cache_mle.RData")
# Loading the function makes the line of codes of all functions called as available in the profvis output 
source("R/build_score_cache_mle_newcache.R")
profvis({
	newBuildScoreCache.mle(data.df = params$data.df,
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
}, prof_output = "scripts/profiler_reports/profvis_buildScoreCache_mle_new.Rprofvis")
#profvis(prof_input = "scripts/profiler_reports/profvis_buildScoreCache_mle_new.Rprofvis")

# Loading the function makes the line of codes of all functions called as available in the profvis output 
source("R/build_score_cache_mle.R")
profvis({
        buildScoreCache.mle(data.df = params$data.df,
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
}, prof_output = "scripts/profiler_reports/profvis_buildScoreCache_mle_orig.Rprofvis")
#profvis(prof_input = "scripts/profiler_reports/profvis_buildScoreCache_mle_orig.Rprofvis")
