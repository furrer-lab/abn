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
                                 control = params$control),
    iterations = 5)
		   )
# Results of the benchmarking:
# A tibble: 2 × 13
#   expression    min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory     time       gc
#   <bch:expr>  <bch> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>     <list>     <list>
# 1 bsc.mle.or… 22.2s  22.2s    0.0449    1.88GB    1.35      1    30      22.2s <named list> <Rprofmem> <bench_tm> <tibble>
# 2 bsc.mle.new 20.9s  20.9s    0.0478    1.84GB    0.717     1    15      20.9s <named list> <Rprofmem> <bench_tm> <tibble>
#
# Comparison after the following changes:
# - implemented the cache creation with matrix preallocation
# - precomupute X and Y (inputs to the Rccp functions)
# - forLoopContent with precomputed X and Y
# - compute  Ymulti in the "multinomial" case by exploiting the data.df.multi
# A tibble: 2 × 13
#  expression        min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result            memory                   time           gc
#  <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>            <list>                   <list>         <list>
#  1 bsc.mle.orig    26.8s    26.9s    0.0372    1.88GB   0.372      5    50      2.24m <named list [18]> <Rprofmem [241,086 × 3]> <bench_tm [5]> <tibble [5 × 3]>
#  2 bsc.mle.new     19.1s    19.2s    0.0521  903.91MB   0.0834     5     8       1.6m <named list [18]> <Rprofmem [84,050 × 3]>  <bench_tm [5]> <tibble [5 × 3]>

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

