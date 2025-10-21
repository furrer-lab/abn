source("scripts/setup_score_cache_mle.R")
source("scripts/cache_creation_refactoring.R")
params <- readRDS(file = "tests/testthat/testdata/input_score_cache_mle.RData")

# Note:
# multinomial adaptation is not included/tested, need to be included in the forLoopNode_maria function defined in cache_creation_refactoring.R

params$max.parents <- unique(params$max.parents)

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
 
## Microbenchmarking
# It also check if the outputs are the same, it returns error otherwise
x <- runif(100)
(lb <- bench::mark(
  computeCache_orig(mleparams$nvars, mleparams$data.df, params$max.parents),
  computeCache_maria(mleparams$nvars, mleparams$data.df, params$max.parents)))

# Results of the benchmark:
# A tibble: 2 × 13
#expression                                           min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result           memory               time           gc
#<bch:expr>                                      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>           <list>               <list>         <list>
#  1 computeCache_orig(nvars, data.df, max.parents)    55.7ms   56.6ms      16.7      47MB    27.9      3     5      179ms <named list [2]> <Rprofmem [701 × 3]> <bench_tm [8]> <tibble [8 × 3]>
#  2 computeCache_maria(nvars, data.df, max.parents)   49.4ms     51ms      19.6    19.5MB     9.81     6     3      306ms <named list [2]> <Rprofmem [869 × 3]> <bench_tm [9]> <tibble [9 × 3]>
#

