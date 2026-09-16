source("scripts/setup_score_cache_mle.R")
source("scripts/cache_creation_refactoring.R")
params <- readRDS(file = "tests/testthat/testdata/input_score_cache_mle.RData")

# Note:
# adjustment is not included/tested in the functions defined in "scripts/cache_creation_refactoring.R"
max_parents <- params$max.parents
params$max.parents <- unique(params$max.parents)

mleparams <- setupScoreCache.mle_orig(data.df = params$data.df,
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

mleparams_mod <- setupScoreCache.mle_mod(data.df = params$data.df,
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

# all.equal(mleparams$data.df.multi, mleparams_mod$data.df.multi, check.attributes = FALSE)
res <- bench::mark(setupScoreCache.mle_orig(data.df = params$data.df,
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
                   setupScoreCache.mle_mod(data.df = params$data.df,
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
                   min_time = 0.5,
                   check = function(x, y) all.equal(x, y, check.attributes = FALSE),
                   iterations = 100)
# Results of the benchmark:
# A tibble: 2 × 13
# expression     min  median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result
# <bch:expr>  <bch:> <bch:t>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>
#  1 setupScore… 9.67ms 10.04ms      96.5    1.63MB     2.98    97     3      1.01s <named list>
#  2 setupScore… 4.64ms  4.81ms     207.     1.92MB     4.22    98     2   473.78ms <named list>
#  # ℹ 3 more variables: memory <list>, time <list>, gc <list>
# The number of iterations not affected by GC is already large, 97%. There is not much improvement on GC calls.
# With filter_gc=FALSE: res_orig$median / res_mod$median = 2x
# The original code is fast enough but the modified one with preallocation is still 2x faster.
# In terms of memory there is an increase of 17%: res_mod$mem_alloc=1.92MB, res_orig$mem_alloc=1.63MB
# because for small matrices the cbind's copying cost is tiny.
# If the size of the matricies are >100×500 the modified version with preallocation will scale better.

## Microbenchmarking cache creation
# It also check if the outputs are the same, it returns error otherwise
# Compare cache creation: pure R using combn VS C implementation also used by the Bayes method
(lb <- bench::mark(
cache_orig = computeCache_orig(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.lv, params$max.parents, params$data.dists),
cache_doCall = computeCache_doCall(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.lv, params$max.parents, params$data.dists),
cache_inForLoop = computeCache_inForLoop(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.lv, params$max.parents, params$data.dists),
cache_prealloc = computeCache_prealloc(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.lv, params$max.parents, params$data.dists),
cache_prealloc_mod = computeCache_prealloc_mod(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.lv, params$max.parents, params$data.dists),
cache_C = computeCache_C(params$adj.vars, mleparams$nvars, mleparams$data.df,
                               mleparams$data.df.multi, mleparams$data.df.lv, max_parents, params$data.dists,
                               params$defn.res, params$dag.banned, params$dag.retained, params$which.nodes),
min_time = 0.5,
iterations = 100))
# Results of the benchmark:
#expression         min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory
#<bch:expr>      <bch:> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>
#  1 cache_orig      69.8ms 70.5ms      14.2    48.2MB   125.      11    97    777.1ms <named list> <Rprofmem>
#  2 cache_doCall    60.1ms 61.7ms      16.2    10.9MB     7.61    68    32      4.21s <named list> <Rprofmem>
#  3 cache_inForLoop 60.1ms 61.6ms      16.2    10.9MB     7.63    68    32      4.19s <named list> <Rprofmem>
#  # ℹ 2 more variables: time <list>, gc <list>
# comparing the original (cbind-in-loop), with the doCall (list preallocation) and the refactored cache in the for loop:
# - cahce_orig has only 11% of iterations not affected by GC calls and intermittent "GC every iteration" warning
# - both doCall and inForLoop have NO "GC every iteration" warning, with about 70% of clean iterations, so decent result
# - using a matrix preallocation instead of a list preallocation might improve.
# Memory allocation is anyhow reduced.

# A tibble: 4 × 13
#expression           min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory
#<bch:expr>      <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>
#  1 cache_orig        75.8ms 84.4ms      11.4   47.88MB    17.6    100   155      8.81s <named list> <Rprofmem>
#  2 cache_doCall      59.7ms 64.6ms      14.1   10.66MB     6.21   100    44      7.08s <named list> <Rprofmem>
#  3 cache_inForLoop   59.6ms 63.3ms      14.7   10.91MB     6.16   100    42      6.82s <named list> <Rprofmem>
#  4 cache_prealloc    37.6ms 39.8ms      23.2    9.38MB     8.83   100    38       4.3s <named list> <Rprofmem>
# In cache_prealloc, the function has been refactored to use matrix preallocation. The calculation of the `tmp`
# matrices has been refactored to reduce duplicated computation.

#expression           min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory
#<bch:expr>      <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>
#  1 cache_orig        77.4ms 84.5ms      11.3   47.88MB    16.7    100   148      8.88s <named list> <Rprofmem>
#  2 cache_doCall      60.8ms 67.2ms      14.2   10.66MB     6.69   100    47      7.03s <named list> <Rprofmem>
#  3 cache_inForLoop   39.5ms 41.4ms      21.1   11.21MB     9.07   100    43      4.74s <named list> <Rprofmem>
#  4 cache_prealloc    38.7ms 40.6ms      21.7    9.38MB     8.01   100    37      4.62s <named list> <Rprofmem>
# Having the cache calculation in the forLoopContent function, doesn't seem performing better than the cache_prealloc
# Performance has to be compared when this logic is integrated in the fit calculation, which would change the length
# of the foreach loop and would avoid passing the cache matrix to the forLoopContent function at every iteration.

#   expression     min  median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
#   <bch:expr> <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
# 1 cache_orig 70.12ms 74.03ms      13.0   47.88MB    20.2    100   156      7.71s
# 2 cache_doC…  56.7ms 60.95ms      15.9   10.66MB     6.99   100    44      6.29s
# 3 cache_inF… 37.07ms 38.14ms      25.0   11.21MB    10.3    100    41      3.99s
# 4 cache_pre… 35.81ms 37.98ms      24.5    9.38MB     8.58   100    35      4.08s
# 5 cache_C     4.05ms  4.24ms     221.     2.04MB    11.1    100     5   452.02ms
# # ℹ 4 more variables: result <list>, memory <list>, time <list>, gc <list>
# Warning message:
# Some expressions had a GC in every iteration; so filtering is disabled.
# The GC warning is triggered by the R implementation in computeCache_orig, the cache_C implementation would solve the issue.

#   expression     min  median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
#   <bch:expr> <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
# 1 cache_orig 68.22ms 68.88ms      14.6   47.88MB   602.       3   124   206.04ms
# 2 cache_doC… 56.83ms  59.9ms      16.8   10.66MB    10.7     61    39      3.63s
# 3 cache_inF… 36.94ms 38.16ms      26.2   11.21MB    11.8     69    31      2.63s
# 4 cache_pre… 35.74ms 37.74ms      26.6    9.38MB    11.9     69    31       2.6s
# 5 cache_pre…  3.09ms  3.13ms     306.     4.27MB    23.0     93     7   303.69ms
# 6 cache_C     4.03ms   4.2ms     233.     2.04MB     9.69    96     4   412.75ms
# # ℹ 4 more variables: result <list>, memory <list>, time <list>, gc <list>
# The cache_prealloc_mod function instead of shifting each element to the temp.list[[i]] matrix is filling the
# node.defn_all.list except for the column -j which is already initialized with 0.

#   expression     min  median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
#   <bch:expr> <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
# 1 cache_orig 70.95ms 74.43ms      12.8   47.88MB    17.1    100   134      7.83s
# 2 cache_doC… 57.06ms 60.61ms      15.8   10.66MB     6.95   100    44      6.33s
# 3 cache_inF… 37.55ms 38.87ms      24.7   11.21MB     8.89   100    36      4.05s
# 4 cache_pre… 36.19ms 38.27ms      24.3    9.38MB     8.03   100    33      4.11s
# 5 cache_pre…  2.09ms  2.14ms     422.     1.74MB    12.7    100     3   237.02ms
# 6 cache_C     4.09ms  4.39ms     215.     2.04MB     8.59   100     4   465.53ms
# New version of the cache_prealloc_mod function: the loop `for (i in 1:(max.parents))` is removed and the `combs <- do.call(rbind, tmp.list)` matrix
# is used to compute `node.defn` and the `children` vector



# To visualize results:
library(ggplot2)
library(dplyr)
library(tidyr)
lb %>%
unnest(c(time, gc)) %>%
filter(gc == "none") %>%
mutate(expression = as.character(expression)) %>%
ggplot(aes(x = mem_alloc, y = time, color = expression)) +
geom_point() +
bench::scale_color_bench_expr(scales::brewer_pal(type = "qual", palette = 3))
