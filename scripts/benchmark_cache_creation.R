source("scripts/setup_score_cache_mle.R")
source("scripts/cache_creation_refactoring.R")
params <- readRDS(file = "tests/testthat/testdata/input_score_cache_mle.RData")

# Note:
# adjustment is not included/tested in the functions defined in "scripts/cache_creation_refactoring.R"

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
(lb <- bench::mark(
cache_orig = computeCache_orig(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.lv, params$max.parents, params$data.dists),
cache_doCall = computeCache_doCall(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.lv, params$max.parents, params$data.dists),
cache_inForLoop = computeCache_inForLoop(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.lv, params$max.parents, params$data.dists),
cache_prealloc = computeCache_prealloc(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.lv, params$max.parents, params$data.dists),
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
