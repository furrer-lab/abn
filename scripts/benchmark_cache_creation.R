source("scripts/setup_score_cache_mle.R")
source("scripts/cache_creation_refactoring.R")
params <- readRDS(file = "tests/testthat/testdata/input_score_cache_mle.RData")

# Note:
# adjustment is not included/tested in the functions defined in "scripts/cache_creation_refactoring.R"

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
res <- bench::mark(setupScoreCache.mle(data.df = params$data.df,
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
