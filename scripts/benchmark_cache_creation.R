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
(lb <- bench::mark(
cache_orig = computeCache_orig(mleparams$nvars, mleparams$data.df, params$max.parents),
cache_mod = computeCache_maria(mleparams$nvars, mleparams$data.df, params$max.parents), min_time = 1))
# Results of the benchmark:
# A tibble: 2 × 13
#   expression     min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory     time       gc
#   <bch:expr> <bch:t> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>     <list>     <list>
# 1 cache_orig    56ms 57.6ms      17.4      47MB     7.90    11     5      633ms <named list> <Rprofmem> <bench_tm> <tibble>
# 2 cache_mod     48ms 49.5ms      20.2    9.78MB     3.57    17     3      841ms <named list> <Rprofmem> <bench_tm> <tibble>


# Results of the the benchmark when the rbind optimization is implemented in bothcode versions:
# A tibble: 2 × 13
#  expression    min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory     time       gc
#   <bch:expr> <bch:> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>     <list>     <list>
# 1 cache_orig   44ms 45.1ms      22.2    9.78MB     5.22    17     4      766ms <named list> <Rprofmem> <bench_tm> <tibble>
# 2 cache_mod  44.1ms   45ms      22.2    9.78MB     3.50    19     3      856ms <named list> <Rprofmem> <bench_tm> <tibble>
#

# To visualize results:
lb %>%
unnest(c(time, gc)) %>%
filter(gc == "none") %>%
mutate(expression = as.character(expression)) %>%
ggplot(aes(x = mem_alloc, y = time, color = expression)) +
geom_point() +
bench::scale_color_bench_expr(scales::brewer_pal(type = "qual", palette = 3))
