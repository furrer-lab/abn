source("scripts/setup_score_cache_mle.R")
source("scripts/cache_creation_refactoring.R")
source("scripts/forLoopContent_refactoring.R")
require(foreach)
params <- readRDS(file = "tests/testthat/testdata/input_score_cache_mle.RData")
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

cache_orig = computeCache_orig(params$adj.vars, mleparams$nvars, mleparams$data.df,
                               mleparams$data.df.lv, params$max.parents, params$data.dists)
cache_precomputeXY = computeCache_precomputeX(params$adj.vars, mleparams$nvars, mleparams$data.df,
                               mleparams$data.df.multi, mleparams$data.df.lv, params$max.parents, params$data.dists)

# Check GC status for the original code
res <- bench::mark(call_forLoopContent_orig(params$adj.vars,mleparams$nvars,mleparams$data.df,
                                            mleparams$data.df.multi, mleparams$data.df.lvl,
                                            params$max.parents,params$data.dists,cache_orig,
                                            params$control,params$verbose))
# Warning message:
#  Some expressions had a GC in every iteration; so filtering is disabled.
# A tibble: 1 × 13
# expression     min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory
# <bch:expr>   <bch> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>
#  1 call_forLoo…   27s    27s    0.0370    1.83GB    0.926     1    25        27s <named list> <Rprofmem>
#  # ℹ 2 more variables: time <list>, gc <list>

res <- bench::mark(call_forLoopContent_orig(params$adj.vars,mleparams$nvars,mleparams$data.df,
                                            mleparams$data.df.multi, mleparams$data.df.lvl,
                                            params$max.parents,params$data.dists,cache_orig,
                                            params$control,params$verbose),
                   call_forLoopContent_foreachrows(params$adj.vars,mleparams$nvars,mleparams$data.df,
                                            mleparams$data.df.multi, mleparams$data.df.lvl,
                                            params$max.parents,params$data.dists,cache_orig,
                                            params$control,params$verbose),
                   call_forLoopContent_noforeach(params$adj.vars,mleparams$nvars,mleparams$data.df,
                                                   mleparams$data.df.multi, mleparams$data.df.lvl,
                                                   params$max.parents,params$data.dists,cache_orig,
                                                   params$control,params$verbose),
                   iterations = 50)

# A tibble: 2 × 13
#expression     min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory
#<bch:expr>   <bch> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>
#  1 call_forLoo… 26.8s  26.8s    0.0373    1.83GB     1.12     1    30      26.8s <named list> <Rprofmem>
#  2 call_forLoo…   25s    25s    0.0400    1.78GB     1.12     1    28        25s

# A tibble: 3 × 13
#expression                                                                                             min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result       memory     time       gc
#<bch:expr>                                                                                           <bch> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>       <list>     <list>     <list>
#  1 call_forLoopContent_orig(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.df.mul… 26.1s  26.4s    0.0377    1.83GB    0.489    50   649      22.1m <named list> <Rprofmem> <bench_tm> <tibble>
#  2 call_forLoopContent_foreachrows(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data… 26.3s  26.4s    0.0361    1.82GB    0.444    50   615      23.1m <named list> <Rprofmem> <bench_tm> <tibble>
#  3 call_forLoopContent_noforeach(params$adj.vars, mleparams$nvars, mleparams$data.df, mleparams$data.d… 24.4s  24.8s    0.0328    1.78GB    0.315    50   480      25.4m <named list> <Rprofmem> <bench_tm> <tibble>


##########
## Test X and Y preallocation

res <- bench::mark( orig = call_forLoopContent_orig(params$adj.vars,mleparams$nvars,mleparams$data.df,
                                                    mleparams$data.df.multi, mleparams$data.df.lvl,
                                                    params$max.parents,params$data.dists,cache_orig,
                                                    params$control,params$verbose),
                    precomputeXY = call_forLoopContent_precomputedXY(params$adj.vars,mleparams$nvars,mleparams$data.df,
                                                                 mleparams$data.df.multi, mleparams$data.df.lvl,
                                                                 params$max.parents,params$data.dists,cache_orig,
                                                                 params$control,params$verbose),

                    iterations = 5)
# A tibble: 2 × 13
# expression        min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result            memory                   time       gc
#<bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>            <list>                   <list>     <list>
#  1 orig              25s      25s    0.0400    1.84GB   0.408      5    51      2.08m <named list [14]> <Rprofmem [245,167 × 3]> <bench_tm> <tibble>
#  2 precomputeXY      19s    19.1s    0.0525  911.79MB   0.0839     5     8      1.59m <named list [14]> <Rprofmem [84,834 × 3]>  <bench_tm> <tibble>



##########
## Test the single function run against the matrix allocation

adj.vars = params$adj.vars
nvars = mleparams$nvar
data.df = mleparams$data.df
data.df.multi = mleparams$data.df.multi
data.df.lvl = mleparams$data.df.lvl
max.parents = params$max.parents
data.dists = params$data.dists
mycache = cache_orig
control = params$control
verbose = params$verbose
row.num <- NULL   # To avoid check comment: 'no visible binding for global variable
out <- list()
rows <- length(mycache[["children"]])
res_test <- bench::mark(
  only_rcpp = {
    for (row.num in 1:rows) {
      forLoopContent(row.num = row.num,
                     mycache = mycache,
                     data.dists = data.dists,
                     data.df.multi = data.df.multi,
                     adj.vars = NULL, #adj.vars,
                     data.df = data.df,
                     data.df.lvl = data.df.lvl,
                     group.var = NULL, #group.var,
                     group.ids = NULL, #group.ids,
                     control = control,
                     n = nvars,
                     verbose = FALSE)
    }
  },
  full_for = {
    res <- matrix(NA_real_, rows, 4)
    for (row.num in 1:rows) {
      res[row.num,] <- forLoopContent(row.num = row.num,
                                      mycache = mycache,
                                      data.dists = data.dists,
                                      data.df.multi = data.df.multi,
                                      adj.vars = NULL, #adj.vars,
                                      data.df = data.df,
                                      data.df.lvl = data.df.lvl,
                                      group.var = NULL, #group.var,
                                      group.ids = NULL, #group.ids,
                                      control = control,
                                      n = nvars,
                                      verbose = FALSE)
    }
  }
)
#Warning message:
#  Some expressions had a GC in every iteration; so filtering is disabled.
#> res_test
## A tibble: 2 × 13
#expression      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
#<bch:expr> <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
#  1 only_rcpp     51.8s  51.8s    0.0193    1.78GB    0.521     1    27      51.8s
# 2 full_for      50.3s  50.3s    0.0199    1.78GB    0.517     1    26      50.3s
## ℹ 4 more variables: result <list>, memory <list>, time <list>, gc <list>

# If full_for is only a bit slower than only_rcpp, then the allocation solution is no longer the main problem; the forLoopContent function is.

