library(abn)
library(nnet, lib.loc = "/usr/lib/R/library")
library(foreach)
library(profvis)

load("data/n_1000_p_16_nodedists_Balanced_s_02_graph1.Rdata")
mydf <- n_1000_p_16_nodedists_Balanced_s_02_graph1$data
mydists <- n_1000_p_16_nodedists_Balanced_s_02_graph1$dists

source("scripts/setup_score_cache.R")
# Uncomment the lines below to run the profiler
#profvis({
params <- setupScoreCache(data.df = mydf,
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
                control = NULL,
                verbose = FALSE,
                debugging = FALSE)
#})

saveRDS(params, file = "tests/testthat/testdata/input_score_cache_mle.RData")


out <-
      buildScoreCache.mle(
        data.df = params$data.df,
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
        control = params$control
      )

saveRDS(out, file = "tests/testthat/testdata/out_build_score_cache_mle.RData")

