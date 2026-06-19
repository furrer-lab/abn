# Performance comparison for the original buildScoreCache() function when multi cores are used.
# Use syrup library for profiling parallel code 
# Run: Rscript scripts/benchmark_orig_ncores.R
# Note: Failing case removed from the computation: child 14, poisson;  parents:     "b1" "g2" "g3" "p1"; function: `/src/irls_poisson_fast.cpp`

library(abn)                     
library(parallel)                
library(doParallel)              
library(purrr)
library(syrup)
library(dplyr)
library(future)
source("scripts/setup_score_cache.R")
load("data/n_1000_p_16_nodedists_Balanced_s_02_graph1.Rdata")
mydf <- n_1000_p_16_nodedists_Balanced_s_02_graph1$data[,-c(14)]
mydists <- n_1000_p_16_nodedists_Balanced_s_02_graph1$dists[-c(14)]

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

grid <- expand.grid(maxp = c(3, 4, 5, 6, 7), ncores = c(2, 3, 4, 1))

run_one <- function(ncores, maxp) {
  print(ncores)
  params <- get_params(mydf, mydists, maxp, ncores)
  if (ncores > 1){
	  future::plan(multicore, workers = ncores) 
  }
  mem <- syrup({
    bsc.mle.orig <- buildScoreCache.mle(data.df = params$data.df,
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
    invisible(bsc.mle.orig)
  }, interval = 0.5, peak = TRUE) #FALSE
   mem = mem %>%
    mutate(ncores = ncores, maxp = maxp)
  saveRDS(mem, paste("bench_memory_syrup_results", ncores, maxp, "peak.rds", sep="_"))
  mem
}

res_syrup <- pmap(grid, run_one) |>
  list_rbind()

res_syrup
saveRDS(res_syrup, "bench_memory_syrup_results_peakTRUE.rds")

worker_ppid <- ps::ps_pid()
worker_ppid


library(ggplot2)
p <- res_syrup %>%
  filter(ppid == worker_ppid | pid == worker_ppid) %>%
  ggplot() +
  aes(x = id, y = rss, group = pid) +
  geom_line() +
  scale_x_continuous(breaks = 1:max(res_syrup$id))
p2 <- res_syrup %>%
  filter(ppid == worker_ppid | pid == worker_ppid) %>%
  mutate(combo = interaction(maxp, ncores, sep = ", ")) %>%
  ggplot(aes(x = id, y = rss, group = pid, colour = combo)) +
  geom_line() +
  scale_x_continuous(breaks = 1:max(res_syrup$id)) +
  labs(colour = "maxp, ncores")

res_syrup2 <- res_syrup %>%
  filter(ppid == worker_ppid | pid == worker_ppid) %>%
  group_by(maxp, ncores) %>%
  mutate(elapsed = time - min(time, na.rm = TRUE),
         combo = interaction(maxp, ncores, sep = "_")) %>%
  ungroup()
p3 <- res_syrup2 %>%
  ggplot(aes(x = elapsed, y = rss, group = pid, colour = combo)) +
  geom_line() +
  labs(colour = "maxp, ncores")
ggsave("plot_syrup_rss_peakTRUE.pdf", p, width = 6, height = 4)
ggsave("plot_syrup_rss_legend_peakTRUE.pdf", p2, width = 6, height = 4)
ggsave("plot_syrup_rss_time_peakTRUE.pdf", p3, width = 6, height = 4)

p <- res_syrup %>%
  filter(ppid == worker_ppid | pid == worker_ppid) %>%
  ggplot() +
  aes(x = id, y = pct_cpu, group = pid) +
  geom_line() +
  scale_x_continuous(breaks = 1:max(res_syrup$id))
p2 <- res_syrup %>%
  filter(ppid == worker_ppid | pid == worker_ppid) %>%
  mutate(combo = interaction(maxp, ncores, sep = ", ")) %>%
  ggplot() +
  aes(x = id, y = pct_cpu, group = pid, colour = combo) +
  geom_line() +
  scale_x_continuous(breaks = 1:max(res_syrup$id)) +
  labs(colour = "maxp, ncores")
p3 <- res_syrup2 %>%
  ggplot(aes(x = elapsed, y = pct_cpu, group = pid, colour = combo)) +
  geom_line() +
  labs(colour = "maxp, ncores")
ggsave("plot_syrup_cpu_peakTRUE.pdf", p, width = 6, height = 4)
ggsave("plot_syrup_cpu_legend_peakTRUE.pdf", p2, width = 6, height = 4)
ggsave("plot_syrup_cpu_time_peakTRUE.pdf", p3, width = 6, height = 4)

