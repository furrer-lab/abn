test_that("searchHeuristic() works.", {
  mydat <- ex1.dag.data ## this data comes with abn see ?ex1.dag.data

  ## setup distribution list for each node
  mydists<-list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                p2="poisson", b3="binomial", g2="gaussian", b4="binomial",
                b5="binomial", g3="gaussian")

  mycache <- buildScoreCache(data.df = mydat, data.dists = mydists, max.parents = 2)

  ## Now peform 10 greedy searches
  expect_no_error({
    heur.res <- searchHeuristic(score.cache = mycache, data.dists = mydists,
                                start.dag = "random", num.searches = 10,
                                max.steps = 50)
  })
})
