test_that("linkStrength() works", {
  # Gaussian
  N <- 1000
  mydists <- list(a="gaussian",
                  b="gaussian",
                  c="gaussian")
  a <- rnorm(n = N, mean = 0, sd = 1)
  b <- 1 + 2*rnorm(n = N, mean = 5, sd = 1)
  c <- 2 + 1*a + 2*b + rnorm(n = N, mean = 2, sd = 1)
  mydf <- data.frame("a" = a,
                     "b" = b,
                     "c" = c)
  mycache.mle <- buildScoreCache(data.df = mydf,
                                 data.dists = mydists,
                                 method = "mle",
                                 max.parents = 2)
  mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
  expect_no_error({
    linkstr <- linkStrength(dag = mydag.mp$dag,
                             data.df = mydf,
                             data.dists = mydists,
                             method = "ls",
                             discretization.method = "sturges")
  })
  expect_equal(dim(linkstr), c(3,3))
})
