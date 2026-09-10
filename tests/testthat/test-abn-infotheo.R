test_that("discretization() works", {
  skip_if_not_installed("entropy")

  mydists <- list(a="gaussian",
                 b="gaussian",
                 c="gaussian")
  a <- rnorm(1000, mean = 1, sd = 0.1)
  b <- a + rnorm(1000, mean = 1, sd = 0.1)
  c <- b + rnorm(1000, mean = 1, sd = 0.1)
  mydf <- data.frame("a" = a,
                     "b" = b,
                     "c" = c)
  mycache.mle <- buildScoreCache(data.df = mydf, data.dists = mydists, method = "mle", max.parents = 1)
  mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
  myfit <- fitAbn(object = mydag.mp, method = "mle")
  out <- invisible(
    simulateAbn(object = myfit, verbose = FALSE)
  )

  y2d.entropy.1 <- entropy::discretize2d(x1=out$a, x2=out$b, numBins1=100, numBins2=100)
  y2d.abn.1 <- abn::discretization(data.df=out[, c(1, 2)], discretization.method=100, data.dists=mydists[c(1, 2)], nb.states=FALSE)

  y2d.entropy.2 <- entropy::discretize2d(x1=out$a, x2=out$c, numBins1=100, numBins2=100)
  y2d.abn.2 <- abn::discretization(data.df=out[, c(1, 3)], discretization.method=100, data.dists=mydists[c(1, 2)], nb.states=FALSE)

  # Not the same dimnames!
  dimnames(y2d.abn.1) <- dimnames(y2d.entropy.1)
  dimnames(y2d.abn.2) <- dimnames(y2d.entropy.2)

  expect_equal(y2d.entropy.1, y2d.abn.1)
  expect_equal(y2d.abn.2, y2d.entropy.2)
})

test_that("miData() works", {
  skip_if_not_installed("entropy")

  mydists <- list(a="binomial",
                  b="binomial",
                  c="binomial")
  a <- rbinom(1000, size = 1, prob = 0.75)
  b <- rbinom(1000, size = 1, prob = 0.75)
  c <- rbinom(1000, size = 1, prob = 0.75)
  mydf <- data.frame("a" = as.factor(a),
                     "b" = as.factor(b),
                     "c" = as.factor(c))
  mycache.mle <- buildScoreCache(data.df = mydf,
                                 data.dists = mydists,
                                 method = "mle",
                                 max.parents = 1,
                                 dag.retained = matrix(c(0,1,0,0,0,1,0,0,0),
                                                       nrow = 3, byrow = TRUE))
  mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
  myfit <- fitAbn(object = mydag.mp, method = "mle")
  out <- invisible(
    simulateAbn(object = myfit, verbose = FALSE)
  )

  x.bc <- table(as.numeric(as.character(out$b)), as.numeric(as.character(out$c)), dnn=c("b", "c"))

  x.ab <- table(as.numeric(out$a), as.numeric(out$b), dnn=c("a", "b"))

  x.ac <- table(as.numeric(out$a), as.numeric(out$c), dnn=c("a", "c"))

  ## MUTUAL INFORMATION

  expect_equal(miData(freqs.table=x.bc, method="mi.raw"), entropy::mi.empirical(y2d=x.bc))
  expect_equal(miData(freqs.table=x.ab, method="mi.raw"), entropy::mi.empirical(y2d=x.ab))
  expect_equal(miData(freqs.table=x.ac, method="mi.raw"), entropy::mi.empirical(y2d=x.ac))
})

test_that("entropyData() works", {
  skip_if_not_installed("entropy")

  mydists <- list(a="gaussian",
                  b="gaussian",
                  c="gaussian")
  a <- rnorm(1000, mean = 1, sd = 0.1)
  b <- a + rnorm(1000, mean = 1, sd = 0.1)
  c <- b + rnorm(1000, mean = 1, sd = 0.1)
  mydf <- data.frame("a" = a,
                     "b" = b,
                     "c" = c)
  mycache.mle <- buildScoreCache(data.df = mydf, data.dists = mydists, method = "mle", max.parents = 1)
  mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
  myfit <- fitAbn(object = mydag.mp, method = "mle")
  out <- invisible(
    simulateAbn(object = myfit, verbose = FALSE)
  )

  y2d.entropy.1 <- entropy::discretize2d(x1=out$a, x2=out$b, numBins1=100, numBins2=100)
  y2d.abn.1 <- abn::discretization(data.df=out[, c(1, 2)], discretization.method=100, data.dists=mydists[c(1, 2)], nb.states=FALSE)

  y2d.entropy.2 <- entropy::discretize2d(x1=out$a, x2=out$c, numBins1=100, numBins2=100)
  y2d.abn.2 <- abn::discretization(data.df=out[, c(1, 3)], discretization.method=100, data.dists=mydists[c(1, 2)], nb.states=FALSE)

  # Not the same dimnames!
  dimnames(y2d.abn.1) <- dimnames(y2d.entropy.1)
  dimnames(y2d.abn.2) <- dimnames(y2d.entropy.2)

  expect_equal(y2d.abn.1, y2d.entropy.1)
  expect_equal(y2d.abn.2, y2d.entropy.2)

  expect_equal(entropyData(freqs.table=y2d.abn.1), entropy::entropy.empirical(y=y2d.entropy.1, unit="log2"))
  expect_equal(entropyData(freqs.table=y2d.abn.2), entropy::entropy.empirical(y=y2d.entropy.2, unit="log2"))
})
