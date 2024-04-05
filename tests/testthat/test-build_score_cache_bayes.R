test_that("General behaviour of buildScoreCache.bayes()", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  df <- airquality[complete.cases(airquality), ]

  # distribution (gaussian)
  dist <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
  names(dist) <- colnames(df)

  ## which.nodes
  expect_warning({
    expect_no_error(buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=6, which.nodes=1:6))
  })

  ## max.parents
  expect_error(buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=6, which.nodes=1:5))
  expect_warning({
    expect_no_error(buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=5, which.nodes=1:5))
  })
  expect_error(buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=list(Ozone=NULL, Solar.R=4, Wind=4, Temp=4, Month=4), which.node=1:5))
  expect_equal(buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=4, which.nodes=1:5),
               buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=list(Ozone=4, Solar.R=4, Wind=4, Temp=4, Month=4), which.node=1:5))
  expect_equal(buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=3, which.nodes=1:5),
               buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=list(Ozone=3, Solar.R=3, Wind=3, Temp=3, Month=3), which.node=1:5))
  expect_no_error(buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=list(Ozone=2, Solar.R=4, Wind=4, Temp=4, Month=4), which.node=1:5))
  expect_error(
    expect_equal(buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=4, which.node=1:5),
                 buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=list(Ozone=2, Solar.R=4, Wind=4, Temp=4, Month=4), which.node=1:5)),
    regexp = "not equal")
})

test_that("buildScoreCache.bayes() with Gaussian nodes", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

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
  expect_no_error({
    mycache <- invisible(buildScoreCache(data.df = mydf,
                                         data.dists = mydists,
                                         method = "bayes",
                                         max.parents = 2,
                                         centre = FALSE))
  })

  # mLik
  skip("Numerical results are not tested for method='bayes'")
  ### ISSUE: numerical testing with method bayes
  # In the "mle" case, we could compare the numerical output from buildScoreCache() with a linear model (see `test-build_score_cache_mle.R`):
  # expect_equal(mycache$mlik[1], as.numeric(logLik(lm(formula=out.sim$a ~ 1))))
  # In the "bayes" case, we need a "Bayesian model". There are different ways to fit such a model. Abn uses INLA or some internal code.
  # We can force abn to use INLA but if INLA fails, abn automatically switches to its internal C code (check mycache$used.INLA) anyway.
  # This makes it hard to create reliable tests.
  # Possible approaches:
  #   - use a test scenario where INLA does not fail an can safely be enforced.
  #   - use internal C code as reference to check (but if there is an issue in the internal C code, we won't catch it...)
  #
  # inla.out <- INLA::inla(formula = a ~ 1, family = "gaussian", data = out.sim)
  # expect_equal(mycache$mlik[1], as.numeric(inla.out[["mlik"]][2]))
  ### EOF ISSUE documentation
  expect_equal(mycache$mlik[2], as.numeric(logLik(lm(formula=mydf$a ~ 1 + mydf$b))))
  expect_equal(mycache$mlik[3], as.numeric(logLik(lm(formula=mydf$b ~ 1))))

  # AIC
  expect_equal(mycache$aic[1], as.numeric(AIC(lm(formula=mydf$a ~ 1))))
  expect_equal(mycache$aic[2], as.numeric(AIC(lm(formula=mydf$a ~ 1 + mydf$b))))
  expect_equal(mycache$aic[3], as.numeric(AIC(lm(formula=mydf$b ~ 1))))

  # BIC
  expect_equal(mycache$bic[1], as.numeric(BIC(lm(formula=mydf$a ~ 1))))
  expect_equal(mycache$bic[2], as.numeric(BIC(lm(formula=mydf$a ~ 1 + mydf$b))))
  expect_equal(mycache$bic[3], as.numeric(BIC(lm(formula=mydf$b ~ 1))))

  ## Gaussian: correlation coefficients

  expect_no_error({
    mycache <- invisible(buildScoreCache(data.df=mydf, data.dists=mydists, method = "bayes", max.parents=2, centre=TRUE))
  })

  expect_equal(mycache$mlik[1], as.numeric(logLik(lm(formula=(mydf$a - mean(mydf$a))/sd(mydf$a) ~ 1))))
})

test_that("buildScoreCache.bayes() with Binomial nodes", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  # Binomial
  mydists <- list(a="binomial",
                  b="binomial")
  a <- rbinom(1000, size = 1, prob = 0.5)
  z <- 1+2*a
  pr <- 1/(1+exp(-z))
  b <- rbinom(1000, size = 1, prob = pr)
  mydf <- data.frame("a" = as.factor(a),
                     "b" = as.factor(b))
  expect_no_error({
    mycache.mle <- buildScoreCache(data.df = mydf,
                                   data.dists = mydists,
                                   method = "bayes",
                                   max.parents = 1)
  })

  skip("Numerical results are not tested for method='bayes'")
  # mLik
  expect_equal(mycache$mlik[1], as.numeric(logLik(glm(formula=out.sim$a ~ 1, family=binomial))), tolerance = 0.00005)
  expect_equal(mycache$mlik[2], as.numeric(logLik(glm(formula=out.sim$a ~ 1 + out.sim$b, family=binomial))), tolerance = 0.00005)
  expect_equal(mycache$mlik[3], as.numeric(logLik(glm(formula=out.sim$b ~ 1, family=binomial))), tolerance = 0.00005)

  # AIC
  expect_equal(mycache$aic[1], as.numeric(AIC(glm(formula=out.sim$a ~ 1, family=binomial))), tolerance = 0.00005)
  expect_equal(mycache$aic[2], as.numeric(AIC(glm(formula=out.sim$a ~ 1 + out.sim$b, family=binomial))), tolerance = 0.00005)
  expect_equal(mycache$aic[3], as.numeric(AIC(glm(formula=out.sim$b ~ 1, family=binomial))), tolerance = 0.00005)

  # BIC
  expect_equal(mycache$bic[1], as.numeric(BIC(glm(formula=out.sim$a ~ 1, family=binomial))), tolerance = 0.00005)
  expect_equal(mycache$bic[2], as.numeric(BIC(glm(formula=out.sim$a ~ 1 + out.sim$b, family=binomial))), tolerance = 0.00005)
  expect_equal(mycache$bic[3], as.numeric(BIC(glm(formula=out.sim$b ~ 1, family=binomial))), tolerance = 0.00005)
})

test_that("buildScoreCache.bayes() simple, historic numeric test", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  load(file="testdata/buildscorecache_ex1.Rdata")
  # load(file='tests/testthat/testdata/buildscorecache_ex1.Rdata')

  invisible(mycache.test <- buildScoreCache(data.df=mydat, data.dists=mydists, method = "bayes", max.parents=max.par))

  expect_equal(mycache.test[1:3], mycache[1:3])
})
