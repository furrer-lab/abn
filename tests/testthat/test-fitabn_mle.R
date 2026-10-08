test_that("Test fitAbn.mle() with real Gaussian nodes", {
  suppressMessages({
    if(.Platform$OS.type == "unix") {
      capture.output({
        # Gaussian
        df <- airquality[complete.cases(airquality), ]

        dist <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
        names(dist) <- colnames(df)

        d <- matrix(data=0, nrow=6, ncol=6)
        d[1, ] <- c(0, 1, 1, 1, 1, 1)
        colnames(d) <- rownames(d) <- names(dist)

        m1 <- fitAbn(method = "mle", dag=d, data.df=df, data.dists=dist, centre=FALSE)
        m2 <- lm(df[, 1] ~ as.matrix(df[, 2:6]))

        expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.lm(object=m2))[, 1])))
        expect_equal(unname(m1$Stderror[[1]]), unname(t(coef(summary.lm(object=m2))[, 2])))

        ## test centre
        m1 <- fitAbn(method = "mle", dag=d, data.df=df, data.dists=dist, centre=TRUE)
        m3 <- fitAbn(method = "mle", dag=d, data.df=df, data.dists=dist)
        d[1, ] <- c(0, 1, 0, 0, 0, 0)
        m2 <- fitAbn(method = "mle", dag=d, data.df=df, data.dists=dist)
        m4 <- cor(df[, 1:6])
      },
      file = "/dev/null")
      expect_equal(m1, m3)
      expect_equal(unname(m2$coef[[1]])[2], m4[1, 2])
    } else {
      skip("fitAbn.mle() is tested mainly on Unix-like systems")
    }
  })
})

test_that("Test fitAbn.mle() with simulated Gaussian nodes", {
  suppressMessages({
    if(.Platform$OS.type == "unix") {
      capture.output({
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
        myfit <- fitAbn(object = mydag.mp, method = "mle", centre = FALSE, verbose = FALSE)
        out.sim <- invisible(
          simulateAbn(object = myfit, n.iter = 10000L, verbose = FALSE)
        )
      },
      file = "/dev/null")
      # with original data
      m1 <- myfit
      m2 <- lm(c ~ a + b, data = mydf)
      expect_equal(unname(m1$coef[[3]]), unname(t(coef(summary(object=m2))[, 1])))

      # with simulated data
      capture.output({
        mycache.mle2 <- buildScoreCache(data.df = out.sim,
                                        data.dists = mydists,
                                        method = "mle",
                                        max.parents = 2)
        mydag.mp2 <- mostProbable(score.cache = mycache.mle2, verbose = FALSE)
        m1 <- fitAbn(object = mydag.mp2, method = "mle", centre=FALSE, verbose = FALSE)
        m2 <- lm(c ~ a + b, data = out.sim)
      },
      file = "/dev/null")
      expect_equal(unname(m1$coef[[3]]), unname(t(coef(summary(object=m2))[, 1])))
      ## new: RF 2022 after big merge: tolerance required.
      expect_equal(unname(m1$Stderror[[3]]), unname(t(coef(summary(object=m2))[, 2])), tolerance=1e-06)
    } else {
      skip("fitAbn.mle() is tested mainly on Unix-like systems")
    }
  })
})

test_that("Test fitAbn.mle() with Binomial nodes", {
  suppressMessages({
    if(.Platform$OS.type == "unix") {
      capture.output({
        # Binomial
        mydists <- list(a="binomial",
                        b="binomial")
        a <- rbinom(1000, size = 1, prob = 0.5)
        z <- 1+2*a
        pr <- 1/(1+exp(-z))
        b <- rbinom(1000, size = 1, prob = pr)
        mydf <- data.frame("a" = as.factor(a),
                           "b" = as.factor(b))
        mycache.mle <- buildScoreCache(data.df = mydf,
                                       data.dists = mydists,
                                       method = "mle",
                                       max.parents = 1)
        mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
        myfit <- fitAbn(object = mydag.mp, method = "mle", centre = TRUE, verbose = FALSE)
        out.sim <- invisible(
          simulateAbn(object = myfit, n.iter = 10000L, verbose = FALSE)
        )
      },
      file = "/dev/null")
      # with original data
      m1 <- myfit
      m2 <- glm(a ~ b, data = mydf, family="binomial")
      expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.glm(object=m2))[, 1])))

      # with simulated data
      capture.output({
        mycache.mle2 <- buildScoreCache(data.df = out.sim,
                                        data.dists = mydists,
                                        method = "mle",
                                        max.parents = 1)
        mydag.mp2 <- mostProbable(score.cache = mycache.mle2, verbose = FALSE)
        m1 <- fitAbn(object = mydag.mp2, method = "mle", centre=TRUE, verbose = FALSE)
        m2 <- glm(a ~ b, data = out.sim, family="binomial")
      },
      file = "/dev/null")
      expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.glm(object=m2))[, 1])))
      ## new: RF 2022 after big merge: tolerance required.
      expect_equal(unname(m1$Stderror[[1]]), unname(t(coef(summary.glm(object=m2))[, 2])), tolerance=1e-06)
    } else {
      skip("fitAbn.mle() is tested mainly on Unix-like systems")
    }
  })
})

test_that("Test fitAbn.mle() with Poisson nodes", {
  skip("BuildScoreCache(method=mle) crashes with Poisson nodes")
  suppressMessages({
    # Poisson
    mydists <- list(a="poisson",
                    b="poisson")
    a <- rpois(1000, lambda = 0.5)
    z <- exp(1+2*a)
    b <- rpois(1000, lambda = z)
    mydf <- data.frame("a" = a,
                       "b" = as.integer(b))
    # "c" = as.factor(c(rep(1, 999), 0)))
    # glm(b~a, data = mydf, family = "poisson")
    mycache.mle <- buildScoreCache(data.df = mydf,
                                   data.dists = mydists,
                                   method = "mle",
                                   max.parents = 1,
                                   # group.var = "c",
                                   # control = list(max.mode.error=100),
                                   verbose = FALSE,
                                   dag.retained = matrix(c(0,0,1,0), nrow = 2, byrow = TRUE))
    mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
    myfit <- fitAbn(object = mydag.mp, method = "mle", centre = FALSE)
    out.sim <- invisible(
      simulateAbn(object = myfit, n.iter = 10000L, verbose = TRUE)
    )

    # with original data
    m1 <- myfit
    m2 <- glm(b ~ a, data = mydf, family="poisson")
    expect_equal(unname(m1$coef[[2]]), unname(t(coef(summary.glm(object=m2))[, 1])))

    # with simulated data
    mycache.mle2 <- buildScoreCache(data.df = out.sim,
                                    data.dists = mydists,
                                    method = "mle",
                                    max.parents = 1)
    mydag.mp2 <- mostProbable(score.cache = mycache.mle2, verbose = FALSE)
    m1 <- fitAbn(object = mydag.mp2, method = "mle", centre=TRUE, verbose = FALSE)
    m2 <- glm(a ~ b, data = out.sim, family="binomial")

    expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.glm(object=m2))[, 1])))
    ## new: RF 2022 after big merge: tolerance required.
    expect_equal(unname(m1$Stderror[[1]]), unname(t(coef(summary.glm(object=m2))[, 2])), tolerance=1e-06)
    # dist <- list(a="poisson", b="poisson")
    #
    # data.param <- matrix(data=c(0, 0.5, 0, 0), nrow=2L, ncol=2L, byrow=TRUE)
    # data.param <- matrix(data=c(0, 1, 0, 0), nrow=2L, ncol=2L, byrow=TRUE)
    #
    # # naming
    # colnames(data.param) <- rownames(data.param) <- names(dist)
    #
    # out.sim <- invisible(simulateAbn(data.dists=dist, n.chains=1, n.adapt=100, n.thin=1, n.iter=100, data.param=data.param,
    #                                  simulate=TRUE, seed=132,verbose=FALSE))

    m1 <- fitAbn(method = "mle", dag=data.param, data.df=out.sim, data.dists=dist, centre=FALSE)
    m2 <- glm(formula=out.sim$a ~ out.sim$b, family="poisson")

    ## pvalues and stderr are computed up to 10e-06 precision!
    expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.glm(object=m2))[, 1])))
    expect_equal(unname(m1$Stderror[[1]]), unname(t(coef(summary.glm(object=m2))[, 2])), tolerance=1e-06)
  })
})

test_that("Test fitAbn.mle() with Multinomial nodes and Gaussians", {
  suppressMessages({
    if(.Platform$OS.type == "unix") {
      capture.output({
        # make data set
        N <- 200
        prob <- c(1,2,3,4)
        coef <- c(1,2,3,4)
        Y <- rmultinom(N, 1, prob=prob)
        Ya <- c(colSums(Y*coef))
        X <- colSums(Y*coef) + rnorm(N,0,sd=1)
        mydf <- data.frame(a=X, b=as.factor(Ya))
        mydists <- list(a="gaussian", b="multinomial")
        mydag <- matrix(c(0,1,0,0), nrow = 2, byrow = TRUE,
                        dimnames = list(c("a", "b"), c("a", "b"))) # b~a

        ###
        # Multinomial parent with one gaussian child
        ###
        mycache.mle <- buildScoreCache(data.df = mydf,data.dists = mydists, method = "mle", max.parents = 1, dag.retained = mydag)
        mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
        myfit <- fitAbn(object = mydag.mp, method = "mle", centre = FALSE)
      },
      file = "/dev/null")
      # with original data
      m1 <- myfit
      m2 <- glm(a ~ -1+ b, data = mydf, family="gaussian")
      expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.glm(object=m2))[, 1])))
      expect_equal(unname(m1$Stderror[[1]]), unname(t(coef(summary.glm(object=m2))[, 2])), tolerance=1e-02)
      # with simulated data
      capture.output({
        out.sim <- invisible(
          simulateAbn(object = myfit, verbose = FALSE, debug = FALSE)
        )
        mycache.mle2 <- buildScoreCache(data.df = out.sim,
                                        data.dists = mydists,
                                        method = "mle",
                                        max.parents = 1,
                                        dag.retained = mydag)
        mydag.mp2 <- mostProbable(score.cache = mycache.mle2, verbose = FALSE)
        m1 <- fitAbn(object = mydag.mp2, method = "mle", centre=FALSE)
        m2 <- glm(a ~ -1+ b, data = out.sim, family="gaussian")
      },
      file = "/dev/null")
      expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.glm(object=m2))[, 1])))
      expect_equal(unname(m1$Stderror[[1]]), unname(t(coef(summary.glm(object=m2))[, 2])), tolerance=1e-3)

      ###
      # Multinomial response with two gaussian parents
      ###
      capture.output({
        # make data set
        mydf3 <- cbind(mydf, c=rnorm(N,0,sd=1))
        mydists3 <- list(a="gaussian", b="multinomial", c="gaussian")
        dag3 <- matrix(c(0,0,0,
                         1,0,1,
                         0,0,0), 3, 3, byrow = T, dimnames = list(names(mydists3), names(mydists3)))
        mycache.mle3 <- buildScoreCache(data.df = mydf3,data.dists = mydists3, method = "mle", max.parents = 2, dag.retained = dag3)
        mydag.mp3 <- mostProbable(score.cache = mycache.mle3, verbose = FALSE)
        myfit3 <- fitAbn(object = mydag.mp3, method = "mle", centre = FALSE)
        # with original data
        m1 <- myfit3
        m2 <- summary(nnet::multinom(b~a+c, data =mydf3, trace = FALSE))
      },
      file = "/dev/null")
      expect_equal(unname(as.vector(m1$coef[[2]])), unname(as.vector(coefficients(m2))))
      expect_equal(unname(as.vector(m1$Stderror[[2]])), unname(as.vector(m2$standard.errors)))
      # with simulated data
      capture.output({
        out.sim3 <- invisible(
          simulateAbn(object = myfit3, verbose = FALSE, debug = FALSE)
        )
        mycache.mle4 <- buildScoreCache(data.df = out.sim3,
                                        data.dists = mydists3,
                                        method = "mle",
                                        max.parents = 2,
                                        dag.retained = dag3)
        mydag.mp4 <- mostProbable(score.cache = mycache.mle4, verbose = FALSE)
        m1 <- fitAbn(object = mydag.mp4, method = "mle", centre=FALSE)
        m2 <- summary(nnet::multinom(b~a+c, data =out.sim3, trace = FALSE))
      },
      file = "/dev/null")
      expect_equal(unname(as.vector(m1$coef[[2]])), unname(as.vector(coefficients(m2))))
      expect_equal(unname(as.vector(m1$Stderror[[2]])), unname(as.vector(m2$standard.errors)))
    } else {
      skip("fitAbn.mle() is tested mainly on Unix-like systems")
    }
  })
})

test_that("Test fitAbn.mle() with Multinomial nodes and Binomials", {
  suppressMessages({
    if(.Platform$OS.type == "unix") {
      capture.output({
        # make data set
        N <- 200
        prob <- c(1,2,3)
        coef <- c(1,2,3)
        Y <- rmultinom(N, 1, prob=prob)
        Ya <- c(colSums(Y*coef))
        X <- rbinom(n = N, size = 1, prob = 0.2) # Bernoulli
        mydf <- data.frame(a=as.factor(X), b=as.factor(Ya))
        mydists <- list(a="binomial", b="multinomial")
        mydag <- matrix(c(0,1,0,0), nrow = 2, byrow = TRUE,
                        dimnames = list(c("a", "b"), c("a", "b"))) # b~a

        ###
        # Multinomial parent with one Binomial child
        ###
        mycache.mle <- buildScoreCache(data.df = mydf,data.dists = mydists, method = "mle", max.parents = 1, dag.retained = mydag)
        mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
        myfit <- fitAbn(object = mydag.mp, method = "mle", centre = FALSE)
        # with original data
        m1 <- myfit
        m2 <- glm(a ~ -1+ b, data = mydf, family="binomial")
      },
      file = "/dev/null")
      expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.glm(object=m2))[, 1])))
      expect_equal(unname(m1$Stderror[[1]]), unname(t(coef(summary.glm(object=m2))[, 2])), tolerance=1e-02)
      # with simulated data
      capture.output({
        out.sim <- invisible(
          simulateAbn(object = myfit, verbose = FALSE, debug = FALSE)
        )
        mycache.mle2 <- buildScoreCache(data.df = out.sim,
                                        data.dists = mydists,
                                        method = "mle",
                                        max.parents = 1,
                                        dag.retained = mydag)
        mydag.mp2 <- mostProbable(score.cache = mycache.mle2, verbose = FALSE)
        m1 <- fitAbn(object = mydag.mp2, method = "mle", centre=FALSE)
        m2 <- glm(a ~ -1+ b, data = out.sim, family="binomial")
      },
      file = "/dev/null")
      expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.glm(object=m2))[, 1])), tolerance=1e-6)
      expect_equal(unname(m1$Stderror[[1]]), unname(t(coef(summary.glm(object=m2))[, 2])), tolerance=1e-3)

      ###
      # Multinomial response with two Binomial parents
      ###
      # make data set
      capture.output({
        mydf3 <- cbind(mydf,
                       c=as.factor(rbinom(n = N, size = 1, prob = 0.8))) # Bernoulli
        mydists3 <- list(a="binomial", b="multinomial", c="binomial")
        dag3 <- matrix(c(0,0,0,
                         1,0,1,
                         0,0,0), 3, 3, byrow = T, dimnames = list(names(mydists3), names(mydists3)))
        mycache.mle3 <- buildScoreCache(data.df = mydf3,data.dists = mydists3, method = "mle", max.parents = 2,
                                        dag.retained = dag3,
                                        dag.banned = matrix(c(1,0,1,
                                                              0,0,0,
                                                              1,0,1), 3, 3, byrow = T, dimnames = list(names(mydists3), names(mydists3))))
        mydag.mp3 <- mostProbable(score.cache = mycache.mle3, verbose = FALSE)
        myfit3 <- fitAbn(object = mydag.mp3, method = "mle", centre = FALSE)
        # with original data
        m1 <- myfit3
        m2 <- summary(nnet::multinom(b~a+c, data =mydf3, trace = FALSE))
      },
      file = "/dev/null")
      expect_equal(unname(as.vector(m1$coef[[2]])), unname(as.vector(coefficients(m2))))
      expect_equal(unname(as.vector(m1$Stderror[[2]])), unname(as.vector(m2$standard.errors)))
      # with simulated data
      capture.output({
        out.sim3 <- invisible(
          simulateAbn(object = myfit3, verbose = FALSE, debug = FALSE)
        )
        mycache.mle4 <- buildScoreCache(data.df = out.sim3,
                                        data.dists = mydists3,
                                        method = "mle",
                                        max.parents = 2,
                                        dag.retained = dag3,
                                        dag.banned = matrix(c(1,0,1,
                                                              0,0,0,
                                                              1,0,1), 3, 3, byrow = T, dimnames = list(names(mydists3), names(mydists3))))
        mydag.mp4 <- mostProbable(score.cache = mycache.mle4, verbose = FALSE)
        m1 <- fitAbn(object = mydag.mp4, method = "mle", centre=FALSE)
        m2 <- summary(nnet::multinom(b~a+c, data =out.sim3, trace = FALSE))
      },
      file = "/dev/null")
      expect_equal(unname(as.vector(m1$coef[[2]])), unname(as.vector(coefficients(m2))))
      expect_equal(unname(as.vector(m1$Stderror[[2]])), unname(as.vector(m2$standard.errors)))
    } else {
      skip("fitAbn.mle() is tested mainly on Unix-like systems")
    }
  })
})

test_that("Test fitAbn.mle() with Multinomial nodes and Poissons", {
  suppressMessages({
    if(.Platform$OS.type == "unix") {
      capture.output({
        # make data set
        N <- 200
        prob <- c(1,2,3)
        coef <- c(1,2,3)
        Y <- rmultinom(N, 1, prob=prob)
        Ya <- c(colSums(Y*coef))
        X <- rpois(n = N, lambda = 0.2) # Poisson
        mydf <- data.frame(a=X, b=as.factor(Ya))
        mydists <- list(a="poisson", b="multinomial")
        mydag <- matrix(c(0,1,0,0), nrow = 2, byrow = TRUE,
                        dimnames = list(c("a", "b"), c("a", "b"))) # b~a

        ###
        # Multinomial parent with one Poisson child
        ###
        mycache.mle <- buildScoreCache(data.df = mydf,data.dists = mydists, method = "mle", max.parents = 1, dag.retained = mydag)
        mydag.mp <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
        suppressWarnings({
          # This is needed to avoid a warning about step halving not being able to find a better fit.
          myfit <- fitAbn(object = mydag.mp, method = "mle", centre = FALSE)
        })
        # with original data
        m1 <- myfit
        m2 <- glm(a ~ -1+ b, data = mydf, family="poisson")
      },
      file = "/dev/null")
      expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.glm(object=m2))[, 1])))
      expect_equal(unname(m1$Stderror[[1]]), unname(t(coef(summary.glm(object=m2))[, 2])), tolerance=1e-02)
      # with simulated data
      capture.output({
        out.sim <- invisible(
          simulateAbn(object = myfit, verbose = FALSE, debug = FALSE)
        )
        mycache.mle2 <- buildScoreCache(data.df = out.sim,
                                        data.dists = mydists,
                                        method = "mle",
                                        max.parents = 1,
                                        dag.retained = mydag)
        mydag.mp2 <- mostProbable(score.cache = mycache.mle2, verbose = FALSE)
        m1 <- fitAbn(object = mydag.mp2, method = "mle", centre=FALSE)
        m2 <- glm(a ~ -1+ b, data = out.sim, family="poisson")
      },
      file = "/dev/null")
      expect_equal(unname(m1$coef[[1]]), unname(t(coef(summary.glm(object=m2))[, 1])), tolerance=1e-6)
      expect_equal(unname(m1$Stderror[[1]]), unname(t(coef(summary.glm(object=m2))[, 2])), tolerance=1e-3)

      ###
      # Multinomial response with two Poisson parents
      ###
      # make data set
      capture.output({
        mydf3 <- cbind(mydf,
                       c=rpois(n = N, lambda = 0.8)) # Poisson
        mydists3 <- list(a="poisson", b="multinomial", c="poisson")
        dag3 <- matrix(c(0,0,0,
                         1,0,1,
                         0,0,0), 3, 3, byrow = T, dimnames = list(names(mydists3), names(mydists3)))
        mycache.mle3 <- buildScoreCache(data.df = mydf3,data.dists = mydists3, method = "mle", max.parents = 2,
                                        dag.retained = dag3,
                                        dag.banned = matrix(c(1,0,1,
                                                              0,0,0,
                                                              1,0,1), 3, 3, byrow = T, dimnames = list(names(mydists3), names(mydists3))))
        suppressWarnings({
          # This is needed to avoid a warning about step halving not being able to find a better fit.
          mydag.mp3 <- mostProbable(score.cache = mycache.mle3, verbose = FALSE)
          myfit3 <- fitAbn(object = mydag.mp3, method = "mle", centre = FALSE)
        })
        # with original data
        m1 <- myfit3
        m2 <- summary(nnet::multinom(b~a+c, data =mydf3, trace = FALSE))
      },
      file = "/dev/null")
      expect_equal(unname(as.vector(m1$coef[[2]])), unname(as.vector(coefficients(m2))))
      expect_equal(unname(as.vector(m1$Stderror[[2]])), unname(as.vector(m2$standard.errors)))
      # with simulated data
      capture.output({
        out.sim3 <- invisible(
          simulateAbn(object = myfit3, verbose = FALSE, debug = FALSE)
        )
        mycache.mle4 <- buildScoreCache(data.df = out.sim3,
                                        data.dists = mydists3,
                                        method = "mle",
                                        max.parents = 2,
                                        dag.retained = dag3,
                                        dag.banned = matrix(c(1,0,1,
                                                              0,0,0,
                                                              1,0,1), 3, 3, byrow = T, dimnames = list(names(mydists3), names(mydists3))))
        suppressWarnings({
          # This is needed to avoid a warning about step halving not being able to find a better fit.
          mydag.mp4 <- mostProbable(score.cache = mycache.mle4, verbose = FALSE)
          m1 <- fitAbn(object = mydag.mp4, method = "mle", centre=FALSE)
        })
        m2 <- summary(nnet::multinom(b~a+c, data =out.sim3, trace = FALSE))
      },
      file = "/dev/null")
      expect_equal(unname(as.vector(m1$coef[[2]])), unname(as.vector(coefficients(m2))))
      expect_equal(unname(as.vector(m1$Stderror[[2]])), unname(as.vector(m2$standard.errors)))
    } else {
      skip("fitAbn.mle() is tested mainly on Unix-like systems")
    }
  })
})

test_that("fitAbn's regressionLoop() works w/o group.var.", {
  # load("tests/testthat/testdata/fitAbn_regressionLoop_data.Rdata")
  load("testdata/fitAbn_regressionLoop_data.Rdata")

  # Add new control option to force glmmTMB for Poisson nodes
  control[["only_glmmTMB_poisson"]] <- FALSE

  verbose <- FALSE

  # Running on one child node (predictor)
  suppressWarnings({
    suppressMessages({
      if(.Platform$OS.type == "unix") {
        capture.output({
          expect_no_error({
            res <- regressionLoop(
              i = 1,
              dag = dag,
              data.df = data.df,
              data.df.multi = data.df.multi,
              data.dists = data.dists,
              group.var = group.var,
              grouped.vars = grouped.vars,
              group.ids = group.ids,
              control = control,
              nvars = nvars,
              nobs = nobs,
              dag.multi = dag.multi,
              verbose = verbose
            )
          })
          expect_equal(names(res), c("mliknode", "mlik", "aicnode", "aic", "bicnode", "bic", "mdlnode", "mdl", "sse", "mse", "df", "mu", "betas", "sigma", "sigma_alpha"))

          # Running on all child nodes (predictor)
          expect_no_error({
            res <- list()
            for (childno in 1:nrow(dag)){
              res[[childno]] <- regressionLoop(
                i = childno,
                dag = dag,
                data.df = data.df,
                data.df.multi = data.df.multi,
                data.dists = data.dists,
                group.var = group.var,
                grouped.vars = grouped.vars,
                group.ids = group.ids,
                control = control,
                nvars = nvars,
                nobs = nobs,
                dag.multi = dag.multi,
                verbose = verbose
              )
            }
          })
          expect_equal(length(res), nrow(dag))

          # Running on all child nodes (predictor) with foreach on single core
          expect_no_error({
            res2 <- foreach(childno = 1:nrow(dag)) %do% {
              regressionLoop(
                i = childno,
                dag = dag,
                data.df = data.df,
                data.df.multi = data.df.multi,
                data.dists = data.dists,
                group.var = group.var,
                grouped.vars = grouped.vars,
                group.ids = group.ids,
                control = control,
                nvars = nvars,
                nobs = nobs,
                dag.multi = dag.multi,
                verbose = verbose
              )
            }
          })
          expect_equal(res, res2)

          # Running on all child nodes (predictor) with foreach in multiple cores
          skip_on_cran() # workaround to not overconsume threads on CRAN. This is related to an issue reported for lme4 (https://github.com/lme4/lme4/issues/627)
          expect_no_error({
            ncores <- 2
            cl <- makeCluster(ncores)
            registerDoParallel(cl)

            res3 <- foreach(childno = 1:nrow(dag),
                            .export = 'regressionLoop') %dopar% {
                              regressionLoop(
                                i = childno,
                                dag = dag,
                                data.df = data.df,
                                data.df.multi = data.df.multi,
                                data.dists = data.dists,
                                group.var = group.var,
                                grouped.vars = grouped.vars,
                                group.ids = group.ids,
                                control = control,
                                nvars = nvars,
                                nobs = nobs,
                                dag.multi = dag.multi,
                                verbose = verbose
                              )
                            }
          })
          expect_equal(res, res3, tolerance = 10e-8)
          expect_equal(res2, res3, tolerance = 10e-8)
        },
        file = "/dev/null")
      } else {
        skip("fitAbn.mle() is tested mainly on Unix-like systems")
      }
    })
  })
})


test_that("fitAbn's regressionLoop() works w/ group.var.", {
  # load("tests/testthat/testdata/fitAbn_regressionLoop_group_data.Rdata")
  load("testdata/fitAbn_regressionLoop_group_data.Rdata")

  # Add new control option to force glmmTMB for Poisson nodes
  control[["only_glmmTMB_poisson"]] <- FALSE

  verbose <- FALSE

  # Running on one child node (predictor)
  suppressWarnings({
    suppressMessages({
      if(.Platform$OS.type == "unix") {
        capture.output({
          expect_no_error({
            res <- regressionLoop(
              i = 1,
              dag = dag,
              data.df = data.df,
              data.df.multi = data.df.multi,
              data.dists = data.dists,
              group.var = group.var,
              grouped.vars = grouped.vars,
              group.ids = group.ids,
              control = control,
              nvars = nvars,
              nobs = nobs,
              dag.multi = dag.multi,
              verbose = verbose
            )
          })
          expect_equal(names(res), c("mliknode", "mlik", "aicnode", "aic", "bicnode", "bic", "mdlnode", "mdl", "sse", "mse", "df", "mu", "betas", "sigma", "sigma_alpha"))

          # Running on all child nodes (predictor)
          expect_no_error({
            res <- list()
            for (childno in 1:nrow(dag)){
              res[[childno]] <- regressionLoop(
                i = childno,
                dag = dag,
                data.df = data.df,
                data.df.multi = data.df.multi,
                data.dists = data.dists,
                group.var = group.var,
                grouped.vars = grouped.vars,
                group.ids = group.ids,
                control = control,
                nvars = nvars,
                nobs = nobs,
                dag.multi = dag.multi,
                verbose = verbose
              )
            }
          })
          expect_equal(length(res), nrow(dag))

          # Running on all child nodes (predictor) with foreach on single core
          expect_no_error({
            res2 <- foreach(childno = 1:nrow(dag)) %do% {
              regressionLoop(
                i = childno,
                dag = dag,
                data.df = data.df,
                data.df.multi = data.df.multi,
                data.dists = data.dists,
                group.var = group.var,
                grouped.vars = grouped.vars,
                group.ids = group.ids,
                control = control,
                nvars = nvars,
                nobs = nobs,
                dag.multi = dag.multi,
                verbose = verbose
              )
            }
          })
          expect_equal(res, res2)

          # Running on all child nodes (predictor) with foreach in multiple cores
          skip_on_cran() # workaround to not overconsume threads on CRAN. This is related to an issue reported for lme4 (https://github.com/lme4/lme4/issues/627)
          expect_no_error({
            ncores <- 2
            cl <- makeCluster(ncores)
            registerDoParallel(cl)

            res3 <- foreach(childno = 1:nrow(dag),
                            .export = 'regressionLoop') %dopar% {
                              regressionLoop(
                                i = childno,
                                dag = dag,
                                data.df = data.df,
                                data.df.multi = data.df.multi,
                                data.dists = data.dists,
                                group.var = group.var,
                                grouped.vars = grouped.vars,
                                group.ids = group.ids,
                                control = control,
                                nvars = nvars,
                                nobs = nobs,
                                dag.multi = dag.multi,
                                verbose = verbose
                              )
                            }
          })
          expect_equal(res, res3, tolerance = 10e-8)
          expect_equal(res2, res3, tolerance = 10e-8)
        },
        file = "/dev/null")
      } else {
        skip("fitAbn.mle() is tested mainly on Unix-like systems")
      }
    })
  })
})

test_that("regressionLoop() prints local model when verbose.", {
  if(.Platform$OS.type == "unix") {
    suppressMessages({
      # Suppress messages that are not related to the test but are printed when verbose
      capture.output({
        # load("tests/testthat/testdata/fitAbn_regressionLoop_group_data.Rdata")
        load("testdata/fitAbn_regressionLoop_group_data.Rdata")

        # Add new control option to force glmmTMB for Poisson nodes
        control[["only_glmmTMB_poisson"]] <- FALSE

        verbose <- TRUE

        # with group.var
        suppressWarnings({
          expect_message({
            res <- regressionLoop(
              i = 3,
              dag = dag,
              data.df = data.df,
              data.df.multi = data.df.multi,
              data.dists = data.dists,
              group.var = group.var,
              grouped.vars = grouped.vars,
              group.ids = group.ids,
              control = control,
              nvars = nvars,
              nobs = nobs,
              dag.multi = dag.multi,
              verbose = verbose
            )
          },
          regex = "using glmer with model")
        })

        # without group.var
        suppressWarnings({
          expect_message({
            res <- regressionLoop(
              i = 3,
              dag = dag,
              data.df = data.df,
              data.df.multi = data.df.multi,
              data.dists = data.dists,
              group.var = NULL, # no group.var
              control = control,
              nvars = nvars,
              nobs = nobs,
              dag.multi = dag.multi,
              verbose = verbose
            )
          },
          regex = "regressing GroupSize on Outdoor")
        })
      },
      file = "/dev/null")
    })
  } else {
    skip("fitAbn.mle() is tested mainly on Unix-like systems")
  }
})


testthat::test_that("Poisson nodes step into calling glmmTMB for the fitabn function.", {
  if(.Platform$OS.type == "unix") {
    suppressMessages({
      # Suppress messages that are not related to the test but are printed when verbose
      suppressWarnings({
        capture.output({
          # load("tests/testthat/testdata/fitAbn_regressionLoop_group_data.Rdata")
          load("testdata/fitAbn_regressionLoop_group_data.Rdata")

          verbose <- TRUE
          # with group.var
          testthat::expect_message(
            regressionLoop(i = 3,
                           dag = dag,
                           data.df = data.df,
                           data.df.multi = data.df.multi,
                           data.dists = data.dists,
                           group.var = group.var,
                           grouped.vars = grouped.vars,
                           group.ids = group.ids,
                           nvars = nvars,
                           nobs = nobs,
                           dag.multi = dag.multi,
                           verbose = verbose,
                           control = build.control(method = "mle", only_glmmTMB_poisson=TRUE)),
            regexp = "trying glmmTMB with model"
          )
        },

        file = "/dev/null")
      })
    })
  } else {
    testthat::skip("`forLoopContent()` is tested mainly on Unix-like systems.")
  }
})

# --- invariants the bayesian-network JSON format relies on ------------------

test_that("binomial nodes model the second factor level with a logit link", {
  spec <- jfx_spec("ex1_mle")
  fit <- jfx_fit("ex1_mle")
  ref <- stats::glm(b3 ~ b1 + g1 + b2, family = stats::binomial(), data = spec$data)
  expect_equal(unname(fit$coef$b3[1, ]), unname(stats::coef(ref)), tolerance = 1e-4)
})

test_that("poisson nodes use a log link", {
  spec <- jfx_spec("ex1_mle")
  fit <- jfx_fit("ex1_mle")
  ref <- stats::glm(p2 ~ b1 + p1, family = stats::poisson(), data = spec$data)
  expect_equal(unname(fit$coef$p2[1, ]), unname(stats::coef(ref)), tolerance = 1e-4)
})

test_that("ungrouped gaussian mse is the unbiased residual variance", {
  spec <- jfx_spec("ex1_mle")
  fit <- jfx_fit("ex1_mle")
  ref <- stats::lm(g2 ~ p1 + g1 + b2, data = spec$data)
  expect_equal(unname(fit$mse[["g2"]]), stats::sigma(ref)^2, tolerance = 1e-6)
})

test_that("a multinomial parent is one-hot encoded over all levels without intercept", {
  spec <- jfx_spec("fcv_mle")
  fit <- jfx_fit("fcv_mle")
  n_levels <- nlevels(spec$data$Sex)
  expect_false(any(grepl("intercept", colnames(fit$coef$Outdoor))))
  expect_equal(ncol(fit$coef$Outdoor), n_levels)
  ref <- stats::glm(Outdoor ~ -1 + Sex, family = stats::binomial(), data = spec$data)
  expect_equal(unname(fit$coef$Outdoor[1, ]), unname(stats::coef(ref)), tolerance = 1e-4)
})

test_that("a multinomial child uses the first level as baseline", {
  spec <- jfx_spec("g2b2c_soft_mle")
  fit <- jfx_fit("g2b2c_soft_mle")
  ref <- nnet::multinom(C ~ B1 + B2, data = spec$data, trace = FALSE)
  non_baseline <- levels(spec$data$C)[-1]
  expect_equal(rownames(stats::coef(ref)), non_baseline)
  intercepts <- fit$coef$C[1, paste0("C|intercept.", non_baseline)]
  expect_equal(unname(intercepts), unname(stats::coef(ref)[, "(Intercept)"]),
               tolerance = 1e-3)
})

test_that("grouped gaussian sigma and sigma_alpha are standard deviations", {
  spec <- jfx_spec("adg_mle_grouped")
  fit <- jfx_fit("adg_mle_grouped")
  ref <- lme4::lmer(adg ~ age + wormCount + (1 | farm), data = spec$data)
  vc <- as.data.frame(lme4::VarCorr(ref))
  expect_equal(unname(fit$sigma_alpha$adg)^2, vc$vcov[vc$grp == "farm"],
               tolerance = 1e-4)
  expect_equal(unname(fit$sigma$adg)^2, vc$vcov[vc$grp == "Residual"],
               tolerance = 1e-4)
})

test_that("grouped multinomial sigma_alpha is a covariance matrix", {
  jfx_skip_unavailable("g2pbcgrp_mle_grouped")
  fit <- jfx_fit("g2pbcgrp_mle_grouped")
  spec <- jfx_spec("g2pbcgrp_mle_grouped")
  k <- nlevels(spec$data$C) - 1
  expect_true(is.matrix(fit$sigma_alpha$C))
  expect_equal(dim(fit$sigma_alpha$C), c(k, k))
  expect_equal(fit$sigma_alpha$C, t(fit$sigma_alpha$C))
})

test_that("gaussian mse uses the correct df with a multinomial parent", {
  spec <- jfx_spec("g2b2c_mle")
  fit <- jfx_fit("g2b2c_mle")
  ref <- stats::lm(G2 ~ G1 + C, data = spec$data)
  expect_equal(unname(fit$df[["G2"]]), stats::df.residual(ref))
  expect_equal(unname(fit$mse[["G2"]]), stats::sigma(ref)^2, tolerance = 1e-6)
})

test_that("multinomial child coefficient names follow the value order", {
  spec <- jfx_spec("g2b2c_soft_mle")
  fit <- jfx_fit("g2b2c_soft_mle")
  expect_equal(levels(spec$data$C), c("a", "b", "c"))
  expect_equal(colnames(fit$coef$C),
               c("C|intercept.b", "C|intercept.c", "B1b", "B1c", "B2b", "B2c"))
  expect_equal(colnames(fit$Stderror$C), colnames(fit$coef$C))

  ref <- nnet::multinom(C ~ B1 + B2, data = spec$data, trace = FALSE)
  ref_coef <- stats::coef(ref)
  for (state in c("b", "c")) {
    expect_equal(unname(fit$coef$C[1, paste0("B1", state)]), unname(ref_coef[state, "B11"]),
                 tolerance = 1e-3)
    expect_equal(unname(fit$coef$C[1, paste0("B2", state)]), unname(ref_coef[state, "B21"]),
                 tolerance = 1e-3)
  }
})

test_that("separation fallback names coefficients correctly", {
  set.seed(1)
  n <- 200
  x <- factor(rbinom(n, 1, 0.5))
  y <- factor(ifelse(as.numeric(x) == 2, sample(c("a", "b"), n, TRUE), "c"))
  dat <- data.frame(x = x, y = y)
  dag <- matrix(0, 2, 2, dimnames = list(c("x", "y"), c("x", "y")))
  dag["y", "x"] <- 1
  expect_warning(
    fit <- suppressMessages(fitAbn(dag = dag, data.df = dat,
                                   data.dists = list(x = "binomial", y = "multinomial"),
                                   method = "mle")),
    "Separation"
  )
  expect_equal(colnames(fit$coef$y), paste0("y|intercept.", c("b", "c")))
  expect_equal(ncol(fit$coef$y), 2L)
})
