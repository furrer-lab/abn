test_that("General behaviour of buildScoreCache.mle()", {
  if(.Platform$OS.type == "unix") {
    capture.output({
      df <- airquality[complete.cases(airquality), ]

      # distribution (gaussian)
      dist <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
      names(dist) <- colnames(df)
    },
    file = "/dev/null")
    ## which.nodes
    expect_warning({
      expect_no_error(buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=6, which.nodes=1:6))
    })

    ## max.parents
    expect_error(buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=6, which.nodes=1:5))
    expect_warning({
      expect_no_error(buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=5, which.nodes=1:5))
    })
    expect_error(buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=list(Ozone=NULL, Solar.R=4, Wind=4, Temp=4, Month=4), which.node=1:5))
    expect_equal(buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=4, which.nodes=1:5),
                 buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=list(Ozone=4, Solar.R=4, Wind=4, Temp=4, Month=4), which.node=1:5))
    expect_equal(buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=3, which.nodes=1:5),
                 buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=list(Ozone=3, Solar.R=3, Wind=3, Temp=3, Month=3), which.node=1:5))
    skip_if(
      expect_error(
        buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=list(Ozone=2, Solar.R=4, Wind=4, Temp=4, Month=4), which.node=1:5),
        regexp = "ISSUE"),
      message = "`max.parents` with node specific values that are not all the same, is not yet implemented.")
  }
})

test_that("buildScoreCache.mle() with Gaussian nodes", {
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
      expect_no_error({
        mycache <- invisible(buildScoreCache(data.df = mydf,
                                             data.dists = mydists,
                                             method = "mle",
                                             max.parents = 2,
                                             centre = FALSE))
      })
    },
    file = "/dev/null")

    # mLik
    expect_equal(mycache$mlik[1], as.numeric(logLik(lm(formula=mydf$a ~ 1))))
    expect_equal(mycache$mlik[2], as.numeric(logLik(lm(formula=mydf$a ~ 1 + mydf$b))))
    expect_equal(mycache$mlik[5], as.numeric(logLik(lm(formula=mydf$b ~ 1))))

    # AIC
    expect_equal(mycache$aic[1], as.numeric(AIC(lm(formula=mydf$a ~ 1))))
    expect_equal(mycache$aic[2], as.numeric(AIC(lm(formula=mydf$a ~ 1 + mydf$b))))
    expect_equal(mycache$aic[5], as.numeric(AIC(lm(formula=mydf$b ~ 1))))

    # BIC
    expect_equal(mycache$bic[1], as.numeric(BIC(lm(formula=mydf$a ~ 1))))
    expect_equal(mycache$bic[2], as.numeric(BIC(lm(formula=mydf$a ~ 1 + mydf$b))))
    expect_equal(mycache$bic[5], as.numeric(BIC(lm(formula=mydf$b ~ 1))))

    ## Gaussian: correlation coefficients
    capture.output({
      expect_no_error({
        mycache <- invisible(buildScoreCache(data.df=mydf, data.dists=mydists, method = "mle", max.parents=2, centre=TRUE))
      })
    },
    file = "/dev/null")
    expect_equal(mycache$mlik[1], as.numeric(logLik(lm(formula=(mydf$a - mean(mydf$a))/sd(mydf$a) ~ 1))))
  }
})

test_that("buildScoreCache.mle() with Binomial nodes", {
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
    },
    file = "/dev/null")
    # mLik
    expect_equal(mycache.mle$mlik[1], as.numeric(logLik(glm(formula=mydf$a ~ 1, family=binomial))), tolerance = 0.00005)
    expect_equal(mycache.mle$mlik[2], as.numeric(logLik(glm(formula=mydf$a ~ 1 + mydf$b, family=binomial))), tolerance = 0.00005)
    expect_equal(mycache.mle$mlik[3], as.numeric(logLik(glm(formula=mydf$b ~ 1, family=binomial))), tolerance = 0.00005)

    # AIC
    expect_equal(mycache.mle$aic[1], as.numeric(AIC(glm(formula=mydf$a ~ 1, family=binomial))), tolerance = 0.00005)
    expect_equal(mycache.mle$aic[2], as.numeric(AIC(glm(formula=mydf$a ~ 1 + mydf$b, family=binomial))), tolerance = 0.00005)
    expect_equal(mycache.mle$aic[3], as.numeric(AIC(glm(formula=mydf$b ~ 1, family=binomial))), tolerance = 0.00005)

    # BIC
    expect_equal(mycache.mle$bic[1], as.numeric(BIC(glm(formula=mydf$a ~ 1, family=binomial))), tolerance = 0.00005)
    expect_equal(mycache.mle$bic[2], as.numeric(BIC(glm(formula=mydf$a ~ 1 + mydf$b, family=binomial))), tolerance = 0.00005)
    expect_equal(mycache.mle$bic[3], as.numeric(BIC(glm(formula=mydf$b ~ 1, family=binomial))), tolerance = 0.00005)
  }
})

test_that("buildScoreCache.mle() with Poisson nodes", {
  if(.Platform$OS.type == "unix") {
    capture.output({
      # Poisson
      mydists <- list(a = "gaussian",
                      b = "poisson")
      a <- cbind(rnorm(1000))
      betasTrue<-c(2)
      etaTrue<-a%*%betasTrue
      b <- rpois(1000,exp(etaTrue))
      retaindag <- matrix(c(0,0,1,0), nrow = 2, byrow = TRUE, dimnames = list(c("a", "b"), c("a", "b")))
      banneddag <- matrix(c(1,1,0,1), nrow = 2, byrow = TRUE, dimnames = list(c("a", "b"), c("a", "b")))
      mydf <- data.frame("a" = a,
                         "b" = b)
      mycache.mle <- buildScoreCache(data.df = mydf,
                                     data.dists = mydists,
                                     method = "mle",
                                     max.parents = 1,
                                     dag.retained = retaindag,
                                     dag.banned = banneddag,
                                     verbose = FALSE)
    },
    file = "/dev/null")

    modglm <- glm(formula=b ~ a, data = mydf, family=poisson)
    # mLik
    expect_equal(mycache.mle$mlik[1], as.numeric(logLik(modglm)), tolerance = 10e0)

    # AIC
    expect_equal(mycache.mle$aic[1], as.numeric(AIC(modglm)), tolerance = 10e0)

    # BIC
    expect_equal(mycache.mle$bic[1], as.numeric(BIC(modglm)), tolerance = 10e0)
  }
})

test_that("buildScoreCache.mle() with Multinomial nodes", {
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
    },
    file = "/dev/null")
    # mLik
    # expect_equal(mycache.mle$mlik[1], as.numeric(logLik(glm(formula=a~b, data = mydf, family=gaussian))))
    # expect_equal(mycache.mle$mlik[3], as.numeric(logLik(glm(formula=a~1, data = mydf, family=gaussian))))

    # AIC
    # expect_equal(mycache.mle$aic[1], as.numeric(AIC(glm(formula=a~b, data = mydf, family=gaussian))))
    # expect_equal(mycache.mle$aic[3], as.numeric(AIC(glm(formula=a~1, data = mydf, family=gaussian))))

    # BIC
    # expect_equal(mycache.mle$bic[1], as.numeric(BIC(glm(formula=a~b, data = mydf, family=gaussian))))
    # expect_equal(mycache.mle$bic[3], as.numeric(BIC(glm(formula=a~1, data = mydf, family=gaussian))))


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
    },
    file = "/dev/null")

    # mLik
    expect_equal(mycache.mle3$mlik[5], as.numeric(logLik(nnet::multinom(b~a+c, data =mydf3, trace = FALSE))))

    # AIC
    # expect_equal(mycache.mle3$mlik[5], as.numeric(AIC(nnet::multinom(b~a+c, data =mydf3, trace = FALSE))))

    # BIC
    # expect_equal(mycache.mle3$mlik[5], as.numeric(BIC(nnet::multinom(b~a+c, data =mydf3, trace = FALSE))))
  }
})

test_that("buildScoreCache.mle() corresponds with simulation results", {
  if(.Platform$OS.type == "unix") {
    capture.output({
      ## simulation data
      n <- 1000
      x1 <- rnorm(n)
      x2 <- rbinom(n, 1, 0.5)
      b0 <- 1
      b1 <- 1.5
      b2 <- 2
      if(requireNamespace("boot", quietly=TRUE)){
        y <- rbinom(n, 1, boot::inv.logit(b0 + b1 * x1 + b2 * x2))
      }

      y <- ifelse(x2 == 1, 1, y)

      dist <- list(a="binomial", b="gaussian", c="binomial")
      dta <- data.frame(y, x1, x2)
      names(dta) <- names(dist)
    },
    file = "/dev/null")

    if(requireNamespace("brglm", quietly=TRUE)){
      # mycache <- invisible(buildScoreCache.mle(data.df=dta, data.dists=dist, max.parents=2, centre=FALSE, dry.run=TRUE))  ##
      # save(mycache, file='tests/testdata/mycache.Rdata')
      load(file='testdata/mycache.Rdata')

      expect_equal(mycache$mlik[1],  suppressWarnings(as.numeric(logLik(brglm::brglm(formula=dta$a ~ 1)))), tolerance=0.01)
      expect_equal(mycache$mlik[2],  as.numeric(logLik(suppressWarnings(brglm::brglm(formula=dta$a ~ dta$b)))), tolerance=0.01)
      expect_equal(mycache$mlik[3],  as.numeric(logLik(suppressWarnings(brglm::brglm(formula=dta$a ~ dta$c)))), tolerance=0.01)
      expect_equal(mycache$mlik[4],  as.numeric(logLik(suppressWarnings(brglm::brglm(formula=dta$a ~ dta$b + dta$c)))), tolerance=0.01)
      expect_equal(mycache$mlik[12], as.numeric(logLik(suppressWarnings(brglm::brglm(formula=dta$c ~ dta$b + dta$a)))), tolerance=0.01)
    }  else {
      cat('Package brglm not available, number of passed tests might be different')
    }
  }
})

test_that("forLoopContent() works as expected.", {
  # load(file = "tests/testthat/testdata/forLoopContent_data.Rdata")
  load(file = "testdata/forLoopContent_data.Rdata")

  # with group.var and binomial, gaussian, poisson, multinomial distributions, max.parents=1
  suppressMessages({
    suppressWarnings({
      if(.Platform$OS.type == "unix") {
        capture.output({
          expect_no_error(
            forLoopContent(row.num = 1,
                           mycache = mycache,
                           data.dists = data.dists,
                           data.df.multi = data.df.multi,
                           adj.vars = adj.vars,
                           data.df = data.df,
                           data.df.lvl = data.df.lvl,
                           group.var = group.var,
                           group.ids = group.ids,
                           control = control,
                           n = nvars,
                           verbose = verbose))

          expect_equal(
            forLoopContent(row.num = 1,
                           mycache = mycache,
                           data.dists = data.dists,
                           data.df.multi = data.df.multi,
                           adj.vars = adj.vars,
                           data.df = data.df,
                           data.df.lvl = data.df.lvl,
                           group.var = group.var,
                           group.ids = group.ids,
                           control = control,
                           n = nvars,
                           verbose = verbose),
            # c(-204.4035, 410.8071, 414.5109, 414.5109), tolerance = 0.01)
            c(-178.9, 361.8, 369.3, 370.9), tolerance = 0.01)
        },
        file = "/dev/null")
      }
    })
  })
})
