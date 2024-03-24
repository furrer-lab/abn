test_that("fit.control() works fine", {
  # basics
  expect_no_error(
    fctrl <- fit.control()
  )
  expect_equal(class(fctrl), class(list()))
  expect_no_error(
    fctrl.bayes <- fit.control(method = "bayes")
  )
  expect_equal(fctrl, fctrl.bayes)

  expect_no_error(
    fctrl.mle <- fit.control(method = "mle")
  )
  expect_equal(class(fctrl), class(list()))
  expect_equal(class(fctrl.bayes), class(list()))
  expect_equal(class(fctrl.mle), class(list()))

  # Arguments
  expect_no_error(
    fit.control(method = "bayes", epsilon = 1e-06) # no error if argument and method don't correspond
  )
})

test_that("fitAbn() wrapper of 'mle' and 'bayes' works", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  # Gaussian
  df <- airquality[complete.cases(airquality), ]

  dist <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
  names(dist) <- colnames(df)

  d <- matrix(data=0, nrow=6, ncol=6)
  d[1, ] <- c(0, 1, 1, 1, 1, 1)
  colnames(d) <- rownames(d) <- names(dist)

  ## test wrapper of method "mle" and "bayes"
  expect_no_error({
    m.0.mle <- abn:::fitAbn.mle(dag=d, data.df=df, data.dists=dist)
  })
  expect_no_error({
    m.0.mle.1 <- fitAbn(dag=d, data.df=df, data.dists=dist, method="mle")
  })
  expect_equal(m.0.mle, unclass(m.0.mle.1))
  expect_s3_class(m.0.mle.1, class = "abnFit")


  expect_no_error({
    m.0.bayes.1 <- fitAbn(dag=d, data.df=df, data.dists=dist, method="bayes")
  })
  expect_s3_class(m.0.bayes.1, class = "abnFit")
})

test_that("fitAbn() works with DAG as formula statement", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  # using 'ex3.dag.data'
  mydists <- list(b1="binomial",  b2="binomial")

  ## test formula statement
  expect_silent(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, method = "bayes"))
  expect_silent(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, method = "mle"))

  ## test formula statement with correlation
  if(!testthat:::on_cran()) {
    if(requireNamespace("INLA", quietly = TRUE)){
      expect_silent(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2,14)], data.dists=mydists, method = "bayes", group.var="group", cor.vars=c("b1","b2")))
    }
  }
  expect_silent(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2,14)], data.dists=mydists, method = "mle", group.var="group", cor.vars=c("b1","b2")))
})

test_that("fitAbn() catches wrong arguments.", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  # using 'ex3.dag.data'
  mydists <- list(b1="binomial",  b2="binomial")

  ## test centre argument
  expect_no_error(
    fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, centre = TRUE))
  expect_no_error(
    fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, centre = FALSE))
  expect_error(
    fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, centre = NULL),
    regexp = "'centre' should be either TRUE or FALSE")

  # TODO: add tests for other arguments here.
  # Test for different method arguments e.g. fitAbn(method="mle", control=fit.control(method="bayes")). -> Consider to remove this from fit.control...
})

test_that("fitAbn() runs in parallel", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  # Prepare some data
  df <- airquality[complete.cases(airquality), ]
  dist <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
  names(dist) <- colnames(df)
  d <- matrix(data=0, nrow=6, ncol=6)
  d[1, ] <- c(0, 1, 1, 1, 1, 1)
  colnames(d) <- rownames(d) <- names(dist)

  # test method="bayes"
  expect_no_error({
    fitAbn(dag=d, data.df=df, data.dists=dist, method="bayes", control = list(ncores = 1))
  })
  expect_no_error({
    fitAbn(dag=d, data.df=df, data.dists=dist, method="bayes", control = list(ncores = 0))
  })
  expect_error({
    fitAbn(dag=d, data.df=df, data.dists=dist, method="bayes", control = list(ncores = -2))
  }, regexp = "Argument 'ncores'")
  skip_on_cran() # workaround to not overconsume threads on CRAN. This is related to an issue reported for lme4
  expect_no_error({
    fitAbn(dag=d, data.df=df, data.dists=dist, method="bayes", control = list(ncores = 2))
  })

  # test method="mle"
  expect_no_error({
    fitAbn(dag=d, data.df=df, data.dists=dist, method="mle", control = list(ncores = 1))
  })
  expect_no_error({
    fitAbn(dag=d, data.df=df, data.dists=dist, method="mle", control = list(ncores = 0))
  })
  expect_error({
    fitAbn(dag=d, data.df=df, data.dists=dist, method="mle", control = list(ncores = -2))
  }, regexp = "Argument 'ncores'")
  skip_on_cran() # workaround to not overconsume threads on CRAN. This is related to an issue reported for lme4 (https://github.com/lme4/lme4/issues/627)
  expect_no_error({
    fitAbn(dag=d, data.df=df, data.dists=dist, method="mle", control = list(ncores = 2))
  })
})

test_that("fitAbn() works with control arguments", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  mydists <- list(b1="binomial",  b2="binomial")

  # Supplying control arguments as direct arguments is deprecated. They should be passed in a list of 'control=list(max.mode.error=0)'
  expect_warning(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, max.mode.error=0))
  expect_no_error(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, method = "bayes", control=list(max.mode.error=0)))
  expect_warning(fitAbn(dag=~b1|b2, data.df=ex3.dag.data[,c(1,2)], data.dists=mydists, method = "mle", control=list(max.mode.error=0)))


  mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")];## take a subset of cols

  ## setup distribution list for each node
  mydists <- list(b1="binomial",
                  b2="binomial",
                  b3="binomial",
                  g1="gaussian",
                  b4="binomial",
                  p2="poisson",
                  p4="poisson"
  );
  ## null model - all independent variables
  mydag.empty <- matrix(data=c(
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0  #
  ), byrow=TRUE,ncol=7);
  colnames(mydag.empty) <- rownames(mydag.empty) <- names(mydat);
  ## now repeat but include some dependencies
  mydag <- mydag.empty;
  mydag["b1","b2"] <- 1; # b1 <- b2
  mydag["b2","p4"] <- 1; # b2 <- p4
  mydag["b2","g1"] <- 1; # b2 <- g1
  mydag["g1","p2"] <- 1; # g1 <- p2
  mydag["b3","g1"] <- 1; # b3 <- g1
  mydag["b4","b1"] <- 1; # b4 <- b1
  mydag["p4","g1"] <- 1; # p4 <- g1

  expect_no_error(fitAbn(dag=mydag,data.df=mydat,data.dists=mydists,centre=TRUE, compute.fixed=TRUE, control=list(n.grid=100)))
  expect_no_error(fitAbn(dag=mydag,data.df=mydat,data.dists=mydists, method = "bayes", centre=TRUE, compute.fixed=TRUE, control=list(n.grid=100)))

  expect_error(fitAbn(dag=mydag,data.df=mydat,data.dists=mydists, method = "mle", centre=TRUE, compute.fixed=TRUE, control=list(n.grid=100)))
  expect_warning(fitAbn(dag=mydag,data.df=mydat,data.dists=mydists, method = "mle", centre=TRUE, compute.fixed=FALSE, control=list(n.grid=100)))
  expect_warning(fitAbn(dag=mydag,data.df=mydat,data.dists=mydists, method = "mle", centre=TRUE, control=list(n.grid=100)))
  expect_no_error(fitAbn(dag=mydag,data.df=mydat,data.dists=mydists, method = "mle", centre=TRUE))

  # TODO: make this tests more explicit.
})

test_that("fitAbn() is backward compatible", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  ## use built-in simulated data set

  mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")];## take a subset of cols

  ## setup distribution list for each node
  mydists <- list(b1="binomial",
                  b2="binomial",
                  b3="binomial",
                  g1="gaussian",
                  b4="binomial",
                  p2="poisson",
                  p4="poisson"
  );
  ## null model - all independent variables
  mydag.empty <- matrix(data=c(
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0  #
  ), byrow=TRUE,ncol=7);
  colnames(mydag.empty) <- rownames(mydag.empty) <- names(mydat);

  ## now fit the model to calculate its goodness of fit
  expect_silent({
    myres.c <- fitAbn(dag=mydag.empty,data.df=mydat,data.dists=mydists,centre=TRUE, compute.fixed=FALSE);
  })
  expect_s3_class(myres.c, class = "abnFit")

  ## now repeat but include some dependencies
  mydag <- mydag.empty;
  mydag["b1","b2"] <- 1; # b1 <- b2
  mydag["b2","p4"] <- 1; # b2 <- p4
  mydag["b2","g1"] <- 1; # b2 <- g1
  mydag["g1","p2"] <- 1; # g1 <- p2
  mydag["b3","g1"] <- 1; # b3 <- g1
  mydag["b4","b1"] <- 1; # b4 <- b1
  mydag["p4","g1"] <- 1; # p4 <- g1

  expect_silent({
    myres.c <- fitAbn(dag=mydag,data.df=mydat,data.dists=mydists,centre=TRUE, compute.fixed=FALSE);
  })

  ## now also plot the graph via graphviz
  expect_silent({
    myres.c <- fitAbn(dag=mydag,data.df=mydat,data.dists=mydists,centre=TRUE, create.graph=TRUE,compute.fixed=FALSE);
  }) # create.graph argument is deprecated and should be ignored.
  plt <- plot(myres.c)
  plt.adjMat <- plt@adjMat
  mydag.adjMat <- mydag
  rownames(mydag.adjMat) <- NULL
  expect_equal(plt.adjMat, mydag.adjMat)


  ## a simple plot of some posterior densities
  ## the algorithm which chooses density points is very simple any may be
  ## rather sparse so also recompute the density over an equally spaced
  ## grid of 100 points between the two end points
  ## which had at f=min.pdf
  expect_silent({
    myres.c <- fitAbn(dag=mydag,data.df=mydat,data.dists=mydists,centre=TRUE, compute.fixed=TRUE, control=list(n.grid=100))
  })
  expect_false(any(myres.c$used.INLA))

  ## repeat but use INLA for the numerics using max.mode.error=100
  ## as using internal code is the default here rather than INLA
  if(!testthat:::on_cran()) {
    if(requireNamespace("INLA", quietly = TRUE)){
      expect_silent({
        myres.inla <- fitAbn(dag=mydag,data.df=mydat,data.dists=mydists,centre=TRUE, compute.fixed=TRUE, control=list(max.mode.error=100))
      })
      expect_true(any(myres.inla$used.INLA))
    }
  }

  ### A very simple mixed model example using built-in data
  ## specify DAG - only two variables using subset of variables from ex3.dag.data
  ## both variables are assumed to need (separate) adjustment for the group variable
  ## i.e. a binomial glmm at each node

  ## model where b1 <- b2
  mydag <- matrix(data=c(
    0,1, # b1
    0,0  # b2
  ), byrow=TRUE,ncol=2);

  colnames(mydag) <- rownames(mydag) <- names(ex3.dag.data[,c(1,2)]);
  ## specific distributions
  mydists <- list(b1="binomial",
                  b2="binomial"
  );

  ## just compute marginal likelihood - use internal code via max.mode.error=0
  ## as using INLA is the default here.
  expect_silent({
    myres.c <- fitAbn(dag=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists, group.var="group",cor.vars=c("b1","b2"),
                      centre=TRUE,compute.fixed=FALSE,control=list(max.mode.error=0))
  })
  expect_false(any(myres.c$used.INLA))

  ## compare with INLA estimate
  if(!testthat:::on_cran()) {
    if(requireNamespace("INLA", quietly = TRUE)){
      expect_silent({
        myres.inla <- fitAbn(dag=mydag,data.df=ex3.dag.data[,c(1,2,14)], data.dists=mydists,group.var="group",cor.vars=c("b1","b2"),
                             centre=TRUE,compute.fixed=FALSE,control=list(max.mode.error=100))
      })
      expect_true(any(myres.inla$used.INLA)) #BUG: this switches back to internal C code eventhough INLA is enforced...
    }
  }

  ## compare log marginal likelihoods for each node and total DAG - should be very similar

  ## now for marginals - INLA is strongly preferable for estimating marginals for nodes
  ## with random effects as it is far faster, but may not be reliable
  ## see www.r-bayesian-networks.org/quality-assurance

  # INLA's estimates of the marginals, using high n.grid=500
  # as this makes the plots smoother - see below.
  # myres.inla <- fitAbn(dag=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
  #                    group.var="group",cor.vars=c("b1","b2"),
  #                    compute.fixed=TRUE,max.mode.error=100,
  #                    n.grid=500,max.hessian.error = 10E-02);

  ## this is NOT recommended - marginal density estimation using fitabn in mixed models
  ## is really just for diagnostic purposes, better to use fitabn.inla() here
  ## but here goes...be patient
  # myres.c <- fitAbn(dag=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
  #                 group.var="group",cor.vars=c("b1","b2"),compute.fixed=TRUE,
  #                 max.mode.error=0,max.hessian.error = 10E-02);


  ### these are very similar although not exactly identical

  ## use internal code but only to compute a single parameter over a specified grid
  ## this can be necessary if the simple auto grid finding functions does a poor job

  # myres.c <- fitAbn(dag=mydag,data.df=ex3.dag.data[,c(1,2,14)],data.dists=mydists,
  #                 group.var="group",
  #                 cor.vars=c("b1","b2"),centre=FALSE,compute.fixed=TRUE,
  #                 marginal.node=1,marginal.param=3,## precision term in node 1
  #                 variate.vec=seq(0.05,1.5,len=25),max.hessian.error = 10E-02);


  ## use built-in simulated data set

  mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")];## take a subset of cols

  ## setup distribution list for each node
  mydists <- list(b1="binomial",
                  b2="binomial",
                  b3="binomial",
                  g1="gaussian",
                  b4="binomial",
                  p2="poisson",
                  p4="poisson"
  );
  ## null model - all independent variables
  mydag.empty <- matrix(data=c(
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0, #
    0,0,0,0,0,0,0  #
  ), byrow=TRUE,ncol=7);
  colnames(mydag.empty) <- rownames(mydag.empty) <- names(mydat);
  ## now repeat but include some dependencies
  mydag <- mydag.empty;
  mydag["b1","b2"] <- 1; # b1 <- b2
  mydag["b2","p4"] <- 1; # b2 <- p4
  mydag["b2","g1"] <- 1; # b2 <- g1
  mydag["g1","p2"] <- 1; # g1 <- p2
  mydag["b3","g1"] <- 1; # b3 <- g1
  mydag["b4","b1"] <- 1; # b4 <- b1
  mydag["p4","g1"] <- 1; # p4 <- g1

  ## now fit the model to calculate its goodness of fit
  expect_silent({
    myres.c <- fitAbn(dag=mydag,data.df=mydat,data.dists=mydists,method="mle",centre=TRUE)
  })
})

test_that("fitabn() works with all distributions, grouping and class abnCache", {
  df <- FCV[, c(11:15)]
  mydists <- list(Pedigree="binomial",
                  Outdoor="binomial",
                  Sex="multinomial",
                  GroupSize="poisson",
                  Age="gaussian")
  mydists <- mydists[-1] # remove grouping variable from distribution list

  ## buildScoreCache -> mostProbable() -> fitAbn()
  expect_error({
    mycache.bayes <- buildScoreCache(data.df = df, data.dists = mydists, method = "bayes",
                                     group.var = "Pedigree", adj.vars = NULL, cor.vars = NULL,
                                     dag.banned = NULL, dag.retained = NULL,
                                     max.parents = 1,
                                     which.nodes = NULL, defn.res = NULL)
  }, regexp = "not yet implemented")
  # expect_no_error({
  #   mp.dag.bayes <- mostProbable(score.cache = mycache.bayes, verbose = FALSE)
  # })
  # myres.bayes <- fitAbn(object = mp.dag.bayes, data.df = df, data.dists = mydists)

  expect_no_error({
    suppressWarnings({
      mycache.mle <- buildScoreCache(data.df = df, data.dists = mydists, method = "mle",
                                     group.var = "Pedigree", adj.vars = NULL, cor.vars = NULL,
                                     dag.banned = NULL, dag.retained = NULL,
                                     max.parents = 1,
                                     which.nodes = NULL, defn.res = NULL)
    }) # ignore non-convergence warnings
  })
  expect_no_error({
    mp.dag.mle <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
  })
  expect_error({
    myres.mle <- fitAbn(object = mp.dag.mle, data.df = df, data.dists = mydists)
  }, regexp = "'data.df' and 'object' provided but can only accept one of them")
  expect_error({
    myres.mle <- fitAbn(object = mp.dag.mle, data.dists = mydists)
  }, regexp = "'data.dists' and 'object' provided but can only accept one of them")
  expect_no_error({
    suppressWarnings({
      myres.mle <- fitAbn(object = mp.dag.mle, method = "mle")
    })
  })


  ## buildScoreCache -> searchHillClimber() -> fitAbn()
  expect_error({
    mycache.bayes <- buildScoreCache(data.df = df, data.dists = mydists, method = "bayes",
                                     group.var = "Pedigree", adj.vars = NULL, cor.vars = NULL,
                                     dag.banned = NULL, dag.retained = NULL,
                                     max.parents = 1,
                                     which.nodes = NULL, defn.res = NULL)
  }, regexp = "not yet implemented")
  # expect_no_error({
  #   hc.dag.bayes <- searchHillClimber(score.cache = mycache.bayes, verbose = FALSE)
  # })
  # myres.bayes <- fitAbn(object = hc.dag.bayes, data.df = df, data.dists = mydists)

  expect_no_error({
    suppressWarnings({
      mycache.mle <- buildScoreCache(data.df = df, data.dists = mydists, method = "mle",
                                     group.var = "Pedigree", adj.vars = NULL, cor.vars = NULL,
                                     dag.banned = NULL, dag.retained = NULL,
                                     max.parents = 1,
                                     which.nodes = NULL, defn.res = NULL)
    }) # ignore non-convergence warnings

  })
  expect_no_error({
    hc.dag.mle <- searchHillClimber(score.cache = mycache.mle, verbose = FALSE)
  })
  expect_error({
    myres.mle <- fitAbn(object = hc.dag.mle, data.df = df, data.dists = mydists)
  }, regexp = "'data.df' and 'object' provided but can only accept one of them")
  expect_error({
    suppressWarnings({
      myres.mle <- fitAbn(object = hc.dag.mle, data.dists = mydists)
    })
  }, regexp = "'data.dists' and 'object' provided but can only accept one of them")
  expect_no_error({
    suppressWarnings({
      myres.mle <- fitAbn(object = hc.dag.mle, method = "mle")
    })
  })
})

