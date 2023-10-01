test_that("Test mostProbable() works", {

  load(file="testdata/buildscorecache_ex1.Rdata")

  invisible(mycache.test <- buildScoreCache(data.df=mydat, data.dists=mydists, method = "bayes", max.parents=1))
  class(mycache.test) <-  c("abnCache")
  expect_silent(mp.dag.test <- mostProbable(score.cache=mycache.test, verbose=FALSE))
})

test_that("mostProbable() is backward compatible with `ex0.dag.data`", {
  ## use built-in simulated data set
  mydat <- ex0.dag.data[,c("b1","b2","g1","g2","p1","p2")];
  ## take a subset of cols

  ## setup distribution list for each node
  mydists <- list(b1="binomial",
                  b2="binomial",
                  g1="gaussian",
                  g2="gaussian",
                  p1="poisson",
                  p2="poisson"
  );

  #use simple banlist with no constraints
  ban <- matrix(c(
    #   1 2 3 4 5 6
    0,0,0,0,0,0, # 1
    0,0,0,0,0,0, # 2
    0,0,0,0,0,0, # 3
    0,0,0,0,0,0, # 4
    0,0,0,0,0,0, # 5
    0,0,0,0,0,0 # 6
  ),byrow=TRUE,ncol=6);

  colnames(ban) <- rownames(ban) <- names(mydat); #names must be set
  ban["b1","b2"] <- 1; # now ban arc from b2 to b1

  retain <- matrix(c(
    #   1 2 3 4 5 6
    0,0,0,0,0,0, # 1
    0,0,0,0,0,0, # 2
    0,0,0,0,0,0, # 3
    0,0,0,0,0,0, # 4
    0,0,0,0,0,0, # 5
    0,0,0,0,0,0 # 6
  ),byrow=TRUE,ncol=6);

  colnames(retain) <- rownames(retain) <- names(mydat); #names must be set
  retain["g1","g2"] <- 1; # always retain arc from g2 to g1
  # parent limits
  max.par <- list("b1"=1,"b2"=1,"g1"=1,"g2"=0,"p1"=1,"p2"=2);
  ## now build cache
  mycache <- buildScoreCache(data.df=mydat,data.dists=mydists,
                             dag.banned=ban, dag.retained=retain,max.parents=max.par);

  #now find the globally best DAG
  expect_no_error({
    mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE);
  })
  expect_equal(mp.dag$dag, retain)

  # get the corresponding best goodness of fit - network score
  expect_message({
    m1 <- fitAbn(dag=mp.dag,data.df=mydat,data.dists=mydists)$mlik;
  })
  expect_error({
    m2 <- fitAbn(object = mp.dag,data.df=mydat,data.dists=mydists)$mlik;
  }, regexp = "`data.df` and `object` provided but can only accept one of them")
  expect_no_message({
    m2 <- fitAbn(object = mp.dag)$mlik;
  })
  expect_equal(m1, m2)
})

test_that("mostProbable() is backward compatible with `ex1.dag.data`", {
  ## Second example ############
  mydat <- ex1.dag.data;## this data comes with abn see ?ex1.dag.data

  ## setup distribution list for each node
  mydists <- list(b1="binomial",
                  p1="poisson",
                  g1="gaussian",
                  b2="binomial",
                  p2="poisson",
                  b3="binomial",
                  g2="gaussian",
                  b4="binomial",
                  b5="binomial",
                  g3="gaussian"
  );

  #use simple banlist with no constraints
  ban <- matrix(c(
    #   1 2 3 4 5 6
    0,0,0,0,0,0,0,0,0,0, # 1
    0,0,0,0,0,0,0,0,0,0, # 2
    0,0,0,0,0,0,0,0,0,0, # 3
    0,0,0,0,0,0,0,0,0,0, # 4
    0,0,0,0,0,0,0,0,0,0, # 5
    0,0,0,0,0,0,0,0,0,0, # 6
    0,0,0,0,0,0,0,0,0,0, # 7
    0,0,0,0,0,0,0,0,0,0, # 8
    0,0,0,0,0,0,0,0,0,0, # 9
    0,0,0,0,0,0,0,0,0,0  # 10
  ),byrow=TRUE,ncol=10);

  colnames(ban) <- rownames(ban) <- names(mydat); #names must be set

  retain <- matrix(c(
    #   1 2 3 4 5 6
    0,0,0,0,0,0,0,0,0,0, # 1
    0,0,0,0,0,0,0,0,0,0, # 2
    0,0,0,0,0,0,0,0,0,0, # 3
    0,0,0,0,0,0,0,0,0,0, # 4
    0,0,0,0,0,0,0,0,0,0, # 5
    0,0,0,0,0,0,0,0,0,0, # 6
    0,0,0,0,0,0,0,0,0,0, # 7
    0,0,0,0,0,0,0,0,0,0, # 8
    0,0,0,0,0,0,0,0,0,0, # 9
    0,0,0,0,0,0,0,0,0,0  # 10
  ),byrow=TRUE,ncol=10);
  colnames(retain) <- rownames(retain) <- names(mydat); #names must be set

  ## parent limits
  max.par <- list("b1"=2,"p1"=2,"g1"=2,"b2"=2,"p2"=2,"b3"=2,"g2"=2,"b4"=2,"b5"=2,"g3"=2);
  ## now build cache
  mycache.exemple1 <- buildScoreCache(data.df=mydat,data.dists=mydists,
                                      dag.banned=ban, dag.retained=retain,max.parents=max.par);

  #now find the globally best DAG
  expect_no_error({
    mp.dag <- mostProbable(score.cache=mycache.exemple1, verbose = FALSE);
  })

  # get the corresponding best goodness of fit - network score
  expect_no_error({
    expect_message({
      m1 <- fitAbn(dag=mp.dag,data.df=mydat,data.dists=mydists)$mlik; # this is ok, because mp.dag is provided to 'dag=', eventhough it actually is an object.
    })
  })
  expect_error({
    m1 <- fitAbn(dag=mp.dag)$mlik;
  }, regexp = "'data.df' is missing but must be provided") # if mp.dag is provided as dag, even it is an object, we need data.df and data.dists.
  expect_error({
    m2 <- fitAbn(object = mp.dag,data.df=mydat,data.dists=mydists)$mlik;
  }, regexp = "`data.df` and `object` provided but can only accept one of them")
  expect_no_message({
    m2 <- fitAbn(object = mp.dag)$mlik;
  })
  expect_equal(m1, m2)

  ## plot the best model
  expect_silent({
    myres <- fitAbn(object=mp.dag,create.graph=TRUE);
    }) # create.graph argument is deprecated and should be ignored.
  plt <- plot(myres)
  plt.adjMat <- plt@adjMat
  mp.dag.adjMat <- mp.dag$dag
  rownames(mp.dag.adjMat) <- NULL
  expect_equal(plt.adjMat, mp.dag.adjMat)
})

test_that("mostProbable() is backward compatible with `ex3.dag.data`", {
  #################################################################
  ## example 3 - models with random effects
  #################################################################

  mydat <- ex3.dag.data[,c(1:4,14)];
  ## this data comes with abn see ?ex3.dag.data

  mydists <- list(b1="binomial",
                  b2="binomial",
                  b3="binomial",
                  b4="binomial"
  );
  max.par <- 1;

  if(requireNamespace("INLA", quietly = TRUE)){
    mycache.mixed <- buildScoreCache(data.df=mydat,data.dists=mydists,
                                     group.var="group",cor.vars=c("b1","b2","b3","b4"),
                                     ## each node uses a random effect adjustment
                                     max.parents=max.par);

    ## find the most probable DAG
    expect_no_error({
      mp.dag <- mostProbable(score.cache=mycache.mixed, verbose = FALSE);
    })

    ## get goodness of fit
    expect_error({
      m <- fitAbn(object = mp.dag,data.df=mydat,data.dists=mydists,group.var="group",cor.vars=c("b1","b2","b3","b4"))$mlik;
    }, regexp = "`data.df` and `object` provided but can only accept one of them")
    expect_no_error({
      m <- fitAbn(object = mp.dag,group.var="group",cor.vars=c("b1","b2","b3","b4"))$mlik;
    })
  }
})

test_that("mostProbable() simple, historic numeric test", {
  load(file="testdata/buildscorecache_ex1.Rdata")
  # load(file='tests/testthat/testdata/buildscorecache_ex1.Rdata')

  invisible(mycache.test <- buildScoreCache(data.df=mydat, data.dists=mydists, method = "bayes", max.parents=max.par))
  class(mycache.test) <-  c("abnCache")
  invisible(mp.dag.test <- mostProbable(score.cache=mycache.test, verbose=FALSE))
  expect_equal(unclass(mp.dag.test[[1]]), (mp.dag))
})
