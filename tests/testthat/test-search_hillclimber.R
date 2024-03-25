test_that("Test searchHillClimber() works", {
  load(file="testdata/buildscorecache_ex1.Rdata")
  expect_error({
    mp.dag.test <- searchHillClimber(score.cache=mycache, verbose=FALSE)
  }, regexp = "should be an object of class 'abnCache'")
  class(mycache) <-  c("abnCache")
  expect_silent({
    mp.dag.test <- searchHillClimber(score.cache=mycache, verbose=FALSE)
  })
  expect_s3_class(mp.dag.test, "abnHillClimber")
  expect_s3_class(mp.dag.test, "abnLearned")

})

test_that("searchHillClimber() is backward compatible with `ex1.dag.data`", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

    ##############################################
    ## example 1: use built-in simulated data set
    ##############################################
    mydat <- ex1.dag.data; ## this data comes with abn see ?ex1.dag.data

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

    ## not run because may take some minutes for buildScoreCache()
    ## parent limits
    max.par <- list("b1"=2,"p1"=2,"g1"=2,"b2"=2,"p2"=2,"b3"=2,"g2"=2,"b4"=2,"b5"=2,"g3"=2);
    ## now build cache
    mycache <- buildScoreCache(data.df=mydat,data.dists=mydists,max.parents=max.par);

    expect_no_error({
      heur.res <- searchHillClimber(score.cache=mycache,
                                    num.searches=100,seed=0,verbose=FALSE,timing.on=FALSE);
    })
    expect_s3_class(heur.res, "abnHillClimber")
    expect_s3_class(heur.res, "abnLearned")
})

test_that("searchHillClimber() is backward compatible with `ex3.dag.data`", {
    ###################################################################################################
    ## example 2 - glmm example - but no difference here as the format of the score cache is identical
    ###################################################################################################
    mydat <- ex3.dag.data[,c(1:5,14)];## this data comes with abn see ?ex1.dag.data

    mydists <- list(b1="binomial",
                    b2="binomial",
                    b3="binomial",
                    b4="binomial",
                    b5="binomial"
    );
    max.par <- 1;

    if(!testthat:::on_cran()) {
      if(requireNamespace("INLA", quietly = TRUE)){
        mycache.mixed <- buildScoreCache(data.df=mydat,data.dists=mydists,group.var="group",
                                         cor.vars=c("b1","b2","b3","b4","b5"),
                                         max.parents=max.par);

        # now peform 1000 greedy searches
        expect_no_error({
          heur.res <- searchHillClimber(score.cache=mycache.mixed,num.searches=100,
                                        seed=0,verbose=FALSE,timing.on=FALSE);
        })
        expect_s3_class(heur.res, "abnHillClimber")
        expect_s3_class(heur.res, "abnLearned")
      }
    } else {
      skip("INLA is not tested on CRAN")
    }
})

test_that("searchHillClimber() simple, historic numeric test", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  load(file="testdata/buildscorecache_ex1.Rdata")
  # load(file='tests/testthat/testdata/buildscorecache_ex1.Rdata')

  invisible(mycache.test <- buildScoreCache(data.df=mydat, data.dists=mydists, method = "bayes", max.parents=max.par))
  class(mycache.test) <-  c("abnCache")

  invisible(heur.res.test <- searchHillClimber(score.cache=mycache.test, num.searches=10, seed=42, verbose=FALSE, timing.on=TRUE))

  expect_equal(heur.res.test[[1]], heur.res[[1]])
  expect_equal(heur.res.test[[2]], heur.res[[2]])
  expect_equal(heur.res.test[[3]], heur.res[[3]])
  expect_equal(heur.res.test[[4]], heur.res[[4]])
  expect_equal(heur.res.test[[5]], heur.res[[5]])
  expect_equal(heur.res.test[[6]], heur.res[[6]])
})
