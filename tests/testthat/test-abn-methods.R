test_that("print.abnDag() works.", {
  mydag <- createAbnDag(dag = ~a+b|a, data.df = data.frame("a"=1, "b"=1))

  expect_output(print(mydag))
})

test_that("summary.abnDag() works.", {
  mydag <- createAbnDag(dag = ~a+b|a, data.df = data.frame("a"=1, "b"=1))

  expect_no_error({
    summary(mydag)
  })
})

test_that("plot.abnDag() works.", {
  mydag <- createAbnDag(dag = ~a+b|a, data.df = data.frame("a"=1, "b"=1))

  if(.Platform$OS.type == "unix") {
    capture.output({
      expect_no_error({
        plot(mydag)
        })
      },
      file = "/dev/null")
  }
})

test_that("print.abnCache() works.", {
  ## Subset of the build-in dataset, see  ?ex0.dag.data
  mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")] ## take a subset of cols

  ## setup distribution list for each node
  mydists <- list(b1="binomial", b2="binomial", g1="gaussian",
                  g2="gaussian", b3="binomial", g3="gaussian")

  # Structural constraints
  # ban arc from b2 to b1
  # always retain arc from g2 to g1

  ## parent limits
  max.par <- list("b1"=2, "b2"=2, "g1"=2, "g2"=2, "b3"=2, "g3"=2)

  ## now build the cache of pre-computed scores accordingly to the structural constraints
  res.c <- buildScoreCache(data.df=mydat, data.dists=mydists,
                           dag.banned= ~b1|b2, dag.retained= ~g1|g2, max.parents=max.par)

  expect_output({
    print(res.c)
  })
})

test_that("print.abnHeuristic() works.", {
  mydat <- ex1.dag.data ## this data comes with abn see ?ex1.dag.data

  ## setup distribution list for each node
  mydists<-list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                p2="poisson", b3="binomial", g2="gaussian", b4="binomial",
                b5="binomial", g3="gaussian")

  mycache <- buildScoreCache(data.df = mydat, data.dists = mydists, max.parents = 2)

  ## Now peform 10 greedy searches
  heur.res <- searchHeuristic(score.cache = mycache, data.dists = mydists,
                              start.dag = "random", num.searches = 10,
                              max.steps = 50)
  expect_output({
    print(heur.res)
  })
})

test_that("plot.abnHeuristic() works.", {
  mydat <- ex1.dag.data ## this data comes with abn see ?ex1.dag.data

  ## setup distribution list for each node
  mydists<-list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                p2="poisson", b3="binomial", g2="gaussian", b4="binomial",
                b5="binomial", g3="gaussian")

  mycache <- buildScoreCache(data.df = mydat, data.dists = mydists, max.parents = 2)

  ## Now peform 10 greedy searches
  heur.res <- searchHeuristic(score.cache = mycache, data.dists = mydists,
                              start.dag = "random", num.searches = 10,
                              max.steps = 50)
  if(.Platform$OS.type == "unix") {
    capture.output({
      expect_no_error({
        plot(heur.res)
      })
    },
    file = "/dev/null")
  }
})

test_that("print.abnHillClimber() works.", {
  ## this data comes with abn see ?ex1.dag.data
  mydat <- ex1.dag.data

  ## setup distribution list for each node
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", b4="binomial",
                  b5="binomial", g3="gaussian")

  ## Build cache may take some minutes for buildScoreCache()
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists,
                             max.parents=2);

  # now peform 10 greedy searches
  heur.res <- searchHillClimber(score.cache=mycache,
                                num.searches=10, timing.on=FALSE)
  expect_output({
    print(heur.res)
  })
})

test_that("plot.abnHillClimber() works.", {
  skip("plotting for abnHillClimber is not up to date with searchHillClimber()'s output.")
  ## this data comes with abn see ?ex1.dag.data
  mydat <- ex1.dag.data

  ## setup distribution list for each node
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", b4="binomial",
                  b5="binomial", g3="gaussian")

  ## Build cache may take some minutes for buildScoreCache()
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists,
                             max.parents=2);

  # now peform 10 greedy searches
  heur.res <- searchHillClimber(score.cache=mycache,
                                num.searches=10, timing.on=FALSE)

  if(.Platform$OS.type == "unix") {
    capture.output({
      expect_no_error({
        plot(heur.res$dag, new=TRUE)
      })
    },
    file = "/dev/null")
  }
})

test_that("print.abnMostprobable() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  expect_output({
    print(mp.dag)
  })
})

test_that("summary.abnMostprobable() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  expect_output({
    summary(mp.dag)
  })
})

test_that("plot.abnMostprobable() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  if(.Platform$OS.type == "unix") {
    capture.output({
      expect_no_error({
        plot(mp.dag)
      })
    },
    file = "/dev/null")
  }
})

test_that("print.abnFit() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  myres <- fitAbn(object=mp.dag,create.graph=TRUE, verbose = FALSE)

  expect_output({
    print(myres)
  })
})

test_that("summary.abnFit() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  myres <- fitAbn(object=mp.dag,create.graph=TRUE, verbose = FALSE)

  expect_output({
    summary(myres)
  })
})

test_that("coef.abnFit() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  myres <- fitAbn(object=mp.dag,create.graph=TRUE, verbose = FALSE)

  expect_output({
    coef(myres)
  })
})

test_that("AIC.abnFit() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  myres <- fitAbn(object=mp.dag,create.graph=TRUE, verbose = FALSE)

  expect_output({
    AIC(myres)
  })
})

test_that("BIC.abnFit() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  myres <- fitAbn(object=mp.dag,create.graph=TRUE, verbose = FALSE)

  expect_output({
    BIC(myres)
  })
})

test_that("logLik.abnFit() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  myres <- fitAbn(object=mp.dag,create.graph=TRUE, verbose = FALSE)

  expect_output({
    logLik(myres)
  })
})

test_that("family.abnFit() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  myres <- fitAbn(object=mp.dag,create.graph=TRUE, verbose = FALSE)

  expect_output({
    family(myres)
  })
})

test_that("nobs.abnFit() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  myres <- fitAbn(object=mp.dag,create.graph=TRUE, verbose = FALSE)

  expect_equal({
    nobs(myres)
  }, 5000)
})

test_that("plot.abnFit() works.", {
  ## This data comes with `abn` see ?ex1.dag.data
  mydat <- ex1.dag.data[1:5000, c(1:7,10)]

  ## Setup distribution list for each node:
  mydists <- list(b1="binomial", p1="poisson", g1="gaussian", b2="binomial",
                  p2="poisson", b3="binomial", g2="gaussian", g3="gaussian")

  ## Parent limits, for speed purposes quite specific here:
  max.par <- list("b1"=0,"p1"=0,"g1"=1,"b2"=1,"p2"=2,"b3"=3,"g2"=3,"g3"=2)
  ## Now build cache (no constraints in ban nor retain)
  mycache <- buildScoreCache(data.df=mydat, data.dists=mydists, max.parents=max.par)

  ## Find the globally best DAG:
  mp.dag <- mostProbable(score.cache=mycache, verbose = FALSE)

  myres <- fitAbn(object=mp.dag,create.graph=TRUE, verbose = FALSE)

  if(.Platform$OS.type == "unix") {
    capture.output({
      expect_no_error({
        plot(myres, new = TRUE)
      })
    },
    file = "/dev/null")
  }
})
