test_that("getmarginals() works", {
  ## Generated test data with:
  # # get data
  # mydat <- ex5.dag.data[,-19] ## get the data - drop group variable
  #
  # # Restrict DAG
  # banned<-matrix(c(
  #   # 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8
  #   0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b1
  #   1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b2
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b3
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b4
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b5
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b6
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g1
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g2
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g3
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g4
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g5
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g6
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g7
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g8
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g9
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g10
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g11
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 # g12
  # ),byrow=TRUE,ncol=18)
  #
  # colnames(banned)<-rownames(banned)<-names(mydat)
  #
  # retain<-matrix(c(
  #   # 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b1
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b2
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b3
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b4
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b5
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # b6
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g1
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g2
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g3
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g4
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g5
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g6
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g7
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g8
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g9
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g10
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # g11
  #   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 # g12
  # ),byrow=TRUE,ncol=18)
  # ## again must set names
  # colnames(retain)<-rownames(retain)<-names(mydat)
  #
  # # set distributions
  # mydists<-list(b1="binomial",
  #               b2="binomial",
  #               b3="binomial",
  #               b4="binomial",
  #               b5="binomial",
  #               b6="binomial",
  #               g1="gaussian",
  #               g2="gaussian",
  #               g3="gaussian",
  #               g4="gaussian",
  #               g5="gaussian",
  #               g6="gaussian",
  #               g7="gaussian",
  #               g8="gaussian",
  #               g9="gaussian",
  #               g10="gaussian",
  #               g11="gaussian",
  #               g12="gaussian"
  # )
  #
  # # Compute score cache
  # mycache.1par <- buildScoreCache(data.df=mydat,data.dists=mydists, max.parents=1,centre=TRUE)
  #
  # # Estimate most probable DAG
  # mp.dag <- mostProbable(score.cache = mycache.1par)
  #
  # # Generate marginal densities
  # marg.f <- fitAbn(object = mp.dag, method = "bayes", compute.fixed = TRUE, n.grid = 1000)

  # Load testdata
  # load("tests/testthat/testdata/getmarginals_1.RData")
  load("../testthat/testdata/getmarginals_1.RData")

  expect_error({
    # if marginal.node is supplied, variate.vec must be supplied as well
    getmarginals(res.list = res.list,
                 data.df = data.df,
                 dag.m = dag,
                 var.types = var.types,
                 max.parents = max.parents,
                 mean = control[["mean"]],
                 prec = control[["prec"]],
                 loggam.shape = control[["loggam.shape"]],
                 loggam.inv.scale = control[["loggam.inv.scale"]],
                 max.iters = control[["max.iters"]],
                 epsabs = control[["epsabs"]],
                 verbose = verbose,
                 error.verbose = control[["error.verbose"]],
                 trace = as.integer(control[["trace"]]),
                 grouped.vars = as.integer(grouped.vars-1),
                 group.ids = as.integer(group.ids),
                 epsabs.inner = control[["epsabs.inner"]],
                 max.iters.inner = control[["max.iters.inner"]],
                 finite.step.size = control[["finite.step.size"]],
                 hessian.params = control[["hessian.params"]],
                 max.iters.hessian = control[["max.iters.hessian"]],
                 min.pdf = control[["min.pdf"]],
                 marginal.node = 1,
                 marginal.param = control[["marginal.param"]],
                 variate.vec = NULL,
                 n.grid = control[["n.grid"]],
                 INLA.marginals = res.list[["used.INLA"]],
                 iter.max = control[["max.grid.iter"]],
                 max.hessian.error = as.double(control[["max.hessian.error"]]),
                 factor.brent = as.double(control[["factor.brent"]]),
                 maxiters.hessian.brent = as.integer(control[["maxiters.hessian.brent"]]),
                 num.intervals.brent = as.double(control[["num.intervals.brent"]])
    )
  })

  expect_no_error({
    res.list <- getmarginals(res.list = res.list, ## rest of arguments as for call to C fitabn
                             data.df = data.df,
                             dag.m = dag,
                             var.types = var.types,
                             max.parents = max.parents,
                             mean = control[["mean"]],
                             prec = control[["prec"]],
                             loggam.shape = control[["loggam.shape"]],
                             loggam.inv.scale = control[["loggam.inv.scale"]],
                             max.iters = control[["max.iters"]],
                             epsabs = control[["epsabs"]],
                             verbose = verbose,
                             error.verbose = control[["error.verbose"]],
                             trace = as.integer(control[["trace"]]),
                             grouped.vars = as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C (changed from earlier fitabn)
                             group.ids = as.integer(group.ids),
                             epsabs.inner = control[["epsabs.inner"]],
                             max.iters.inner = control[["max.iters.inner"]],
                             finite.step.size = control[["finite.step.size"]],
                             hessian.params = control[["hessian.params"]],
                             max.iters.hessian = control[["max.iters.hessian"]],
                             min.pdf = control[["min.pdf"]],
                             marginal.node = control[["marginal.node"]],
                             marginal.param = control[["marginal.param"]],
                             variate.vec = control[["variate.vec"]],
                             n.grid = control[["n.grid"]],
                             INLA.marginals = res.list[["used.INLA"]],
                             iter.max = control[["max.grid.iter"]],
                             max.hessian.error = as.double(control[["max.hessian.error"]]),
                             factor.brent = as.double(control[["factor.brent"]]),
                             maxiters.hessian.brent = as.integer(control[["maxiters.hessian.brent"]]),
                             num.intervals.brent = as.double(control[["num.intervals.brent"]])
    )
  })

  # Check that the result is a list
  expect_type(res.list, "list")

  # Check that the result has the correct number of elements
  expect_equal(length(res.list), 10)

  # Check that the result has the correct names
  expect_equal(names(res.list), c("modes", "error.code", "hessian.accuracy", "error.code.desc", "mliknode", "mlik", "mse", "coef", "used.INLA", "marginals"))

  # Check that the result has the correct class
  expect_equal(class(res.list[["modes"]]), "list")
  expect_equal(class(res.list[["error.code"]]), "numeric")
  expect_equal(class(res.list[["hessian.accuracy"]]), "numeric")
  expect_equal(class(res.list[["error.code.desc"]]), "character")
  expect_equal(class(res.list[["mliknode"]]), "numeric")
  expect_equal(class(res.list[["mlik"]]), "numeric")
  expect_equal(class(res.list[["mse"]]), "numeric")
  expect_equal(class(res.list[["coef"]]), "list")
  expect_equal(class(res.list[["used.INLA"]]), "logical")
  expect_equal(class(res.list[["marginals"]]), "list")

  # Check that the result has the correct dimensions
  expect_equal(length(res.list[["modes"]]), 18)
  expect_equal(length(res.list[["error.code"]]), 18)
  expect_equal(length(res.list[["hessian.accuracy"]]), 18)
  expect_equal(length(res.list[["error.code.desc"]]), 18)
  expect_equal(length(res.list[["mliknode"]]), 18)
  expect_equal(length(res.list[["mlik"]]), 1)
  expect_equal(length(res.list[["mse"]]), 12) # only for gaussian nodes
  expect_equal(length(res.list[["coef"]]), 18)
  expect_equal(length(res.list[["used.INLA"]]), 18)
  expect_equal(length(res.list[["marginals"]]), 18)

  # Check that there were no errors in the computation
  expect_equal(unname(res.list[["error.code"]]), rep(0, 18))
})
