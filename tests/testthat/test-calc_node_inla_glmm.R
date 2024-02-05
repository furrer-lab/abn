test_that("calc.node.inla.glmm() works", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  # We can use the same data as for calc.node.inla.glm() for this rough test.
  load(file = "testdata/calc.node.inla.glm_1.RData")
  # load(file = 'tests/testthat/testdata/calc.node.inla.glm_1.RData')

  expect_no_error({
    res <- calc.node.inla.glmm(
      child,
      dag,
      data.df,
      data.dists,
      rep(1, dim(data.df)[1]),
      ## ntrials
      rep(1, dim(data.df)[1]),
      ## exposure
      TRUE,
      mean.intercept,
      prec.intercept,
      control[["mean"]],
      control[["prec"]],
      control[["loggam.shape"]],
      control[["loggam.inv.scale"]],
      verbose,
      nthreads = control[["ncores"]] # single threaded
    )
  })

  expect_no_error({
    calc.node.inla.glmm(
      child,
      dag,
      data.df,
      data.dists,
      rep(1, dim(data.df)[1]),
      ## ntrials
      rep(1, dim(data.df)[1]),
      ## exposure
      TRUE,
      mean.intercept,
      prec.intercept,
      control[["mean"]],
      control[["prec"]],
      control[["loggam.shape"]],
      control[["loggam.inv.scale"]],
      verbose.loc = FALSE,
      nthreads = 2
    )
  })
})
