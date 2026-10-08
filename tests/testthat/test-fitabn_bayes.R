test_that("Test fitAbn.bayes()", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package


  load(file="testdata/fitabn_ex0.Rdata")
  # load(file='tests/testthat/testdata/fitabn_ex0.Rdata')

  myres.c.test <- fitAbn(dag=mydag, data.df=mydat, data.dists=mydists, method = "bayes")

  expect_equal(unclass(myres.c.test)[["mliknode"]][[1]], myres.c[[1]])
  expect_equal(unclass(myres.c.test)[["mliknode"]][[2]], myres.c[[2]])
  expect_equal(unclass(myres.c.test)[["mliknode"]][[3]], myres.c[[3]])
  expect_equal(unclass(myres.c.test)[["mliknode"]][[4]], myres.c[[4]])
  expect_equal(unclass(myres.c.test)[["mliknode"]][[5]], myres.c[[5]])
  expect_equal(unclass(myres.c.test)[["mliknode"]][[6]], myres.c[[6]])
  expect_equal(unclass(myres.c.test)[["mliknode"]][[7]], myres.c[[7]])
  expect_equal(unclass(myres.c.test)[["modes"]], myres.c[[8]])
  expect_equal(unclass(myres.c.test)[["error.code"]], myres.c[[9]])
  expect_equal(unclass(myres.c.test)[["hessian.accuracy"]], myres.c[[10]])
  expect_equal(unname(unclass(myres.c.test)[["error.code.desc"]]), unname(myres.c[[11]])) # reference predates the restored names
  expect_equal(unclass(myres.c.test)[["mlik"]], myres.c[[12]])
  expect_equal(unname(unclass(myres.c.test)[["used.INLA"]]), myres.c[[13]]) # historical reasons. Can be updated in the future.
})

# --- invariants the bayesian-network JSON format relies on ------------------

test_that("bayes modes carry gaussian precisions", {
  jfx_skip_if_no_bayes()
  fit <- jfx_fit("ex1_bayes")
  expect_true("g1|precision" %in% names(fit$modes$g1))
  expect_gt(fit$modes$g1[["g1|precision"]], 0)
})

test_that("coef of grouped gaussian bayes nodes contains no precision", {
  jfx_skip_if_no_bayes()
  fit <- jfx_fit("adg_bayes_grouped")
  expect_false(any(grepl("precision", colnames(fit$coef$adg))))
  expect_equal(ncol(fit$coef$adg), 2L)
  expect_true(all(c("adg|precision", "adg|group.precision") %in% names(fit$modes$adg)))
})

test_that("bayes error.code.desc is named by node", {
  jfx_skip_if_no_bayes()
  fit <- jfx_fit("ex1_bayes")
  expect_equal(names(fit$error.code.desc), names(fit$error.code))
  expect_equal(names(fit$error.code.desc), colnames(fit$abnDag$dag))
})
