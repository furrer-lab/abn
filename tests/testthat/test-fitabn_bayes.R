test_that("Test fitAbn.bayes()", {

  load(file="testdata/fitabn_ex0.Rdata")
  # load(file='tests/testthat/testdata/fitabn_ex0.Rdata')

  myres.c.test <- fitAbn(dag=mydag, data.df=mydat, data.dists=mydists, method = "bayes")

  expect_equal(unclass(myres.c.test)[[1]], myres.c[[1]])
  expect_equal(unclass(myres.c.test)[[2]], myres.c[[2]])
  expect_equal(unclass(myres.c.test)[[3]], myres.c[[3]])
  expect_equal(unclass(myres.c.test)[[4]], myres.c[[4]])
  expect_equal(unclass(myres.c.test)[[5]], myres.c[[5]])
  expect_equal(unclass(myres.c.test)[[6]], myres.c[[6]])
  expect_equal(unclass(myres.c.test)[[7]], myres.c[[7]])
  expect_equal(unclass(myres.c.test)[[8]], myres.c[[8]])
  expect_equal(unclass(myres.c.test)[[9]], myres.c[[9]])
  expect_equal(unclass(myres.c.test)[[10]], myres.c[[10]])
  expect_equal(unclass(myres.c.test)[[11]], myres.c[[11]])
  expect_equal(unclass(myres.c.test)[[12]], myres.c[[12]])
  expect_equal(unclass(myres.c.test)[[15]], myres.c[[13]]) # historical reasons. Can be updated in the future.
})
