test_that("Markov blanket is correctly retrieved", {
  m2 <- matrix( 0, 4,4, dimnames=list(c("a","b","c","d"), c("a","b","c","d")))
  m2[c(2,4,12,8)] <- 1

  expect_equal(sort(mb(dag = m2, node = "c")), sort(c("a","b","d")))
  expect_error(mb(dag = m2))

  dists <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="binomial", f="binomial")
  data.param <- matrix(data=c(0, 0.2, 0.5, 0, 0.01, 0, 0, 0, 0.3, 0.1, 0, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0), nrow=6L, ncol=6L, byrow=TRUE)
  colnames(data.param) <- rownames(data.param) <- names(dists)
  a <- mb(dag=data.param, node="b", data.dists=dists)
  b <- mb(dag=data.param, node="e", data.dists=dists)
  c <- mb(dag=data.param, node=c("b", "e"), data.dists=dists)

  expect_equal(a, c("a", "c", "d", "f", "e"))
  expect_equal(b, c("a", "f", "b", "c"))
  expect_equal(c, c("a", "c", "d", "f", "e", "b"))
})
