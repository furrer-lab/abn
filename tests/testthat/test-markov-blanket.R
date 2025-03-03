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

test_that("Markov blanket is correctly retrieved 2", {
  dag <- matrix(c(0,1,0,1,1,0,0,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),nrow=5,byrow=TRUE)
  colnames(dag) =rownames(dag) <- c("g1","g2","b1","b2","b3")

  data.dists_wrongOrder <- list("b1" = "binomial", "b2" = "binomial", b3 = "binomial", g1 = "gaussian", g2 = "gaussian")
  data.dists_correctOrder <- list("g1" = "gaussian", "g2" = "gaussian", "b1" = "binomial", "b2" = "binomial", "b3" = "binomial")


  # TODO: Check order of data.dists must be the same as in colnames(dag)
  expect_equal(abn::mb(dag,node="b1",data.dists_correctOrder), c("g2", "b3", "b2"))
  expect_error(abn::mb(dag,node="b1",data.dists_wrongOrder))
})
