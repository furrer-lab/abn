test_that("Markov blanket is correctly retrieved", {
  m2 <- matrix( 0, 4,4, dimnames=list(c("a","b","c","d"), c("a","b","c","d")))
  m2[c(2,4,12,8)] <- 1
  dists <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian")

  expect_equal(sort(mb(dag = m2, node = "c", data.dists = dists)), sort(c("a","b","d")))
  expect_error(mb(dag = m2))
})

test_that("Markov blanket is correctly retrieved 2", {
  dag <- matrix(c(0,1,0,1,1,0,0,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),nrow=5,byrow=TRUE)
  colnames(dag) =rownames(dag) <- c("g1","g2","b1","b2","b3")

  data.dists_wrongOrder <- list("b1" = "binomial", "b2" = "binomial", b3 = "binomial", g1 = "gaussian", g2 = "gaussian")
  data.dists_correctOrder <- list("g1" = "gaussian", "g2" = "gaussian", "b1" = "binomial", "b2" = "binomial", "b3" = "binomial")

  # Check order of data.dists must be the same as in colnames(dag) or will be sorted with a warning
  expect_equal(abn::mb(dag,node="b1", data.dists_correctOrder), c("g2", "b3", "b2"))
  expect_warning(abn::mb(dag,node="b1", data.dists_wrongOrder))
})

test_that("Markov blanket works with formula statement", {
  ## Defining distribution, dag and (random) data
  dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian",
               e="binomial", f="binomial")
  data.df <- data.frame(a=rnorm(100), b=rnorm(100), c=rnorm(100),
                        d=rnorm(100), e=rbinom(100, 1, 0.5), f=rbinom(100, 1, 0.5))
  dag <- matrix(c(0,1,1,0,1,0,
                  0,0,1,1,0,1,
                  0,0,0,0,0,0,
                  0,0,0,0,0,0,
                  0,0,0,0,0,1,
                  0,0,0,0,0,0), nrow = 6L, ncol = 6L, byrow = TRUE)
  colnames(dag) <- rownames(dag) <- names(dist)

  # test with adjacency matrix and one node
  expect_equal(mb(dag, node = "b", data.dists = dist), c("a","c","d","f","e"))

  # test with adjacency matrix and multiple nodes
  expect_equal(mb(dag, node = c("b","e"), data.dists = dist), c("a","c","d","f", "e", "b"))

  # test with formula statement
  expect_equal(mb(~a|b:c:e+b|c:d:f+e|f, node = "e", data.dists = dist, data.df = data.df), c("a","f", "b", "c"))
})
