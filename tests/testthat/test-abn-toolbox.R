test_that("logit() works", {
  expect_equal(abn::logit(x=0.678), boot::logit(p=0.678))
  expect_equal(abn::logit(x=0.783491741), boot::logit(p=0.783491741))

  x.mat <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
  expect_equal(abn::logit(abn::expit(x=x.mat)), x.mat)
})

test_that("expit() works", {
  expect_equal(abn::expit(x=0.678), boot::inv.logit(x=0.678))
  expect_equal(abn::expit(x=-0.783492343421741), boot::inv.logit(x=-0.783492343421741))

  x.mat <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
  expect_equal(abn::logit(abn::expit(x=x.mat)), x.mat)
})

test_that("or() works", {
  # Binomial
  N <- 100000
  mydists <- list(a="binomial",
                  b="binomial",
                  c="binomial")
  a <- rbinom(N, size = 1, prob = 0.75)
  z <- 1+2*a
  pr <- 1/(1+exp(-z))
  b <- rbinom(N, size = 1, prob = pr)
  z2 <- 1+b
  pr2 <- 1/(1+exp(-z2))
  c <- rbinom(N, size = 1, prob = pr2)
  mydf <- data.frame("a" = as.factor(a),
                     "b" = as.factor(b),
                     "c" = as.factor(c))

  x.bc <- table(as.numeric(as.character(mydf$b)), as.numeric(as.character(mydf$c)), dnn=c("b", "c"))

  x.ab <- table(as.numeric(mydf$a), as.numeric(mydf$b), dnn=c("a", "b"))

  x.ac <- table(as.numeric(mydf$a), as.numeric(mydf$c), dnn=c("a", "c"))

  ## OR()
  expect_equal(round(abn::or(x=x.bc)), 3)
  expect_equal(round(abn::or(x=x.ab)), 7, tolerance = 1)
  expect_equal(round(abn::or(x=x.ac)), 1)

  expect_error(abn::or(x=matrix(c(-1, 1, 2, 3), nrow = 2, byrow = TRUE)))
  expect_error(abn::or(x=matrix(c(1,2,3,4,5,6), nrow = 3)))
})

test_that("odds() works", {
  expect_equal({odds(0.5)}, 1)
  expect_equal({odds(0)}, 0)
  expect_equal({odds(1)}, Inf)

  expect_error({odds(1.1)})
  expect_error({odds(-0.1)})
  expect_error({odds(2)})
})

test_that("compareDag() works", {
  a <- matrix(data=0, nrow=3, ncol=3)
  colnames(a) <- rownames(a) <- c("a", "b", "c")

  a1 <- matrix(data=c(0, 0, 0, 1, 0, 0, 1, 0, 0), nrow=3, ncol=3)
  colnames(a1) <- rownames(a1) <- c("a", "b", "c")
  a2 <- matrix(data=c(0, 0, 0, 1, 0, 0, 1, 1, 0), nrow=3, ncol=3)
  colnames(a2) <- rownames(a2) <- c("a", "b", "c")
  b <- matrix(data=c(0, 0, 0, 1, 0, 1, 1, 0, 0), nrow=3, ncol=3)
  colnames(b) <- rownames(b) <- c("a", "b", "c")

  expect_no_error(compareDag(ref=a1, test=b, checkDAG = TRUE))
  expect_equal(suppressWarnings(compareDag(ref=a, test=b))$`Hamming-distance`, expected=3)
  expect_equal(compareDag(ref=a1, test=b)$`Hamming-distance`, expected=1)
  expect_equal(compareDag(ref=a1, test=b)$TPR, expected=1)
  expect_equal(compareDag(ref=a1, test=b)$PPV, expected=2/3)
  expect_equal(compareDag(ref=a2, test=b)$`Hamming-distance`, expected=1)
  expect_equal(compareDag(ref=a2, test=b)$PPV, expected=2/3)
  expect_equal(compareDag(ref=a2, test=b)$FDR, expected=1/3)
  expect_equal(compareDag(ref=a2, test=b)$TPR, expected=2/3)
})

test_that("infoDag() works", {
  dag <- matrix(data=0, nrow=6, ncol=6)
  dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="gaussian", f="gaussian")
  colnames(dag) <- rownames(dag) <- names(dist)

  infoDag.out1 <- infoDag(dag, node.names=names(dist))

  dag <- matrix(data=c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0), nrow=6, ncol=6)
  colnames(dag) <- rownames(dag) <- names(dist)

  infoDag.out2 <- infoDag(dag, node.names=names(dist))

  dag <- matrix(data=c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0), nrow=6, ncol=6)
  colnames(dag) <- rownames(dag) <- names(dist)

  infoDag.out3 <- infoDag(dag, node.names=names(dist))

  expect_equal(infoDag.out1$n.nodes, 6)
  expect_equal(infoDag.out1$n.arcs, 0)
  expect_equal(infoDag.out1$mb.average, 0)
  expect_equal(infoDag.out1$nh.average, 0)
  expect_equal(infoDag.out1$parent.average, 0)
  expect_equal(infoDag.out1$children.average, 0)

  expect_equal(infoDag.out2$n.nodes, 6)
  expect_equal(infoDag.out2$n.arcs, 2)
  expect_equal(infoDag.out2$mb.average, 1)
  expect_equal(infoDag.out2$nh.average, 2/3)
  expect_equal(infoDag.out2$parent.average, 1/3)
  expect_equal(infoDag.out2$children.average, 1/3)

  expect_equal(infoDag.out3$n.nodes, 6)
  expect_equal(infoDag.out3$n.arcs, 3)
  expect_equal(infoDag.out3$mb.average, 2)
  expect_equal(infoDag.out3$nh.average, 1)
  expect_equal(infoDag.out3$parent.average, 0.5)
  expect_equal(infoDag.out3$children.average, 0.5)
})

test_that("simulateDag() works", {
    expect_error(simulateDag(node.name = NULL), regexp = "number of nodes in the DAG")
    expect_error(simulateDag(node.name = c("a", "b"), edge.density = 1.1),
                 regexp = "should be a real number")
    expect_error(simulateDag(node.name = c("a", "b"), edge.density = -0.1),
                 regexp = "should be a real number")
    expect_error(
      simulateDag(node.name = c("a", "b"), edge.density = 0.5, data.dists = c("gaussian", "binomial")),
      regexp = "Node distribution has to be a named object")
    expect_error(
      simulateDag(node.name = c("a", "b", "c"), edge.density = 0.5, data.dists = c(a = "gaussian", b = "binomial")),
      regexp = "The number of distributions provided does not match the number of nodes")
    expect_error(
      simulateDag(node.name = c("a", "b", "c"), edge.density = 0.5, data.dists = c(a = "gaussian", b = "binomial", c = "uniform")),
      regexp = "poisson, binomial, gaussian, multinomial")

    expect_no_error({
      simdag <- simulateDag(node.name = c("a", "b", "c", "d"), edge.density = 0.5, data.dists = list(a = "gaussian", b= "binomial", c= "poisson", d= "multinomial"))
    })
    expect_s3_class(simdag, "abnDag")
    expect_equal(simdag$data.dists, list(a = "gaussian", b= "binomial", c= "poisson", d= "multinomial"))
    expect_equal(simdag$data.df, NULL)
    expect_equal(class(simdag$dag), c("matrix", "array"))
})

test_that("skewness() works", {
  skip_if_not_installed("moments")

  data <- c(19.09, 19.55, 17.89, 17.73, 25.15, 27.27, 25.24, 21.05, 21.65, 20.92, 22.61, 15.71, 22.04, 22.6, 24.25)
  expect_equal(abn:::skewness(x=data), moments::skewness(x=data))
})

test_that("essentialGraph() works", {
  dist.test <- list(a="gaussian", b="gaussian", c="gaussian")

  ## essentialGraph()
  expect_equal(object=essentialGraph(  dag=~a | b + c | b, node.names=names(dist.test)),
               expected=essentialGraph(dag=~b | a + c | b, node.names=names(dist.test)))

  expect_equal(object=  essentialGraph(dag=~a | b + b | c, node.names=names(dist.test)),
               expected=essentialGraph(dag=~b | a + c | b, node.names=names(dist.test)))

  # more complex
  dist.test <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="gaussian", f="gaussian")
  # examples from 'Bayesian Networks in Bioinformatics, Kyu-Baek Hwang'

  minimal.dag <- matrix(data=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                               1, 0, 0, 0, 1, 0), nrow=6, byrow=TRUE)
  colnames(minimal.dag) <- rownames(minimal.dag) <- names(dist.test)
  completed.dag <- matrix(data=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                 1, 0, 0, 0, 1, 0), nrow=6, byrow=TRUE)
  colnames(completed.dag) <- rownames(completed.dag) <- names(dist.test)

  expect_equal(essentialGraph(dag=~f | e:a + e | c + c | a:b + d | b, node.names=names(dist.test),
                              PDAG="minimal"), minimal.dag)
  expect_equal(essentialGraph(dag=~f | e:a + e | c + c | a:b + d | b, node.names=names(dist.test),
                              PDAG="completed"),  completed.dag)
})

test_that("scoreContribution() works", {
  suppressMessages({
    suppressWarnings({
      mydat <- ex1.dag.data[,c("b1","g1","p1")]
      ## take a subset of cols

      ## setup distribution list for each node
      mydists <- list(b1="binomial",
                      g1="gaussian",
                      p1="poisson"
      )

      ## now build cache
      mycache <- buildScoreCache(data.df=mydat,data.dists=mydists,max.parents=1, method="mle")

      #now find the globally best DAG
      mp.dag <- mostProbable(score.cache=mycache, score="bic", verbose=FALSE)

      out <- scoreContribution(object=mp.dag)

      out.fit <- fitAbn(object=mp.dag, method="mle")

      expect_equal(unname(colSums(out$mlik))/1000, unname(unlist(out.fit$mliknode))/1000,tolerance=1e-2)
    })
  })
})
