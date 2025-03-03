test_that("plotAbn() works", {

  #Define distribution list
  dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="binomial", f="multinomial")

  #Plot from a formula and markov blanket for multinomial node
  expect_no_error(plotAbn(dag=~a|b:c:e+b|c:d:f+e|f, markov.blanket.node="b", data.dists =dist))

  # Gaussian
  N <- 1000
  mydists <- list(a="gaussian",
                  b="gaussian")
  a <- rnorm(n = N, mean = 0, sd = 1)
  b <- 2 + a + rnorm(n = N, mean = 5, sd = 1)
  mydf <- data.frame("a" = a,
                     "b" = b)
  mycache.mle <- buildScoreCache(data.df = mydf,
                                 data.dists = mydists,
                                 method = "mle",
                                 max.parents = 1)
  dag <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
  expect_silent({
    plt_dag <- plotAbn(dag=dag, data.dists = mydists)
    })
  expect_s4_class(plt_dag, class = "graph")
  expect_s4_class(plt_dag, class = "graphAM")
  expect_s4_class(plt_dag, class = "graphBase")

  data.param.adjMat <- matrix(data=c(0, 0, 1, 0), nrow=2L, byrow=TRUE, dimnames = list(NULL, c("a", "b")))
  # expect_equal(plt_dag@adjMat, data.param.adjMat) # TODO: For some reason, this does not work on CRAN
  expect_equal(plt_dag@graphData$edgemode, "directed")

  # Having an output
  expect_snapshot_output(plotAbn(dag=dag, plot = TRUE))
})

# construct example
dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="binomial", f="binomial")
edge.strength <- matrix(c(0,0.5,0.5,0.7,0.1,0,   #Define a matrix formulation
                          0,0,0.3,0.1,0,0.8,
                          0,0,0,0.35,0.66,0,
                          0,0,0,0,0.9,0,
                          0,0,0,0,0,0.8,
                          0,0,0,0,0,0),nrow = 6L, ncol = 6L, byrow = TRUE)
colnames(edge.strength) <- rownames(edge.strength) <- names(dist)  #Naming of the matrix
data <- data.frame(a=rnorm(100), b=rnorm(100), c=rnorm(100), d=rnorm(100), e=rbinom(100, 1, 0.5), f=rbinom(100, 1, 0.5))  #Generate random data

test_that("Plot from a formula works", {
  dag.mat <- formula_abn(f = ~a|b:c:e+b|c:d:f+e|f, name = letters[1:6])
  rownames(dag.mat) <- NULL
  plt <- plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dist = dist, node.fillcolor.list= "e")
  expect_equal(plt@adjMat, dag.mat)
  expect_equal(plt@renderInfo@nodes$fill[["e"]], "brown3")

  # fill color
  plt <- plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dist = dist, node.fillcolor.list= "e", markov.blanket.node = "b")
  expect_equal(plt@renderInfo@nodes$fill[["b"]], "brown3")

  # shapes
  plt <- plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dist = dist, node.shape='diamond')
  expect_equal(unname(plt@renderInfo@nodes$shape), rep("diamond", 6))

  plt <- plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dist = dist, node.shape=c('diamond','box'))
  shapes_by_dists <- c()
  for (i in 1:length(unlist(dist))) {
    if(unlist(dist)[i] == "gaussian"){
      shapes_by_dists[i] <- "diamond"
    } else {
      shapes_by_dists[i] <- "box"
    }
  }
  expect_equal(unname(plt@renderInfo@nodes$shape), shapes_by_dists)

  # Markov blanket
  plt <- plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, markov.blanket.node = "e", data.df = data)
  plt_fill <- plt@renderInfo@nodes$fill
  true_mb <- mb(dag= ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, node = "e", data.df = data)
  expect_equal(length(unique(plt_fill[true_mb])), 1) # check for one color among those nodes from true_mb
  idx <- which(plt_fill %in% plt_fill[true_mb]) # get indexes from true_mb nodes
  expect_false(any(names(plt_fill[-idx]) %in% true_mb)) # None of the other nodes (that are in true_mb) should be in true_mb

  plt <- plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, markov.blanket.node = "c", data.df = data)
  plt_fill <- plt@renderInfo@nodes$fill
  true_mb <- mb(dag= ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, node = "c", data.df = data)
  expect_equal(length(unique(plt_fill[true_mb])), 1) # check for one color among those nodes from true_mb
  idx <- which(plt_fill %in% plt_fill[true_mb]) # get indexes from true_mb nodes
  expect_false(any(names(plt_fill[-idx]) %in% true_mb)) # None of the other nodes (that are in true_mb) should be in true_mb

  plt <- plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, markov.blanket.node = c("d"), data.df = data)
  plt_fill <- plt@renderInfo@nodes$fill
  true_mb <- mb(dag= ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, node = "d", data.df = data)
  expect_equal(length(unique(plt_fill[true_mb])), 1) # check for one color among those nodes from true_mb
  idx <- which(plt_fill %in% plt_fill[true_mb]) # get indexes from true_mb nodes
  expect_false(any(names(plt_fill[-idx]) %in% true_mb)) # None of the other nodes (that are in true_mb) should be in true_mb

  plt <- plotAbn(dag = ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, markov.blanket.node = c("d","f"), data.df = data)
  plt_fill <- plt@renderInfo@nodes$fill
  true_mb <- mb(dag= ~a|b:c:e+b|c:d:f+e|f, data.dists = dist, node = c("d","f"), data.df = data)
  expect_equal(length(unique(plt_fill[true_mb])), 2) # check for one color among those nodes from true_mb
  idx <- which(plt_fill %in% plt_fill[true_mb]) # get indexes from true_mb nodes
  expect_false(any(names(plt_fill[-idx]) %in% true_mb)) # None of the other nodes (that are in true_mb) should be in true_mb

  # edge.strength
  expect_no_error({
    plt <- plotAbn(edge.strength, edge.strength=edge.strength, data.dist = dist) # TODO: think about a better test for edgestrength
    })
  vals <- c(t(edge.strength))
  expect_error(plotAbn(edge.strength, edge.strength=edge.strength, data.dist = dist,
                       fitted.values=vals[vals>0]), "argument is of length zero") # TODO: consider to make this error more meaningful.
  tmp <- edge.strength
  tmp[1,6] <- 1
  expect_error(plotAbn(tmp, edge.strength=-edge.strength-4, data.dist = dist),
               "'edge.strength' should be positive")
  expect_no_error(plotAbn(tmp, edge.strength=edge.strength, data.dist = dist)) # zeros allowed
  expect_error(plotAbn(edge.strength, edge.strength=tmp, data.dist = dist),"'edge.strength' does not match dag")

  # classic view
  skip("Plot classic view is deprecated.")
  tmp <- plotAbn(dag = edge.strength, edge.strength=edge.strength, data.dist = dist)
  slot(tmp, "renderInfo", check=FALSE) <- NULL
  plot(tmp)
  # renderGraph(tmp)
})

test_that("Plot from fitAbn() works", {
  skip_on_cran() # Skipped on CRAN because it requires the INLA package

  mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")]
  mydists <- list(b1="binomial", b2="binomial", b3="binomial", g1="gaussian",
                  b4="binomial", p2="poisson", p4="poisson")
  mydag.empty <- matrix(0, nrow=7, ncol=7)
  colnames(mydag.empty) <- rownames(mydag.empty) <- names(mydat)

  myres <- fitAbn(dag = ~b1|b2+b2|p4+g1+g1|p2+b3|g1+b4|b1+p4|g1, data.df = mydat, data.dists = mydists)
  expect_no_error({
    g <- plotAbn(myres$abnDag$dag, fitted.values = myres$modes, data.dists = mydists, edge.direction = 'pc')
  })

  # check if all edges have labels
  expect_false(any(unique(g@renderInfo@edges$label) %in% c(" "))) # TODO: consider to relax to only some edges have labels.

  myres <- fitAbn(dag = ~b1|b2+b2|p4+g1+g1|p2+b3|g1+b4|b1:g1+p4|g1:b3:p2, data.df = mydat, data.dists = mydists)
  expect_no_error({
    g <- plotAbn(myres$abnDag$dag, fitted.values = myres$modes, data.dists = mydists, edge.direction = 'pc') # TODO: Think about more specific test
  })

  # check multinomial with "mle"
  mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")]
  mydat1 <- cbind(m1=as.factor(as.numeric(mydat[,1])*2-as.numeric(mydat[,2])),
                  mydat[3:7])
  mydists1 <- list(m1="multinomial", b3="binomial", g1="gaussian",
                   b4="binomial", p2="poisson", p4="poisson")
  myres1 <- fitAbn(dag = ~m1|b3:g1:b4+b3|p2:p4+p2|p4, data.df = mydat1, data.dists = mydists1, method='mle')
  expect_no_error(plotAbn(dag = ~m1|b3:g1:b4+b3|p2:p4+p2|p4, data.dist=mydists1, node.shape=rep('box',4)))
})
