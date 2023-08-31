test_that("toGraphviz() is backwards compatible.", {
  ## Subset of a build-in dataset
  mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")]

  ## setup distribution list for each node
  mydists <- list(b1="binomial", b2="binomial", b3="binomial",
                  g1="gaussian", b4="binomial", p2="poisson",
                  p4="poisson")
  ## specify DAG model
  mydag <- matrix(c(   0,1,0,0,1,0,0, #
                       0,0,0,0,0,0,0, #
                       0,1,0,0,1,0,0, #
                       1,0,0,0,0,0,1, #
                       0,0,0,0,0,0,0, #
                       0,0,0,1,0,0,0, #
                       0,0,0,0,1,0,0  #
  ), byrow=TRUE, ncol=7)

  colnames(mydag) <- rownames(mydag) <- names(mydat)

  ## create file for processing with graphviz
  # outfile <- "graph1.dot"
  outfile <- tempfile(fileext = ".dot")
  expect_no_error({
    toGraphviz(dag=mydag, data.df=mydat, data.dists=mydists, outfile=outfile)
  })
  expect_snapshot_file(path = outfile, name = "graph1.dot")

  ## check if works without a specific outfile
  expect_no_error({
    toGraphviz(dag=mydag, data.df=mydat, data.dists=mydists, outfile=NULL)
  })

  ## Example using data with a group variable  where b1<-b2
  mydag <- matrix(c(0,1, 0,0), byrow=TRUE, ncol=2)

  colnames(mydag) <- rownames(mydag) <- names(ex3.dag.data[,c(1,2)])
  ## specific distributions
  mydists <- list(b1="binomial", b2="binomial")

  ## create file for processing with graphviz
  outfile <- tempfile(fileext = ".dot")
  # outfile <- "graph_group.dot"
  expect_no_error({
    toGraphviz(dag=mydag, data.df=ex3.dag.data[,c(1,2,14)], data.dists=mydists,
               group.var="group",
               outfile=outfile, directed=FALSE)
  })
  expect_snapshot_file(outfile, "graph_group.dot")
})

test_that("toGraphviz() works with all distributions.", {
  capture.output({
    suppressMessages({

      df <- FCV[, c(11:15)]
      mydists <- list(Pedigree="binomial",
                      Outdoor="binomial",
                      Sex="multinomial",
                      GroupSize="poisson",
                      Age="gaussian")
      mydists <- mydists[-1] # remove grouping variable from distribution list

      ## buildScoreCache -> mostProbable() -> fitAbn()
      expect_no_error({
        suppressWarnings({
          mycache.mle <- buildScoreCache(data.df = df, data.dists = mydists, method = "mle",
                                         group.var = "Pedigree",
                                         adj.vars = NULL, cor.vars = NULL,
                                         dag.banned = NULL, dag.retained = NULL,
                                         max.parents = 3,
                                         which.nodes = NULL, defn.res = NULL)
        }) # ignore non-convergence warnings
      })
      expect_no_error({
        mp.dag.mle <- mostProbable(score.cache = mycache.mle, verbose = FALSE)
      })
      expect_no_error({
        suppressWarnings({
          myres.mle <- fitAbn(object = mp.dag.mle, method = "mle")
        })
      })
    })
  })

  # Check if works as it should
  # ...with dag as "abnFit"
  expect_no_error({
    toGraphviz(dag=myres.mle,
               # data.df=df,
               # data.dists=mydists,
               # group.var = "Pedigree",
               outfile=NULL)
  })
  # ...with dag as matrix
  expect_no_error({
    toGraphviz(dag=myres.mle$abnDag$dag,
               data.df=df,
               data.dists=myres.mle$abnDag$data.dists,
               group.var = myres.mle$group.var,
               outfile=NULL)
  })
  # ...with dag as "abnLearned"
  expect_no_error({
    toGraphviz(dag=mp.dag.mle,
               data.df=df,
               data.dists=mydists,
               group.var = "Pedigree",
               outfile=NULL)
  })

  # Check if error if missing group.var
  # ...with dag as matrix
  expect_error({
    expect_warning({
      toGraphviz(dag=myres.mle$abnDag$dag,
                 data.df=df,
                 data.dists=myres.mle$abnDag$data.dists,
                 group.var = NULL,
                 outfile=NULL)
    })
  })
  # ...with dag as "abnLearned"
  expect_error({
    toGraphviz(dag=mp.dag.mle,
               data.df=df,
               data.dists=mydists,
               group.var = NULL,
               outfile=NULL)
  })

  # Check if warning if dag as "abnFit" for variables that are being ignored in this case
  expect_warning({
    toGraphviz(dag=myres.mle,
               data.df=df,
               data.dists=mydists,
               group.var = "Pedigree",
               outfile=NULL)
  })

  # Check if error if unknown dists provided
  wrong_dists <- myres.mle$abnDag$data.dists
  wrong_dists[[1]] <- "foo"
  # ...with dag as matrix
  expect_error({
    toGraphviz(dag=myres.mle$abnDag$dag,
               data.df=df,
               data.dists=wrong_dists,
               group.var = myres.mle$group.var,
               outfile=NULL)
  })
  # ...with dag as "abnLearned"
  expect_error({
    toGraphviz(dag=mp.dag.mle,
               data.df=df,
               data.dists=wrong_dists,
               group.var = "Pedigree",
               outfile=NULL)
  })

  # check if works with a specific outfile
  outfile <- tempfile(fileext = ".dot")
  # ...with dag as "abnFit"
  expect_no_error({
    toGraphviz(dag=myres.mle,
               # data.df=df,
               # data.dists=mydists,
               # group.var = "Pedigree",
               outfile=outfile)
  })
  expect_snapshot_file(outfile, "graph.dot")
  # ...with dag as matrix
  outfile <- tempfile(fileext = ".dot")
  expect_no_error({
    toGraphviz(dag=myres.mle$abnDag$dag,
               data.df=df,
               data.dists=myres.mle$abnDag$data.dists,
               group.var = myres.mle$group.var,
               outfile=outfile)
  })
  expect_snapshot_file(outfile, "graph.dot")
  # ...with dag as "abnLearned"
  outfile <- tempfile(fileext = ".dot")
  expect_no_error({
    toGraphviz(dag=mp.dag.mle,
               data.df=df,
               data.dists=mydists,
               group.var = "Pedigree",
               outfile=outfile)
  })
  expect_snapshot_file(outfile, "graph.dot")
})
