test_that("build.control() works fine", {
  # basics
  expect_no_error(
    bctrl <- build.control()
  )
  expect_equal(class(bctrl), class(list()))
  expect_no_error(
    bctrl.bayes <- build.control(method = "bayes")
  )
  expect_equal(bctrl, bctrl.bayes)

  expect_no_error(
    bctrl.mle <- build.control(method = "mle")
  )
  expect_equal(class(bctrl), class(list()))
  expect_equal(class(bctrl.bayes), class(list()))
  expect_equal(class(bctrl.mle), class(list()))

  # Arguments
  expect_no_error(
    build.control(method = "bayes", epsilon = 1e-06) # no error if argument and method don't correspond
  )
})

test_that("buildScoreCache() checks work fine", {
  df <- airquality[complete.cases(airquality), ]

  # distribution (gaussian)
  dists <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
  names(dists) <- colnames(df)

  # Check arguments
  suppressWarnings({
    suppressMessages({

      expect_error(buildScoreCache(data.df=NULL, data.dists=dists))
      expect_error(buildScoreCache(data.df=df, data.dists=NULL))
      expect_error(buildScoreCache(data.df=df, data.dists=dists[1:3]))
      expect_error(buildScoreCache(data.df=df, data.dists=dists, group.var = NULL, cor.vars = c("Ozone")))

      expect_no_error(buildScoreCache(data.df=df, data.dists=dists))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "Bayes"))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "bayes"))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "MLE"))
      expect_error(buildScoreCache(data.df=df, data.dists=dists, method = "foo"), regexp = "unknown")
      expect_error(buildScoreCache(data.df=df, data.dists=dists, method = NA), regexp = "is NA")
      expect_error(buildScoreCache(data.df=df, data.dists=dists, method = NULL), regexp = "not provided")

      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, dag.banned = ~ Ozone | Solar.R))

      ban <- matrix(0, nrow = dim(df)[2], ncol = dim(df)[2], dimnames = list(names(df), names(df)))
      ban[1,2] <- ban[2,1] <- 1
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "bayes", dag.banned = ban))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "mle", dag.banned = ban))

      reta <- matrix(0, nrow = dim(df)[2], ncol = dim(df)[2], dimnames = list(names(df), names(df)))
      # ban and retain can not be on the same edge
      reta[1,2] <- reta[2,3] <- 1
      expect_error(buildScoreCache(data.df=df, data.dists=dists, method = "bayes", dag.banned = ban, dag.retained = reta))
      expect_error(buildScoreCache(data.df=df, data.dists=dists, method = "mle", dag.banned = ban, dag.retained = reta))
      reta[1,2] <- 0
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "bayes", dag.banned = ban, dag.retained = reta))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "mle", dag.banned = ban, dag.retained = reta))

      max.parents.faulty <- list(Ozone=6, Solar.R=6, Wind=6, Temp=6, Month=6, Day=6)
      expect_warning(buildScoreCache(data.df=df, data.dists=dists, method = "bayes", max.parents = max.parents.faulty), regexp = "555555")
      expect_warning(buildScoreCache(data.df=df, data.dists=dists, method = "mle", max.parents = max.parents.faulty), regexp = "555555")
      max.parents <- list(Ozone=5, Solar.R=5, Wind=5, Temp=5, Month=5, Day=5)
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "bayes", max.parents = max.parents))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "mle", max.parents = max.parents))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "bayes", max.parents = 5))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "mle", max.parents = 5))

      # Test compatability of max.parents and dag.retain
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "bayes", max.parents = max.parents, dag.banned = ban, dag.retained = reta))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "mle", max.parents = max.parents, dag.banned = ban, dag.retained = reta))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "bayes", max.parents = 5, dag.banned = ban, dag.retained = reta))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "mle", max.parents = 5, dag.banned = ban, dag.retained = reta))
      reta[2,4] <- 1
      max.parents.nope <- list(Ozone=1, Solar.R=1, Wind=1, Temp=1, Month=1, Day=1)
      expect_error(buildScoreCache(data.df=df, data.dists=dists, method = "bayes", max.parents = max.parents.nope, dag.banned = ban, dag.retained = reta)) # 2 nodes restricted but only 1 max parent allowed should throw error.
      expect_error(buildScoreCache(data.df=df, data.dists=dists, method = "mle", max.parents = max.parents.nope, dag.banned = ban, dag.retained = reta))

      # check group.var argument
      group <- sample(as.factor(c(1,2)), size = nrow(df), prob = c(0.3, 0.7), replace = TRUE)
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "bayes", group.var = NULL))
      expect_no_error(buildScoreCache(data.df=df, data.dists=dists, method = "mle", group.var = NULL))

      if(!testthat:::on_cran()) {
        if(requireNamespace("INLA", quietly = TRUE)){
          expect_no_error(buildScoreCache(data.df=data.frame(df[,1:2], group=group), data.dists=dists[1:2], method = "bayes", group.var = "group")) # taking only the first two columns to increase performance
        }
      }
      expect_no_error(buildScoreCache(data.df=data.frame(df[,1:2], group=group), data.dists=dists[1:2], method = "mle", group.var = "group"))

      expect_error(buildScoreCache(data.df=data.frame(df[,1:2], group1=group, group2=group), data.dists=dists[1:2], method = "bayes", group.var = c("group1", "group2"))) # only one grouping variable is implemented
      expect_error(buildScoreCache(data.df=data.frame(df[,1:2], group1=group, group2=group), data.dists=dists[1:2], method = "mle", group.var = c("group1", "group2")))

      # check control argument
      expect_error(buildScoreCache(data.df=data.frame(df[,1:2], group=group), data.dists=dists[1:2], method = "bayes", group.var = "group",
                                   control= list(max.mode.error=101)))
      if(!testthat:::on_cran()) {
        if(requireNamespace("INLA", quietly = TRUE)){
          expect_no_error(buildScoreCache(data.df=data.frame(df[,1:2], group=group), data.dists=dists[1:2], method = "bayes", group.var = "group",
                                          control= list(max.mode.error=100)))
        }
      }
      # increase convergence tolerances for GLMMs
      if(!testthat:::on_cran()) {
        if(requireNamespace("INLA", quietly = TRUE)){
          expect_no_error(buildScoreCache(data.df=data.frame(df[,1:2], group=group), data.dists=dists[1:2], method = "bayes", group.var = "group",
                                          control= list(xtol_abs=1e-8, ftol_abs=1e-8, epsilon=1e-10)))
        }
      }
      expect_no_error(buildScoreCache(data.df=data.frame(df[,1:2], group=group), data.dists=dists[1:2], method = "mle", group.var = "group",
                                      control= list(xtol_abs=1e-8, ftol_abs=1e-8, epsilon=1e-10)))
    })
  })
})

test_that("buildScoreCache()'s methods `bayes` and `mle` behave similarly", {
  df <- airquality[complete.cases(airquality), ]

  # distribution (gaussian)
  dist <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
  names(dist) <- colnames(df)



  expect_warning({
    mycache.mle <- buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=6)
  })
  expect_warning({
    mycache.bayes <- buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=6, dry.run=TRUE)
  })

  mycache.mle.1 <- buildScoreCache(data.df=df, data.dists=dist, method = "mle", max.parents=3)
  mycache.bayes.1 <- buildScoreCache(data.df=df, data.dists=dist, method = "bayes", max.parents=3, dry.run=TRUE)

  ## test

  ## dag retain
  mycache.mle.2 <- buildScoreCache(data.df=df, data.dists=dist, method = "mle",
                                   dag.banned=NULL, dag.retained=~Ozone|Solar.R, max.parents=3, dry.run=TRUE)
  mycache.bayes.2 <- buildScoreCache(data.df=df, data.dists=dist, method = "bayes",
                                     max.parents=3, dry.run=TRUE, dag.retained=~Ozone|Solar.R, dag.banned=NULL)

  expect_warning({
    mycache.mle.3 <- buildScoreCache(data.df=df, data.dists=dist, method = "mle",
                                     dag.banned=NULL, dag.retained=~Wind |., max.parents=6, dry.run=TRUE)
  })
  expect_warning({
    mycache.bayes.3 <- buildScoreCache(data.df=df, data.dists=dist, method = "bayes",
                                       max.parents=6, dry.run=TRUE, dag.retained=~Wind|., dag.banned=NULL)
  })


  ## dag ban
  mycache.mle.4 <- buildScoreCache(data.df=df, data.dists=dist, method = "mle",
                                   dag.banned=~Ozone | Solar.R:Wind, max.parents=3, dry.run=TRUE)
  mycache.bayes.4 <- buildScoreCache(data.df=df, data.dists=dist, method = "bayes",
                                     max.parents=3, dry.run=TRUE, dag.banned=~Ozone|Solar.R:Wind)

  ## test cache
  expect_equal(mycache.mle$children, mycache.bayes$children)
  expect_equal(mycache.mle$node.defn, mycache.bayes$node.defn)
  expect_equal(mycache.mle.1$children, mycache.bayes.1$children)
  expect_equal(mycache.mle.1$node.defn, mycache.bayes.1$node.defn)
  expect_equal(mycache.mle.2$children, mycache.bayes.2$children)
  expect_equal(mycache.mle.2$node.defn, mycache.bayes.2$node.defn)

  expect_equal(mycache.mle.3$children, mycache.bayes.3$children)
  expect_equal(mycache.mle.3$node.defn, mycache.bayes.3$node.defn)
  expect_equal(mycache.mle.4$children, mycache.bayes.4$children)
  expect_equal(mycache.mle.4$node.defn, mycache.bayes.4$node.defn)
})

test_that("buildScoreCache() is backward compatible", {
  suppressMessages({
    # test_unit <- function(){
    expect_no_error({
      mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")];## take a subset of cols

      ## setup distribution list for each node
      mydists <- list(b1="binomial",
                      b2="binomial",
                      g1="gaussian",
                      g2="gaussian",
                      b3="binomial",
                      g3="gaussian"
      );

      ban <- matrix(rep(0,dim(mydat)[2]^2),ncol=dim(mydat)[2]);# ban nothing
      colnames(ban) <- rownames(ban) <- names(mydat); #names must be set
      ban["b1","b2"] <- 1; # now ban arc from b2 to b1
      retain <- matrix(rep(0,dim(mydat)[2]^2),ncol=dim(mydat)[2]);# retain nothing
      colnames(retain) <- rownames(retain) <- names(mydat); #names must be set
      retain["g1","g3"] <- 1; # always retain arc from g3 to g1
      # parent limits
      max.par <- list("b1"=2,"b2"=2,"g1"=2,"g2"=0,"b3"=2,"g3"=3);

      ## now build cache of scores (goodness of fits for each node)

      res.c <- buildScoreCache(data.df=mydat,data.dists=mydists,
                               dag.banned=ban, dag.retained=retain,max.parents=max.par
      );

      ## repeat but using R-INLA. The mlik's should be virtually identical.
      ## now build cache

      if(!testthat:::on_cran()) {
        if(requireNamespace("INLA", quietly = TRUE)){
          res.inla <- buildScoreCache(data.df=mydat,data.dists=mydists,
                                      dag.banned=ban, dag.retained=retain,max.parents=max.par,
                                      control=list(max.mode.error=100));
        }
      }
    })

    expect_no_error({
      #################################################################
      ## Example 2 - much bigger problem using glm - make take a while
      #################################################################

      mydat <- ex2.dag.data;## this data comes with abn see ?ex2.dag.data

      ## setup distribution list for each node
      mydists <- list(b1="binomial",
                      g1="gaussian",
                      p1="poisson",
                      b2="binomial",
                      g2="gaussian",
                      p2="poisson",
                      b3="binomial",
                      g3="gaussian",
                      p3="poisson",
                      b4="binomial",
                      g4="gaussian",
                      p4="poisson",
                      b5="binomial",
                      g5="gaussian",
                      p5="poisson",
                      b6="binomial",
                      g6="gaussian",
                      p6="poisson"
      );

      ## parent limits
      max.par <- list("b1"=2,"g1"=2,"p1"=2,"b2"=2,"g2"=2,"p2"=2,"b3"=2,
                      "g3"=2,"p3"=2,"b4"=2,"g4"=2,
                      "p4"=2,"b5"=2,"g5"=2,"p5"=2,"b6"=2,"g6"=2,"p6"=2);

      ## no explicit ban or retain restrictions set so dont need to supply ban
      ##  or retain matrices

      ## now build cache using internal code just for nodes 1,2 and 3
      ## e.g. "b1", "p1" and "g1"
      mycache.c <- buildScoreCache(data.df=mydat,data.dists=mydists,
                                   max.parents=max.par[1:3], which.nodes=c(1:3));
    })

    expect_no_error({
      ###################################################################
      ## Example 3 - grouped data - random effects example e.g. glmm
      ###################################################################

      mydat <- ex3.dag.data;## this data comes with abn see ?ex3.dag.data

      mydists <- list(b1="binomial",
                      b2="binomial",
                      b3="binomial",
                      b4="binomial",
                      b5="binomial",
                      b6="binomial",
                      b7="binomial",
                      b8="binomial",
                      b9="binomial",
                      b10="binomial",
                      b11="binomial",
                      b12="binomial",
                      b13="binomial"
      );
      max.par <- 1;

      ## in this example INLA is used as default since these are glmm nodes
      ## when running this at node-parent combination 71 the default accuracy check on the
      ## INLA modes is exceeded (default is a max. of 10 percent difference from
      ## modes estimated using internal code) and a message is given that internal code
      ## will be used in place of INLA's results.

      # which.nodes doesn't work yet in combination with group.var
      # mycache <- buildScoreCache(data.df=mydat,data.dists=mydists,group.var="group",
      #                          cor.vars=c("b1","b2"),
      #                          max.parents=max.par, which.nodes=c(1));
      if(!testthat:::on_cran()) {
        if(requireNamespace("INLA", quietly = TRUE)){
          mycache <- buildScoreCache(data.df=mydat,data.dists=mydists,group.var="group",
                                     cor.vars=c("b1","b2"), method = "bayes",
                                     max.parents=max.par);
        }
      }
      # mle ignores group.var
      mycache <- buildScoreCache(data.df=mydat,data.dists=mydists,group.var="group",
                                 cor.vars=c("b1","b2"), method = "mle",
                                 max.parents=max.par);
    })

    expect_no_error({
      mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")] ## take a subset of cols

      ## setup distribution list for each node
      mydists <- list(b1="binomial",
                      b2="binomial",
                      g1="gaussian",
                      g2="gaussian",
                      b3="binomial",
                      g3="gaussian")

      ## now build cache of scores (goodness of fits for each node)
      res.mle <- buildScoreCache(data.df=mydat,data.dists=mydists,max.parents=3,method="mle")
      res.abn <- buildScoreCache(data.df=mydat,data.dists=mydists,max.parents=3,method="bayes")
    })
  })
})

test_that("buildScoreCache() works with all distributions", {
  df <- FCV[, c(11:15)]
  mydists <- list(Pedigree="binomial",
                  Outdoor="binomial",
                  Sex="multinomial",
                  GroupSize="poisson",
                  Age="gaussian")

  # Simple no grouping
  suppressMessages({
    if(.Platform$OS.type == "unix") {
      capture.output({
        expect_error({
          mycache.bayes <- buildScoreCache(data.df = df, data.dists = mydists, method = "bayes",
                                           group.var = NULL, adj.vars = NULL, cor.vars = NULL,
                                           dag.banned = NULL, dag.retained = NULL,
                                           max.parents = 1,
                                           which.nodes = NULL, defn.res = NULL)
        }, regexp = "not yet implemented") # Multinomial nodes are not yet implemented for method 'bayes'. Try with method='mle'.

        expect_no_error({
          mycache.mle <- buildScoreCache(data.df = df, data.dists = mydists, method = "mle",
                                         group.var = NULL, adj.vars = NULL, cor.vars = NULL,
                                         dag.banned = NULL, dag.retained = NULL,
                                         max.parents = 1,
                                         which.nodes = NULL, defn.res = NULL)
        })

        # With grouping
        mydists <- mydists[-1] # remove grouping variable from distribution list
        expect_error({
          mycache.bayes <- buildScoreCache(data.df = df, data.dists = mydists, method = "bayes",
                                           group.var = "Pedigree", adj.vars = NULL, cor.vars = NULL,
                                           dag.banned = NULL, dag.retained = NULL,
                                           max.parents = 1,
                                           which.nodes = NULL, defn.res = NULL)
        }, regexp = "not yet implemented") # Multinomial nodes are not yet implemented for method 'bayes'. Try with method='mle'.

        expect_no_error({
          suppressWarnings({
            mycache.mle <- buildScoreCache(data.df = df, data.dists = mydists, method = "mle",
                                           group.var = "Pedigree", adj.vars = NULL, cor.vars = NULL,
                                           dag.banned = NULL, dag.retained = NULL,
                                           max.parents = 1,
                                           which.nodes = NULL, defn.res = NULL)
          }) # ignore non-convergence warnings
        })
      },
      file = "/dev/null")
    }
  })
})

test_that("buildScoreCache() computes in parallel", {
  skip_on_cran() # workaround to not overconsume threads on CRAN. This is related to an issue reported for lme4 (https://github.com/lme4/lme4/issues/627)
  df <- FCV[, c(11:15)]
  mydists <- list(Pedigree="binomial",
                  Outdoor="binomial",
                  Sex="multinomial",
                  GroupSize="poisson",
                  Age="gaussian")
  ncores <- 2

  # Simple no grouping
  suppressWarnings({
    suppressMessages({
      if(.Platform$OS.type == "unix") {
        capture.output({
          expect_error({
            expect_warning({
              expect_warning({
                mycache.bayes <- buildScoreCache(data.df = df, data.dists = mydists, method = "bayes",
                                                 group.var = NULL, adj.vars = NULL, cor.vars = NULL,
                                                 dag.banned = NULL, dag.retained = NULL,
                                                 max.parents = 1,
                                                 which.nodes = NULL, defn.res = NULL,
                                                 control = list("ncores" = ncores))
              }, regexp = "Control parameters provided that are not used with method bayes are ignored")
            }, regexp = "Multithreading is currently only implemented for method")
          }, regexp = "Multinomial nodes are not yet implemented for method")

          expect_no_error({
            mycache.mle <- buildScoreCache(data.df = df, data.dists = mydists, method = "mle",
                                           group.var = NULL, adj.vars = NULL, cor.vars = NULL,
                                           dag.banned = NULL, dag.retained = NULL,
                                           max.parents = 1,
                                           which.nodes = NULL, defn.res = NULL,
                                           control = list("ncores" = ncores))
          })
        },
        file = "/dev/null")
      }
    })
  })
})
