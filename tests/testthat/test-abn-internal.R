test_that("Parsing DAG as formula works", {
  expected_output <- matrix(c(0,0,0,0,0,0,0,1,0), nrow = 3, dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  expect_equal(formula_abn(~ A + B | C, name = c("A", "B", "C")), expected_output)


  dist <- list(a="gaussian", b="gaussian", c="gaussian", d="gaussian", e="gaussian", f="gaussian")
  df.1 <- data.frame(a=1, b=1, c=1, d=1, e=1, f=1)
  m.formula.1 <- createAbnDag(dag=~a | b:c + b | c:d + a | e:f, data.dists=dist, data.df = df.1)$dag

  m.formula.2 <- createAbnDag(dag=~a | ., data.dists=dist, data.df = df.1)$dag
  m.true.1 <- matrix(data=c(0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0), nrow=6, ncol=6, byrow=TRUE)
  m.true.2 <- matrix(data=c(0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0), nrow=6, ncol=6, byrow=TRUE)
  colnames(m.true.1) <- rownames(m.true.1) <- colnames(m.true.2) <- rownames(m.true.2) <- names(dist)
  expect_equal(m.formula.1, m.true.1)
  expect_equal(m.formula.2, m.true.2)

  ## formula with real data
  df <- airquality[complete.cases(airquality), ]

  # distribution (gaussian)
  dist <- list(Ozone="gaussian", Solar.R="gaussian", Wind="gaussian", Temp="gaussian", Month="gaussian", Day="gaussian")
  names(dist) <- colnames(df)
  m.formula.1 <- createAbnDag(dag=~Ozone | Solar.R, data.dists=dist, data.df = df)$dag
  m.formula.2 <- createAbnDag(dag=~Solar.R | ., data.dists=dist, data.df = df)$dag
  m.true.1 <- matrix(data=c(0, 1, 0, 0, 0, 0, rep(0, 30)), nrow=6, ncol=6, byrow=TRUE)
  m.true.2 <- matrix(data=c(0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, rep(0, 24)), nrow=6, ncol=6, byrow=TRUE)
  colnames(m.true.1) <- rownames(m.true.1) <- colnames(m.true.2) <- rownames(m.true.2) <- names(dist)

  expect_equal(m.formula.1, m.true.1)
  expect_equal(m.formula.2, m.true.2)
})

test_that("Data is checked properly", {
  set.seed(1234)
  df <- data.frame(A=as.factor(rbinom(20, 1, 0.6)),
                   B=rpois(20, 0.5),
                   C=rnorm(20),
                   D=as.factor(sample(c("x", "y", "z"), 20, replace = T)))

  dists <- list(A="binomial", B="poisson", C="gaussian", D="multinomial")
  names(dists) <- colnames(df)

  expect_equal(check.valid.data(data.df = df, data.dists = dists, group.var = NULL),
               list(gaus=3, bin=1, pois=2, mult=4))

  # test empty factor levels
  df$D <- factor(df$D, levels=c("x", "y", "z", "w"))
  expect_error(check.valid.data(data.df = df, data.dists = dists, group.var = NULL))
})

test_that("DAG is checked properly", {
  set.seed(1234)
  df <- data.frame(A=as.factor(rbinom(20, 1, 0.6)),
                   B=rpois(20, 0.5),
                   C=rnorm(20),
                   D=as.factor(sample(c("x", "y", "z"), 20, replace = T)))

  d <- matrix(data=0, nrow=4, ncol=4)
  d[1, ] <- c(0, 1, 1, 1)
  colnames(d) <- rownames(d) <- names(df)

  expect_equal(check.valid.dag(dag = d, data.df = df, group.var = NULL), d)

  expect_error(check.valid.dag(dag = ~A))   # data.df missing
  expect_error(check.valid.dag(dag = ~a, data.df = c(a=1)))
  expect_equal(check.valid.dag(dag = ~a, data.df = data.frame(a=1)), matrix(0,1,1, dimnames=list("a","a")))

  expect_error(check.valid.dag(dag = ~a|b+b|a, data.df = data.frame(a = c(0,1,4,5,6), b=c(1,2,3,4,5))))  # 1-diag(2), cyclic

  m1 <- matrix( 0, 2,2, dimnames=list(c("a","b"), c("a","b")))
  expect_error(check.valid.dag(dag = m1[2:1,])) # rownames in reverse order to columnnames.
})

test_that("Parent nodes limit is checked properly", {
  set.seed(1234)
  df <- data.frame(A=as.factor(rbinom(20, 1, 0.6)),
                   B=rpois(20, 0.5),
                   C=rnorm(20),
                   D=as.factor(sample(c("x", "y", "z"), 20, replace = T)))

  expect_error(check.valid.parents(data.df = NULL, max.parents = NULL, group.var = NULL))
  expect_error(check.valid.parents(data.df = df, max.parents = list(A=NULL, B=1, C=2, D=0), group.var = NULL))
  expect_error(check.valid.parents(data.df = df, max.parents = dim(df)[2]+1, group.var = NULL))
  expect_equal(check.valid.parents(data.df = df, max.parents = NULL, group.var = NULL), rep(dim(df)[2], dim(df)[2]))
  expect_equal(check.valid.parents(data.df = df, max.parents = c(1,1,1,1)), c(1,1,1,1))
  expect_equal(check.valid.parents(data.df = df, max.parents = 0, group.var = NULL), c(0,0,0,0))
  expect_equal(check.valid.parents(data.df = df, max.parents = 1, group.var = NULL), check.valid.parents(data.df = df, max.parents = 1))
})

test_that("Node validitiy is assessed properly", {
  set.seed(1234)
  df <- data.frame(A=as.factor(rbinom(20, 1, 0.6)),
                   B=rpois(20, 0.5),
                   C=rnorm(20),
                   D=as.factor(sample(c("x", "y", "z"), 20, replace = T)))

  expect_error(check.which.valid.nodes(data.df = df, which.nodes = "A"))
  expect_error(check.which.valid.nodes(data.df = df, which.nodes = seq(1,dim(df)[2]+1))) # error if more nodes possible than variables in df
  expect_error(check.which.valid.nodes(data.df = df, which.nodes = -1))
  expect_no_error(check.which.valid.nodes(data.df = df, which.nodes = NULL)) # no error if group.var is not passed
  expect_no_error(check.which.valid.nodes(data.df = df, which.nodes = c(1,2,3,4)))

  expect_equal(check.which.valid.nodes(data.df = df, which.nodes = 1), 1)
})

test_that("Grouping variables are checked properly", {
  set.seed(1234)
  df <- data.frame(A=as.factor(rbinom(20, 1, 0.6)),
                   B=rpois(20, 0.5),
                   C=rnorm(20),
                   D=as.factor(sample(c("x", "y", "z"), 20, replace = T)),
                   G=as.factor(c(rep("A", 10), rep("B", 10))))

  expect_error(check.valid.groups(group.var = NULL, data.df = NULL, cor.vars = NULL))
  expect_error(check.valid.groups(group.var = 5, data.df = df, cor.vars = c("A", "B", "C", "D")))
  expect_error(check.valid.groups(group.var = NULL, data.df = df, cor.vars = c("A", "B", "C", "D")))
  expect_warning(check.valid.groups(group.var = "G", data.df = df, cor.vars = NULL, verbose = TRUE))

  checkout <- check.valid.groups(group.var = "G", data.df = df, cor.vars = c("A", "B", "C", "D"))
  expect_equal(checkout$data.df,df[, which(names(df) != "G")])
  expect_equal(checkout$grouped.vars, which(colnames(df[, which(names(df) != "G")]) %in% colnames(df)))
  expect_equal(checkout$group.ids, as.integer(df$G))
})

test_that("Distribution of variables is properly assessed.", {
  dists <- list(A="binomial", B="poisson", C="gaussian", D="multinomial")

  expect_equal(get.var.types(data.dists = dists), c(1,3,2,4))

  wrongdists <- dists
  wrongdists$E <- "thisisnodistribution"
  expect_error(get.var.types(wrongdists))
})

test_that("check.valid.buildControls() works properly", {
  verbose <- TRUE
  for (method in c("bayes", "mle")) {
    expect_error({
      check.valid.buildControls(control = "foo", method = method, verbose = verbose)
    }, regexp = "Control arguments must be provided as named list")
    expect_no_error({
      check.valid.buildControls(control = NULL, method = method, verbose = verbose)
    })

    expect_error({
      check.valid.buildControls(control = list(foo = "bar"), method = method, verbose = verbose)
    }, regexp = "Unknown control")

    expect_no_error({
      ctrl <- check.valid.buildControls(control = list(), method = method, verbose = verbose)
    })
    expect_equal(ctrl, build.control(method = method))

    expect_no_error({
      ctrl <- check.valid.buildControls(control = list(seed = 42L), method = method, verbose = verbose)
    })
    expect_equal(ctrl[["seed"]], 42L)
    expect_error({
      check.valid.buildControls(control = list(seed = -1), method = method, verbose = verbose)
    }, regexp = "must be a non-negative integer")

    expect_error({
      check.valid.buildControls(control = list(ncores = -10), method = method, verbose = verbose)
    })
    expect_message({
      check.valid.buildControls(control = list(ncores = 1e6), method = method, verbose = verbose)
    }, regexp = "Running in parallel with")
    expect_message({
      check.valid.buildControls(control = list(ncores = 0), method = method, verbose = verbose)
    }, regexp = "Running in single core mode")
    expect_message({
      check.valid.buildControls(control = list(ncores = -1), method = method, verbose = verbose)
    })

    expect_error({
      suppressWarnings(check.valid.buildControls(control = list(max.mode.error = -10), method = method, verbose = verbose))
    })
    expect_error({
      suppressWarnings(check.valid.buildControls(control = list(max.mode.error = 101), method = method, verbose = verbose))
    })
    expect_no_error({
      suppressWarnings(check.valid.buildControls(control = list(max.mode.error = 0), method = method, verbose = verbose))
      suppressWarnings(check.valid.buildControls(control = list(max.mode.error = 100), method = method, verbose = verbose))
    })
  }

  expect_error({
    check.valid.buildControls(control = list(), method = "foo")
  }, regexp = "unknown")

  expect_warning({
    ctrl <- check.valid.buildControls(control = list(epsilon = 1e-08), method = "bayes") # epsilon parameter is only used in mle
  }, regexp = "Control parameters provided that are not used with method")
  expect_equal(length(ctrl), length(build.control(method = "bayes")))
})

test_that("check.valid.fitControls() works properly", {
  verbose <- TRUE
  for (method in c("bayes", "mle")) {
    expect_error({
      check.valid.fitControls(control = "foo", method = method, verbose = verbose)
    }, regexp = "Control arguments must be provided as named list")
    expect_no_error({
      check.valid.fitControls(control = NULL, method = method, verbose = verbose)
    })

    expect_error({
      check.valid.fitControls(control = list(foo = "bar"), method = method, verbose = verbose)
    }, regexp = "Unknown control")

    expect_no_error({
      ctrl <- check.valid.fitControls(control = list(), method = method, verbose = verbose)
    })
    expect_equal(ctrl, fit.control(method = method))

    expect_no_error({
      ctrl <- check.valid.fitControls(control = list(seed = 42L), method = method, verbose = verbose)
    })
    expect_equal(ctrl[["seed"]], 42L)
    expect_error({
      check.valid.fitControls(control = list(seed = -1), method = method, verbose = verbose)
    }, regexp = "must be a non-negative integer")

    expect_error({
      check.valid.fitControls(control = list(ncores = -10), method = method, verbose = verbose)
    })
    expect_message({
      check.valid.fitControls(control = list(ncores = 1e6), method = method, verbose = verbose)
    }, regexp = "Running in parallel with")
    expect_message({
      check.valid.fitControls(control = list(ncores = 0), method = method, verbose = verbose)
    }, regexp = "Running in single core mode")
    expect_message({
      check.valid.fitControls(control = list(ncores = -1), method = method, verbose = verbose)
    })

    expect_error({
      suppressWarnings(check.valid.fitControls(control = list(max.mode.error = -10), method = method, verbose = verbose))
    })
    expect_error({
      suppressWarnings(check.valid.fitControls(control = list(max.mode.error = 101), method = method, verbose = verbose))
    })
    expect_no_error({
      suppressWarnings(check.valid.fitControls(control = list(max.mode.error = 0), method = method, verbose = verbose))
      suppressWarnings(check.valid.fitControls(control = list(max.mode.error = 100), method = method, verbose = verbose))
    })
  }

  expect_error({
    check.valid.fitControls(control = list(), method = "foo")
  }, regexp = "unknown")

  expect_warning({
    ctrl <- check.valid.fitControls(control = list(epsilon = 1e-08), method = "bayes") # epsilon parameter is only used in mle
  }, regexp = "Control parameters provided that are not used with method")
  expect_equal(length(ctrl), length(fit.control(method = "bayes")))
})
