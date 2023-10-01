test_that("AbnDag class works", {
  set.seed(1234)
  df <- data.frame(A=as.factor(rbinom(20, 1, 0.6)),
                   B=rpois(20, 0.5),
                   C=rnorm(20),
                   D=as.factor(sample(c("x", "y", "z"), 20, replace = T)))
  dists <- list(A="binomial", B="poisson", C="gaussian", D="multinomial")
  d <- matrix(data=0, nrow=4, ncol=4)
  d[1, ] <- c(0, 1, 1, 1)
  colnames(d) <- rownames(d) <- names(df)

  expect_error(createAbnDag(dag = NULL, data.df = df, data.dists = NULL))
  expect_no_error(createAbnDag(dag = d, data.df = df, data.dists = NULL))
  expect_no_error(createAbnDag(dag = d, data.df = NULL, data.dists = NULL))
  expect_no_error({
    abn_d <- createAbnDag(dag = d, data.df = df, data.dists = dists)
  })
  expect_true(inherits(abn_d, "abnDag"))
  expect_equal(abn_d[["dag"]], d)
  expect_equal(abn_d[["data.df"]], df)
  expect_equal(abn_d[["data.dists"]], dists)
})

test_that("AbnDag class validates correctly", {
  set.seed(1234)
  df <- data.frame(A=as.factor(rbinom(20, 1, 0.6)),
                   B=rpois(20, 0.5),
                   C=rnorm(20),
                   D=as.factor(sample(c("x", "y", "z"), 20, replace = T)))
  dists <- list(A="binomial", B="poisson", C="gaussian", D="multinomial")
  d <- matrix(data=0, nrow=4, ncol=4)
  d[1, ] <- c(0, 1, 1, 1)
  colnames(d) <- rownames(d) <- names(df)

  abn_d <- createAbnDag(dag = d, data.df = df, data.dists = dists)

  expect_error(validate_abnDag(dag = NULL, data.df = NULL, returnDag = TRUE))
  expect_no_error(validate_abnDag(dag = abn_d, data.df = NULL, returnDag = TRUE))
  expect_no_error({
    validabn_d <- validate_abnDag(dag = abn_d, data.df = df, returnDag = TRUE)
  })
  expect_true(validate_abnDag(dag = abn_d, data.df = df, returnDag = FALSE))
  expect_true(validate_abnDag(dag =  ~A|B+A|C+A|D, data.df = df, returnDag = FALSE))
  expect_equal(validate_abnDag(dag =  ~A|B+A|C+A|D, data.df = df, returnDag = TRUE),
               formula_abn(f = ~A|B+A|C+A|D, name = c("A", "B", "C", "D"))
  )
  expect_equal(validabn_d, d)
})

test_that("Distributions are validated correctly", {
  dists <- list(A="binomial", B="poisson", C="gaussian", D="multinomial")

  expect_equal(validate_dists(data.dists = dists, returnDists = TRUE), dists)
  expect_true(validate_dists(data.dists = dists, returnDists = FALSE))
  expect_error(validate_dists(data.dists = NULL))
  expect_error(validate_dists(data.dists = list(A="foo"), returnDists = TRUE))
})
