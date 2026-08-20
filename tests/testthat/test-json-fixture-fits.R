test_that("JSON stage 1: g2b2c_data constructs compact mixed MLE fit", {
  suppressMessages({
    suppressWarnings({
      expect_no_error(fit <- json_fixture_fit_g2b2c_mle())
    })
  })

  dists <- list(
    G1 = "gaussian",
    B1 = "binomial",
    B2 = "binomial",
    C = "multinomial",
    G2 = "gaussian"
  )
  json_fixture_assert_abnfit(fit, method = "mle", dists = dists,
                             multinomial = TRUE)
})

test_that("JSON stage 1: g2pbcgrp constructs grouped mixed MLE fit", {
  suppressMessages({
    suppressWarnings({
      expect_no_error(fit <- json_fixture_fit_g2pbcgrp_mle_grouped())
    })
  })

  dists <- list(
    G1 = "gaussian",
    P = "poisson",
    B = "binomial",
    C = "multinomial",
    G2 = "gaussian"
  )
  json_fixture_assert_abnfit(fit, method = "mle", dists = dists,
                             grouped = TRUE, multinomial = TRUE)
})

test_that("JSON stage 1: ex1.dag.data constructs classic mixed MLE fit", {
  suppressMessages({
    suppressWarnings({
      expect_no_error(fit <- json_fixture_fit_ex1_mle())
    })
  })

  dists <- list(
    b1 = "binomial",
    p1 = "poisson",
    g1 = "gaussian",
    b2 = "binomial",
    p2 = "poisson",
    b3 = "binomial",
    g2 = "gaussian"
  )
  json_fixture_assert_abnfit(fit, method = "mle", dists = dists)
})

test_that("JSON stage 1: ex1.dag.data constructs supported Bayesian fit", {
  suppressMessages({
    suppressWarnings({
      expect_no_error(fit <- json_fixture_fit_ex1_bayes())
    })
  })

  dists <- list(
    b1 = "binomial",
    p1 = "poisson",
    g1 = "gaussian",
    p2 = "poisson"
  )
  json_fixture_assert_abnfit(fit, method = "bayes", dists = dists)
})

test_that("JSON stage 1: ex3.dag.data constructs grouped binary MLE fit", {
  suppressMessages({
    suppressWarnings({
      expect_no_error(fit <- json_fixture_fit_ex3_mle_grouped())
    })
  })

  dists <- list(b1 = "binomial", b2 = "binomial")
  json_fixture_assert_abnfit(fit, method = "mle", dists = dists,
                             grouped = TRUE)
})

test_that("JSON stage 1: adg constructs real grouped mixed MLE fit", {
  suppressMessages({
    suppressWarnings({
      expect_no_error(fit <- json_fixture_fit_adg_mle_grouped())
    })
  })

  dists <- list(
    AR = "binomial",
    eggs = "binomial",
    wormCount = "poisson",
    age = "gaussian",
    adg = "gaussian"
  )
  json_fixture_assert_abnfit(fit, method = "mle", dists = dists,
                             grouped = TRUE)
})

test_that("JSON stage 1: FCV constructs real mixed multinomial MLE fit", {
  suppressMessages({
    suppressWarnings({
      expect_no_error(fit <- json_fixture_fit_fcv_mle())
    })
  })

  dists <- list(
    Outdoor = "binomial",
    Sex = "multinomial",
    GroupSize = "poisson",
    Age = "gaussian"
  )
  json_fixture_assert_abnfit(fit, method = "mle", dists = dists,
                             multinomial = TRUE)
})

test_that("JSON stage 1: Bayesian multinomial remains intentionally unsupported", {
  dists <- list(G1 = "gaussian", C = "multinomial")
  data.df <- droplevels(g2b2c_data[, names(dists)])
  dag <- json_fixture_empty_dag(dists)
  dag["C", "G1"] <- 1

  expect_error(
    fitAbn(dag = dag, data.df = data.df, data.dists = dists, method = "bayes"),
    regexp = "Multinomial nodes are not yet implemented"
  )
})
