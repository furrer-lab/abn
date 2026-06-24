test_that("JSON stage 2: g2b2c_data mixed MLE export is complete", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_g2b2c_mle()
      dists <- list(
        G1 = "gaussian",
        B1 = "binomial",
        B2 = "binomial",
        C = "multinomial",
        G2 = "gaussian"
      )
      json_fixture_assert_export(fit, dists = dists, multinomial = TRUE)
    })
  })
})

test_that("JSON stage 2: g2pbcgrp grouped mixed MLE export is complete", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_g2pbcgrp_mle_grouped()
      dists <- list(
        G1 = "gaussian",
        P = "poisson",
        B = "binomial",
        C = "multinomial",
        G2 = "gaussian"
      )
      json_fixture_assert_export(fit, dists = dists, grouped = TRUE,
                                 multinomial = TRUE)
    })
  })
})

test_that("JSON stage 2: ex1.dag.data classic mixed MLE export is complete", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_mle()
      dists <- list(
        b1 = "binomial",
        p1 = "poisson",
        g1 = "gaussian",
        b2 = "binomial",
        p2 = "poisson",
        b3 = "binomial",
        g2 = "gaussian"
      )
      json_fixture_assert_export(fit, dists = dists)
    })
  })
})

test_that("JSON stage 2: ex3.dag.data grouped binary MLE export is complete", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex3_mle_grouped()
      dists <- list(b1 = "binomial", b2 = "binomial")
      json_fixture_assert_export(fit, dists = dists, grouped = TRUE)
    })
  })
})

test_that("JSON stage 2: adg real grouped mixed MLE export is complete", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_adg_mle_grouped()
      dists <- list(
        AR = "binomial",
        eggs = "binomial",
        wormCount = "poisson",
        age = "gaussian",
        adg = "gaussian"
      )
      json_fixture_assert_export(fit, dists = dists, grouped = TRUE)
    })
  })
})

test_that("JSON stage 2: FCV real mixed multinomial MLE export is complete", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_fcv_mle()
      dists <- list(
        Outdoor = "binomial",
        Sex = "multinomial",
        GroupSize = "poisson",
        Age = "gaussian"
      )
      json_fixture_assert_export(fit, dists = dists, multinomial = TRUE)
    })
  })
})

test_that("JSON stage 2: ex1.dag.data Bayesian export is complete", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes()
      dists <- list(
        b1 = "binomial",
        p1 = "poisson",
        g1 = "gaussian",
        p2 = "poisson"
      )
      parsed <- json_fixture_assert_export(fit, dists = dists)
    })
  })

  expect_equal(parsed$method, "bayes")
  expect_true("original_model" %in% names(parsed))
  expect_gt(length(json_fixture_collect_rows(parsed$variables)), 0)
  expect_gt(length(json_fixture_collect_rows(parsed$parameters)), 0)
  expect_gt(length(json_fixture_collect_rows(parsed$arcs)), 0)
})
