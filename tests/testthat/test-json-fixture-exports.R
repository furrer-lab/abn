test_that("generic export covers mixed MLE networks", {
  fit <- json_fixture_fit_g2b2c_mle()
  dists <- list(G1 = "gaussian", B1 = "binomial", B2 = "binomial",
                C = "multinomial", G2 = "gaussian")

  json_fixture_assert_export(fit, dists = dists, multinomial = TRUE)
})

test_that("generic export covers grouped MLE networks", {
  fit <- json_fixture_fit_g2pbcgrp_mle_grouped()
  dists <- list(G1 = "gaussian", P = "poisson", B = "binomial",
                C = "multinomial", G2 = "gaussian")

  json_fixture_assert_export(fit, dists = dists, grouped = TRUE,
                             multinomial = TRUE)
})

test_that("generic export covers standard MLE networks", {
  fit <- json_fixture_fit_ex1_mle()
  dists <- list(b1 = "binomial", p1 = "poisson", g1 = "gaussian",
                b2 = "binomial", p2 = "poisson", b3 = "binomial",
                g2 = "gaussian")

  json_fixture_assert_export(fit, dists = dists)
})

test_that("generic export covers grouped binary MLE networks", {
  fit <- json_fixture_fit_ex3_mle_grouped()
  json_fixture_assert_export(
    fit,
    dists = list(b1 = "binomial", b2 = "binomial"),
    grouped = TRUE
  )
})

test_that("generic export covers grouped real MLE networks", {
  fit <- json_fixture_fit_adg_mle_grouped()
  dists <- list(AR = "binomial", eggs = "binomial", wormCount = "poisson",
                age = "gaussian", adg = "gaussian")

  json_fixture_assert_export(fit, dists = dists, grouped = TRUE)
})

test_that("generic export covers multinomial MLE networks", {
  fit <- json_fixture_fit_fcv_mle()
  dists <- list(Outdoor = "binomial", Sex = "multinomial",
                GroupSize = "poisson", Age = "gaussian")

  json_fixture_assert_export(fit, dists = dists, multinomial = TRUE)
})

test_that("Bayesian export uses generic inference and ABN extensions", {
  fit <- json_fixture_fit_ex1_bayes()
  dists <- list(b1 = "binomial", p1 = "poisson", g1 = "gaussian",
                p2 = "poisson")
  parsed <- json_fixture_assert_export(fit, dists = dists)

  expect_equal(parsed$inference$type, "bayesian")
  expect_true("abn" %in% names(parsed$metadata$extensions))
  expect_false("original_model" %in% names(parsed))
})

test_that("Bayesian posterior data uses generic inference fields", {
  fit <- json_fixture_fit_ex1_bayes_fixed()
  parsed <- json_fixture_parse_export(fit)

  expect_equal(parsed$inference$type, "bayesian")
  expect_true("posterior" %in% names(parsed$inference))
  expect_true("marginals" %in% names(parsed$inference$posterior))
  expect_true("quantiles" %in% names(parsed$inference$posterior))
  expect_false("marginals" %in% names(parsed))
  expect_false("marginal_quantiles" %in% names(parsed))
})

test_that("generic export is deterministic", {
  fit <- json_fixture_fit_g2b2c_mle()

  expect_identical(export_abnFit(fit), export_abnFit(fit))
})
