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

test_that("JSON stage 2: Bayesian parameter values match fitted coefficients", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes()
      parsed <- json_fixture_parse_export(fit)
    })
  })

  variables <- json_fixture_collect_rows(parsed$variables)
  variable_names <- stats::setNames(
    vapply(variables, json_fixture_field, character(1), "attribute_name"),
    vapply(variables, json_fixture_field, character(1), "variable_id")
  )

  for (parameter in json_fixture_collect_rows(parsed$parameters)) {
    source <- parameter$source
    if (is.list(source) && length(source) == 1 && is.list(source[[1]])) {
      source <- source[[1]]
    }
    node_name <- variable_names[[json_fixture_field(source, "variable_id")]]
    coefficients <- json_fixture_collect_rows(parameter$coefficients)
    for (coefficient in coefficients) {
      condition_type <- json_fixture_field(coefficient, "condition_type")
      if (identical(condition_type, "intercept")) {
        expected_name <- paste0(node_name, "|intercept")
      } else if (identical(condition_type, "linear_term")) {
        conditions <- json_fixture_collect_rows(coefficient$conditions)
        expect_length(conditions, 1)
        parent_id <- json_fixture_field(conditions[[1]], "parent_variable_id")
        expected_name <- variable_names[[parent_id]]
      } else {
        next
      }

      expect_true(expected_name %in% colnames(fit$coef[[node_name]]),
                  info = paste("Missing Bayesian coefficient", expected_name,
                               "for node", node_name))
      expect_equal(as.numeric(coefficient$value),
                   unname(as.numeric(fit$coef[[node_name]][1, expected_name])),
                   tolerance = 1e-12)
    }
  }
})

test_that("JSON stage 2: Bayesian original_model preserves fit metadata", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes()
      parsed <- json_fixture_parse_export(fit)
    })
  })

  expect_true("original_model" %in% names(parsed))
  expect_equal(parsed$original_model$mlik, fit$mlik)
  expect_equal(unlist(parsed$original_model$mliknode), unname(fit$mliknode))
  expect_equal(unlist(parsed$original_model$used_INLA), unname(fit$used.INLA))
  expect_equal(unlist(parsed$original_model$error_code), unname(fit$error.code))
  expect_equal(unlist(parsed$original_model$error_code_desc), unname(fit$error.code.desc))
  expect_equal(unlist(parsed$original_model$hessian_accuracy),
               unname(fit$hessian.accuracy))
})

test_that("JSON stage 2: Bayesian fixed marginals and quantiles are exported", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes_fixed()
      parsed <- json_fixture_parse_export(fit)
    })
  })

  expect_true("marginals" %in% names(fit))
  expect_true("marginal.quantiles" %in% names(fit))
  expect_true("original_model" %in% names(parsed))
  expect_true("marginals" %in% names(parsed$original_model))
  expect_true("marginal_quantiles" %in% names(parsed$original_model))
  expect_equal(length(parsed$original_model$marginals), length(fit$marginals))
  expect_equal(length(parsed$original_model$marginal_quantiles),
               length(fit$marginal.quantiles))
})

test_that("JSON stage 2: Bayesian marginal export is JSON-native", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes_fixed()
      json <- export_abnFit(fit)
      parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
    })
  })

  expect_true(jsonlite::validate(json))
  expect_false(inherits(parsed$original_model$marginals, "abnFit"))
  expect_false(inherits(parsed$original_model$marginal_quantiles, "abnFit"))
  expect_type(parsed$original_model$marginals, "list")
  expect_type(parsed$original_model$marginal_quantiles, "list")
})

test_that("JSON stage 2: JSON sanitizer handles recursive abnFit objects", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes()
    })
  })

  unsafe <- list(model = fit, nested = list(model = fit))
  safe <- export_json_safe(unsafe)
  json <- jsonlite::toJSON(safe, auto_unbox = TRUE, null = "null")

  expect_true(jsonlite::validate(json))
  expect_equal(safe$model$method, "bayes")
  expect_true("node_names" %in% names(safe$model))
  expect_false(inherits(safe$model, "abnFit"))
})

test_that("JSON stage 2: JSON sanitizer strips classed atomic abnFit values", {
  x <- structure(c(a = 1, b = 2), class = "abnFit")
  safe <- export_json_safe(x)
  json <- jsonlite::toJSON(safe, auto_unbox = TRUE, null = "null")

  expect_true(jsonlite::validate(json))
  expect_equal(safe, c(a = 1, b = 2))
  expect_false(inherits(safe, "abnFit"))
})
