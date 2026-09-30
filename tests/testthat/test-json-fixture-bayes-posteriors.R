test_that("JSON stage 4: Bayesian fixed marginal names round-trip through JSON", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes_fixed()
      parsed <- jsonlite::fromJSON(export_abnFit(fit), simplifyVector = FALSE)
    })
  })

  expect_true("posterior" %in% names(parsed$inference))
  expect_true("marginals" %in% names(parsed$inference$posterior))
  expect_true("quantiles" %in% names(parsed$inference$posterior))
  expect_setequal(names(parsed$inference$posterior$marginals), names(fit$marginals))
  expect_setequal(names(parsed$inference$posterior$quantiles),
                  names(fit$marginal.quantiles))
})

test_that("JSON stage 4: Bayesian fixed marginal numeric contents are exported", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes_fixed()
      parsed <- jsonlite::fromJSON(export_abnFit(fit), simplifyVector = FALSE)
    })
  })

  first_node <- names(fit$marginals)[[1]]
  first_param <- names(fit$marginals[[first_node]])[[1]]
  exported <- parsed$inference$posterior$marginals[[first_node]][[first_param]]
  original <- fit$marginals[[first_node]][[first_param]]

  expect_true("values" %in% names(exported))
  expect_equal(unlist(exported$column_names), colnames(original))
  expect_equal(length(exported$values), nrow(original))
  expect_equal(unname(unlist(exported$values[[1]])), unname(as.numeric(original[1, ])),
               tolerance = .Machine$double.eps^0.5)
})

test_that("JSON stage 4: Bayesian fixed quantile numeric contents are exported", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes_fixed()
      parsed <- jsonlite::fromJSON(export_abnFit(fit), simplifyVector = FALSE)
    })
  })

  first_node <- names(fit$marginal.quantiles)[[1]]
  exported <- parsed$inference$posterior$quantiles[[first_node]]
  original <- json_fixture_json_roundtrip_value(
    export_json_safe(fit$marginal.quantiles[[first_node]])
  )

  expect_equal(names(exported), names(original))
  expect_equal(exported, original, tolerance = .Machine$double.eps^0.5)
})
