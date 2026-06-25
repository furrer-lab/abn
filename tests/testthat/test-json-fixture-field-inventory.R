json_fixture_inventory_cases <- list(
  g2b2c_mle = list(
    fit = json_fixture_fit_g2b2c_mle,
    must_preserve = c("abnDag", "coef", "Stderror", "method"),
    intentionally_dropped = c("pvalue", "aicnode", "aic", "bicnode", "bic",
                              "mdlnode", "mdl", "df", "sse", "mse", "mliknode", "mlik"),
    derived = character(),
    optional_if_present = c("multinomial.states", "pvalue", "call"),
    unsupported_yet = character()
  ),
  g2pbcgrp_grouped_mle = list(
    fit = json_fixture_fit_g2pbcgrp_mle_grouped,
    must_preserve = c("abnDag", "mu", "betas", "sigma", "sigma_alpha", "method",
                      "group.var", "group.ids", "grouped.vars"),
    intentionally_dropped = c("aicnode", "aic", "bicnode", "bic", "mdlnode", "mdl",
                              "df", "sse", "mse", "mliknode", "mlik"),
    derived = character(),
    optional_if_present = c("multinomial.states", "call"),
    unsupported_yet = character()
  ),
  ex1_bayes = list(
    fit = json_fixture_fit_ex1_bayes,
    must_preserve = c("abnDag", "coef", "method", "modes", "mliknode", "mlik",
                      "mse", "used.INLA", "error.code", "error.code.desc",
                      "hessian.accuracy"),
    intentionally_dropped = character(),
    derived = character(),
    optional_if_present = c("call"),
    unsupported_yet = character()
  ),
  ex1_bayes_fixed = list(
    fit = json_fixture_fit_ex1_bayes_fixed,
    must_preserve = c("abnDag", "coef", "method", "modes", "mliknode", "mlik",
                      "mse", "used.INLA", "error.code", "error.code.desc",
                      "hessian.accuracy", "marginals", "marginal.quantiles"),
    intentionally_dropped = character(),
    derived = character(),
    optional_if_present = c("call"),
    unsupported_yet = character()
  )
)

test_that("JSON field inventory classifies every fitted field", {
  for (case_name in names(json_fixture_inventory_cases)) {
    case <- json_fixture_inventory_cases[[case_name]]
    suppressMessages({
      suppressWarnings({
        fit <- case$fit()
      })
    })

    classified <- unique(c(case$must_preserve, case$intentionally_dropped,
                           case$derived, case$optional_if_present,
                           case$unsupported_yet))
    expect_true(all(names(fit) %in% classified))
    expect_true(all(case$must_preserve %in% names(fit)))
  }
})

test_that("JSON field inventory must-preserve fields survive the current contract", {
  for (case_name in names(json_fixture_inventory_cases)) {
    case <- json_fixture_inventory_cases[[case_name]]
    suppressMessages({
      suppressWarnings({
        fit <- case$fit()
        imported <- import_abnFit(json = export_abnFit(fit))
      })
    })

    for (field in case$must_preserve) {
      if (field %in% c("abnDag", "coef", "Stderror", "method", "mu", "betas",
                      "sigma", "sigma_alpha", "multinomial.states", "group.var")) {
        expect_true(field %in% names(imported), info = paste(case_name, field))
      } else if (field %in% c("modes", "mliknode", "mlik", "mse", "used.INLA",
                             "error.code", "error.code.desc", "hessian.accuracy",
                             "marginals", "marginal.quantiles")) {
        json_name <- switch(field,
                            "used.INLA" = "used_INLA",
                            "error.code" = "error_code",
                            "error.code.desc" = "error_code_desc",
                            "hessian.accuracy" = "hessian_accuracy",
                            "marginal.quantiles" = "marginal_quantiles",
                            field)
        expect_true("original_model" %in% names(imported), info = case_name)
        expect_true(json_name %in% names(imported$original_model),
                    info = paste(case_name, field))
      }
    }
  }
})
