test_that("JSON stage 3: g2b2c_data mixed MLE export imports as abnFit", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_g2b2c_mle()
      json <- export_abnFit(fit)
      parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
      expect_no_error(imported <- import_abnFit(json = json))
    })
  })

  dists <- list(
    G1 = "gaussian",
    B1 = "binomial",
    B2 = "binomial",
    C = "multinomial",
    G2 = "gaussian"
  )
  json_fixture_assert_import(imported, method = "mle", dists = dists,
                             multinomial = TRUE, parsed = parsed)
})

test_that("JSON stage 3: g2pbcgrp grouped mixed MLE export imports as abnFit", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_g2pbcgrp_mle_grouped()
      json <- export_abnFit(fit)
      parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
      expect_no_error(imported <- import_abnFit(json = json))
    })
  })

  dists <- list(
    G1 = "gaussian",
    P = "poisson",
    B = "binomial",
    C = "multinomial",
    G2 = "gaussian"
  )
  json_fixture_assert_import(imported, method = "mle", dists = dists,
                             grouped = TRUE, multinomial = TRUE, parsed = parsed)
})

test_that("JSON stage 3: ex1.dag.data classic mixed MLE export imports as abnFit", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_mle()
      json <- export_abnFit(fit)
      parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
      expect_no_error(imported <- import_abnFit(json = json))
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
  json_fixture_assert_import(imported, method = "mle", dists = dists,
                             parsed = parsed)
})

test_that("JSON stage 3: ex3.dag.data grouped binary MLE export imports as abnFit", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex3_mle_grouped()
      json <- export_abnFit(fit)
      parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
      expect_no_error(imported <- import_abnFit(json = json))
    })
  })

  dists <- list(b1 = "binomial", b2 = "binomial")
  json_fixture_assert_import(imported, method = "mle", dists = dists,
                             grouped = TRUE, parsed = parsed)
})

test_that("JSON stage 3: adg real grouped mixed MLE export imports as abnFit", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_adg_mle_grouped()
      json <- export_abnFit(fit)
      parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
      expect_no_error(imported <- import_abnFit(json = json))
    })
  })

  dists <- list(
    AR = "binomial",
    eggs = "binomial",
    wormCount = "poisson",
    age = "gaussian",
    adg = "gaussian"
  )
  json_fixture_assert_import(imported, method = "mle", dists = dists,
                             grouped = TRUE, parsed = parsed)
})

test_that("JSON stage 3: FCV real mixed multinomial MLE export imports as abnFit", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_fcv_mle()
      json <- export_abnFit(fit)
      parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
      expect_no_error(imported <- import_abnFit(json = json))
    })
  })

  dists <- list(
    Outdoor = "binomial",
    Sex = "multinomial",
    GroupSize = "poisson",
    Age = "gaussian"
  )
  json_fixture_assert_import(imported, method = "mle", dists = dists,
                             multinomial = TRUE, parsed = parsed)
})

test_that("JSON stage 3: ex1.dag.data Bayesian export imports as abnFit", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes()
      json <- export_abnFit(fit)
      parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
      expect_no_error(imported <- import_abnFit(json = json))
    })
  })

  dists <- list(
    b1 = "binomial",
    p1 = "poisson",
    g1 = "gaussian",
    p2 = "poisson"
  )
  json_fixture_assert_import(imported, method = "bayes", dists = dists,
                             original_model = TRUE, parsed = parsed)
})

test_that("JSON stage 3: Bayesian fixed marginal export imports with metadata", {
  suppressMessages({
    suppressWarnings({
      fit <- json_fixture_fit_ex1_bayes_fixed()
      json <- export_abnFit(fit)
      parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
      expect_no_error(imported <- import_abnFit(json = json))
    })
  })

  expect_equal(imported$method, "bayes")
  expect_true("original_model" %in% names(imported))
  expect_true("marginals" %in% names(imported$original_model))
  expect_true("marginal_quantiles" %in% names(imported$original_model))
  expect_equal(length(imported$original_model$marginals),
               length(parsed$original_model$marginals))
  expect_equal(length(imported$original_model$marginal_quantiles),
               length(parsed$original_model$marginal_quantiles))
})

test_that("JSON stage 3: imported MLE coefficient names follow abn conventions", {
  suppressMessages({
    suppressWarnings({
      ex1_fit <- json_fixture_fit_ex1_mle()
      ex1_imported <- import_abnFit(json = export_abnFit(ex1_fit))
      fcv_fit <- json_fixture_fit_fcv_mle()
      fcv_imported <- import_abnFit(json = export_abnFit(fcv_fit))
    })
  })

  expect_true("b1" %in% colnames(ex1_imported$coef[["p2"]]))
  expect_true("p1" %in% colnames(ex1_imported$coef[["p2"]]))
  expect_false("p2|b1" %in% colnames(ex1_imported$coef[["p2"]]))
  expect_false("p2|p1" %in% colnames(ex1_imported$coef[["p2"]]))

  expect_true(all(c("Sexm", "Sexmc", "Sexf", "Sexfc") %in%
                    colnames(fcv_imported$coef[["Outdoor"]])))
  expect_false(any(grepl("^Outdoor\\|Sex", colnames(fcv_imported$coef[["Outdoor"]]))))
})

test_that("JSON stage 3: grouped multinomial fields preserve abn names", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_g2pbcgrp_mle_grouped()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  expect_setequal(names(imported$betas$G2), names(original$betas$G2))
  expect_setequal(rownames(imported$betas$C), rownames(original$betas$C))
  expect_equal(dim(imported$sigma_alpha$C), dim(original$sigma_alpha$C))
  expect_setequal(rownames(imported$sigma_alpha$C), rownames(original$sigma_alpha$C))
  expect_setequal(colnames(imported$sigma_alpha$C), colnames(original$sigma_alpha$C))
})
