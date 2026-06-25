test_that("JSON stage 4: g2b2c_data mixed MLE round-trips", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_g2b2c_mle()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  json_fixture_assert_roundtrip(original, imported, multinomial = TRUE)
})

test_that("JSON stage 4: g2pbcgrp grouped mixed MLE round-trips", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_g2pbcgrp_mle_grouped()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  json_fixture_assert_roundtrip(original, imported, grouped = TRUE,
                                multinomial = TRUE)
})

test_that("JSON stage 4: ex1.dag.data classic mixed MLE round-trips", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_ex1_mle()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  json_fixture_assert_roundtrip(original, imported)
})

test_that("JSON stage 4: ex3.dag.data grouped binary MLE round-trips", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_ex3_mle_grouped()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  json_fixture_assert_roundtrip(original, imported, grouped = TRUE)
})

test_that("JSON stage 4: adg real grouped mixed MLE round-trips", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_adg_mle_grouped()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  json_fixture_assert_roundtrip(original, imported, grouped = TRUE)
})

test_that("JSON stage 4: FCV real mixed multinomial MLE round-trips", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_fcv_mle()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  json_fixture_assert_roundtrip(original, imported, multinomial = TRUE)
})

test_that("JSON stage 4: ex1.dag.data Bayesian round-trips", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_ex1_bayes()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  json_fixture_assert_roundtrip(original, imported, bayes = TRUE)
})

test_that("JSON stage 4: Bayesian fixed marginal metadata round-trips", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_ex1_bayes_fixed()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  json_fixture_assert_roundtrip(original, imported, bayes = TRUE, fixed = TRUE)
})
