test_that("JSON stage 4: mixed MLE canonical export is idempotent", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_g2b2c_mle()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  expect_equal(json_fixture_canonical_export(imported),
               json_fixture_json_roundtrip_value(json_fixture_canonical_export(original)))
})

test_that("JSON stage 4: grouped MLE canonical export is idempotent", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_g2pbcgrp_mle_grouped()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  expect_equal(json_fixture_canonical_export(imported),
               json_fixture_json_roundtrip_value(json_fixture_canonical_export(original)))
})

test_that("JSON stage 4: Bayesian canonical export is idempotent", {
  suppressMessages({
    suppressWarnings({
      original <- json_fixture_fit_ex1_bayes()
      imported <- import_abnFit(json = export_abnFit(original))
    })
  })

  expect_equal(json_fixture_canonical_export(imported),
               json_fixture_json_roundtrip_value(json_fixture_canonical_export(original)))
})
