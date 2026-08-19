test_that("MLE export-import preserves all relevant native fields", {
  original <- json_fixture_fit_g2b2c_mle()
  imported <- import_abnFit(json = export_abnFit(original))

  json_fixture_assert_roundtrip(original, imported, multinomial = TRUE)
})

test_that("grouped MLE export-import preserves all relevant native fields", {
  original <- json_fixture_fit_g2pbcgrp_mle_grouped()
  imported <- import_abnFit(json = export_abnFit(original))

  json_fixture_assert_roundtrip(original, imported, grouped = TRUE,
                                multinomial = TRUE)
})

test_that("Bayesian export-import preserves all relevant native fields", {
  original <- json_fixture_fit_ex1_bayes()
  imported <- import_abnFit(json = export_abnFit(original))

  json_fixture_assert_roundtrip(original, imported, bayes = TRUE)
})

test_that("Bayesian posterior export-import preserves all relevant native fields", {
  original <- json_fixture_fit_ex1_bayes_fixed()
  imported <- import_abnFit(json = export_abnFit(original))

  json_fixture_assert_roundtrip(original, imported, bayes = TRUE, fixed = TRUE)
})

test_that("metadata is optional for generic fit import", {
  original <- json_fixture_fit_ex1_mle()
  document <- jsonlite::fromJSON(export_abnFit(original), simplifyVector = FALSE)
  document$metadata <- NULL

  imported <- import_abnFit(json = jsonlite::toJSON(
    document, auto_unbox = TRUE, null = "null", digits = NA
  ))

  expect_s3_class(imported, "abnFit")
  expect_s3_class(imported$abnDag, "abnDag")
  expect_equal(imported$method, original$method)
  expect_equal(imported$abnDag$dag, original$abnDag$dag)
  expect_equal(imported$abnDag$data.dists, original$abnDag$data.dists)
  json_fixture_assert_native_fields(original, imported)
})
