test_that("import-export preserves generic MLE values", {
  original <- json_fixture_fit_g2b2c_mle()
  first <- export_abnFit(original)
  imported <- import_abnFit(json = first)
  second <- export_abnFit(imported)

  json_fixture_assert_generic_roundtrip(first, second)
})

test_that("import-export preserves generic grouped MLE values", {
  original <- json_fixture_fit_g2pbcgrp_mle_grouped()
  first <- export_abnFit(original)
  imported <- import_abnFit(json = first)
  second <- export_abnFit(imported)

  json_fixture_assert_generic_roundtrip(first, second)
})

test_that("import-export preserves generic Bayesian values", {
  original <- json_fixture_fit_ex1_bayes_fixed()
  first <- export_abnFit(original)
  imported <- import_abnFit(json = first)
  second <- export_abnFit(imported)

  json_fixture_assert_generic_roundtrip(first, second)
})

test_that("metadata may differ without changing genuine JSON values", {
  original <- json_fixture_fit_ex1_mle()
  first <- jsonlite::fromJSON(export_abnFit(original), simplifyVector = FALSE)
  imported <- import_abnFit(json = jsonlite::toJSON(
    first, auto_unbox = TRUE, null = "null", digits = NA
  ))
  second <- export_abnFit(imported)

  expect_false(identical(first$metadata, jsonlite::fromJSON(
    second, simplifyVector = FALSE
  )$metadata))
  json_fixture_assert_generic_roundtrip(first, second)
})
