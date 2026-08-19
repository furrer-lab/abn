json_fixture_inventory_cases <- list(
  g2b2c_mle = json_fixture_fit_g2b2c_mle,
  g2pbcgrp_grouped_mle = json_fixture_fit_g2pbcgrp_mle_grouped,
  ex1_bayes = json_fixture_fit_ex1_bayes,
  ex1_bayes_fixed = json_fixture_fit_ex1_bayes_fixed
)

test_that("every fitted field has an explicit round-trip contract", {
  for (case_name in names(json_fixture_inventory_cases)) {
    fit <- json_fixture_inventory_cases[[case_name]]()
    relevant <- setdiff(names(fit), "call")

    expect_gt(length(relevant), 0, info = case_name)
    expect_true("abnDag" %in% relevant, info = case_name)
    expect_true("method" %in% relevant, info = case_name)
  }
})

test_that("every relevant fitted field survives export-import", {
  for (case_name in names(json_fixture_inventory_cases)) {
    fit <- json_fixture_inventory_cases[[case_name]]()
    imported <- import_abnFit(json = export_abnFit(fit))

    json_fixture_assert_native_fields(fit, imported)
  }
})
