library(testthat)
library(abn)

test_that("Grouped MLE: grouping column name survives JSON round-trip", {
  test_file <- test_path("testdata", "abnfit_mle_groups.Rdata")
  if (!file.exists(test_file)) skip("Grouped model test data not found")

  load(test_file)
  direct_model <- abn_fit
  group_var <- direct_model$group.var

  expect_true(is.character(group_var))
  expect_length(group_var, 1)

  json_str <- export_abnFit(direct_model)
  parsed <- jsonlite::fromJSON(json_str)
  expect_equal(parsed$metadata$issuer, "abn::export_abnFit")
  expect_identical(parsed$metadata$extensions$abn$configs$group_var, group_var)

  imported_model <- import_abnFit(json = json_str)
  expect_identical(imported_model$group.var, group_var)
})
