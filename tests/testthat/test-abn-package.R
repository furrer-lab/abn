# Package exports are declared explicitly (no exportPattern).

test_that("NAMESPACE declares exports explicitly", {
  ns_file <- system.file("NAMESPACE", package = "abn")
  skip_if(!nzchar(ns_file), "installed NAMESPACE not found")
  expect_false(any(grepl("^exportPattern", readLines(ns_file))))
})

test_that("the JSON API is exported", {
  exports <- getNamespaceExports("abn")
  expect_true(all(c("export_abnFit", "import_abnFit", "export_abnData",
                    "import_abnData") %in% exports))
})

test_that("JSON internals are not exported", {
  exports <- getNamespaceExports("abn")
  internals <- c("abn_json_property_spec", "export_json_safe", "import_json_safe",
                 "reconstruct_abnfit_mle", "reconstruct_abnfit_bayes",
                 "validate_json_structure", "normalize_abn_network_document",
                 "export_to_json", "compute_data_json_summary")
  expect_equal(intersect(internals, exports), character(0))
})
