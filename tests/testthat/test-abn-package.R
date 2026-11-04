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

test_that("JSON privacy preserves the exact pre-change non-JSON export set", {
  # Captured with load_all(export_all = FALSE) before replacing exportPattern.
  # This pins all legacy names, including helpers, S3 methods and C bindings.
  legacy <- readLines(test_path("testdata", "non-json-exports.txt"))
  json_api <- c("export_abnFit", "import_abnFit", "export_abnData", "import_abnData",
                "validate_data_json_summary", "repair_data_json_summary")
  expect_setequal(getNamespaceExports("abn"), c(legacy, json_api))
})

test_that("all current JSON helper families stay private", {
  exports <- getNamespaceExports("abn")
  json_helpers <- c("compare_data_json_summary", "compute_data_json_summary",
                    "parse_data_json_input", "validate_abn_data_inputs",
                    "validate_abn_data_json_object", "validate_data_json_document_semantics",
                    "validate_data_json_shape", "validate_data_json_type",
                    "validate_data_json_values")
  expect_false(any(grepl("^(abn_json_|data_json_)", exports)))
  expect_length(intersect(json_helpers, exports), 0L)
})
