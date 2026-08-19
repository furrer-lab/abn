library(testthat)
library(abn)

test_that("import_abnFit accepts generic JSON files", {
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  writeLines(json_import_string(json_import_generic_document()), file)

  imported <- import_abnFit(file = file)

  expect_s3_class(imported, "abnFit")
  expect_equal(imported$method, "mle")
  expect_equal(imported$abnDag$dag["height", "status"], 1)
})

test_that("file and string imports produce equivalent native fits", {
  json <- json_import_string(json_import_generic_document())
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  writeLines(json, file)

  from_file <- import_abnFit(file = file)
  from_string <- import_abnFit(json = json)

  expect_true(abnfit_objects_equal(from_file, from_string))
})

test_that("imported generic models are valid abnFit objects", {
  imported <- import_abnFit(json = json_import_string(json_import_generic_document()))

  expect_s3_class(imported, "abnFit")
  expect_s3_class(imported$abnDag, "abnDag")
  expect_true(is.matrix(imported$abnDag$dag))
  expect_true(imported$method %in% c("mle", "bayes"))
})

test_that("generic JSON round-trip preserves genuine values", {
  first <- json_import_string(json_import_generic_document())
  imported <- import_abnFit(json = first)
  second <- export_abnFit(imported)

  json_fixture_assert_generic_roundtrip(first, second)
})

test_that("import validates required generic structure", {
  expect_error(
    import_abnFit(json = '{"invalid": "structure"}'),
    "Invalid JSON structure"
  )
  expect_error(
    import_abnFit(json = '{"metadata": {}, "parameters": []}'),
    "Invalid JSON structure"
  )
})

test_that("import reports missing files", {
  expect_error(
    import_abnFit(file = "nonexistent_file.json"),
    "File 'nonexistent_file.json' does not exist"
  )
})
