library(testthat)
library(abn)

get_test_data <- function(filename) {
  path <- test_path("testdata", filename)
  if (!file.exists(path)) {
    alt_path <- system.file("testdata", filename, package = "abn")
    if (alt_path != "" && file.exists(alt_path)) return(alt_path)
    return("")
  }
  return(path)
}

test_that("import_abnFit correctly reads test JSON file (embedded format)", {
  test_file <- get_test_data("test_model.json")
  if (test_file == "") skip("Required test_model.json not found")

  imported_model <- import_abnFit(file = test_file)

  expect_s3_class(imported_model, "abnFit")
  expect_equal(imported_model$method, "mle")
  expect_equal(ncol(imported_model$abnDag$data.df), 5)
})

test_that("import_abnFit correctly reads test JSON file (normalized format)", {
  test_file <- get_test_data("test_model_normalized.json")
  if (test_file == "") skip("Required test_model_normalized.json not found")

  imported_model <- import_abnFit(file = test_file)

  expect_s3_class(imported_model, "abnFit")
  expect_equal(imported_model$method, "mle")
  expect_equal(ncol(imported_model$abnDag$data.df), 5)

  expect_true(all(c("b1", "p1", "g1", "c1", "b2") %in% colnames(imported_model$abnDag$dag)))
  expect_true(all(c("b1", "p1", "g1", "c1", "b2") %in% names(imported_model$coef)))
})

test_that("imported model is valid abnFit object", {
  test_file <- get_test_data("test_model.json")
  if (test_file == "") skip("Test data file not found")

  imported_model <- import_abnFit(file = test_file)

  expect_true("abnDag" %in% names(imported_model))
  expect_true("method" %in% names(imported_model))
  expect_true(imported_model$method %in% c("mle", "bayes"))

  expect_true(is.matrix(imported_model$abnDag$dag))
  expect_equal(nrow(imported_model$abnDag$dag), 5)
})

test_that("round-trip export-import-export produces an equivalent abnFit", {
  test_file <- get_test_data("test_model.json")
  if (test_file == "") skip("Test data file not found")

  first_model <- import_abnFit(file = test_file)
  exported_json <- export_abnFit(first_model)

  expect_type(exported_json, "character")
  expect_true(jsonlite::validate(exported_json))

  second_model <- import_abnFit(json = exported_json)

  expect_equal(first_model$abnDag$dag, second_model$abnDag$dag)
  expect_equal(first_model$abnDag$data.dists, second_model$abnDag$data.dists)
  expect_equal(names(first_model$coef), names(second_model$coef))
  for (nm in names(first_model$coef)) {
    expect_equal(colnames(first_model$coef[[nm]]), colnames(second_model$coef[[nm]]),
                 info = paste("coef column names differ for node:", nm))
    expect_equal(colnames(first_model$Stderror[[nm]]), colnames(second_model$Stderror[[nm]]),
                 info = paste("Stderror column names differ for node:", nm))
  }
  expect_equal(first_model$scenario_id, second_model$scenario_id)
  expect_equal(first_model$label, second_model$label)
})

test_that("import_abnFit works with JSON string directly", {
  test_file <- get_test_data("test_model.json")
  if (test_file == "") skip("Test data file not found")
  json_string <- paste(readLines(test_file, warn = FALSE), collapse = "\n")

  imported_model <- import_abnFit(json = json_string)

  expect_s3_class(imported_model, "abnFit")
  expect_equal(imported_model$method, "mle")
})

test_that("import_abnFit validates JSON structure", {
  invalid_json <- '{"invalid": "structure"}'
  expect_error(import_abnFit(json = invalid_json), "Invalid JSON structure")
})

test_that("import_abnFit handles missing required fields", {
  incomplete_json <- '{"variables": [], "parameters": []}'
  expect_error(import_abnFit(json = incomplete_json), "Invalid JSON structure")
})

test_that("import_abnFit handles file not found error", {
  expect_error(
    import_abnFit(file = "nonexistent_file.json"),
    "File 'nonexistent_file.json' does not exist"
  )
})

test_that("embedded and normalized JSON produce equivalent abnFit objects", {
  embedded_file <- get_test_data("test_model.json")
  normalized_file <- get_test_data("test_model_normalized.json")
  if (embedded_file == "") skip("test_model.json not found")
  if (normalized_file == "") skip("test_model_normalized.json not found")

  embedded_model <- import_abnFit(file = embedded_file)
  normalized_model <- import_abnFit(file = normalized_file)

  expect_equal(colnames(embedded_model$abnDag$dag), colnames(normalized_model$abnDag$dag))
  expect_equal(rownames(embedded_model$abnDag$dag), rownames(normalized_model$abnDag$dag))
  expect_equal(embedded_model$abnDag$dag, normalized_model$abnDag$dag)
  expect_equal(embedded_model$abnDag$data.dists, normalized_model$abnDag$data.dists)
  expect_equal(names(embedded_model$coef), names(normalized_model$coef))
  expect_equal(names(embedded_model$Stderror), names(normalized_model$Stderror))

  for (nm in names(embedded_model$coef)) {
    expect_equal(colnames(embedded_model$coef[[nm]]), colnames(normalized_model$coef[[nm]]),
                 info = paste("coef column names differ for node:", nm))
    expect_equal(colnames(embedded_model$Stderror[[nm]]), colnames(normalized_model$Stderror[[nm]]),
                 info = paste("Stderror column names differ for node:", nm))
  }
})

test_that("import_abnFit preserves scenario_id and label", {
  test_file <- get_test_data("test_model.json")
  if (test_file == "") skip("test_model.json not found")

  imported_model <- import_abnFit(file = test_file)

  expect_equal(imported_model$scenario_id, "test_scenario_1")
  expect_equal(imported_model$label, "Test ABN Model")
})

test_that("normalized format round-trip preserves multinomial states", {
  test_file <- get_test_data("test_model_normalized.json")
  if (test_file == "") skip("test_model_normalized.json not found")

  first_model <- import_abnFit(file = test_file)
  exported_json <- export_abnFit(first_model)

  second_model <- import_abnFit(json = exported_json)
  parsed <- jsonlite::fromJSON(exported_json)

  c1_states <- parsed$variables$states[[which(parsed$variables$variable_id == "c1")]]
  expect_equal(nrow(c1_states), 3)
  expect_equal(c1_states$value_name, c("Low", "Medium", "High"))
})
