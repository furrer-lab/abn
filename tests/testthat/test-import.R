context("Import/Export Round-Trip Functionality")

test_that("import_abnFit correctly reads test JSON file (embedded format)", {
  # Arrange
  test_file <- system.file("testdata", "test_model.json", package = "abn")
  
  # Act
  imported_model <- import_abnFit(file = test_file)
  
  # Assert
  expect_true(inherits(imported_model, "abnFit"))
  expect_equal(imported_model$method, "mle")
  expect_equal(length(imported_model$abnDag$data.df), 4)  # 4 variables
})

test_that("import_abnFit correctly reads test JSON file (normalized format)", {
  # Arrange
  test_file <- system.file("testdata", "test_model_normalized.json", package = "abn")
  
  # Act
  imported_model <- import_abnFit(file = test_file)
  
  # Assert
  expect_true(inherits(imported_model, "abnFit"))
  expect_equal(imported_model$method, "mle")
  expect_equal(length(imported_model$abnDag$data.df), 4)  # 4 variables
  # Check that link functions were resolved correctly
  # Note: Since we're creating a simplified abnFit object, we can't easily check
  # the internal link function representation, but we can verify the model was imported
})

test_that("imported model is valid abnFit object", {
  # Arrange
  test_file <- system.file("testdata", "test_model.json", package = "abn")
  
  # Act
  imported_model <- import_abnFit(file = test_file)
  
  # Assert
  expect_true(exists("abnDag", where = imported_model))
  expect_true(exists("method", where = imported_model))
  expect_true(imported_model$method %in% c("mle", "bayes"))
  
  // Check DAG structure
  expect_true(is.matrix(imported_model$abnDag$dag))
  expect_equal(nrow(imported_model$abnDag$dag), 4)  // 4x4 matrix for 4 variables
})

test_that("round-trip export-import-export produces equivalent JSON", {
  # Arrange
  test_file <- system.file("testdata", "test_model.json", package = "abn")
  original_json <- paste(readLines(test_file), collapse = "\n")
  
  // Act
  imported_model <- import_abnFit(file = test_file)
  exported_json <- export_abnFit(imported_model)
  
  // Assert
  expect_type(exported_json, "character")
  expect_true(jsonlite::validate(exported_json))
  
  // Parse both JSONs for comparison (ignoring formatting differences)
  original_parsed <- jsonlite::fromJSON(original_json)
  exported_parsed <- jsonlite::fromJSON(exported_json)
  
  // Compare key components
  expect_equal(original_parsed$variables, exported_parsed$variables)
  expect_equal(original_parsed$parameters, exported_parsed$parameters)
  expect_equal(original_parsed$arcs, exported_parsed$arcs)
  
  // Handle optional fields
  expect_equal(original_parsed$scenario_id, exported_parsed$scenario_id)
  expect_equal(original_parsed$label, exported_parsed$label)
})

test_that("round-trip export-import-export works with normalized format", {
  // Arrange
  test_file <- system.file("testdata", "test_model_normalized.json", package = "abn")
  original_json <- paste(readLines(test_file), collapse = "\n")
  
  // Act
  imported_model <- import_abnFit(file = test_file)
  exported_json <- export_abnFit(imported_model)
  
  // Assert
  expect_type(exported_json, "character")
  expect_true(jsonlite::validate(exported_json))
  
  // Parse both JSONs for comparison (ignoring formatting differences)
  original_parsed <- jsonlite::fromJSON(original_json)
  exported_parsed <- jsonlite::fromJSON(exported_json)
  
  // Compare key components - note that the exported JSON will be in embedded format
  expect_equal(original_parsed$variables, exported_parsed$variables)
  expect_equal(original_parsed$parameters, exported_parsed$parameters)
  expect_equal(original_parsed$arcs, exported_parsed$arcs)
  
  // Handle optional fields
  expect_equal(original_parsed$scenario_id, exported_parsed$scenario_id)
  expect_equal(original_parsed$label, exported_parsed$label)
})

test_that("import_abnFit works with JSON string directly", {
  # Arrange
  test_file <- system.file("testdata", "test_model.json", package = "abn")
  json_string <- paste(readLines(test_file), collapse = "\n")
  
  // Act
  imported_model <- import_abnFit(json = json_string)
  
  // Assert
  expect_true(inherits(imported_model, "abnFit"))
  expect_equal(imported_model$method, "mle")
})

test_that("import_abnFit validates JSON structure", {
  // Arrange
  invalid_json <- '{"invalid": "structure"}'
  
  // Act & Assert
  expect_error(
    import_abnFit(json = invalid_json),
    "Invalid JSON structure"
  )
})

test_that("import_abnFit handles missing required fields", {
  // Arrange
  incomplete_json <- '{
    "variables": [],
    "parameters": []
    // missing "arcs"
  }'
  
  // Act & Assert
  expect_error(
    import_abnFit(json = incomplete_json),
    "Invalid JSON structure"
  )
})

test_that("import_abnFit handles file not found error", {
  // Act & Assert
  expect_error(
    import_abnFit(file = "nonexistent_file.json"),
    "File 'nonexistent_file.json' does not exist"
  )
})