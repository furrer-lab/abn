test_that("export_abn_json works for MLE method without group.var", {
  # ARRANGE
  abn_fit <- create_test_abnfit_mle()

  # ACT
  result <- export_abnFit(abn_fit)

  # ASSERT
  expect_type(result, "character")
  expect_true(jsonlite::validate(result))
})

test_that("export_abn_json writes to file when 'file' argument is provided", {
  # ARRANGE
  abn_fit <- create_test_abnfit_mle()
  temp_file <- tempfile(fileext = ".json")
  on.exit(unlink(temp_file), add = TRUE)

  # ACT
  result <- export_abnFit(abn_fit, file = temp_file)

  # ASSERT
  expect_null(result) # Function should return NULL invisibly
  expect_true(file.exists(temp_file))
  expect_true(file.info(temp_file)$size > 0)
})

test_that("exported JSON has correct top-level structure", {
  # ARRANGE
  abn_fit <- create_test_abnfit_mle()

  # ACT
  json_str <- export_abnFit(abn_fit)
  parsed <- jsonlite::fromJSON(json_str)

  # ASSERT
  expect_true("graph" %in% names(parsed))
  expect_true("nodes" %in% names(parsed))
  expect_true("arcs" %in% names(parsed))
})

test_that("graph metadata is correctly exported for MLE", {
  # ARRANGE
  abn_fit <- create_test_abnfit_mle()

  # ACT
  parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

  # ASSERT
  expect_equal(parsed$graph$method, "mle")
  expect_equal(parsed$graph$n_nodes, 4)
  expect_equal(parsed$graph$n_observations, 100)
  expect_null(parsed$graph$group_var)
})
