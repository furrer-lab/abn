make_valid_data_json <- function() {
  list(
    metadata = list(
      schema_version = "bn-data-v1",
      orientation = "columnar",
      column_order = c("G1", "B1", "C"),
      missing_values = list(
        missing = NULL,
        not_a_number = "NaN",
        positive_infinity = "Infinity",
        negative_infinity = "-Infinity"
      ),
      columns = list(
        list(
          name = "G1",
          role = "node",
          semantic_type = "continuous",
          model_distribution = "gaussian",
          value_encoding = "number",
          nullable = FALSE
        ),
        list(
          name = "B1",
          role = "node",
          semantic_type = "binary",
          model_distribution = "binomial",
          value_encoding = "labels",
          categories = c("0", "1"),
          ordered = FALSE,
          nullable = TRUE
        ),
        list(
          name = "C",
          role = "node",
          semantic_type = "categorical",
          model_distribution = "multinomial",
          value_encoding = "labels",
          categories = c("a", "b"),
          ordered = FALSE,
          nullable = FALSE
        )
      ),
      summary = list(
        row_count = 3,
        column_count = 3,
        missing_count = 1,
        columns = list(
          G1 = list(missing_count = 0, unique_count = 3),
          B1 = list(missing_count = 1, unique_count = 2),
          C = list(missing_count = 0, unique_count = 2)
        )
      )
    ),
    data = list(
      G1 = list(0.1, -0.4, 1.2),
      B1 = list("0", NULL, "1"),
      C = list("a", "b", "a")
    )
  )
}

make_invalid_summary_data_json <- function() {
  x <- make_valid_data_json()
  x$metadata$summary$row_count <- 99
  x$metadata$summary$column_count <- 2
  x$metadata$summary$missing_count <- 0
  x$metadata$summary$columns$B1$missing_count <- 0
  x$metadata$summary$columns$C$unique_count <- 99
  x
}

test_that("validate_data_json_summary accepts matching summary metadata", {
  result <- validate_data_json_summary(make_valid_data_json(), on_error = "return")

  expect_true(result$valid)
  expect_false(result$repaired)
  expect_length(result$errors, 0)
  expect_equal(result$summary_expected$row_count, 3)
  expect_equal(result$summary_expected$column_count, 3)
  expect_equal(result$summary_expected$missing_count, 1)
  expect_equal(result$summary_expected$columns$B1$missing_count, 1)
  expect_equal(result$summary_expected$columns$B1$unique_count, 2)
})

test_that("validate_data_json_summary errors by default on stale summary metadata", {
  expect_error(
    validate_data_json_summary(make_invalid_summary_data_json()),
    regexp = "summary"
  )
})

test_that("validate_data_json_summary can return summary diagnostics silently", {
  result <- validate_data_json_summary(
    make_invalid_summary_data_json(),
    on_error = "return"
  )

  expect_false(result$valid)
  expect_false(result$repaired)
  expect_gt(length(result$errors), 0)
  expect_equal(result$summary_found$row_count, 99)
  expect_equal(result$summary_expected$row_count, 3)
  expect_equal(result$summary_found$columns$C$unique_count, 99)
  expect_equal(result$summary_expected$columns$C$unique_count, 2)
})

test_that("validate_data_json_summary can warn and return summary diagnostics", {
  expect_warning(
    result <- validate_data_json_summary(
      make_invalid_summary_data_json(),
      on_error = "warning"
    ),
    regexp = "summary"
  )

  expect_false(result$valid)
  expect_false(result$repaired)
  expect_gt(length(result$warnings), 0)
  expect_equal(result$summary_expected$row_count, 3)
})

test_that("validate_data_json_summary can repair stale summary metadata", {
  result <- validate_data_json_summary(
    make_invalid_summary_data_json(),
    on_error = "return",
    repair = TRUE
  )

  expect_false(result$valid)
  expect_true(result$repaired)
  expect_gt(length(result$errors), 0)
  expect_equal(result$object$metadata$summary, result$summary_expected)

  repaired <- validate_data_json_summary(result$object, on_error = "return")
  expect_true(repaired$valid)
  expect_false(repaired$repaired)
})

test_that("validate_data_json_summary can add a missing summary when repairing", {
  x <- make_valid_data_json()
  x$metadata$summary <- NULL

  result <- validate_data_json_summary(x, on_error = "return", repair = TRUE)

  expect_true(result$valid)
  expect_true(result$repaired)
  expect_equal(result$object$metadata$summary$row_count, 3)
  expect_equal(result$object$metadata$summary$column_count, 3)
  expect_equal(result$object$metadata$summary$columns$B1$missing_count, 1)
})

test_that("validate_data_json_summary rejects inconsistent data column lengths", {
  x <- make_valid_data_json()
  x$data$C <- list("a", "b")

  expect_error(
    validate_data_json_summary(x, repair = TRUE),
    regexp = "length"
  )
})
