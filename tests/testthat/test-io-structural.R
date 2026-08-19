library(testthat)
library(abn)

.parse_json <- function(json_str) {
  jsonlite::fromJSON(json_str, simplifyVector = FALSE)
}

test_that("generic parameter IDs and references are structural", {
  fit <- create_test_abnfit_mle()
  parsed <- .parse_json(export_abnFit(fit))
  variable_ids <- vapply(parsed$variables, function(x) as.character(x$`_id`), character(1))
  parameter_ids <- vapply(parsed$parameters, function(x) as.character(x$`_id`), character(1))

  expect_equal(anyDuplicated(variable_ids), 0L)
  expect_equal(anyDuplicated(parameter_ids), 0L)
  for (parameter in parsed$parameters) {
    expect_true(parameter$target %in% variable_ids)
    if (!is.null(parameter$parents)) {
      expect_true(all(unlist(parameter$parents) %in% variable_ids))
    }
  }
})

test_that("re-export preserves genuine values after import", {
  fit <- create_test_abnfit_mle()
  first <- export_abnFit(fit)
  imported <- import_abnFit(json = first)
  second <- export_abnFit(imported)

  json_fixture_assert_generic_roundtrip(first, second)
})

test_that("parameter display labels do not control generic import", {
  fit <- create_test_abnfit_mle()
  parsed <- .parse_json(export_abnFit(fit))
  base_import <- import_abnFit(json = export_abnFit(fit))

  parsed$parameters <- lapply(parsed$parameters, function(parameter) {
    parameter$label <- "irrelevant display label"
    parameter
  })
  mutated_json <- jsonlite::toJSON(parsed, auto_unbox = TRUE,
                                   pretty = FALSE, digits = NA, null = "null")
  mutated_import <- import_abnFit(json = mutated_json)

  expect_true(abnfit_objects_equal(base_import, mutated_import))
})
