json_export_rows <- function(value) {
  if (is.null(value)) return(list())
  if (is.data.frame(value)) return(split(value, seq_len(nrow(value))))
  value
}

json_export_ids <- function(rows) {
  vapply(rows, function(row) as.character(row$`_id`), character(1))
}

test_that("export_abnFit emits the generic network document", {
  fit <- create_test_abnfit_mle()
  document <- jsonlite::fromJSON(
    export_abnFit(fit),
    simplifyVector = FALSE
  )

  expect_named(document, c("metadata", "variables", "arcs", "parameters",
                           "inference"), ignore.order = TRUE)
  expect_equal(document$metadata$schema_version, "bayesian-network-v1")
  expect_equal(document$metadata$issuer, "abn::export_abnFit")
  expect_true(is.list(document$metadata$configs))
  expect_true(is.list(document$metadata$extensions))
  expect_equal(document$inference$type, "maximum_likelihood")

  expect_length(json_export_rows(document$variables), length(fit$coef))
  expect_gt(length(json_export_rows(document$parameters)), 0)
  expect_true(is.list(document$inference$estimates))
  expect_true(is.list(document$inference$uncertainty))
  expect_true(is.list(document$inference$diagnostics))
})

test_that("network objects use stable generic IDs and references", {
  fit <- create_test_abnfit_mle()
  document <- jsonlite::fromJSON(
    export_abnFit(fit),
    simplifyVector = FALSE
  )
  variables <- json_export_rows(document$variables)
  variable_ids <- json_export_ids(variables)
  parameter_ids <- json_export_ids(json_export_rows(document$parameters))

  expect_equal(anyDuplicated(variable_ids), 0L)
  expect_equal(anyDuplicated(parameter_ids), 0L)
  expect_equal(variable_ids, as.character(seq_along(variable_ids)))

  for (variable in variables) {
    expect_true(all(c("_id", "name", "type") %in% names(variable)))
    optional <- setdiff(names(variable), c("_id", "name", "type"))
    expect_setequal(optional,
                    if (identical(variable$type, "categorical")) "states" else character())
  }

  for (arc in json_export_rows(document$arcs)) {
    expect_named(arc, c("source", "target"), ignore.order = TRUE)
    expect_true(as.character(arc$source) %in% variable_ids)
    expect_true(as.character(arc$target) %in% variable_ids)
    expect_false(identical(arc$source, arc$target))
  }

  for (parameter in json_export_rows(document$parameters)) {
    expect_true(parameter$`_id` %in% as.integer(parameter_ids))
    expect_true(as.character(parameter$target) %in% variable_ids)
    if (!is.null(parameter$parents)) {
      expect_true(all(unlist(parameter$parents) %in% variable_ids))
    }
    expect_true(parameter$kind %in% c("intercept", "coefficient", "variance",
                                      "random_variance", "random_covariance"))
  }
})

test_that("Bayesian ABN details are namespaced as extensions", {
  fit <- json_fixture_fit_ex1_bayes()
  document <- jsonlite::fromJSON(
    export_abnFit(fit),
    simplifyVector = FALSE
  )

  expect_equal(document$inference$type, "bayesian")
  expect_true("abn" %in% names(document$metadata$extensions))
  expect_false("method" %in% names(document))
  expect_false("original_model" %in% names(document))
  expect_false("used.INLA" %in% names(document))
  expect_true(is.list(document$inference$estimates))
  expect_true(is.list(document$inference$uncertainty))
  expect_true(is.list(document$inference$diagnostics))
})

test_that("include_network includes or omits the DAG", {
  fit <- create_test_abnfit_mle()
  included <- jsonlite::fromJSON(
    export_abnFit(fit, include_network = TRUE),
    simplifyVector = FALSE
  )
  omitted <- jsonlite::fromJSON(
    export_abnFit(fit, include_network = FALSE),
    simplifyVector = FALSE
  )

  expect_true("variables" %in% names(included))
  expect_true("arcs" %in% names(included))
  expect_false("variables" %in% names(omitted))
  expect_false("arcs" %in% names(omitted))
  expect_true("parameters" %in% names(omitted))
  expect_true("inference" %in% names(omitted))
})

test_that("export_abnFit writes the generic document to a file", {
  fit <- create_test_abnfit_mle()
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)

  result <- export_abnFit(fit, file = file)
  document <- jsonlite::fromJSON(paste(readLines(file), collapse = "\n"),
                                 simplifyVector = FALSE)

  expect_identical(unname(result), file)
  expect_equal(document$metadata$schema_version, "bayesian-network-v1")
})

test_that("raw data export identifies the generic data schema", {
  data(g2b2c_data)
  json <- export_abnData(
    g2b2c_data,
    list(G1 = "gaussian", B1 = "binomial", B2 = "binomial",
         C = "multinomial", G2 = "gaussian")
  )
  document <- jsonlite::fromJSON(json, simplifyVector = FALSE)

  expect_equal(document$metadata$schema_version, "bn-data-v1")
  expect_equal(document$metadata$issuer, "abn::export_abnData")
})
