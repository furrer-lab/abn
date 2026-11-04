structure_data_document <- function() {
  jdoc_parse(export_abnData(data.frame(C = c("a", "b"), X = c(1.5, 2.5)),
                            list(C = "multinomial", X = "gaussian")))
}

structure_bad_networks <- function() {
  base <- jdoc_foreign_document()
  collection <- base
  names(collection$variables) <- c("first", "second")
  missing_value <- base
  missing_value$parameters[[1]]$value <- NULL
  scalar_array <- base
  scalar_array$parameters[[1]]$value <- list(1.2)
  bad_boolean <- base
  bad_boolean$variables[[2]]$states[[2]]$baseline <- "false"
  string_transform <- base
  string_transform$variables[[1]]$transform <- list(center = "0", scale = 1)
  string_id <- base
  string_id$variables[[1]]$`_id` <- "1"
  list(collection = collection, missing_value = missing_value,
       scalar_array = scalar_array, bad_boolean = bad_boolean,
       string_transform = string_transform, string_id = string_id)
}

structure_bad_data <- function() {
  base <- structure_data_document()
  scalar_order <- base
  scalar_order$metadata$column_order <- "C"
  scalar_categories <- base
  scalar_categories$metadata$columns[[1]]$categories <- "a"
  numeric_categories <- base
  numeric_categories$metadata$columns[[1]]$categories <- list(1, 2)
  missing_nullable <- base
  missing_nullable$metadata$columns[[1]]$nullable <- NULL
  scalar_column <- base
  scalar_column$data$X <- 1.5
  nested_scalar <- base
  nested_scalar$metadata$columns[[2]]$nullable <- TRUE
  nested_scalar$data$X[[1]] <- list(nested = 1)
  list(scalar_order = scalar_order, scalar_categories = scalar_categories,
       numeric_categories = numeric_categories, missing_nullable = missing_nullable,
       scalar_column = scalar_column, nested_scalar = nested_scalar)
}

test_that("network import rejects malformed shapes before reconstructing values", {
  for (doc in structure_bad_networks()) {
    expect_error(import_abnFit(json = jdoc_serialize(doc)))
  }
})

test_that("the network schema also rejects those malformed structural documents", {
  for (doc in structure_bad_networks()) jdoc_expect_invalid(jdoc_serialize(doc))
})

test_that("data import rejects malformed shapes before coercing columns", {
  for (doc in structure_bad_data()) {
    expect_error(import_abnData(json = jdoc_serialize(doc)))
  }
})

test_that("the data schema also rejects those malformed structural documents", {
  for (doc in structure_bad_data()) {
    jdoc_expect_invalid(jdoc_serialize(doc), "bn-data.schema.json")
  }
})

test_that("explicit null parameter estimates differ from absent required keys", {
  doc <- jdoc_foreign_document()
  doc$parameters[[1]]["value"] <- list(NULL)
  fit <- import_abnFit(json = jdoc_serialize(doc))
  expect_true(is.na(fit$coef$height[1, 1]))
})

test_that("recognized extension flags must be actual booleans", {
  doc <- jdoc_foreign_document()
  doc$metadata$extensions <- list(abn = list(nodes = list(
    list(variable = 1, used_inla = "true"))))
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "used_inla.*boolean")
})

test_that("unknown extension namespaces remain ignorable", {
  doc <- jdoc_foreign_document()
  doc$metadata$extensions <- list(other_tool = list(anything = list(1, "two")))
  expect_s3_class(import_abnFit(json = jdoc_serialize(doc)), "abnFit")
})

test_that("non-finite in-memory parameter values cannot be exported as JSON numbers", {
  for (value in c(Inf, -Inf, NaN)) {
    doc <- jdoc_foreign_document()
    doc$parameters[[1]]$value <- value
    expect_error(abn:::abn_json_check_document(doc), "scalar number")
  }
})
