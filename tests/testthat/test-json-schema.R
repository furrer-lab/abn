# JSON Schema files and validation of exported documents.

test_that("schema files are installed and are valid JSON", {
  for (name in c("bayesian-network.schema.json", "bn-data.schema.json")) {
    path <- jdoc_schema_file(name)
    expect_true(nzchar(path), info = name)
    schema <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    expect_equal(schema$`$schema`, "http://json-schema.org/draft-07/schema#")
  }
})

test_that("the network schema pins the schema version", {
  schema <- jsonlite::fromJSON(jdoc_schema_file(), simplifyVector = FALSE)
  expect_setequal(unlist(schema$required), JDOC_BLOCKS)
  expect_equal(schema$definitions$metadata$properties$schema_version$const,
               JDOC_SCHEMA_VERSION)
})

test_that("exports of all fixtures validate against the schema", {
  for (name in names(jfx_registry())) {
    jfx_skip_bayes(name)
    jdoc_expect_valid(jfx_json(name))
  }
})

test_that("the foreign reference document validates against the schema", {
  jdoc_expect_valid(jdoc_serialize(jdoc_foreign_document()))
})

test_that("the schema rejects structurally invalid documents", {
  base <- jdoc_foreign_document()

  wrong_version <- base
  wrong_version[["metadata"]][["schema_version"]] <- "not-a-known-format"
  jdoc_expect_invalid(jdoc_serialize(wrong_version))

  missing_block <- base
  missing_block[["groups"]] <- NULL
  jdoc_expect_invalid(jdoc_serialize(missing_block))

  bad_kind <- base
  bad_kind[["parameters"]][[1]][["kind"]] <- "slope"
  jdoc_expect_invalid(jdoc_serialize(bad_kind))

  bad_type <- base
  bad_type[["variables"]][[1]][["type"]] <- "numeric"
  jdoc_expect_invalid(jdoc_serialize(bad_type))

  bad_scale <- base
  bad_scale[["parameters"]][[2]][["scale"]] <- "sd"
  jdoc_expect_invalid(jdoc_serialize(bad_scale))

  missing_states <- base
  missing_states[["variables"]][[2]][["states"]] <- NULL
  jdoc_expect_invalid(jdoc_serialize(missing_states))
})

test_that("exported data documents validate against the data schema", {
  spec <- jfx_spec("ex1_mle")
  json <- export_abnData(spec$data, spec$dists)
  jdoc_expect_valid(json, "bn-data.schema.json")
})

test_that("a parameter value may be null (not estimable)", {
  doc <- jdoc_foreign_document()
  doc[["parameters"]][[4]]["value"] <- list(NULL)
  jdoc_expect_valid(jdoc_serialize(doc))
  data <- data.frame(height = c(1, 2), status = factor(c("no", "yes")))
  # only the intercept contributes
  expect_equal(jdoc_linear_predictor(doc, data, "status"), c(-0.42, -0.42))
})
