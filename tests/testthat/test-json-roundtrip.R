# Round-trip: export -> import must reproduce the fit; export is idempotent.

# Fields that are not stored in the network document (see property spec).
JRT_NOT_STORED <- c("call", "group.ids")

jrt_expect_same_fit <- function(original, imported, info) {
  expect_s3_class(imported, "abnFit")
  expect_s3_class(imported$abnDag, "abnDag")
  expect_setequal(setdiff(names(imported), JRT_NOT_STORED),
                  setdiff(names(original), JRT_NOT_STORED))
  for (field in setdiff(names(original), c(JRT_NOT_STORED, "abnDag"))) {
    expect_equal(imported[[field]], original[[field]],
                 tolerance = .Machine$double.eps^0.5,
                 info = paste(info, "field:", field))
  }
  expect_equal(imported$abnDag$dag, original$abnDag$dag, info = info)
  expect_equal(imported$abnDag$data.dists, original$abnDag$data.dists, info = info)
}

jrt_strip_created <- function(json) {
  doc <- jdoc_parse(json)
  doc[["metadata"]][["created"]] <- NULL
  doc
}

test_that("every fixture round-trips all stored fields", {
  for (name in jfx_names()) {
    jfx_skip_bayes(name)
    original <- jfx_fit(name)
    imported <- import_abnFit(json = jfx_json(name))
    jrt_expect_same_fit(original, imported, info = name)
  }
})

test_that("every fixture round-trips completely when the data is supplied", {
  for (name in jfx_names()) {
    jfx_skip_bayes(name)
    original <- jfx_fit(name)
    imported <- import_abnFit(json = jfx_json(name), data = jfx_spec(name)$data)
    jrt_expect_same_fit(original, imported, info = name)
    expect_equal(imported$abnDag$data.df, original$abnDag$data.df, info = name)
    expect_equal(imported$group.ids, original$group.ids, info = name)
  }
})

test_that("export -> import -> export is idempotent", {
  for (name in jfx_names()) {
    jfx_skip_bayes(name)
    first <- jfx_json(name)
    second <- export_abnFit(import_abnFit(json = first))
    expect_equal(jrt_strip_created(second), jrt_strip_created(first), info = name)
  }
})

test_that("the core alone reproduces all model values", {
  model_fields <- c("coef", "Stderror", "mu", "betas", "sigma", "sigma_alpha", "modes",
                    "mse", "mliknode", "mlik", "aicnode", "aic", "bicnode", "bic",
                    "mdlnode", "df", "sse", "marginals", "marginal.quantiles", "priors",
                    "centre", "group.var", "grouped.vars", "method")
  for (name in jfx_names()) {
    jfx_skip_bayes(name)
    original <- jfx_fit(name)
    core <- jdoc_drop_extensions(jfx_doc(name))
    imported <- import_abnFit(json = jdoc_serialize(core))
    for (field in intersect(model_fields, names(original))) {
      o <- original[[field]]
      i <- imported[[field]]
      # native names may differ without the extension; values must not
      expect_equal(unname(unlist(i)), unname(unlist(o)),
                   tolerance = .Machine$double.eps^0.5,
                   info = paste(name, "field:", field))
    }
    expect_equal(imported$abnDag$dag, original$abnDag$dag, info = name)
  }
})

test_that("metadata passes through import and re-export", {
  ref <- list(uri = "data.json", schema_version = "bn-data", sha256 = "abc123")
  json <- export_abnFit(jfx_fit("ex1_mle"), label = "L", scenario_id = "S",
                        data_reference = ref)
  again <- jdoc_parse(export_abnFit(import_abnFit(json = json)))
  expect_equal(again[["metadata"]][["label"]], "L")
  expect_equal(again[["metadata"]][["scenario_id"]], "S")
  expect_equal(again[["metadata"]][["data_reference"]], ref)
})

test_that("a foreign document round-trips its core", {
  doc <- jdoc_foreign_document()
  again <- jdoc_parse(export_abnFit(import_abnFit(json = jdoc_serialize(doc))))
  for (block in c("variables", "arcs", "groups")) {
    expect_equal(again[[block]], jdoc_parse(jdoc_serialize(doc))[[block]], info = block)
  }
  expect_setequal(jdoc_values(again[["parameters"]]), jdoc_values(doc[["parameters"]]))
})

test_that("round-trip through a file works", {
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path))
  export_abnFit(jfx_fit("g2b2c_mle"), file = path)
  jrt_expect_same_fit(jfx_fit("g2b2c_mle"), import_abnFit(file = path), info = "file")
})
