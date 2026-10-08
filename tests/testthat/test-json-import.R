# Import: documents from other issuers, validation and attaching data.

test_that("import accepts a JSON string or a file", {
  json <- jfx_json("ex1_mle")
  from_string <- import_abnFit(json = json)
  expect_s3_class(from_string, "abnFit")

  path <- tempfile(fileext = ".json")
  on.exit(unlink(path))
  writeLines(json, path)
  from_file <- import_abnFit(file = path)
  expect_equal(from_file, from_string)
})

test_that("import requires exactly one source", {
  expect_error(import_abnFit(), "file|json")
  expect_error(import_abnFit(file = "x.json", json = "{}"), "file|json")
  expect_error(import_abnFit(file = tempfile(fileext = ".json")), "exist")
})

test_that("removed pre-release arguments are rejected", {
  expect_error(import_abnFit(json = jfx_json("ex1_mle"), conflict = "prefer_core"))
})

test_that("a foreign MLE document is imported from the core only", {
  fit <- import_abnFit(json = jdoc_serialize(jdoc_foreign_document()))
  expect_s3_class(fit, "abnFit")
  expect_s3_class(fit$abnDag, "abnDag")
  expect_equal(fit$method, "mle")
  expect_equal(fit$abnDag$data.dists, list(height = "gaussian", status = "binomial"))
  expect_equal(unname(fit$abnDag$dag["status", "height"]), 1)
  expect_equal(unname(fit$abnDag$dag["height", "status"]), 0)
  expect_equal(unname(fit$coef$height[1, ]), 1.2)
  expect_equal(unname(fit$coef$status[1, ]), c(-0.42, 0.87))
  expect_equal(unname(fit$Stderror$status[1, ]), c(0.2, 0.3))
  expect_equal(unname(fit$mse[["height"]]), 0.25)
  expect_equal(levels(fit$abnDag$data.df$status), c("no", "yes"))
  expect_equal(nrow(fit$abnDag$data.df), 0L)
})

test_that("a foreign document with non-first baseline is imported with that baseline", {
  doc <- jdoc_foreign_document()
  doc[["variables"]][[2]][["states"]] <- list(list(`_id` = 1, label = "no", baseline = FALSE),
                                    list(`_id` = 2, label = "yes", baseline = TRUE))
  fit <- import_abnFit(json = jdoc_serialize(doc))
  # abn models the second level, so the baseline must become the first level
  expect_equal(levels(fit$abnDag$data.df$status), c("yes", "no"))
})

test_that("unknown extensions are ignored", {
  doc <- jdoc_foreign_document()
  doc[["metadata"]][["extensions"]] <- list(othertool = list(anything = 1))
  expect_silent(fit <- import_abnFit(json = jdoc_serialize(doc)))
  expect_s3_class(fit, "abnFit")
})

test_that("invalid documents are rejected with informative errors", {
  base <- jdoc_foreign_document()
  check <- function(doc, pattern) {
    expect_error(import_abnFit(json = jdoc_serialize(doc)), pattern)
  }

  wrong_version <- base
  wrong_version[["metadata"]][["schema_version"]] <- "not-a-known-format"
  check(wrong_version, "schema_version")

  dup <- base
  dup[["variables"]][[2]]$`_id` <- 1
  check(dup, "_id")

  dangling_arc <- base
  dangling_arc[["arcs"]][[1]][["source"]] <- 99
  check(dangling_arc, "99")

  dangling_param <- base
  dangling_param[["parameters"]][[4]][["parent"]] <- 99
  check(dangling_param, "99")

  not_an_arc <- base
  not_an_arc[["arcs"]] <- list()
  check(not_an_arc, "arc")

  cyclic <- base
  cyclic[["arcs"]] <- list(list(source = 1, target = 2), list(source = 2, target = 1))
  check(cyclic, "cycl")

  no_baseline <- base
  no_baseline[["variables"]][[2]][["states"]][[1]][["baseline"]] <- FALSE
  check(no_baseline, "baseline")

  bad_state <- base
  bad_state[["parameters"]][[4]][["parent_state"]] <- 7
  check(bad_state, "state")
})

test_that("import attaches a data frame and restores the transformed data", {
  spec <- jfx_spec("ex1_mle_centred")
  fit <- jfx_fit("ex1_mle_centred")
  imported <- import_abnFit(json = jfx_json("ex1_mle_centred"), data = spec$data)
  expect_equal(imported$abnDag$data.df, fit$abnDag$data.df)
})

test_that("import attaches a data document", {
  spec <- jfx_spec("ex1_mle")
  fit <- jfx_fit("ex1_mle")
  data_file <- tempfile(fileext = ".json")
  on.exit(unlink(data_file))
  export_abnData(spec$data, spec$dists, file = data_file)
  imported <- import_abnFit(json = jfx_json("ex1_mle"), data = data_file)
  expect_equal(imported$abnDag$data.df, fit$abnDag$data.df)
})

test_that("import restores group.ids from the data", {
  spec <- jfx_spec("adg_mle_grouped")
  fit <- jfx_fit("adg_mle_grouped")
  without <- import_abnFit(json = jfx_json("adg_mle_grouped"))
  expect_null(without$group.ids)
  with <- import_abnFit(json = jfx_json("adg_mle_grouped"), data = spec$data)
  expect_equal(with$group.ids, fit$group.ids)
})

test_that("import rejects data that does not match the variables", {
  spec <- jfx_spec("ex1_mle")
  json <- jfx_json("ex1_mle")

  missing_col <- spec$data[, setdiff(names(spec$data), "g1")]
  expect_error(import_abnFit(json = json, data = missing_col), "g1")

  bad_levels <- spec$data
  bad_levels$b1 <- factor(ifelse(bad_levels$b1 == levels(bad_levels$b1)[1], "u", "v"))
  expect_error(import_abnFit(json = json, data = bad_levels), "b1")
})

test_that("import accepts the result of import_abnData()", {
  spec <- jfx_spec("ex1_mle")
  doc <- export_abnData(spec$data, spec$dists)
  via_list <- import_abnFit(json = jfx_json("ex1_mle"),
                            data = import_abnData(json = doc))
  via_frame <- import_abnFit(json = jfx_json("ex1_mle"), data = spec$data)
  expect_equal(via_list$abnDag$data.df, via_frame$abnDag$data.df)
})

test_that("import rejects data documents with mismatching distributions", {
  spec <- jfx_spec("ex1_mle")
  obj <- jsonlite::fromJSON(export_abnData(spec$data, spec$dists,
                                           include_summary = FALSE),
                             simplifyVector = FALSE)
  obj$metadata$data_dists$g1 <- "poisson"
  data_file <- tempfile(fileext = ".json")
  on.exit(unlink(data_file))
  writeLines(jsonlite::toJSON(obj, auto_unbox = TRUE, null = "null", digits = NA),
             data_file)
  expect_error(import_abnFit(json = jfx_json("ex1_mle"), data = data_file), "g1")
})
