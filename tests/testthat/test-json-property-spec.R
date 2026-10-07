# abn_json_property_spec() is the single mapping between fitAbn fields and the
# JSON document. One row per (field, fit_type):
#   field    : name of the top-level fitAbn element
#   fit_type : one of "mle", "mle_grouped", "bayes"
#   location : JSON path of the one place the field is stored, NA if not stored
#   reason   : why a field is not stored (NA if stored): "derived", "data", "excluded"

JSPEC_FIT_TYPES <- c("mle", "mle_grouped", "bayes")

jspec_fit_type <- function(name) {
  spec <- jfx_spec(name)
  if (identical(spec$method, "bayes")) "bayes" else if (spec$grouped) "mle_grouped" else "mle"
}

jspec_row <- function(field, fit_type) {
  spec <- abn_json_property_spec()
  spec[spec$field == field & spec$fit_type == fit_type, , drop = FALSE]
}

test_that("the property specification has the expected shape", {
  spec <- abn_json_property_spec()
  expect_s3_class(spec, "data.frame")
  expect_true(all(c("field", "fit_type", "location", "reason") %in% names(spec)))
  expect_true(all(spec$fit_type %in% JSPEC_FIT_TYPES))
  expect_equal(anyDuplicated(paste(spec$field, spec$fit_type)), 0L)
  expect_type(spec$location, "character")
})

test_that("every field is stored exactly once per fit type or has a reason not to be", {
  spec <- abn_json_property_spec()
  stored <- !is.na(spec$location)
  expect_true(all(is.na(spec$reason[stored])))
  expect_true(all(spec$reason[!stored] %in% c("derived", "data", "excluded")))
  for (type in JSPEC_FIT_TYPES) {
    locations <- spec$location[stored & spec$fit_type == type]
    expect_equal(anyDuplicated(locations), 0L, info = type)
  }
})

test_that("locations point into a top-level block of the document", {
  spec <- abn_json_property_spec()
  roots <- sub("[.\\[].*$", "", stats::na.omit(spec$location))
  expect_true(all(roots %in% JDOC_BLOCKS))
})

test_that("only abn-internal fields are stored in the abn extension", {
  spec <- abn_json_property_spec()
  in_extension <- spec$field[!is.na(spec$location) &
                               startsWith(spec$location, "metadata.extensions.abn")]
  expect_true(all(in_extension %in%
                    c("used.INLA", "error.code", "error.code.desc", "hessian.accuracy")))
})

test_that("known field decisions are encoded", {
  expect_equal(jspec_row("call", "mle")$reason, "excluded")
  expect_equal(jspec_row("group.ids", "mle_grouped")$reason, "data")
  expect_equal(jspec_row("method", "bayes")$location, "inference.type")
  expect_equal(jspec_row("priors", "bayes")$location, "inference.priors")
  expect_equal(jspec_row("centre", "mle")$location, "variables[].transform")
  expect_equal(jspec_row("marginals", "bayes")$location, "inference.posterior.marginals")
  expect_equal(jspec_row("modes", "bayes")$location, "parameters[].value")
  expect_equal(jspec_row("coef", "mle")$location, "parameters[].value")
  expect_equal(jspec_row("coef", "bayes")$reason, "derived")
  expect_equal(jspec_row("mse", "bayes")$reason, "derived")
  expect_equal(jspec_row("abnDag", "mle")$location, "arcs")
})

test_that("the specification covers every field of every fixture fit", {
  spec <- abn_json_property_spec()
  for (name in names(jfx_registry())) {
    jfx_skip_bayes(name)
    fit <- jfx_fit(name)
    type <- jspec_fit_type(name)
    missing <- setdiff(names(fit), spec$field[spec$fit_type == type])
    expect_equal(missing, character(0), info = paste("fixture:", name))
  }
})
