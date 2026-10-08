# Export: metadata.extensions.abn contains only abn-internal information.

test_that("the abn extension only contains the allowed keys", {
  for (name in jfx_names()) {
    jfx_skip_bayes(name)
    ext <- jfx_doc(name)$metadata[["extensions"]][["abn"]]
    expect_true(all(names(ext) %in% JDOC_EXTENSION_KEYS), info = name)
    expect_null(ext[["native_fields"]])
    expect_null(ext[["native_presence"]])
  }
})

test_that("the extension maps every parameter to its native abn names", {
  for (name in jfx_names()) {
    jfx_skip_bayes(name)
    doc <- jfx_doc(name)
    names_map <- doc[["metadata"]][["extensions"]][["abn"]][["parameter_names"]]
    param_ids <- jdoc_ids(doc[["parameters"]])
    mapped <- vapply(names_map, function(x) jdoc_chr(x[["parameter"]]), character(1))
    expect_true(all(mapped %in% param_ids), info = name)
    expect_setequal(unique(mapped), param_ids)
    fields <- vapply(names_map, function(x) x[["field"]], character(1))
    expect_true(all(fields %in% c("coef", "Stderror", "mse", "mu", "betas", "sigma",
                                  "sigma_alpha", "modes")), info = name)
    for (x in names_map) expect_false(any(c("value", "values") %in% names(x)))
  }
})

test_that("native names in the extension agree with the fit", {
  fit <- jfx_fit("g2b2c_mle")
  doc <- jfx_doc("g2b2c_mle")
  names_map <- doc[["metadata"]][["extensions"]][["abn"]][["parameter_names"]]
  coef_names <- Filter(function(x) identical(x[["field"]], "coef"), names_map)
  for (x in coef_names) {
    p <- Filter(function(p) identical(jdoc_chr(p$`_id`), jdoc_chr(x[["parameter"]])),
                doc[["parameters"]])[[1]]
    target <- jdoc_variable_by_id(doc, p[["target"]])$name
    expect_equal(p[["value"]], unname(fit$coef[[target]][1, x[["name"]]]), info = x[["name"]])
  }
})

test_that("bayes node flags are stored per variable in the extension", {
  jfx_skip_if_no_bayes()
  fit <- jfx_fit("ex1_bayes")
  doc <- jfx_doc("ex1_bayes")
  nodes <- doc[["metadata"]][["extensions"]][["abn"]][["nodes"]]
  expect_length(nodes, length(fit$modes))
  na_to_null <- function(x) if (length(x) == 1 && is.na(x)) NULL else unname(x)
  for (n in nodes) {
    node <- jdoc_variable_by_id(doc, n[["variable"]])$name
    expect_equal(n[["used_inla"]], na_to_null(fit$used.INLA[[node]]))
    expect_equal(n[["error_code"]], na_to_null(fit$error.code[[node]]))
    expect_equal(n[["error_code_desc"]], na_to_null(fit$error.code.desc[[node]]))
  }
})

test_that("MLE exports carry no bayes node flags", {
  expect_null(jfx_doc("ex1_mle")$metadata[["extensions"]][["abn"]][["nodes"]])
})
