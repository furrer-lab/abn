test_that("generic MLE export imports as an equivalent abnFit", {
  original <- json_fixture_fit_ex1_mle()
  imported <- import_abnFit(json = export_abnFit(original))

  expect_true(abnfit_objects_equal(original, imported))
  expect_equal(imported$method, "mle")
  expect_equal(imported$abnDag$dag, original$abnDag$dag)
})

test_that("generic grouped MLE export imports ABN grouped fields", {
  load(test_path("testdata", "abnfit_mle_groups.Rdata"))
  original <- abn_fit
  imported <- import_abnFit(json = export_abnFit(original))

  expect_true(abnfit_objects_equal(original, imported))
  expect_identical(imported$group.var, original$group.var)
  for (field in c("mu", "betas", "sigma", "sigma_alpha")) {
    expect_true(field %in% names(imported))
  }
})

test_that("generic Bayesian export imports ABN inference extensions", {
  original <- json_fixture_fit_ex1_bayes()
  imported <- import_abnFit(json = export_abnFit(original))

  expect_equal(imported$method, "bayes")
  expect_equal(imported$mlik, original$mlik)
  expect_equal(unname(imported$mliknode), unname(original$mliknode))
  expect_equal(unname(imported$used.INLA), unname(original$used.INLA))
  expect_equal(unname(imported$error.code), unname(original$error.code))
  expect_equal(unname(imported$error.code.desc), unname(original$error.code.desc))
  expect_equal(unname(imported$hessian.accuracy),
               unname(original$hessian.accuracy))
})

test_that("generic Bayesian posterior extensions import into native fields", {
  original <- json_fixture_fit_ex1_bayes_fixed()
  imported <- import_abnFit(json = export_abnFit(original))

  expect_equal(length(imported$marginals), length(original$marginals))
  expect_equal(length(imported$marginal.quantiles),
               length(original$marginal.quantiles))
})

test_that("import ignores non-ABN extension properties", {
  document <- json_import_generic_document()
  document$metadata$extensions <- list(
    other_tool = list(
      variables = list(`2` = list(display_color = "red")),
      parameters = list(`3` = list(display_precision = 4)),
      inference = list(private_trace = list(1, 2, 3))
    )
  )

  imported <- import_abnFit(json = json_import_string(document))

  expect_s3_class(imported, "abnFit")
  expect_false("display_color" %in% names(imported))
  expect_false("display_precision" %in% names(imported))
  expect_false("private_trace" %in% names(imported))
})

test_that("import warns but ignores extension references to missing IDs", {
  document <- json_import_generic_document()
  document$metadata$extensions <- list(
    other_tool = list(
      variables = list(missing_variable = list(display_color = "red")),
      parameters = list(missing_parameter = list(display_precision = 4))
    )
  )

  expect_warning(
    imported <- import_abnFit(json = json_import_string(document)),
    "unknown|missing|not found"
  )
  expect_s3_class(imported, "abnFit")
})

test_that("import rejects a network document without network structure", {
  document <- json_import_generic_document()
  document$variables <- NULL
  document$arcs <- NULL

  expect_error(
    import_abnFit(json = json_import_string(document)),
    "variables|network|DAG"
  )
})
