json_import_generic_document <- function() {
  list(
    metadata = list(
      schema_version = "bayesian-network-v1",
      issuer = "independent-bn-tool::export",
      configs = list(),
      extensions = list()
    ),
    variables = list(
      list(id = "height", name = "height", type = "continuous"),
      list(
        id = "status",
        name = "status",
        type = "binary",
        states = list(
          list(id = "no", label = "no"),
          list(id = "yes", label = "yes")
        )
      )
    ),
    arcs = list(list(source = "height", target = "status")),
    parameters = list(
      list(
        id = "height-intercept",
        target = "height",
        kind = "intercept",
        link = "identity",
        value = 1.2,
        uncertainty = list(standard_error = 0.1)
      ),
      list(
        id = "status-intercept",
        target = "status",
        kind = "intercept",
        link = "logit",
        value = -0.42,
        uncertainty = list(standard_error = 0.2)
      ),
      list(
        id = "status-height",
        target = "status",
        parents = list("height"),
        kind = "coefficient",
        link = "logit",
        value = 0.87,
        uncertainty = list(standard_error = 0.3)
      )
    ),
    inference = list(
      type = "maximum_likelihood",
      estimates = list(),
      uncertainty = list(),
      diagnostics = list()
    )
  )
}

json_import_string <- function(document) {
  jsonlite::toJSON(document, auto_unbox = TRUE, null = "null", digits = NA)
}

test_that("import_abnFit imports a genuine generic network document", {
  imported <- import_abnFit(json = json_import_string(json_import_generic_document()))

  expect_s3_class(imported, "abnFit")
  expect_s3_class(imported$abnDag, "abnDag")
  expect_equal(rownames(imported$abnDag$dag), c("height", "status"))
  expect_equal(colnames(imported$abnDag$dag), c("height", "status"))
  expect_equal(imported$abnDag$dag["height", "status"], 1)
  expect_equal(unname(unlist(imported$abnDag$data.dists)),
               c("gaussian", "binomial"))
  expect_equal(imported$method, "mle")
  expect_equal(unname(imported$coef$height[1, "height|intercept"]), 1.2)
  expect_equal(unname(imported$coef$status[1, "height"]), 0.87)
  expect_equal(unname(imported$Stderror$status[1, "height"]), 0.3)
})

test_that("import_abnFit ingests ABN-specific fit properties", {
  fit <- json_fixture_fit_ex1_bayes_fixed()
  document <- json_import_generic_document()
  document$metadata$issuer <- "abn::export_abnFit"
  document$metadata$extensions <- list(
    abn = list(
      configs = list(
        method = "bayes",
        compute_fixed = TRUE,
        used_inla = unname(fit$used.INLA)
      ),
      inference = list(
        mlik = fit$mlik,
        mliknode = fit$mliknode,
        error_code = fit$error.code,
        error_code_desc = fit$error.code.desc,
        hessian_accuracy = fit$hessian.accuracy,
        marginals = fit$marginals,
        marginal_quantiles = fit$marginal.quantiles
      )
    )
  )
  document$inference$type <- "bayesian"

  imported <- import_abnFit(json = json_import_string(document))

  expect_equal(imported$method, "bayes")
  expect_equal(imported$mlik, fit$mlik)
  expect_equal(unname(imported$mliknode), unname(fit$mliknode))
  expect_equal(unname(imported$used.INLA), unname(fit$used.INLA))
  expect_equal(unname(imported$error.code), unname(fit$error.code))
  expect_equal(unname(imported$error.code.desc), unname(fit$error.code.desc))
  expect_equal(unname(imported$hessian.accuracy), unname(fit$hessian.accuracy))
  expect_length(imported$marginals, length(fit$marginals))
  expect_length(imported$marginal.quantiles, length(fit$marginal.quantiles))
})

test_that("ABN grouped properties are read from the ABN extension", {
  load(test_path("testdata", "abnfit_mle_groups.Rdata"))
  fit <- abn_fit
  document <- jsonlite::fromJSON(export_abnFit(fit), simplifyVector = FALSE)

  document$metadata$issuer <- "abn::export_abnFit"
  document$metadata$extensions$abn$configs$group_var <- fit$group.var
  document$metadata$extensions$abn$configs$group_ids <- fit$group.ids
  document$metadata$extensions$abn$configs$grouped_vars <- fit$grouped.vars

  imported <- import_abnFit(json = json_import_string(document))

  expect_identical(imported$group.var, fit$group.var)
  expect_equal(unname(imported$group.ids), unname(fit$group.ids))
  expect_equal(unname(imported$grouped.vars), unname(fit$grouped.vars))
  expect_true(all(c("mu", "betas", "sigma", "sigma_alpha") %in%
                    names(imported)))
})

test_that("non-ABN extensions are ignored without changing generic import", {
  document <- json_import_generic_document()
  document$metadata$extensions <- list(
    another_bn_tool = list(
      configs = list(private_engine_setting = "ignore-me"),
      variables = list(status = list(private_color = "red")),
      parameters = list(`status-height` = list(private_precision = 4)),
      inference = list(private_trace = list(1, 2, 3))
    )
  )

  imported <- import_abnFit(json = json_import_string(document))

  expect_s3_class(imported, "abnFit")
  expect_equal(imported$method, "mle")
  expect_false("private_engine_setting" %in% names(imported))
  expect_false("private_color" %in% names(imported))
  expect_false("private_precision" %in% names(imported))
  expect_false("private_trace" %in% names(imported))
})

test_that("missing extension IDs are ignored with warnings", {
  document <- json_import_generic_document()
  document$metadata$extensions <- list(
    another_bn_tool = list(
      variables = list(does_not_exist = list(private_color = "red")),
      parameters = list(missing_parameter = list(private_precision = 4))
    )
  )

  expect_warning(
    imported <- import_abnFit(json = json_import_string(document)),
    "unknown|missing|not found"
  )
  expect_s3_class(imported, "abnFit")
})

test_that("network-less generic documents are rejected by the abn adapter", {
  document <- json_import_generic_document()
  document$variables <- NULL
  document$arcs <- NULL

  expect_error(
    import_abnFit(json = json_import_string(document)),
    "variables|network|DAG"
  )
})
