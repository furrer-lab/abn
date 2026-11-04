# Foreign documents: these tests exercise import without fitting or INLA.
contract_bayes_document <- function() {
  doc <- jdoc_foreign_document()
  doc$inference$type <- "bayesian"
  doc$parameters[[2]]$scale <- "precision"
  doc$parameters[[2]]$value <- 4
  doc
}

test_that("posterior quantiles import independently of marginal densities", {
  for (densities in c(FALSE, TRUE)) {
    for (quantiles in c(FALSE, TRUE)) {
      doc <- contract_bayes_document()
      if (densities) {
        doc$inference$posterior <- list(marginals = list(
          list(parameter = 1, x = list(0, 1, 2), density = list(0.1, 0.8, 0.1))))
      }
      if (quantiles) {
        doc$parameters[[1]]$uncertainty$posterior_quantiles <- list(
          list(probability = 0.25, value = 0.7),
          list(probability = 0.75, value = 1.4))
      }
      fit <- import_abnFit(json = jdoc_serialize(doc))
      name <- "height|(Intercept)"
      expect_identical(!is.null(fit$marginals), densities)
      expect_identical(!is.null(fit$marginal.quantiles), quantiles)
      if (densities) {
        expect_equal(unname(fit$marginals$height[[name]]),
                     matrix(c(0, 1, 2, 0.1, 0.8, 0.1), ncol = 2))
      }
      if (quantiles) {
        expect_equal(unname(fit$marginal.quantiles$height[[name]]),
                     matrix(c(0.25, 0.75, 0.7, 1.4), ncol = 2))
        expect_length(fit$marginal.quantiles$status, 0L)
      }
      back <- import_abnFit(json = export_abnFit(fit))
      expect_equal(back$marginals, fit$marginals)
      expect_equal(back$marginal.quantiles, fit$marginal.quantiles)
    }
  }
})

test_that("the adapter rejects cross-scale estimates rather than reinterpreting them", {
  mle <- jdoc_foreign_document()
  mle$parameters[[2]]$scale <- "precision"
  mle$parameters[[2]]$value <- 4
  expect_error(import_abnFit(json = jdoc_serialize(mle)), "Parameter 2.*cross-scale")
  expect_error(import_abnFit(json = jdoc_serialize(mle), validate = FALSE), "cross-scale")

  bayes <- contract_bayes_document()
  bayes$parameters[[2]]$scale <- "variance"
  bayes$parameters[[2]]$value <- 0.25
  expect_error(import_abnFit(json = jdoc_serialize(bayes)), "Parameter 2.*precision")

  native <- import_abnFit(json = jdoc_serialize(contract_bayes_document()))
  expect_equal(unname(native$modes$height[["height|precision"]]), 4)
  expect_equal(unname(native$mse[["height"]]), 0.25)
})

test_that("variance and precision parameters obey their scalar domains", {
  doc <- jdoc_foreign_document()
  doc$parameters[[2]]$value <- -0.5
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "non-negative variance")
  doc <- contract_bayes_document()
  doc$parameters[[2]]$value <- 0
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "strictly positive precision")
})

test_that("multiple groups are rejected as an adapter limitation", {
  doc <- jdoc_foreign_document()
  doc$groups <- list(list(`_id` = 1, name = "farm", variables = list(1, 2)),
                     list(`_id` = 2, name = "year", variables = list(1, 2)))
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "adapter supports at most one group")
  expect_error(import_abnFit(json = jdoc_serialize(doc), validate = FALSE),
               "adapter supports at most one group")
})

test_that("random parameters must target a member of their referenced group", {
  doc <- jdoc_foreign_document()
  doc$groups <- list(list(`_id` = 1, name = "farm", variables = list(2)))
  doc$parameters[[5]] <- list(`_id` = 5, target = 1, group = 1,
                             kind = "random_variance", value = 0.5, scale = "variance")
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "outside its group")
})

test_that("missing value fields and nearly integral ids are not silently accepted", {
  doc <- jdoc_foreign_document()
  doc$parameters[[1]]$value <- NULL
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "missing required field.*value")
  doc <- jdoc_foreign_document()
  doc$variables[[1]]$`_id` <- 1 + 1e-10
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "invalid.*_id")
})

test_that("singleton posterior grids and quantile collections remain JSON arrays", {
  doc <- contract_bayes_document()
  doc$inference$posterior <- list(marginals = list(
    list(parameter = 1, x = list(1.2), density = list(1))))
  doc$parameters[[1]]$uncertainty$posterior_quantiles <- list(
    list(probability = 0.5, value = 1.2))
  fit <- import_abnFit(json = jdoc_serialize(doc))
  exported <- jdoc_parse(export_abnFit(fit))
  grid <- exported$inference$posterior$marginals[[1]]
  expect_type(grid$x, "list")
  expect_length(grid$x, 1L)
  expect_type(grid$density, "list")
  quantiles <- exported$parameters[[1]]$uncertainty$posterior_quantiles
  expect_type(quantiles, "list")
  expect_length(quantiles, 1L)
  jdoc_expect_valid(jdoc_serialize(exported))
})
