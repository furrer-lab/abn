# Serialization fixtures, not fitted estimates. Keep these tests independent of
# the known separated g2pbcgrp regression and its upstream fitting failure.
json_grouped_multinomial_fixture <- function() {
  dat <- data.frame(
    X = c(-1, 0, 1, 2, 0.5, -0.5),
    C = factor(c("a", "b", "c", "a", "b", "c"), levels = c("a", "b", "c")),
    group = factor(rep(c("g1", "g2"), each = 3))
  )
  dists <- list(X = "gaussian", C = "multinomial")
  dag <- jfx_empty_dag(dists)
  dag["C", "X"] <- 1
  covariance <- matrix(c(0.7, -0.12, -0.12, 0.4), 2, 2,
                        dimnames = list(c("C.b~1", "C.c~1"), c("C.b~1", "C.c~1")))
  fit <- structure(list(
    method = "mle",
    abnDag = structure(list(dag = dag, data.df = dat[c("X", "C")],
                            data.dists = dists), class = "abnDag"),
    levels = list(C = levels(dat$C)),
    group.var = "group",
    grouped.vars = c(1L, 2L),
    group.ids = as.integer(dat$group),
    mu = list(X = 0.1, C = c(C.b = 0.3, C.c = -0.2)),
    betas = list(X = stats::setNames(numeric(0), character(0)),
                 C = matrix(c(0.6, -0.4), 2, 1,
                            dimnames = list(c("b", "c"), "X"))),
    sigma = list(X = 0.5, C = NA),
    sigma_alpha = list(X = 0.2, C = covariance)
  ), class = "abnFit")
  # No fabricated likelihoods: scores are intentionally unavailable for this
  # hand-constructed fixture. Only the parameter and metadata contract is tested.
  list(fit = fit, data = dat)
}

test_that("grouped multinomial JSON maps each fixed effect to its child state", {
  fixture <- json_grouped_multinomial_fixture()
  doc <- jdoc_parse(export_abnFit(fixture$fit))
  for (state in c("b", "c")) {
    intercept <- jdoc_parameter(doc, "C", "intercept", target_state = state)
    slope <- jdoc_parameter(doc, "C", "coefficient", parent = "X", target_state = state)
    expect_equal(intercept$value, unname(fixture$fit$mu$C[[paste0("C.", state)]]))
    expect_equal(slope$value, unname(fixture$fit$betas$C[state, "X"]))
    expected <- intercept$value + slope$value * fixture$data$X
    expect_equal(jdoc_linear_predictor(doc, fixture$data, "C", state), expected)
  }
  expect_length(jdoc_parameters(doc, "C", "intercept"), 2L)
  expect_length(jdoc_parameters(doc, "C", "coefficient"), 2L)
})

test_that("grouped multinomial JSON preserves covariance signs and native shapes", {
  fixture <- json_grouped_multinomial_fixture()
  json <- export_abnFit(fixture$fit)
  doc <- jdoc_parse(json)
  covariance <- jdoc_parameters(doc, "C", "random_covariance")
  expect_length(covariance, 1L)
  expect_equal(covariance[[1]]$value, -0.12)
  expect_equal(covariance[[1]]$scale, "variance")
  for (state in c("b", "c")) {
    variance <- jdoc_parameter(doc, "C", "random_variance", target_state = state)
    expect_equal(variance$value, unname(fixture$fit$sigma_alpha$C[
      paste0("C.", state, "~1"), paste0("C.", state, "~1")]))
  }
  back <- import_abnFit(json = json, data = fixture$data)
  for (field in c("mu", "betas", "sigma", "sigma_alpha", "levels",
                  "group.var", "grouped.vars", "group.ids")) {
    expect_equal(back[[field]], fixture$fit[[field]], info = field)
  }
  expect_equal(back$abnDag, fixture$fit$abnDag)
})

test_that("grouped multinomial covariance and effects import from the generic core", {
  fixture <- json_grouped_multinomial_fixture()
  doc <- jdoc_drop_extensions(jdoc_parse(export_abnFit(fixture$fit)))
  back <- import_abnFit(json = jdoc_serialize(doc))
  expect_equal(unname(back$mu$C), unname(fixture$fit$mu$C))
  expect_equal(back$betas$C, fixture$fit$betas$C)
  expect_equal(unname(back$sigma_alpha$C), unname(fixture$fit$sigma_alpha$C))
  expect_equal(dim(back$sigma_alpha$C), c(2L, 2L))
  expect_equal(back$group.var, "group")
  expect_equal(back$abnDag$dag, fixture$fit$abnDag$dag)
})

test_that("constructed grouped multinomial exports conform to the network schema", {
  fixture <- json_grouped_multinomial_fixture()
  jdoc_expect_valid(export_abnFit(fixture$fit))
  core <- jdoc_drop_extensions(jdoc_parse(export_abnFit(fixture$fit)))
  jdoc_expect_valid(jdoc_serialize(core))
})
