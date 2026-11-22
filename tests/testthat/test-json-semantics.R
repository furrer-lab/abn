# Semantics: the generic core alone (extensions removed) must reproduce the
# model. Linear predictors computed from the JSON parameters are compared with
# independent reference fits (glm, lm, nnet, lme4) on the raw data.

jsem_family <- function(dist) {
  switch(dist, gaussian = stats::gaussian(), binomial = stats::binomial(),
         poisson = stats::poisson())
}

jsem_formula <- function(node, parents, random = NULL) {
  rhs <- if (length(parents) == 0) "1" else paste(parents, collapse = " + ")
  if (!is.null(random)) rhs <- paste(rhs, "+ (1 |", random, ")")
  stats::as.formula(paste(node, "~", rhs))
}

jsem_core <- function(name) jdoc_drop_extensions(jfx_doc(name))

jsem_standardise <- function(data, dists) {
  for (node in names(dists)) {
    if (identical(dists[[node]], "gaussian")) {
      data[[node]] <- (data[[node]] - mean(data[[node]])) / stats::sd(data[[node]])
    }
  }
  data
}

# Compare every non-multinomial node with a glm fitted on `ref_data`.
jsem_expect_glm_nodes <- function(name, ref_data = jfx_spec(name)$data) {
  spec <- jfx_spec(name)
  fit <- jfx_fit(name)
  core <- jsem_core(name)
  for (node in names(spec$dists)) {
    dist <- spec$dists[[node]]
    if (dist == "multinomial") next
    parents <- jdoc_parents(fit, node)
    ref <- stats::glm(jsem_formula(node, parents), family = jsem_family(dist),
                      data = ref_data)
    expect_equal(jdoc_linear_predictor(core, spec$data, node),
                 unname(stats::predict(ref, type = "link")),
                 tolerance = 1e-4, info = paste(name, node))
    if (dist == "gaussian") {
      expect_equal(jdoc_parameter(core, node, "residual_variance")$value,
                   stats::sigma(ref)^2, tolerance = 1e-6, info = paste(name, node))
    }
  }
}

test_that("core reproduces glm linear predictors (binary and continuous parents)", {
  jsem_expect_glm_nodes("ex1_mle")
})

test_that("core reproduces glm linear predictors for centred fits via transform", {
  spec <- jfx_spec("ex1_mle_centred")
  jsem_expect_glm_nodes("ex1_mle_centred",
                        ref_data = jsem_standardise(spec$data, spec$dists))
})

test_that("core reproduces glm linear predictors with multinomial parents", {
  jsem_expect_glm_nodes("fcv_mle")
  jsem_expect_glm_nodes("g2b2c_mle")
})

test_that("core reproduces baseline-category logits of a multinomial child", {
  spec <- jfx_spec("g2b2c_soft_mle")
  core <- jsem_core("g2b2c_soft_mle")
  ref <- nnet::multinom(C ~ B1 + B2, data = spec$data, trace = FALSE)
  probs <- stats::predict(ref, type = "probs")
  for (state in c("b", "c")) {
    expect_equal(jdoc_linear_predictor(core, spec$data, "C", state),
                 unname(log(probs[, state] / probs[, "a"])),
                 tolerance = 1e-3, info = state)
  }
})

test_that("core reproduces lme4 fixed effects and variance components", {
  spec <- jfx_spec("adg_mle_grouped")
  fit <- jfx_fit("adg_mle_grouped")
  core <- jsem_core("adg_mle_grouped")
  for (node in names(spec$dists)) {
    dist <- spec$dists[[node]]
    f <- jsem_formula(node, jdoc_parents(fit, node), random = "farm")
    ref <- if (dist == "gaussian") {
      lme4::lmer(f, data = spec$data)
    } else {
      lme4::glmer(f, data = spec$data, family = jsem_family(dist))
    }
    info <- paste("adg", node)
    expect_equal(jdoc_linear_predictor(core, spec$data, node),
                 unname(stats::predict(ref, re.form = NA, type = "link")),
                 tolerance = 1e-3, info = info)
    vc <- as.data.frame(lme4::VarCorr(ref))
    expect_equal(jdoc_parameter(core, node, "random_variance")$value,
                 vc$vcov[vc$grp == "farm"], tolerance = 1e-3, info = info)
    if (dist == "gaussian") {
      expect_equal(jdoc_parameter(core, node, "residual_variance")$value,
                   vc$vcov[vc$grp == "Residual"], tolerance = 1e-3, info = info)
    }
  }
})

test_that("core-only bayes documents convert precisions consistently", {
  jfx_skip_if_no_bayes()
  fit <- jfx_fit("adg_bayes_grouped")
  core <- jsem_core("adg_bayes_grouped")
  resid <- jdoc_parameter(core, "adg", "residual_variance")
  random <- jdoc_parameter(core, "adg", "random_variance")
  expect_equal(resid[["scale"]], "precision")
  expect_equal(random[["scale"]], "precision")
  expect_equal(resid[["value"]], unname(fit$modes$adg[["adg|precision"]]))
  expect_equal(random[["value"]], unname(fit$modes$adg[["adg|group.precision"]]))
})
