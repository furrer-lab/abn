# Exercise the INLA-produced JSON path explicitly, not just namespace loading
# or a fit which silently fell back to the internal Bayesian implementation.
test_that("INLA-produced Bayesian fits preserve engine flags and posterior fields", {
  jfx_skip_if_no_bayes()
  spec <- jfx_spec("adg_bayes_grouped")
  spec$extra <- list(
    compute.fixed = TRUE,
    control = fit.control(method = "bayes", max.mode.error = 100, n.grid = 100)
  )
  fit <- jfx_fit_spec(spec)
  expect_true(all(fit$used.INLA %in% TRUE),
              info = "This test requires actual INLA results, not internal-code fallback")
  expect_length(fit$used.INLA, length(spec$dists))
  expect_true(length(fit$marginals) > 0L)

  json <- export_abnFit(fit)
  doc <- jdoc_parse(json)
  for (node in doc$metadata$extensions$abn$nodes) {
    expect_true(node$used_inla)
  }
  expect_true(length(doc$inference$posterior$marginals) > 0L)
  back <- import_abnFit(json = json)
  for (field in c("used.INLA", "modes", "marginals", "marginal.quantiles",
                  "priors", "mliknode", "mlik", "group.var", "grouped.vars")) {
    expect_equal(back[[field]], fit[[field]], info = field)
  }
  core <- jdoc_drop_extensions(doc)
  from_core <- import_abnFit(json = jdoc_serialize(core))
  for (field in c("modes", "marginals", "marginal.quantiles", "mlik")) {
    expect_equal(from_core[[field]], fit[[field]], info = field)
  }
})
