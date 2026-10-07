# Phase 0: model semantics the JSON format relies on, and required fixes in the
# fitting code (U1-U6, see json_plan.md). These tests check abn's own fit objects
# against independent reference implementations and do not touch JSON.

# --- Established semantics ---------------------------------------------------

test_that("binomial nodes model the second factor level with a logit link", {
  spec <- jfx_spec("ex1_mle")
  fit <- jfx_fit("ex1_mle")
  ref <- stats::glm(b3 ~ b1 + g1 + b2, family = stats::binomial(), data = spec$data)
  expect_equal(unname(fit$coef$b3[1, ]), unname(stats::coef(ref)), tolerance = 1e-4)
})

test_that("poisson nodes use a log link", {
  spec <- jfx_spec("ex1_mle")
  fit <- jfx_fit("ex1_mle")
  ref <- stats::glm(p2 ~ b1 + p1, family = stats::poisson(), data = spec$data)
  expect_equal(unname(fit$coef$p2[1, ]), unname(stats::coef(ref)), tolerance = 1e-4)
})

test_that("ungrouped gaussian mse is the unbiased residual variance", {
  spec <- jfx_spec("ex1_mle")
  fit <- jfx_fit("ex1_mle")
  ref <- stats::lm(g2 ~ p1 + g1 + b2, data = spec$data)
  expect_equal(unname(fit$mse[["g2"]]), stats::sigma(ref)^2, tolerance = 1e-6)
})

test_that("a multinomial parent is one-hot encoded over all levels without intercept", {
  spec <- jfx_spec("fcv_mle")
  fit <- jfx_fit("fcv_mle")
  n_levels <- nlevels(spec$data$Sex)
  expect_false(any(grepl("intercept", colnames(fit$coef$Outdoor))))
  expect_equal(ncol(fit$coef$Outdoor), n_levels)
  ref <- stats::glm(Outdoor ~ -1 + Sex, family = stats::binomial(), data = spec$data)
  expect_equal(unname(fit$coef$Outdoor[1, ]), unname(stats::coef(ref)), tolerance = 1e-4)
})

test_that("a multinomial child uses the first level as baseline", {
  spec <- jfx_spec("g2b2c_mle")
  fit <- jfx_fit("g2b2c_mle")
  ref <- nnet::multinom(C ~ B1 + B2, data = spec$data, trace = FALSE)
  non_baseline <- levels(spec$data$C)[-1]
  expect_equal(rownames(stats::coef(ref)), non_baseline)
  intercepts <- fit$coef$C[1, paste0("C|intercept.", non_baseline)]
  expect_equal(unname(intercepts), unname(stats::coef(ref)[, "(Intercept)"]),
               tolerance = 1e-3)
})

test_that("grouped gaussian sigma and sigma_alpha are standard deviations", {
  spec <- jfx_spec("adg_mle_grouped")
  fit <- jfx_fit("adg_mle_grouped")
  ref <- lme4::lmer(adg ~ age + wormCount + (1 | farm), data = spec$data)
  vc <- as.data.frame(lme4::VarCorr(ref))
  expect_equal(unname(fit$sigma_alpha$adg)^2, vc$vcov[vc$grp == "farm"],
               tolerance = 1e-4)
  expect_equal(unname(fit$sigma$adg)^2, vc$vcov[vc$grp == "Residual"],
               tolerance = 1e-4)
})

test_that("grouped multinomial sigma_alpha is a covariance matrix", {
  fit <- jfx_fit("g2pbcgrp_mle_grouped")
  spec <- jfx_spec("g2pbcgrp_mle_grouped")
  k <- nlevels(spec$data$C) - 1
  expect_true(is.matrix(fit$sigma_alpha$C))
  expect_equal(dim(fit$sigma_alpha$C), c(k, k))
  expect_equal(fit$sigma_alpha$C, t(fit$sigma_alpha$C))
})

test_that("bayes modes carry gaussian precisions", {
  jfx_skip_if_no_bayes()
  fit <- jfx_fit("ex1_bayes")
  expect_true("g1|precision" %in% names(fit$modes$g1))
  expect_gt(fit$modes$g1[["g1|precision"]], 0)
})

# --- U1: centring is recorded -----------------------------------------------

test_that("U1: centred MLE fits record centre and scale of gaussian nodes", {
  spec <- jfx_spec("ex1_mle_centred")
  fit <- jfx_fit("ex1_mle_centred")
  expect_setequal(names(fit$centre), c("g1", "g2"))
  for (node in c("g1", "g2")) {
    expect_equal(fit$centre[[node]][["center"]], mean(spec$data[[node]]))
    expect_equal(fit$centre[[node]][["scale"]], stats::sd(spec$data[[node]]))
  }
})

test_that("U1: uncentred fits record no centring", {
  expect_null(jfx_fit("ex1_mle")$centre)
})

# --- U2: centre is honoured by Bayes fits ------------------------------------

test_that("U2: bayes fits honour centre = FALSE", {
  jfx_skip_if_no_bayes()
  spec <- jfx_spec("ex1_bayes")
  uncentred <- jfx_fit("ex1_bayes")
  centred <- jfx_fit("ex1_bayes_centred")
  expect_equal(uncentred$modes$g1[["g1|(Intercept)"]], mean(spec$data$g1),
               tolerance = 1e-2)
  expect_equal(centred$modes$g1[["g1|(Intercept)"]], 0, tolerance = 1e-2)
  expect_null(uncentred$centre)
  expect_setequal(names(centred$centre), "g1")
})

# --- U3: Bayes fits store the grouping --------------------------------------

test_that("U3: grouped bayes fits store group.var, group.ids and grouped.vars", {
  jfx_skip_if_no_bayes()
  bayes <- jfx_fit("ex3_bayes_grouped")
  mle <- jfx_fit("ex3_mle_grouped")
  expect_equal(bayes$group.var, "group")
  expect_equal(bayes$group.ids, mle$group.ids)
  expect_equal(bayes$grouped.vars, mle$grouped.vars)
})

# --- U4: Bayes fits store the priors ------------------------------------------

test_that("U4: bayes fits store the priors used", {
  jfx_skip_if_no_bayes()
  expect_equal(jfx_fit("ex1_bayes")$priors,
               list(mean = 0, prec = 0.001, loggam.shape = 1, loggam.inv.scale = 5e-05))
  expect_equal(jfx_fit("ex1_bayes_priors")$priors,
               list(mean = 0.5, prec = 0.01, loggam.shape = 2, loggam.inv.scale = 1e-3))
})

test_that("U4: mle fits store no priors", {
  expect_null(jfx_fit("ex1_mle")$priors)
})

# --- U5: modes2coefs drops all precisions -----------------------------------

test_that("U5: coef of grouped gaussian bayes nodes contains no precision", {
  jfx_skip_if_no_bayes()
  fit <- jfx_fit("adg_bayes_grouped")
  expect_false(any(grepl("precision", colnames(fit$coef$adg))))
  expect_equal(ncol(fit$coef$adg), 2L)
  expect_true(all(c("adg|precision", "adg|group.precision") %in% names(fit$modes$adg)))
})

# --- U7: residual variance with multinomial parents -------------------------

test_that("U7: gaussian mse uses the correct df with a multinomial parent", {
  spec <- jfx_spec("g2b2c_mle")
  fit <- jfx_fit("g2b2c_mle")
  ref <- stats::lm(G2 ~ G1 + C, data = spec$data)
  expect_equal(unname(fit$df[["G2"]]), stats::df.residual(ref))
  expect_equal(unname(fit$mse[["G2"]]), stats::sigma(ref)^2, tolerance = 1e-6)
})

# --- U6: multinomial coefficient names match their values -------------------

test_that("U6: multinomial child coefficient names follow the value order", {
  spec <- jfx_spec("g2b2c_mle")
  fit <- jfx_fit("g2b2c_mle")
  expect_equal(levels(spec$data$C), c("a", "b", "c"))
  expect_equal(colnames(fit$coef$C),
               c("C|intercept.b", "C|intercept.c", "B1b", "B1c", "B2b", "B2c"))
  expect_equal(colnames(fit$Stderror$C), colnames(fit$coef$C))

  ref <- nnet::multinom(C ~ B1 + B2, data = spec$data, trace = FALSE)
  ref_coef <- stats::coef(ref)
  for (state in c("b", "c")) {
    expect_equal(unname(fit$coef$C[1, paste0("B1", state)]), unname(ref_coef[state, "B11"]),
                 tolerance = 1e-3)
    expect_equal(unname(fit$coef$C[1, paste0("B2", state)]), unname(ref_coef[state, "B21"]),
                 tolerance = 1e-3)
  }
})

# --- U8: original factor levels are recorded --------------------------------

test_that("U8: fits record the factor levels of binomial and multinomial nodes", {
  # fcv_mle: Sex has non-alphabetical factor levels (m, mc, f, fc)
  for (name in c("ex1_mle", "fcv_mle", "g2pbcgrp_mle_grouped", "ex1_bayes")) {
    jfx_skip_bayes(name)
    spec <- jfx_spec(name)
    fit <- jfx_fit(name)
    categorical <- names(spec$dists)[unlist(spec$dists) %in% c("binomial", "multinomial")]
    expect_setequal(names(fit$levels), categorical)
    for (node in categorical) {
      x <- spec$data[[node]]
      # ungrouped MLE multinomial nodes: nnet re-derives levels alphabetically
      expected <- if (spec$method == "mle" && !spec$grouped &&
                        spec$dists[[node]] == "multinomial") {
        levels(factor(as.character(x)))
      } else {
        levels(factor(x))
      }
      expect_equal(fit$levels[[node]], expected, info = paste(name, node))
    }
  }
})

# --- U9: error.code.desc keeps node names -----------------------------------

test_that("U9: bayes error.code.desc is named by node", {
  jfx_skip_if_no_bayes()
  fit <- jfx_fit("ex1_bayes")
  expect_equal(names(fit$error.code.desc), names(fit$error.code))
  expect_equal(names(fit$error.code.desc), colnames(fit$abnDag$dag))
})

# --- #272: ungrouped bayes cache via the object= path -----------------------

test_that("#272: ungrouped bayes fit from cache has no grouping fields", {
  jfx_skip_if_no_bayes()
  df <- FCV[, c(12, 14:15)]
  mydists <- list(Outdoor = "binomial", GroupSize = "poisson", Age = "gaussian")
  cache <- suppressWarnings(buildScoreCache(data.df = df, data.dists = mydists,
                                            method = "bayes", max.parents = 1))
  expect_true("group.var" %in% names(cache))
  expect_null(cache[["group.var"]])
  mp <- mostProbable(score.cache = cache, verbose = FALSE)
  fit <- suppressWarnings(fitAbn(object = mp, method = "bayes", centre = FALSE))
  expect_null(fit$group.var)
  expect_null(fit$group.ids)
  expect_null(fit$grouped.vars)
})
