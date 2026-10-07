# Export: the generic core of the document, compared with the fit's values.

# Levels in the order abn uses them (see U8): ungrouped MLE multinomial nodes are
# re-derived alphabetically by nnet, everything else keeps the factor order.
jcore_levels <- function(spec, node) {
  x <- spec$data[[node]]
  if (spec$method == "mle" && !spec$grouped && spec$dists[[node]] == "multinomial") {
    return(levels(factor(as.character(x))))
  }
  levels(factor(x))
}

# --- metadata ---------------------------------------------------------------

test_that("metadata identifies format, issuer and version", {
  doc <- jfx_doc("ex1_mle")
  expect_equal(doc[["metadata"]][["schema_version"]], JDOC_SCHEMA_VERSION)
  expect_equal(doc[["metadata"]][["issuer"]], JDOC_ISSUER)
  expect_equal(doc[["metadata"]][["issuer_version"]], as.character(utils::packageVersion("abn")))
  created <- as.POSIXct(doc[["metadata"]][["created"]], format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
  expect_false(is.na(created))
  expect_null(doc[["metadata"]][["label"]])
  expect_null(doc[["metadata"]][["scenario_id"]])
  expect_null(doc[["metadata"]][["data_reference"]])
})

test_that("label, scenario_id and data_reference are written to metadata", {
  ref <- list(uri = "data.json", schema_version = "bn-data", sha256 = "abc123")
  doc <- jdoc_parse(export_abnFit(jfx_fit("ex1_mle"), label = "my model",
                                  scenario_id = "s-42", data_reference = ref))
  expect_equal(doc[["metadata"]][["label"]], "my model")
  expect_equal(doc[["metadata"]][["scenario_id"]], "s-42")
  expect_equal(doc[["metadata"]][["data_reference"]], ref)
})

test_that("removed pre-release arguments are rejected", {
  fit <- jfx_fit("ex1_mle")
  expect_error(export_abnFit(fit, include_network = FALSE))
  expect_error(export_abnFit(fit, format = "json"))
})

test_that("export writes to a file and returns the path invisibly", {
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path))
  expect_invisible(res <- export_abnFit(jfx_fit("ex1_mle"), file = path))
  expect_equal(res, path)
  expect_equal(jdoc_parse(paste(readLines(path), collapse = "\n"))$metadata[["schema_version"]],
               JDOC_SCHEMA_VERSION)
})

# --- structure for every fixture --------------------------------------------

test_that("all blocks are present and all references resolve", {
  for (name in names(jfx_registry())) {
    jfx_skip_bayes(name)
    doc <- jfx_doc(name)
    info <- paste("fixture:", name)
    expect_true(all(JDOC_BLOCKS %in% names(doc)), info = info)

    var_ids <- jdoc_ids(doc[["variables"]])
    expect_equal(anyDuplicated(var_ids), 0L, info = info)
    param_ids <- jdoc_ids(doc[["parameters"]])
    expect_equal(anyDuplicated(param_ids), 0L, info = info)
    group_ids <- jdoc_ids(doc[["groups"]])

    for (v in doc[["variables"]]) {
      if (!is.null(v[["states"]])) {
        expect_equal(anyDuplicated(jdoc_ids(v[["states"]])), 0L, info = info)
      }
    }
    arc_keys <- vapply(doc[["arcs"]], function(a) {
      paste(jdoc_chr(a[["source"]]), jdoc_chr(a[["target"]]))
    }, character(1))
    for (a in doc[["arcs"]]) {
      expect_true(jdoc_chr(a[["source"]]) %in% var_ids, info = info)
      expect_true(jdoc_chr(a[["target"]]) %in% var_ids, info = info)
    }
    for (p in doc[["parameters"]]) {
      expect_true(p[["kind"]] %in% JDOC_KINDS, info = info)
      expect_true(jdoc_chr(p[["target"]]) %in% var_ids, info = info)
      target <- jdoc_variable_by_id(doc, p[["target"]])
      if (!is.null(p[["target_state"]])) {
        expect_true(jdoc_chr(p[["target_state"]]) %in% jdoc_ids(target[["states"]]), info = info)
      }
      if (!is.null(p[["parent"]])) {
        expect_true(paste(jdoc_chr(p[["parent"]]), jdoc_chr(p[["target"]])) %in% arc_keys,
                    info = info)
        parent <- jdoc_variable_by_id(doc, p[["parent"]])
        if (!is.null(p[["parent_state"]])) {
          expect_true(jdoc_chr(p[["parent_state"]]) %in% jdoc_ids(parent[["states"]]), info = info)
        }
      }
      if (!is.null(p[["group"]])) {
        expect_true(jdoc_chr(p[["group"]]) %in% group_ids, info = info)
      }
      if (p[["kind"]] %in% JDOC_VARIANCE_KINDS) {
        expect_true(p[["scale"]] %in% c("variance", "precision"), info = info)
      } else {
        expect_null(p[["scale"]])
      }
      expect_true(is.numeric(p[["value"]]), info = info)
    }
  }
})

test_that("variables describe distribution, link and states", {
  for (name in names(jfx_registry())) {
    jfx_skip_bayes(name)
    spec <- jfx_spec(name)
    doc <- jfx_doc(name)
    info <- paste("fixture:", name)
    fit <- jfx_fit(name)
    expect_equal(vapply(doc[["variables"]], function(v) v[["name"]], character(1)),
                 colnames(fit$abnDag$dag), info = info)
    for (node in names(spec$dists)) {
      v <- jdoc_variable(doc, node)
      expected <- JDOC_DIST_MAP[[spec$dists[[node]]]]
      expect_equal(v[["type"]], expected$type, info = info)
      expect_equal(v[["distribution"]], expected$distribution, info = info)
      expect_equal(v[["link"]], expected$link, info = info)
      if (spec$dists[[node]] %in% c("binomial", "multinomial")) {
        labels <- vapply(v[["states"]], function(s) as.character(s[["label"]]), character(1))
        expect_equal(labels, jcore_levels(spec, node), info = info)
        expect_equal(jdoc_baseline_label(v), jcore_levels(spec, node)[1], info = info)
      } else {
        expect_null(v[["states"]])
      }
    }
  }
})

test_that("gaussian transforms are written only for centred fits", {
  doc <- jfx_doc("ex1_mle")
  for (v in doc[["variables"]]) expect_null(v[["transform"]])

  spec <- jfx_spec("ex1_mle_centred")
  doc <- jfx_doc("ex1_mle_centred")
  for (node in c("g1", "g2")) {
    tr <- jdoc_variable(doc, node)$transform
    expect_equal(tr[["center"]], mean(spec$data[[node]]))
    expect_equal(tr[["scale"]], stats::sd(spec$data[[node]]))
  }
  for (node in c("b1", "p1")) expect_null(jdoc_variable(doc, node)$transform)
})

test_that("arcs reproduce the DAG", {
  for (name in names(jfx_registry())) {
    jfx_skip_bayes(name)
    fit <- jfx_fit(name)
    doc <- jfx_doc(name)
    dag <- fit$abnDag$dag
    rebuilt <- dag
    rebuilt[] <- 0
    for (a in doc[["arcs"]]) {
      rebuilt[jdoc_variable_by_id(doc, a[["target"]])$name,
              jdoc_variable_by_id(doc, a[["source"]])$name] <- 1
    }
    expect_equal(rebuilt, dag, info = paste("fixture:", name))
  }
})

test_that("groups describe the random effect structure", {
  expect_equal(length(jfx_doc("ex1_mle")$groups), 0L)

  for (name in jfx_names(grouped = TRUE)) {
    jfx_skip_bayes(name)
    fit <- jfx_fit(name)
    doc <- jfx_doc(name)
    info <- paste("fixture:", name)
    expect_equal(length(doc[["groups"]]), 1L, info = info)
    group <- jdoc_group(doc, fit$group.var)
    grouped_names <- colnames(fit$abnDag$dag)[fit$grouped.vars]
    member_names <- vapply(group[["variables"]], function(id) {
      jdoc_variable_by_id(doc, id)$name
    }, character(1))
    expect_setequal(member_names, grouped_names)
  }
})

# --- parameters: ungrouped MLE ------------------------------------------------

test_that("ungrouped MLE: one parameter per coef column with its standard error", {
  for (name in jfx_names(method = "mle", grouped = FALSE)) {
    fit <- jfx_fit(name)
    doc <- jfx_doc(name)
    for (node in names(fit$coef)) {
      info <- paste(name, node)
      params <- jdoc_parameters(doc, node, c("intercept", "coefficient"))
      expect_equal(length(params), ncol(fit$coef[[node]]), info = info)
      expect_equal(sort(jdoc_values(params)), sort(unname(fit$coef[[node]][1, ])),
                   info = info)
      se <- vapply(params, function(p) as.numeric(p[["uncertainty"]][["standard_error"]]),
                   numeric(1))
      expect_equal(sort(se), sort(unname(fit$Stderror[[node]][1, ])), info = info)
    }
  }
})

test_that("ungrouped MLE: binary and continuous parents are mapped by name", {
  spec <- jfx_spec("ex1_mle")
  fit <- jfx_fit("ex1_mle")
  doc <- jfx_doc("ex1_mle")
  b1_success <- jcore_levels(spec, "b1")[2]
  b2_success <- jcore_levels(spec, "b2")[2]
  coef <- fit$coef$b3[1, ]
  expect_equal(jdoc_parameter(doc, "b3", "intercept")$value, coef[["b3|intercept"]])
  expect_equal(jdoc_parameter(doc, "b3", "coefficient", "b1", b1_success)$value,
               coef[["b1"]])
  expect_equal(jdoc_parameter(doc, "b3", "coefficient", "g1")$value, coef[["g1"]])
  expect_equal(jdoc_parameter(doc, "b3", "coefficient", "b2", b2_success)$value,
               coef[["b2"]])
})

test_that("ungrouped MLE: a multinomial parent yields one coefficient per level", {
  spec <- jfx_spec("fcv_mle")
  fit <- jfx_fit("fcv_mle")
  doc <- jfx_doc("fcv_mle")
  expect_equal(length(jdoc_parameters(doc, "Outdoor", "intercept")), 0L)
  for (level in jcore_levels(spec, "Sex")) {
    expect_equal(jdoc_parameter(doc, "Outdoor", "coefficient", "Sex", level)$value,
                 unname(fit$coef$Outdoor[1, paste0("Sex", level)]))
  }
})

test_that("ungrouped MLE: a multinomial child has parameters per non-baseline state", {
  fit <- jfx_fit("g2b2c_mle")
  doc <- jfx_doc("g2b2c_mle")
  coef <- fit$coef$C[1, ]
  for (state in c("b", "c")) {
    expect_equal(jdoc_parameter(doc, "C", "intercept", target_state = state)$value,
                 coef[[paste0("C|intercept.", state)]])
    expect_equal(jdoc_parameter(doc, "C", "coefficient", "B1", "1", state)$value,
                 coef[[paste0("B1", state)]])
    expect_equal(jdoc_parameter(doc, "C", "coefficient", "B2", "1", state)$value,
                 coef[[paste0("B2", state)]])
  }
  baseline_params <- Filter(function(p) {
    identical(jdoc_chr(p[["target_state"]]),
              jdoc_state_id(jdoc_variable(doc, "C"), "a"))
  }, doc[["parameters"]])
  expect_equal(length(baseline_params), 0L)
})

test_that("ungrouped MLE: gaussian nodes have a residual variance equal to mse", {
  fit <- jfx_fit("ex1_mle")
  doc <- jfx_doc("ex1_mle")
  for (node in c("g1", "g2")) {
    p <- jdoc_parameter(doc, node, "residual_variance")
    expect_equal(p[["scale"]], "variance")
    expect_equal(p[["value"]], unname(fit$mse[[node]]))
  }
  for (node in c("b1", "p1", "b2", "p2", "b3")) {
    expect_equal(length(jdoc_parameters(doc, node, JDOC_VARIANCE_KINDS)), 0L)
  }
})

# --- parameters: grouped MLE --------------------------------------------------

test_that("grouped MLE: gaussian node has fixed effects and variances", {
  fit <- jfx_fit("adg_mle_grouped")
  doc <- jfx_doc("adg_mle_grouped")
  group_id <- jdoc_chr(jdoc_group(doc, "farm")$`_id`)

  expect_equal(jdoc_parameter(doc, "adg", "intercept")$value, unname(fit$mu$adg))
  expect_equal(jdoc_parameter(doc, "adg", "coefficient", "age")$value,
               unname(fit$betas$adg[["age"]]))
  expect_equal(jdoc_parameter(doc, "adg", "coefficient", "wormCount")$value,
               unname(fit$betas$adg[["wormCount"]]))

  resid <- jdoc_parameter(doc, "adg", "residual_variance")
  expect_equal(resid[["scale"]], "variance")
  expect_equal(resid[["value"]], unname(fit$sigma$adg)^2)

  random <- jdoc_parameter(doc, "adg", "random_variance")
  expect_equal(random[["scale"]], "variance")
  expect_equal(jdoc_chr(random[["group"]]), group_id)
  expect_equal(random[["value"]], unname(fit$sigma_alpha$adg)^2)
})

test_that("grouped MLE: binary and count nodes have no residual variance", {
  fit <- jfx_fit("adg_mle_grouped")
  doc <- jfx_doc("adg_mle_grouped")
  for (node in c("AR", "eggs", "wormCount")) {
    expect_equal(length(jdoc_parameters(doc, node, "residual_variance")), 0L)
    expect_equal(jdoc_parameter(doc, node, "random_variance")$value,
                 unname(fit$sigma_alpha[[node]])^2)
    expect_equal(jdoc_parameter(doc, node, "intercept")$value, unname(fit$mu[[node]]))
  }
})

test_that("grouped MLE: fixed effects have no standard errors", {
  doc <- jfx_doc("adg_mle_grouped")
  for (p in doc[["parameters"]]) expect_null(p[["uncertainty"]][["standard_error"]])
})

test_that("grouped MLE: multinomial child exports its full covariance matrix", {
  spec <- jfx_spec("g2pbcgrp_mle_grouped")
  fit <- jfx_fit("g2pbcgrp_mle_grouped")
  doc <- jfx_doc("g2pbcgrp_mle_grouped")
  states <- jcore_levels(spec, "C")[-1]
  k <- length(states)

  intercepts <- jdoc_parameters(doc, "C", "intercept")
  expect_equal(length(intercepts), k)
  expect_equal(sort(jdoc_values(intercepts)), sort(unname(fit$mu$C)))

  coefs <- jdoc_parameters(doc, "C", "coefficient")
  expect_equal(length(coefs), length(fit$betas$C))
  expect_equal(sort(jdoc_values(coefs)), sort(as.numeric(fit$betas$C)))
  for (p in coefs) expect_false(is.null(p[["target_state"]]))

  variances <- jdoc_parameters(doc, "C", "random_variance")
  expect_equal(length(variances), k)
  expect_equal(sort(jdoc_values(variances)), sort(unname(diag(fit$sigma_alpha$C))))
  for (p in variances) expect_equal(p[["scale"]], "variance")

  covariances <- jdoc_parameters(doc, "C", "random_covariance")
  expect_equal(length(covariances), k * (k - 1) / 2)
  sa <- fit$sigma_alpha$C
  expect_equal(sort(jdoc_values(covariances)), sort(unname(sa[upper.tri(sa)])))
  for (p in covariances) {
    expect_length(p[["states"]], 2)
    expect_false(identical(p[["states"]][[1]], p[["states"]][[2]]))
  }
  expect_equal(length(jdoc_parameters(doc, "C", "residual_variance")), 0L)
})

# --- parameters: Bayes --------------------------------------------------------

jcore_bayes_parameter <- function(doc, spec, node, native) {
  term <- sub(paste0("^", node, "\\|"), "", native)
  if (term == "(Intercept)") return(jdoc_parameter(doc, node, "intercept"))
  if (term == "precision") return(jdoc_parameter(doc, node, "residual_variance"))
  if (term == "group.precision") return(jdoc_parameter(doc, node, "random_variance"))
  parent_state <- if (spec$dists[[term]] == "binomial") jcore_levels(spec, term)[2] else NULL
  jdoc_parameter(doc, node, "coefficient", term, parent_state)
}

test_that("bayes: parameter values are the posterior modes", {
  jfx_skip_if_no_bayes()
  for (name in jfx_names(method = "bayes")) {
    spec <- jfx_spec(name)
    fit <- jfx_fit(name)
    doc <- jfx_doc(name)
    for (node in names(fit$modes)) {
      info <- paste(name, node)
      expect_equal(length(jdoc_parameters(doc, node)), length(fit$modes[[node]]),
                   info = info)
      for (native in names(fit$modes[[node]])) {
        p <- jcore_bayes_parameter(doc, spec, node, native)
        expect_equal(p[["value"]], unname(fit$modes[[node]][[native]]), info = native)
        expect_null(p[["uncertainty"]][["standard_error"]])
        if (grepl("precision", native)) expect_equal(p[["scale"]], "precision")
      }
    }
  }
})

test_that("bayes: group precisions reference the group", {
  jfx_skip_if_no_bayes()
  doc <- jfx_doc("adg_bayes_grouped")
  p <- jdoc_parameter(doc, "adg", "random_variance")
  expect_equal(jdoc_chr(p[["group"]]), jdoc_chr(jdoc_group(doc, "farm")$`_id`))
})

test_that("bayes: posterior quantiles and marginals are attached to parameters", {
  jfx_skip_if_no_bayes()
  spec <- jfx_spec("ex1_bayes_fixed")
  fit <- jfx_fit("ex1_bayes_fixed")
  doc <- jfx_doc("ex1_bayes_fixed")
  marginals <- doc[["inference"]][["posterior"]][["marginals"]]
  for (node in names(fit$marginals)) {
    for (native in names(fit$marginals[[node]])) {
      p <- jcore_bayes_parameter(doc, spec, node, native)
      q <- fit$marginal.quantiles[[node]][[native]]
      pq <- p[["uncertainty"]][["posterior_quantiles"]]
      expect_equal(vapply(pq, function(x) x[["probability"]], numeric(1)),
                   unname(q[, "P(X<=x)"]), info = native)
      expect_equal(vapply(pq, function(x) x[["value"]], numeric(1)), unname(q[, "x"]),
                   info = native)

      m <- Filter(function(x) identical(jdoc_chr(x[["parameter"]]), jdoc_chr(p$`_id`)),
                  marginals)
      expect_length(m, 1)
      expect_equal(unlist(m[[1]][["x"]]), unname(fit$marginals[[node]][[native]][, "x"]))
      expect_equal(unlist(m[[1]][["density"]]),
                   unname(fit$marginals[[node]][[native]][, "f(x)"]))
    }
  }
})

test_that("bayes without compute.fixed has no posterior block", {
  jfx_skip_if_no_bayes()
  doc <- jfx_doc("ex1_bayes")
  expect_null(doc[["inference"]][["posterior"]])
  for (p in doc[["parameters"]]) expect_null(p[["uncertainty"]][["posterior_quantiles"]])
})

# --- inference ----------------------------------------------------------------

test_that("inference type follows the fit method", {
  expect_equal(jfx_doc("ex1_mle")$inference[["type"]], "maximum_likelihood")
  jfx_skip_if_no_bayes()
  expect_equal(jfx_doc("ex1_bayes")$inference[["type"]], "bayesian")
})

test_that("bayes priors are exported, MLE has none", {
  expect_null(jfx_doc("ex1_mle")$inference[["priors"]])
  jfx_skip_if_no_bayes()
  priors <- jfx_doc("ex1_bayes_priors")$inference[["priors"]]
  fixed <- Filter(function(x) identical(x[["applies_to"]], "fixed_effects"), priors)
  precisions <- Filter(function(x) identical(x[["applies_to"]], "precisions"), priors)
  expect_length(fixed, 1)
  expect_length(precisions, 1)
  expect_equal(fixed[[1]][["family"]], "normal")
  expect_equal(fixed[[1]][["mean"]], 0.5)
  expect_equal(fixed[[1]][["precision"]], 0.01)
  expect_equal(precisions[[1]][["family"]], "gamma")
  expect_equal(precisions[[1]][["shape"]], 2)
  expect_equal(precisions[[1]][["rate"]], 1e-3)
})

test_that("MLE diagnostics are exported once, per network and per node", {
  fit <- jfx_fit("ex1_mle")
  doc <- jfx_doc("ex1_mle")
  diag <- doc[["inference"]][["diagnostics"]]
  expect_equal(diag[["log_marginal_likelihood"]], unname(fit$mlik))
  expect_equal(diag[["aic"]], unname(fit$aic))
  expect_equal(diag[["bic"]], unname(fit$bic))
  expect_length(diag[["nodes"]], length(fit$mliknode))
  for (node in names(fit$mliknode)) {
    d <- jdoc_node_diagnostics(doc, node)
    expect_equal(d[["log_marginal_likelihood"]], unname(fit$mliknode[[node]]))
    expect_equal(d[["aic"]], unname(fit$aicnode[[node]]))
    expect_equal(d[["bic"]], unname(fit$bicnode[[node]]))
    expect_equal(d[["mdl"]], unname(fit$mdlnode[[node]]))
    expect_equal(d[["df"]], unname(fit$df[[node]]))
    expect_equal(d[["sse"]], unname(fit$sse[[node]]))
    # the gaussian mse is the residual_variance parameter and not repeated here
    if (identical(jfx_spec("ex1_mle")$dists[[node]], "gaussian")) {
      expect_null(d[["mse"]])
    } else {
      expect_equal(d[["mse"]], unname(fit$mse[[node]]))
    }
  }
})

test_that("bayes diagnostics contain the marginal likelihoods only", {
  jfx_skip_if_no_bayes()
  fit <- jfx_fit("ex1_bayes")
  doc <- jfx_doc("ex1_bayes")
  diag <- doc[["inference"]][["diagnostics"]]
  expect_equal(diag[["log_marginal_likelihood"]], unname(fit$mlik))
  expect_null(diag[["aic"]])
  for (node in names(fit$mliknode)) {
    d <- jdoc_node_diagnostics(doc, node)
    expect_equal(d[["log_marginal_likelihood"]], unname(fit$mliknode[[node]]))
    expect_null(d[["mse"]])
  }
})
