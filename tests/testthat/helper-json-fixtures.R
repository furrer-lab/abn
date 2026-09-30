json_fixture_empty_dag <- function(dists) {
  dag <- matrix(0, nrow = length(dists), ncol = length(dists))
  colnames(dag) <- rownames(dag) <- names(dists)
  dag
}

json_import_generic_document <- function() {
  list(
    metadata = list(schema_version = "bayesian-network-v1",
                    issuer = "independent-bn-tool::export"),
    variables = list(
      list(`_id` = 1, name = "height", type = "continuous"),
      list(`_id` = 2, name = "status", type = "binary",
           states = list(list(`_id` = 1, label = "no"),
                         list(`_id` = 2, label = "yes")))
    ),
    arcs = list(list(source = 1, target = 2)),
    parameters = list(
      list(`_id` = 1, target = 1, kind = "intercept", link = "identity",
           value = 1.2, uncertainty = list(standard_error = 0.1)),
      list(`_id` = 2, target = 2, kind = "intercept", link = "logit",
           value = -0.42, uncertainty = list(standard_error = 0.2)),
      list(`_id` = 3, target = 2, parents = list(1), kind = "coefficient",
           link = "logit", value = 0.87,
           uncertainty = list(standard_error = 0.3))
    ),
    inference = list(type = "maximum_likelihood", estimates = list(),
                     uncertainty = list(), diagnostics = list())
  )
}

json_import_string <- function(document) {
  jsonlite::toJSON(document, auto_unbox = TRUE, null = "null", digits = NA)
}

json_fixture_assert_abnfit <- function(fit, method, dists, grouped = FALSE,
                                       multinomial = FALSE) {
  expect_s3_class(fit, "abnFit")
  expect_s3_class(fit$abnDag, "abnDag")
  expect_equal(fit$method, method)
  expect_named(fit$abnDag$data.dists, names(dists))
  expect_equal(unname(unlist(fit$abnDag$data.dists)), unname(unlist(dists)))
  expect_equal(dim(fit$abnDag$dag), c(length(dists), length(dists)))
  expect_setequal(colnames(fit$abnDag$dag), names(dists))
  expect_setequal(rownames(fit$abnDag$dag), names(dists))

  if (!grouped) {
    expect_equal(length(fit$coef), length(dists))
    expect_named(fit$coef, names(dists))
    if (identical(method, "mle")) {
      expect_named(fit$Stderror, names(dists))
    }
  } else {
    if (!is.null(fit$coef)) expect_named(fit$coef, names(dists))
    if (!is.null(fit$Stderror)) expect_named(fit$Stderror, names(dists))
  }

  if (multinomial) {
    expect_true("multinomial" %in% unlist(dists))
    expect_true("multinomial" %in% unlist(fit$abnDag$data.dists))
    expect_true(!is.null(fit$multinomial.states) ||
                  any(unlist(fit$abnDag$data.dists) == "multinomial"))
  }

  if (grouped) {
    expect_true(!is.null(fit$group.var) || !is.null(fit$mu))
    expect_true(!is.null(fit$mu) || !is.null(fit$betas) ||
                  !is.null(fit$sigma) || !is.null(fit$sigma_alpha))
  }

  invisible(fit)
}

json_fixture_parse_export <- function(fit) {
  json <- export_abnFit(fit)
  expect_type(json, "character")
  expect_length(json, 1)
  expect_true(jsonlite::validate(json))
  jsonlite::fromJSON(json, simplifyVector = FALSE)
}

json_fixture_collect_rows <- function(x) {
  if (is.null(x)) return(list())
  if (is.data.frame(x)) return(split(x, seq_len(nrow(x))))
  if (is.list(x) && length(x) > 0 && !is.null(names(x[[1]]))) return(x)
  list()
}

json_fixture_field <- function(x, name) {
  value <- x[[name]]
  if (is.null(value) || length(value) == 0) return(NA_character_)
  as.character(value[[1]])
}

json_fixture_expected_arc_count <- function(fit) {
  sum(fit$abnDag$dag == 1)
}

json_fixture_assert_export <- function(fit, dists, grouped = FALSE,
                                       multinomial = FALSE,
                                       require_parameters = TRUE) {
  parsed <- json_fixture_parse_export(fit)
  expect_true("metadata" %in% names(parsed))
  expect_equal(parsed$metadata$schema_version, "bayesian-network-v1")
  expect_equal(parsed$metadata$issuer, "abn::export_abnFit")
  expect_true("configs" %in% names(parsed$metadata))
  expect_true("extensions" %in% names(parsed$metadata))
  expect_true("inference" %in% names(parsed))
  expected_inference <- if (identical(fit$method, "bayes")) {
    "bayesian"
  } else {
    "maximum_likelihood"
  }
  expect_equal(parsed$inference$type, expected_inference)
  expect_true("variables" %in% names(parsed))
  expect_true("parameters" %in% names(parsed))
  expect_true("arcs" %in% names(parsed))

  variables <- json_fixture_collect_rows(parsed$variables)
  parameters <- json_fixture_collect_rows(parsed$parameters)
  arcs <- json_fixture_collect_rows(parsed$arcs)
  expect_equal(length(variables), length(dists))
  if (require_parameters) expect_gt(length(parameters), 0)
  expect_equal(length(arcs), json_fixture_expected_arc_count(fit))

  variable_ids <- vapply(variables, function(variable) {
    expect_true("_id" %in% names(variable))
    expect_true("name" %in% names(variable))
    expect_true("type" %in% names(variable))
    json_fixture_field(variable, "_id")
  }, character(1))
  variable_names <- vapply(variables, json_fixture_field, character(1), "name")
  variable_types <- vapply(variables, json_fixture_field, character(1), "type")
  expect_equal(anyDuplicated(variable_ids), 0L)
  expect_setequal(variable_names, names(dists))
  expected_types <- unname(unlist(stats::setNames(lapply(dists, function(x) {
    switch(x, gaussian = "continuous", binomial = "binary",
           poisson = "count", multinomial = "categorical")
  }), names(dists))))
  expect_equal(unname(variable_types[match(names(dists), variable_names)]),
               expected_types)

  for (arc in arcs) {
    expect_true("source" %in% names(arc))
    expect_true("target" %in% names(arc))
    expect_true(as.character(arc$source) %in% variable_ids)
    expect_true(as.character(arc$target) %in% variable_ids)
  }

  multinomial_variable_ids <- character(0)
  for (variable in variables) {
    if (identical(variable$type, "categorical")) {
      multinomial_variable_ids <- c(multinomial_variable_ids,
                                     as.character(variable$`_id`))
      expect_true("states" %in% names(variable))
      states <- json_fixture_collect_rows(variable$states)
      expect_gt(length(states), 0)
      state_ids <- vapply(states, function(state) {
        expect_true("_id" %in% names(state))
        expect_true("label" %in% names(state))
        json_fixture_field(state, "_id")
      }, character(1))
      expect_equal(anyDuplicated(state_ids), 0L)
    }
  }
  if (multinomial) expect_gt(length(multinomial_variable_ids), 0)

  condition_types <- character(0)
  for (parameter in parameters) {
    expect_true("_id" %in% names(parameter))
    expect_true("target" %in% names(parameter))
    expect_true("kind" %in% names(parameter))
    expect_true("value" %in% names(parameter))
    expect_true(parameter$target %in% variable_ids)
    condition_types <- c(condition_types, as.character(parameter$kind))
    if (!is.null(parameter$parents)) {
      expect_true(all(unlist(parameter$parents) %in% variable_ids))
    }
  }

  if (grouped) {
    expect_true("abn" %in% names(parsed$metadata$extensions))
    expect_true(any(condition_types %in% c("variance", "random_variance",
                                          "random_covariance")))
  }

  invisible(parsed)
}

json_fixture_assert_import <- function(imported, method, dists, grouped = FALSE,
                                       multinomial = FALSE,
                                       original_model = FALSE, parsed = NULL) {
  json_fixture_assert_abnfit(imported, method = method, dists = dists,
                             grouped = grouped, multinomial = multinomial)

  expect_equal(imported$method, method)
  expect_equal(length(imported$coef), length(dists))
  expect_named(imported$coef, names(dists))
  expect_named(imported$Stderror, names(dists))

  expect_true(is.matrix(imported$abnDag$dag))
  expect_equal(dim(imported$abnDag$dag), c(length(dists), length(dists)))
  expect_setequal(rownames(imported$abnDag$dag), names(dists))
  expect_setequal(colnames(imported$abnDag$dag), names(dists))

  for (node in names(imported$coef)) {
    expect_false(any(is.na(colnames(imported$coef[[node]]))))
    expect_equal(anyDuplicated(colnames(imported$coef[[node]])), 0L)
  }

  if (grouped) {
    expect_true(!is.null(imported$mu))
    expect_true(!is.null(imported$betas))
    expect_true(!is.null(imported$sigma))
    expect_true(!is.null(imported$sigma_alpha))
  }

  if (original_model) {
    expect_true("metadata" %in% names(parsed))
    expect_equal(parsed$metadata$issuer, "abn::export_abnFit")
    expect_true("abn" %in% names(parsed$metadata$extensions))
  }

  if (!is.null(parsed)) {
    variables <- json_fixture_collect_rows(parsed$variables)
    parameters <- json_fixture_collect_rows(parsed$parameters)
    arcs <- json_fixture_collect_rows(parsed$arcs)

    id_to_name <- stats::setNames(
      vapply(variables, json_fixture_field, character(1), "name"),
      vapply(variables, json_fixture_field, character(1), "_id")
    )
    expect_equal(length(variables), length(dists))
    expect_equal(length(arcs), sum(imported$abnDag$dag == 1))

    expected_dag <- matrix(0, nrow = length(dists), ncol = length(dists),
                           dimnames = list(names(dists), names(dists)))
    for (arc in arcs) {
      source <- id_to_name[[json_fixture_field(arc, "source")]]
      target <- id_to_name[[json_fixture_field(arc, "target")]]
      expect_true(source %in% names(dists))
      expect_true(target %in% names(dists))
      expected_dag[source, target] <- 1
    }
    expect_equal(imported$abnDag$dag[names(dists), names(dists)], expected_dag)

    exported_condition_types <- character(0)
    for (parameter in parameters) {
      exported_condition_types <- c(exported_condition_types,
                                    json_fixture_field(parameter, "kind"))
    }
    expected_coef_count <- sum(exported_condition_types %in% c("intercept", "linear_term"))
    imported_coef_count <- sum(vapply(imported$coef, ncol, integer(1)))
    expect_equal(imported_coef_count, expected_coef_count)

    if (grouped) {
      expect_true(any(exported_condition_types %in% c("variance", "random_variance",
                                                      "random_covariance")))
    }

    if (original_model) {
      expect_true("abn" %in% names(parsed$metadata$extensions))
      expect_true("inference" %in% names(parsed$metadata$extensions$abn))
    }
  }

  invisible(imported)
}

json_fixture_compare_numeric <- function(x, y, tolerance = .Machine$double.eps^0.5) {
  if (is.null(x) && is.null(y)) return(invisible(TRUE))
  expect_false(xor(is.null(x), is.null(y)))
  expect_equal(unname(as.numeric(x)), unname(as.numeric(y)), tolerance = tolerance)
}

json_fixture_absent_grouped_value <- function(x) {
  length(x) == 0 || (length(x) == 1 && is.na(x))
}

json_fixture_as_named_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (is.null(x) || length(x) == 0) return(matrix(numeric(0), nrow = 0, ncol = 0))
  values <- as.numeric(x)
  nms <- names(x)
  out <- matrix(values, nrow = 1)
  if (!is.null(nms)) colnames(out) <- nms
  out
}

json_fixture_compare_matrix <- function(original, imported,
                                        tolerance = .Machine$double.eps^0.5,
                                        allow_imported_all_na = FALSE) {
  original <- json_fixture_as_named_matrix(original)
  imported <- json_fixture_as_named_matrix(imported)

  if (allow_imported_all_na && length(original) == 0 &&
      length(imported) > 0 && all(is.na(imported))) {
    return(invisible(TRUE))
  }

  expect_equal(dim(imported), dim(original))

  original_names <- colnames(original)
  imported_names <- colnames(imported)
  if (is.null(original_names) || is.null(imported_names)) {
    expect_equal(imported_names, original_names)
    imported_reordered <- imported
  } else {
    expect_setequal(imported_names, original_names)
    imported_reordered <- imported[, original_names, drop = FALSE]
  }

  if (allow_imported_all_na && all(is.na(imported_reordered))) {
    return(invisible(TRUE))
  }

  expect_equal(unname(as.numeric(imported_reordered)),
               unname(as.numeric(original)), tolerance = tolerance)
  invisible(TRUE)
}

json_fixture_compare_named_list <- function(original, imported,
                                            tolerance = .Machine$double.eps^0.5) {
  expect_equal(length(imported), length(original))
  expect_setequal(names(imported), names(original))
  for (node in names(original)) {
    o <- original[[node]]
    i <- imported[[node]]
    if (json_fixture_absent_grouped_value(o) && json_fixture_absent_grouped_value(i)) {
      next
    }
    if (is.matrix(o) || is.matrix(i)) {
      json_fixture_compare_matrix(o, i, tolerance = tolerance)
    } else if (length(o) == 1 && is.na(o) && length(i) == 1 && is.na(i)) {
      next
    } else {
      if (!is.null(names(o)) && !is.null(names(i))) {
        expect_setequal(names(i), names(o))
        i <- i[names(o)]
      } else if (is.null(names(o)) && is.null(names(i))) {
        # Compare unnamed scalars/vectors positionally.
      } else if (length(o) <= 1 && length(i) <= 1) {
        # Allow scalar grouped fields to differ only in incidental naming.
        names(o) <- NULL
        names(i) <- NULL
      } else {
        stop("Grouped field lost names for a multi-value component", call. = FALSE)
      }
      json_fixture_compare_numeric(o, i, tolerance = tolerance)
    }
  }
  invisible(TRUE)
}

json_fixture_assert_roundtrip <- function(original, imported, grouped = FALSE,
                                          multinomial = FALSE, bayes = FALSE,
                                          fixed = FALSE,
                                          tolerance = .Machine$double.eps^0.5) {
  expect_s3_class(imported, "abnFit")
  expect_equal(imported$method, original$method)
  expect_s3_class(imported$abnDag, "abnDag")

  expect_setequal(colnames(imported$abnDag$dag), colnames(original$abnDag$dag))
  expect_setequal(rownames(imported$abnDag$dag), rownames(original$abnDag$dag))
  node_order <- colnames(original$abnDag$dag)
  expect_equal(imported$abnDag$dag[node_order, node_order],
               original$abnDag$dag[node_order, node_order])
  expect_equal(unlist(imported$abnDag$data.dists[node_order]),
               unlist(original$abnDag$data.dists[node_order]))

  if (!is.null(names(original$coef)) && !is.null(names(imported$coef))) {
    expect_setequal(names(imported$coef), names(original$coef))
    for (node in names(original$coef)) {
      json_fixture_compare_matrix(original$coef[[node]], imported$coef[[node]],
                                  tolerance = tolerance)
      json_fixture_compare_matrix(original$Stderror[[node]], imported$Stderror[[node]],
                                  tolerance = tolerance,
                                  allow_imported_all_na = bayes)
    }
  } else if (!grouped) {
    expect_setequal(names(imported$coef), names(original$coef))
  }

  if (multinomial) {
    expect_true(!is.null(imported$multinomial.states))
    for (node in names(original$coef)) {
      expect_equal(anyDuplicated(colnames(imported$coef[[node]])), 0L)
    }
  }

  if (grouped) {
    for (field in c("mu", "betas", "sigma", "sigma_alpha")) {
      expect_true(field %in% names(imported))
      json_fixture_compare_named_list(original[[field]], imported[[field]],
                                      tolerance = tolerance)
    }
  }

  if (bayes) {
    expect_equal(imported$method, "bayes")
    expect_equal(imported$mlik, original$mlik, tolerance = tolerance)
    expect_equal(unlist(imported$mliknode), unname(original$mliknode),
                 tolerance = tolerance)
    expect_equal(unlist(imported$used.INLA), unname(original$used.INLA))
    expect_equal(unlist(imported$error.code), unname(original$error.code))
    expect_equal(unlist(imported$error.code.desc), unname(original$error.code.desc))
    expect_equal(unlist(imported$hessian.accuracy),
                 unname(original$hessian.accuracy), tolerance = tolerance)
  }

  if (fixed) {
    expect_true("marginals" %in% names(imported))
    expect_true("marginal.quantiles" %in% names(imported))
    expect_equal(length(imported$marginals), length(original$marginals))
    expect_equal(length(imported$marginal.quantiles),
                 length(original$marginal.quantiles))
  }

  json_fixture_assert_native_fields(original, imported, tolerance = tolerance)

  invisible(TRUE)
}

json_fixture_assert_native_fields <- function(original, imported,
                                              tolerance = .Machine$double.eps^0.5) {
  expect_setequal(
    setdiff(names(imported), c("call", "abnDag")),
    setdiff(names(original), c("call", "abnDag"))
  )
  for (field in setdiff(names(original), c("call", "abnDag"))) {
    expect_equal(imported[[field]], original[[field]], tolerance = tolerance,
                 info = paste("native field:", field))
  }
  invisible(TRUE)
}

json_fixture_generic_document <- function(json) {
  document <- if (is.character(json)) {
    jsonlite::fromJSON(json, simplifyVector = FALSE)
  } else {
    json
  }
  document$metadata <- NULL
  document <- json_fixture_resolve_internal_ids(document)
  document
}

json_fixture_resolve_internal_ids <- function(document) {
  variables <- json_fixture_collect_rows(document$variables)
  if (length(variables) == 0) return(document)
  ids <- vapply(variables, function(x) as.character(x$`_id`), character(1))
  names_by_id <- stats::setNames(vapply(variables, function(x) x$name, character(1)), ids)
  document$variables <- lapply(variables, function(x) {
    x$`_id` <- NULL
    x
  })
  document$arcs <- lapply(json_fixture_collect_rows(document$arcs), function(x) {
    x$source <- names_by_id[[as.character(x$source)]]
    x$target <- names_by_id[[as.character(x$target)]]
    x
  })
  document$parameters <- lapply(json_fixture_collect_rows(document$parameters), function(x) {
    x$`_id` <- NULL
    x$target <- names_by_id[[as.character(x$target)]]
    if (!is.null(x$parents)) {
      x$parents <- lapply(x$parents, function(parent) names_by_id[[as.character(parent)]])
    }
    x
  })
  document
}

json_fixture_assert_generic_roundtrip <- function(first, second,
                                                  tolerance = .Machine$double.eps^0.5) {
  expect_equal(
    json_fixture_canonical_json_value(json_fixture_generic_document(first)),
    json_fixture_canonical_json_value(json_fixture_generic_document(second)),
    tolerance = tolerance
  )
  invisible(TRUE)
}

json_fixture_canonical_json_value <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    out <- lapply(x, json_fixture_canonical_json_value)
    if (!is.null(names(out))) {
      out <- out[sort(names(out))]
    }
    return(out)
  }
  if (is.numeric(x)) return(as.numeric(x))
  if (is.logical(x)) return(as.logical(x))
  if (is.character(x)) return(as.character(x))
  x
}

json_fixture_sort_json_rows <- function(rows, key_fun) {
  rows <- json_fixture_collect_rows(rows)
  if (length(rows) == 0) return(list())
  rows <- lapply(rows, json_fixture_canonical_json_value)
  rows[order(vapply(rows, key_fun, character(1)))]
}

json_fixture_canonical_abn_json <- function(parsed) {
  parsed <- json_fixture_canonical_json_value(parsed)

  if (!is.null(parsed$variables)) {
    parsed$variables <- json_fixture_sort_json_rows(parsed$variables, function(variable) {
      json_fixture_field(variable, "_id")
    })
    parsed$variables <- lapply(parsed$variables, function(variable) {
      if (!is.null(variable$states)) {
        variable$states <- json_fixture_sort_json_rows(variable$states, function(state) {
          json_fixture_field(state, "_id")
        })
      }
      variable
    })
  }

  if (!is.null(parsed$parameters)) {
    parsed$parameters <- json_fixture_sort_json_rows(parsed$parameters, function(parameter) {
      json_fixture_field(parameter, "_id")
    })
  }

  if (!is.null(parsed$arcs)) {
    parsed$arcs <- json_fixture_sort_json_rows(parsed$arcs, function(arc) {
      paste(json_fixture_field(arc, "source"),
            json_fixture_field(arc, "target"), sep = "->")
    })
  }

  parsed
}

json_fixture_canonical_export <- function(fit) {
  json_fixture_canonical_abn_json(
    jsonlite::fromJSON(export_abnFit(fit), simplifyVector = FALSE)
  )
}

json_fixture_json_roundtrip_value <- function(x) {
  jsonlite::fromJSON(
    jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA),
    simplifyVector = FALSE
  )
}

json_fixture_fit_g2b2c_mle <- function() {
  dists <- list(
    G1 = "gaussian",
    B1 = "binomial",
    B2 = "binomial",
    C = "multinomial",
    G2 = "gaussian"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["B1", "G1"] <- 1
  dag["B2", "B1"] <- 1
  dag["C", c("B1", "B2")] <- 1
  dag["G2", c("G1", "C")] <- 1

  suppressWarnings(fitAbn(dag = dag, data.df = g2b2c_data, data.dists = dists,
                           method = "mle", centre = FALSE))
}

json_fixture_fit_g2pbcgrp_mle_grouped <- function() {
  dists <- list(
    G1 = "gaussian",
    P = "poisson",
    B = "binomial",
    C = "multinomial",
    G2 = "gaussian"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["P", "G1"] <- 1
  dag["B", "G1"] <- 1
  dag["C", c("P", "B")] <- 1
  dag["G2", c("G1", "C")] <- 1

  suppressWarnings(fitAbn(dag = dag, data.df = droplevels(g2pbcgrp),
                           data.dists = dists, method = "mle",
                           group.var = "group", centre = FALSE))
}

json_fixture_fit_ex1_mle <- function() {
  data.df <- droplevels(ex1.dag.data[seq_len(1000),
                                    c("b1", "p1", "g1", "b2", "p2", "b3", "g2")])
  dists <- list(
    b1 = "binomial",
    p1 = "poisson",
    g1 = "gaussian",
    b2 = "binomial",
    p2 = "poisson",
    b3 = "binomial",
    g2 = "gaussian"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["p2", c("b1", "p1")] <- 1
  dag["b3", c("b1", "g1", "b2")] <- 1
  dag["g2", c("p1", "g1", "b2")] <- 1

  suppressWarnings(fitAbn(dag = dag, data.df = data.df, data.dists = dists,
                           method = "mle", centre = FALSE))
}

json_fixture_fit_ex1_bayes <- function() {
  data.df <- droplevels(ex1.dag.data[seq_len(1000), c("b1", "p1", "g1", "p2")])
  dists <- list(
    b1 = "binomial",
    p1 = "poisson",
    g1 = "gaussian",
    p2 = "poisson"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["p2", c("b1", "p1")] <- 1

  suppressWarnings(fitAbn(dag = dag, data.df = data.df, data.dists = dists,
                           method = "bayes", centre = FALSE))
}

json_fixture_fit_ex1_bayes_fixed <- function() {
  data.df <- droplevels(ex1.dag.data[seq_len(1000), c("b1", "p1", "g1", "p2")])
  dists <- list(
    b1 = "binomial",
    p1 = "poisson",
    g1 = "gaussian",
    p2 = "poisson"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["p2", c("b1", "p1")] <- 1

  suppressWarnings(fitAbn(dag = dag, data.df = data.df, data.dists = dists,
                           method = "bayes", centre = FALSE, compute.fixed = TRUE,
                           control = fit.control(method = "bayes", n.grid = 100)))
}

json_fixture_fit_ex3_mle_grouped <- function() {
  dists <- list(b1 = "binomial", b2 = "binomial")

  suppressWarnings(fitAbn(dag = ~ b1 | b2, data.df = ex3.dag.data[, c(1, 2, 14)],
                           data.dists = dists, method = "mle", group.var = "group",
                           cor.vars = c("b1", "b2"), centre = FALSE))
}

json_fixture_fit_adg_mle_grouped <- function() {
  data.df <- droplevels(adg[, c("AR", "eggs", "wormCount", "age", "adg", "farm")])
  data.df$AR <- factor(data.df$AR)
  data.df$eggs <- factor(data.df$eggs)
  data.df$farm <- factor(data.df$farm)

  dists <- list(
    AR = "binomial",
    eggs = "binomial",
    wormCount = "poisson",
    age = "gaussian",
    adg = "gaussian"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["eggs", "AR"] <- 1
  dag["wormCount", "eggs"] <- 1
  dag["adg", c("age", "wormCount")] <- 1

  suppressWarnings(fitAbn(dag = dag, data.df = data.df, data.dists = dists,
                           method = "mle", group.var = "farm", centre = FALSE))
}

json_fixture_fit_fcv_mle <- function() {
  data.df <- droplevels(FCV[, c("Outdoor", "Sex", "GroupSize", "Age")])
  dists <- list(
    Outdoor = "binomial",
    Sex = "multinomial",
    GroupSize = "poisson",
    Age = "gaussian"
  )
  dag <- json_fixture_empty_dag(dists)
  dag["Outdoor", "Sex"] <- 1
  dag["GroupSize", "Outdoor"] <- 1
  dag["Age", c("Sex", "GroupSize")] <- 1

  suppressWarnings(fitAbn(dag = dag, data.df = data.df, data.dists = dists,
                           method = "mle", centre = FALSE))
}
