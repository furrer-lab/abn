json_fixture_empty_dag <- function(dists) {
  dag <- matrix(0, nrow = length(dists), ncol = length(dists))
  colnames(dag) <- rownames(dag) <- names(dists)
  dag
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
    expect_named(fit$Stderror, names(dists))
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
  expect_true("variables" %in% names(parsed))
  expect_true("parameters" %in% names(parsed))
  expect_true("arcs" %in% names(parsed))
  if (!is.null(parsed$method)) expect_equal(parsed$method, fit$method)

  variables <- json_fixture_collect_rows(parsed$variables)
  parameters <- json_fixture_collect_rows(parsed$parameters)
  arcs <- json_fixture_collect_rows(parsed$arcs)
  expect_equal(length(variables), length(dists))
  if (require_parameters) expect_gt(length(parameters), 0)
  expect_equal(length(arcs), json_fixture_expected_arc_count(fit))

  variable_ids <- vapply(variables, function(variable) {
    expect_true("variable_id" %in% names(variable))
    expect_true("attribute_name" %in% names(variable))
    expect_true("model_type" %in% names(variable))
    json_fixture_field(variable, "variable_id")
  }, character(1))
  variable_names <- vapply(variables, json_fixture_field, character(1), "attribute_name")
  variable_types <- vapply(variables, json_fixture_field, character(1), "model_type")
  expect_equal(anyDuplicated(variable_ids), 0L)
  expect_setequal(variable_names, names(dists))
  expect_equal(unname(variable_types[match(names(dists), variable_names)]),
               unname(unlist(dists)))

  for (arc in arcs) {
    expect_true("source_variable_id" %in% names(arc))
    expect_true("target_variable_id" %in% names(arc))
    expect_true(as.character(arc$source_variable_id) %in% variable_ids)
    expect_true(as.character(arc$target_variable_id) %in% variable_ids)
  }

  multinomial_variable_ids <- character(0)
  for (variable in variables) {
    if (identical(variable$model_type, "multinomial")) {
      multinomial_variable_ids <- c(multinomial_variable_ids,
                                    as.character(variable$variable_id))
      expect_true("states" %in% names(variable))
      states <- json_fixture_collect_rows(variable$states)
      expect_gt(length(states), 0)
      state_ids <- vapply(states, function(state) {
        expect_true("state_id" %in% names(state))
        expect_true("value_name" %in% names(state))
        expect_true("is_baseline" %in% names(state))
        json_fixture_field(state, "state_id")
      }, character(1))
      expect_equal(anyDuplicated(state_ids), 0L)
    }
  }
  if (multinomial) expect_gt(length(multinomial_variable_ids), 0)

  condition_types <- character(0)
  for (parameter in parameters) {
    expect_true("parameter_id" %in% names(parameter))
    expect_true("name" %in% names(parameter))
    expect_true("link_function_name" %in% names(parameter))
    expect_true("source" %in% names(parameter))
    expect_true("coefficients" %in% names(parameter))

    source <- parameter$source
    if (is.list(source) && length(source) == 1 && is.list(source[[1]])) {
      source <- source[[1]]
    }
    expect_true("variable_id" %in% names(source))
    expect_true(json_fixture_field(source, "variable_id") %in% variable_ids)

    coefficients <- json_fixture_collect_rows(parameter$coefficients)
    expect_gt(length(coefficients), 0)
    for (coefficient in coefficients) {
      expect_true("value" %in% names(coefficient))
      expect_true("stderr" %in% names(coefficient))
      expect_true("condition_type" %in% names(coefficient))
      expect_true("conditions" %in% names(coefficient))
      condition_types <- c(condition_types, as.character(coefficient$condition_type))

      conditions <- json_fixture_collect_rows(coefficient$conditions)
      for (condition in conditions) {
        expect_true("parent_variable_id" %in% names(condition))
        expect_true(json_fixture_field(condition, "parent_variable_id") %in% variable_ids)
        if (!is.null(condition$parent_state_id)) {
          expect_true(json_fixture_field(condition, "parent_variable_id") %in% variable_ids)
        }
      }
    }
  }

  if (grouped) {
    expect_true(any(condition_types %in% c("variance", "random_variance",
                                          "random_covariance")))
  }

  invisible(parsed)
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

  fitAbn(dag = dag, data.df = g2b2c_data, data.dists = dists,
         method = "mle", centre = FALSE)
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

  fitAbn(dag = dag, data.df = droplevels(g2pbcgrp), data.dists = dists,
         method = "mle", group.var = "group", centre = FALSE)
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

  fitAbn(dag = dag, data.df = data.df, data.dists = dists,
         method = "mle", centre = FALSE)
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

  fitAbn(dag = dag, data.df = data.df, data.dists = dists,
         method = "bayes", centre = FALSE)
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

  fitAbn(dag = dag, data.df = data.df, data.dists = dists,
         method = "bayes", centre = FALSE, compute.fixed = TRUE,
         control = fit.control(method = "bayes", n.grid = 100))
}

json_fixture_fit_ex3_mle_grouped <- function() {
  dists <- list(b1 = "binomial", b2 = "binomial")

  fitAbn(dag = ~ b1 | b2, data.df = ex3.dag.data[, c(1, 2, 14)],
         data.dists = dists, method = "mle", group.var = "group",
         cor.vars = c("b1", "b2"), centre = FALSE)
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

  fitAbn(dag = dag, data.df = data.df, data.dists = dists,
         method = "mle", group.var = "farm", centre = FALSE)
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

  fitAbn(dag = dag, data.df = data.df, data.dists = dists,
         method = "mle", centre = FALSE)
}
