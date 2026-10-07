# Test helpers for the `bayesian-network` JSON format.
#
# jfx_*  : fixture specifications and cached fits
# jdoc_* : accessors for parsed JSON documents
#
# See json_plan.md for the format definition these tests pin down.

JDOC_SCHEMA_VERSION <- "bayesian-network"
JDOC_ISSUER <- "abn::export_abnFit"
JDOC_BLOCKS <- c("metadata", "variables", "groups", "arcs", "parameters", "inference")
JDOC_KINDS <- c("intercept", "coefficient", "residual_variance", "random_variance",
                "random_covariance")
JDOC_VARIANCE_KINDS <- c("residual_variance", "random_variance", "random_covariance")
JDOC_EXTENSION_KEYS <- c("parameter_names", "nodes")

# Map from abn distribution to the generic variable description.
JDOC_DIST_MAP <- list(
  gaussian = list(type = "continuous", distribution = "gaussian", link = "identity"),
  binomial = list(type = "binary", distribution = "bernoulli", link = "logit"),
  poisson = list(type = "count", distribution = "poisson", link = "log"),
  multinomial = list(type = "categorical", distribution = "categorical",
                     link = "baseline_logit")
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

jfx_empty_dag <- function(dists) {
  dag <- matrix(0, nrow = length(dists), ncol = length(dists))
  colnames(dag) <- rownames(dag) <- names(dists)
  dag
}

jfx_spec_g2b2c <- function() {
  dists <- list(G1 = "gaussian", B1 = "binomial", B2 = "binomial",
                C = "multinomial", G2 = "gaussian")
  dag <- jfx_empty_dag(dists)
  dag["B1", "G1"] <- 1
  dag["B2", "B1"] <- 1
  dag["C", c("B1", "B2")] <- 1
  dag["G2", c("G1", "C")] <- 1
  list(data = g2b2c_data[, names(dists)], dists = dists, dag = dag)
}

jfx_spec_fcv <- function() {
  dists <- list(Outdoor = "binomial", Sex = "multinomial", GroupSize = "poisson",
                Age = "gaussian")
  dag <- jfx_empty_dag(dists)
  dag["Outdoor", "Sex"] <- 1
  dag["GroupSize", "Outdoor"] <- 1
  dag["Age", c("Sex", "GroupSize")] <- 1
  list(data = droplevels(FCV[, names(dists)]), dists = dists, dag = dag)
}

jfx_spec_ex1 <- function() {
  dists <- list(b1 = "binomial", p1 = "poisson", g1 = "gaussian", b2 = "binomial",
                p2 = "poisson", b3 = "binomial", g2 = "gaussian")
  dag <- jfx_empty_dag(dists)
  dag["p2", c("b1", "p1")] <- 1
  dag["b3", c("b1", "g1", "b2")] <- 1
  dag["g2", c("p1", "g1", "b2")] <- 1
  list(data = droplevels(ex1.dag.data[seq_len(1000), names(dists)]),
       dists = dists, dag = dag)
}

jfx_spec_ex1_small <- function() {
  dists <- list(b1 = "binomial", p1 = "poisson", g1 = "gaussian", p2 = "poisson")
  dag <- jfx_empty_dag(dists)
  dag["p2", c("b1", "p1")] <- 1
  list(data = droplevels(ex1.dag.data[seq_len(1000), names(dists)]),
       dists = dists, dag = dag)
}

jfx_spec_g2pbcgrp <- function() {
  dists <- list(G1 = "gaussian", P = "poisson", B = "binomial", C = "multinomial",
                G2 = "gaussian")
  dag <- jfx_empty_dag(dists)
  dag["P", "G1"] <- 1
  dag["B", "G1"] <- 1
  dag["C", c("P", "B")] <- 1
  dag["G2", c("G1", "C")] <- 1
  list(data = droplevels(g2pbcgrp[, c(names(dists), "group")]), dists = dists,
       dag = dag, group.var = "group")
}

jfx_spec_ex3 <- function() {
  dists <- list(b1 = "binomial", b2 = "binomial")
  dag <- jfx_empty_dag(dists)
  dag["b1", "b2"] <- 1
  list(data = ex3.dag.data[, c("b1", "b2", "group")], dists = dists, dag = dag,
       group.var = "group", extra = list(cor.vars = c("b1", "b2")))
}

jfx_spec_adg <- function() {
  data.df <- droplevels(adg[, c("AR", "eggs", "wormCount", "age", "adg", "farm")])
  data.df$AR <- factor(data.df$AR)
  data.df$eggs <- factor(data.df$eggs)
  data.df$farm <- factor(data.df$farm)
  dists <- list(AR = "binomial", eggs = "binomial", wormCount = "poisson",
                age = "gaussian", adg = "gaussian")
  dag <- jfx_empty_dag(dists)
  dag["eggs", "AR"] <- 1
  dag["wormCount", "eggs"] <- 1
  dag["adg", c("age", "wormCount")] <- 1
  list(data = data.df, dists = dists, dag = dag, group.var = "farm")
}

jfx_spec_adg_gaussian <- function() {
  data.df <- droplevels(adg[, c("age", "adg", "farm")])
  data.df$farm <- factor(data.df$farm)
  dists <- list(age = "gaussian", adg = "gaussian")
  dag <- jfx_empty_dag(dists)
  dag["adg", "age"] <- 1
  list(data = data.df, dists = dists, dag = dag, group.var = "farm")
}

# Registry of all fixtures. Flags describe what each fixture exercises.
jfx_registry <- function() {
  mk <- function(base, method, centre = FALSE, grouped = FALSE, multinomial = FALSE,
                 fixed = FALSE, extra = list()) {
    spec <- base()
    spec$method <- method
    spec$centre <- centre
    spec$grouped <- grouped
    spec$multinomial <- multinomial
    spec$fixed <- fixed
    spec$extra <- c(spec$extra, extra)
    spec
  }
  list(
    g2b2c_mle = function() mk(jfx_spec_g2b2c, "mle", multinomial = TRUE),
    fcv_mle = function() mk(jfx_spec_fcv, "mle", multinomial = TRUE),
    ex1_mle = function() mk(jfx_spec_ex1, "mle"),
    ex1_mle_centred = function() mk(jfx_spec_ex1, "mle", centre = TRUE),
    g2pbcgrp_mle_grouped = function() mk(jfx_spec_g2pbcgrp, "mle", grouped = TRUE,
                                         multinomial = TRUE),
    ex3_mle_grouped = function() mk(jfx_spec_ex3, "mle", grouped = TRUE),
    adg_mle_grouped = function() mk(jfx_spec_adg, "mle", grouped = TRUE),
    ex1_bayes = function() mk(jfx_spec_ex1_small, "bayes"),
    ex1_bayes_centred = function() mk(jfx_spec_ex1_small, "bayes", centre = TRUE),
    ex1_bayes_fixed = function() mk(
      jfx_spec_ex1_small, "bayes", fixed = TRUE,
      extra = list(compute.fixed = TRUE,
                   control = fit.control(method = "bayes", n.grid = 100))
    ),
    ex1_bayes_priors = function() mk(
      jfx_spec_ex1_small, "bayes",
      extra = list(control = fit.control(method = "bayes", mean = 0.5, prec = 0.01,
                                         loggam.shape = 2, loggam.inv.scale = 1e-3))
    ),
    ex3_bayes_grouped = function() mk(jfx_spec_ex3, "bayes", grouped = TRUE),
    adg_bayes_grouped = function() mk(jfx_spec_adg_gaussian, "bayes", grouped = TRUE)
  )
}

jfx_names <- function(method = NULL, grouped = NULL, multinomial = NULL, fixed = NULL) {
  reg <- jfx_registry()
  keep <- vapply(names(reg), function(name) {
    spec <- jfx_spec(name)
    (is.null(method) || identical(spec$method, method)) &&
      (is.null(grouped) || identical(spec$grouped, grouped)) &&
      (is.null(multinomial) || identical(spec$multinomial, multinomial)) &&
      (is.null(fixed) || identical(spec$fixed, fixed))
  }, logical(1))
  names(reg)[keep]
}

jfx_cache <- new.env(parent = emptyenv())

jfx_spec <- function(name) {
  key <- paste0("spec:", name)
  if (is.null(jfx_cache[[key]])) {
    reg <- jfx_registry()
    if (!name %in% names(reg)) stop("Unknown JSON fixture: ", name, call. = FALSE)
    jfx_cache[[key]] <- reg[[name]]()
  }
  jfx_cache[[key]]
}

jfx_fit_spec <- function(spec) {
  args <- list(dag = spec$dag, data.df = spec$data, data.dists = spec$dists,
               method = spec$method, centre = spec$centre)
  if (!is.null(spec$group.var)) args$group.var <- spec$group.var
  args <- c(args, spec$extra)
  fit <- NULL
  # capture.output silences optimiser traces (e.g. mblogit iterations)
  utils::capture.output(
    fit <- suppressMessages(suppressWarnings(do.call(fitAbn, args)))
  )
  fit
}

# Fitted object for a fixture, cached for the whole test run.
jfx_fit <- function(name) {
  key <- paste0("fit:", name)
  if (is.null(jfx_cache[[key]])) jfx_cache[[key]] <- jfx_fit_spec(jfx_spec(name))
  jfx_cache[[key]]
}

# Exported JSON string for a fixture, cached.
jfx_json <- function(name) {
  key <- paste0("json:", name)
  if (is.null(jfx_cache[[key]])) jfx_cache[[key]] <- export_abnFit(jfx_fit(name))
  jfx_cache[[key]]
}

jfx_doc <- function(name) jdoc_parse(jfx_json(name))

# Bayes fits need a loadable INLA (which in turn needs e.g. 'sf').
jfx_bayes_available <- function() {
  if (is.null(jfx_cache[["bayes_available"]])) {
    jfx_cache[["bayes_available"]] <- isTRUE(tryCatch({
      suppressMessages(suppressWarnings(loadNamespace("INLA")))
      TRUE
    }, error = function(e) FALSE))
  }
  jfx_cache[["bayes_available"]]
}

jfx_skip_if_no_bayes <- function() {
  skip_on_cran()
  skip_if_not(jfx_bayes_available(), "INLA (and its dependencies) not loadable")
}

jfx_skip_bayes <- function(name) {
  if (identical(jfx_spec(name)$method, "bayes")) jfx_skip_if_no_bayes()
}

# ---------------------------------------------------------------------------
# Document accessors
# ---------------------------------------------------------------------------

jdoc_parse <- function(json) jsonlite::fromJSON(json, simplifyVector = FALSE)

jdoc_serialize <- function(doc) {
  jsonlite::toJSON(doc, auto_unbox = TRUE, null = "null", digits = NA)
}

jdoc_chr <- function(x) if (is.null(x)) NA_character_ else as.character(x)

jdoc_ids <- function(rows) vapply(rows, function(r) jdoc_chr(r$`_id`), character(1))

jdoc_variable <- function(doc, name) {
  hit <- Filter(function(v) identical(v[["name"]], name), doc[["variables"]])
  if (length(hit) != 1) stop("Expected exactly one variable named ", name, call. = FALSE)
  hit[[1]]
}

jdoc_variable_by_id <- function(doc, id) {
  hit <- Filter(function(v) identical(jdoc_chr(v$`_id`), jdoc_chr(id)), doc[["variables"]])
  if (length(hit) != 1) stop("Expected exactly one variable with _id ", id, call. = FALSE)
  hit[[1]]
}

jdoc_variable_id <- function(doc, name) jdoc_chr(jdoc_variable(doc, name)$`_id`)

jdoc_state_id <- function(variable, label) {
  hit <- Filter(function(s) identical(as.character(s[["label"]]), as.character(label)),
                variable[["states"]])
  if (length(hit) != 1) {
    stop("Expected exactly one state ", label, " in ", variable[["name"]], call. = FALSE)
  }
  jdoc_chr(hit[[1]]$`_id`)
}

jdoc_state_label <- function(variable, id) {
  hit <- Filter(function(s) identical(jdoc_chr(s$`_id`), jdoc_chr(id)), variable[["states"]])
  if (length(hit) != 1) {
    stop("Expected exactly one state id ", id, " in ", variable[["name"]], call. = FALSE)
  }
  as.character(hit[[1]][["label"]])
}

jdoc_baseline_label <- function(variable) {
  hit <- Filter(function(s) isTRUE(s[["baseline"]]), variable[["states"]])
  if (length(hit) != 1) stop("No unique baseline in ", variable[["name"]], call. = FALSE)
  as.character(hit[[1]][["label"]])
}

jdoc_group <- function(doc, name) {
  hit <- Filter(function(g) identical(g[["name"]], name), doc[["groups"]])
  if (length(hit) != 1) stop("Expected exactly one group named ", name, call. = FALSE)
  hit[[1]]
}

# All parameters of a target (by variable name), optionally filtered by kind.
jdoc_parameters <- function(doc, target = NULL, kind = NULL) {
  target_id <- if (is.null(target)) NULL else jdoc_variable_id(doc, target)
  Filter(function(p) {
    (is.null(target_id) || identical(jdoc_chr(p[["target"]]), target_id)) &&
      (is.null(kind) || p[["kind"]] %in% kind)
  }, doc[["parameters"]])
}

# Exactly one parameter identified by names/labels (not ids).
jdoc_parameter <- function(doc, target, kind, parent = NULL, parent_state = NULL,
                           target_state = NULL) {
  target_var <- jdoc_variable(doc, target)
  want_parent <- if (is.null(parent)) NA_character_ else jdoc_variable_id(doc, parent)
  want_parent_state <- if (is.null(parent_state)) {
    NA_character_
  } else {
    jdoc_state_id(jdoc_variable(doc, parent), parent_state)
  }
  want_target_state <- if (is.null(target_state)) {
    NA_character_
  } else {
    jdoc_state_id(target_var, target_state)
  }
  hit <- Filter(function(p) {
    identical(p[["kind"]], kind) &&
      identical(jdoc_chr(p[["parent"]]), want_parent) &&
      identical(jdoc_chr(p[["parent_state"]]), want_parent_state) &&
      identical(jdoc_chr(p[["target_state"]]), want_target_state)
  }, jdoc_parameters(doc, target))
  if (length(hit) != 1) {
    stop(sprintf("Expected one %s parameter for %s (parent=%s, parent_state=%s, state=%s), got %d",
                 kind, target, want_parent, want_parent_state, want_target_state,
                 length(hit)), call. = FALSE)
  }
  hit[[1]]
}

jdoc_values <- function(parameters) {
  vapply(parameters, function(p) as.numeric(p[["value"]]), numeric(1))
}

jdoc_node_diagnostics <- function(doc, name) {
  id <- jdoc_variable_id(doc, name)
  hit <- Filter(function(n) identical(jdoc_chr(n[["variable"]]), id),
                doc[["inference"]][["diagnostics"]][["nodes"]])
  if (length(hit) != 1) stop("Expected one diagnostics entry for ", name, call. = FALSE)
  hit[[1]]
}

jdoc_drop_extensions <- function(doc) {
  doc[["metadata"]][["extensions"]] <- NULL
  doc
}

jdoc_parents <- function(fit, node) {
  dag <- fit$abnDag$dag
  colnames(dag)[dag[node, ] == 1]
}

# ---------------------------------------------------------------------------
# Core-only evaluation (independent of abn internals)
# ---------------------------------------------------------------------------

# Fixed-effect linear predictor of `target` (and `target_state` for categorical
# targets) for every row of `data`, computed from the generic core only.
jdoc_linear_predictor <- function(doc, data, target, target_state = NULL) {
  target_var <- jdoc_variable(doc, target)
  want_state <- if (is.null(target_state)) {
    NA_character_
  } else {
    jdoc_state_id(target_var, target_state)
  }
  params <- Filter(function(p) {
    p[["kind"]] %in% c("intercept", "coefficient") &&
      identical(jdoc_chr(p[["target_state"]]), want_state)
  }, jdoc_parameters(doc, target))
  eta <- numeric(nrow(data))
  for (p in params) {
    # value null = not estimable, contributes nothing
    if (is.null(p[["value"]])) next
    if (identical(p[["kind"]], "intercept")) {
      eta <- eta + as.numeric(p[["value"]])
      next
    }
    parent_var <- jdoc_variable_by_id(doc, p[["parent"]])
    x <- data[[parent_var[["name"]]]]
    if (!is.null(p[["parent_state"]])) {
      label <- jdoc_state_label(parent_var, p[["parent_state"]])
      eta <- eta + as.numeric(p[["value"]]) * as.numeric(as.character(x) == label)
    } else {
      x <- as.numeric(x)
      if (!is.null(parent_var[["transform"]])) {
        x <- (x - as.numeric(parent_var[["transform"]][["center"]])) /
          as.numeric(parent_var[["transform"]][["scale"]])
      }
      eta <- eta + as.numeric(p[["value"]]) * x
    }
  }
  eta
}

# ---------------------------------------------------------------------------
# Schema validation
# ---------------------------------------------------------------------------

jdoc_schema_file <- function(name = "bayesian-network.schema.json") {
  system.file("schemas", name, package = "abn")
}

jdoc_validator <- function(name = "bayesian-network.schema.json") {
  skip_if_not_installed("jsonvalidate")
  path <- jdoc_schema_file(name)
  if (!nzchar(path)) stop("Schema file not installed: ", name, call. = FALSE)
  jsonvalidate::json_validator(path, engine = "ajv")
}

jdoc_expect_valid <- function(json, name = "bayesian-network.schema.json") {
  validator <- jdoc_validator(name)
  ok <- validator(json, verbose = TRUE, greedy = TRUE)
  errors <- attr(ok, "errors")
  expect_true(isTRUE(as.logical(ok)),
              info = paste(utils::capture.output(print(errors)), collapse = "\n"))
}

jdoc_expect_invalid <- function(json, name = "bayesian-network.schema.json") {
  validator <- jdoc_validator(name)
  expect_false(isTRUE(as.logical(validator(json))))
}

# ---------------------------------------------------------------------------
# A minimal document from a foreign issuer
# ---------------------------------------------------------------------------

# height (continuous) -> status (binary); status is the child.
jdoc_foreign_document <- function() {
  list(
    metadata = list(schema_version = JDOC_SCHEMA_VERSION,
                    issuer = "independent-bn-tool::export"),
    variables = list(
      list(`_id` = 1, name = "height", type = "continuous", distribution = "gaussian",
           link = "identity"),
      list(`_id` = 2, name = "status", type = "binary", distribution = "bernoulli",
           link = "logit",
           states = list(list(`_id` = 1, label = "no", baseline = TRUE),
                         list(`_id` = 2, label = "yes", baseline = FALSE)))
    ),
    groups = list(),
    arcs = list(list(source = 1, target = 2)),
    parameters = list(
      list(`_id` = 1, target = 1, kind = "intercept", value = 1.2,
           uncertainty = list(standard_error = 0.1)),
      list(`_id` = 2, target = 1, kind = "residual_variance", scale = "variance",
           value = 0.25, uncertainty = list(standard_error = NULL)),
      list(`_id` = 3, target = 2, kind = "intercept", value = -0.42,
           uncertainty = list(standard_error = 0.2)),
      list(`_id` = 4, target = 2, kind = "coefficient", parent = 1, value = 0.87,
           uncertainty = list(standard_error = 0.3))
    ),
    inference = list(type = "maximum_likelihood", diagnostics = list(nodes = list()))
  )
}
