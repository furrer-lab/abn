# Internal null-coalescing operator shared by the JSON import/export code.
# (Base R provides `%||%` from 4.4.0; the package still supports R >= 4.0.)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Export a fitted abn model to JSON (bayesian-network)
#'
#' Writes an \code{abnFit} object as a \code{bayesian-network} JSON document.
#' The generic core (variables, groups, arcs, parameters, inference) fully
#' describes the model; abn-internal details are stored in
#' \code{metadata.extensions.abn}. The data are not included; use
#' \code{\link{export_abnData}} for a separate data document.
#'
#' @param object An object of class \code{abnFit}.
#' @param file Optional path. If given, the JSON is written there and the path is
#'   returned invisibly.
#' @param pretty Logical, pretty-print the JSON.
#' @param label Optional descriptive label (\code{metadata.label}).
#' @param scenario_id Optional scenario identifier (\code{metadata.scenario_id}).
#' @param data_reference Optional list \code{(uri, schema_version, sha256)}
#'   pointing to the data document.
#' @return A JSON string, or invisibly \code{file}.
#' @details
#' The document format is specified by \code{inst/schemas/bayesian-network.schema.json}
#' and described in the vignettes \code{vignette("fitabn-json-specification", package = "abn")}
#' and \code{vignette("json-format", package = "abn")}.
#' Re-export reuses stored document identifiers and metadata when available.
#' Each export writes a new \code{metadata.created} timestamp, so exact byte
#' identity is not guaranteed even for an otherwise unchanged model document.
#' @examples
#' mydists <- list(b1 = "binomial", p1 = "poisson", g1 = "gaussian", b3 = "binomial")
#' mydag <- matrix(0, 4, 4, dimnames = list(names(mydists), names(mydists)))
#' mydag["b3", c("b1", "g1")] <- 1
#' fit <- fitAbn(dag = mydag, data.df = ex1.dag.data[, names(mydists)],
#'               data.dists = mydists, method = "mle")
#' json <- export_abnFit(fit, label = "tiny example")
#' stopifnot(isTRUE(all.equal(import_abnFit(json = json)$coef, fit$coef)))
#' @seealso \code{\link{import_abnFit}}, \code{\link{export_abnData}}
#' @export
export_abnFit <- function(object, file = NULL, pretty = TRUE, label = NULL,
                          scenario_id = NULL, data_reference = NULL) {
  if (!inherits(object, "abnFit")) stop("'object' must be of class 'abnFit'.", call. = FALSE)
  if (!object$method %in% c("mle", "bayes")) {
    stop("Unsupported fit method: ", object$method, call. = FALSE)
  }
  stored <- attr(object, "abn_json") %||% list()
  meta <- list(label = label %||% stored$label,
               scenario_id = scenario_id %||% stored$scenario_id,
               data_reference = data_reference %||% stored$data_reference)

  doc <- abn_json_build_document(object, meta)
  abn_json_check_document(doc)

  json <- jsonlite::toJSON(doc, auto_unbox = TRUE, pretty = pretty, null = "null",
                           digits = NA)
  if (!is.null(file)) {
    writeLines(json, con = file)
    return(invisible(file))
  }
  json
}

# ---------------------------------------------------------------------------
# Document assembly
# ---------------------------------------------------------------------------

abn_json_build_document <- function(fit, meta) {
  ctx <- abn_json_context(fit)
  params <- abn_json_build_parameters(fit, ctx)
  list(
    metadata = abn_json_build_metadata(fit, ctx, meta, params),
    variables = abn_json_build_variables(fit, ctx),
    groups = abn_json_build_groups(fit, ctx),
    arcs = abn_json_build_arcs(fit, ctx),
    parameters = lapply(params, function(p) p$json),
    inference = abn_json_build_inference(fit, ctx, params)
  )
}

# Shared lookups: node order, ids, distributions, levels, state ids.
abn_json_context <- function(fit) {
  dag <- fit$abnDag$dag
  nodes <- colnames(dag)
  dists <- unlist(fit$abnDag$data.dists)[nodes]
  levels <- fit$levels %||% list()
  for (node in nodes[dists %in% c("binomial", "multinomial")]) {
    if (is.null(levels[[node]])) {
      stop("Fit has no recorded levels for node '", node,
           "'. Refit with the current version of abn.", call. = FALSE)
    }
  }
  stored <- attr(fit, "abn_json")$ids
  var_id <- stats::setNames(seq_along(nodes), nodes)
  if (!is.null(stored$variables)) {
    hit <- stored$variables[nodes]
    if (!anyNA(hit)) var_id <- hit
  }
  group_id <- NULL
  if (!is.null(fit$group.var)) {
    group_id <- if (!is.null(stored$groups)) stored$groups[[fit$group.var]] else NULL
    if (is.null(group_id)) group_id <- 1L
  }
  list(
    nodes = nodes,
    dists = dists,
    var_id = var_id,
    levels = levels,
    parents = stats::setNames(lapply(nodes, function(n) nodes[dag[n, ] == 1]), nodes),
    group_id = group_id,
    ids = stored
  )
}

abn_json_state_id <- function(ctx, node, label) {
  stored <- ctx$ids$states[[node]][[as.character(label)]]
  if (!is.null(stored)) return(stored)
  id <- match(as.character(label), ctx$levels[[node]])
  if (is.na(id)) stop("Unknown level '", label, "' of node '", node, "'.", call. = FALSE)
  id
}

abn_json_build_metadata <- function(fit, ctx, meta, params) {
  out <- list(
    schema_version = "bayesian-network",
    issuer = "abn::export_abnFit",
    issuer_version = as.character(utils::packageVersion("abn")),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
  if (!is.null(meta$label)) out$label <- meta$label
  if (!is.null(meta$scenario_id)) out$scenario_id <- meta$scenario_id
  if (!is.null(meta$data_reference)) out$data_reference <- meta$data_reference
  out$extensions <- list(abn = abn_json_build_extension(fit, ctx, params))
  out
}

abn_json_dist_map <- function(dist) {
  switch(dist,
         gaussian = list(type = "continuous", distribution = "gaussian", link = "identity"),
         binomial = list(type = "binary", distribution = "bernoulli", link = "logit"),
         poisson = list(type = "count", distribution = "poisson", link = "log"),
         multinomial = list(type = "categorical", distribution = "categorical",
                            link = "baseline_logit"),
         stop("Unsupported distribution: ", dist, call. = FALSE))
}

abn_json_build_variables <- function(fit, ctx) {
  lapply(ctx$nodes, function(node) {
    v <- c(list(`_id` = ctx$var_id[[node]], name = node),
           abn_json_dist_map(ctx$dists[[node]]))
    lv <- ctx$levels[[node]]
    if (!is.null(lv)) {
      v$states <- lapply(seq_along(lv), function(i) {
        list(`_id` = abn_json_state_id(ctx, node, lv[i]), label = lv[i],
             baseline = (i == 1L))
      })
    }
    tr <- fit$centre[[node]]
    if (!is.null(tr)) {
      v$transform <- list(center = unname(tr[["center"]]), scale = unname(tr[["scale"]]))
    }
    v
  })
}

abn_json_build_groups <- function(fit, ctx) {
  if (is.null(fit$group.var)) return(list())
  members <- names(fit$abnDag$data.dists)[fit$grouped.vars]
  list(list(`_id` = ctx$group_id, name = fit$group.var,
            variables = as.list(unname(ctx$var_id[members]))))
}

abn_json_build_arcs <- function(fit, ctx) {
  arcs <- list()
  for (child in ctx$nodes) {
    for (parent in ctx$parents[[child]]) {
      arcs[[length(arcs) + 1L]] <- list(source = ctx$var_id[[parent]],
                                        target = ctx$var_id[[child]])
    }
  }
  arcs
}

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
# Internally each parameter is list(json = <JSON object without _id>, node,
# field, name) where field/name identify the native abn element (extension and
# Bayes marginals lookup). Ids are assigned at the end.

abn_json_param <- function(ctx, node, kind, value, field, name, se = NULL,
                           parent = NULL, parent_state = NULL, target_state = NULL,
                           states = NULL, scale = NULL, group = FALSE) {
  value <- unname(as.numeric(value))
  j <- list(target = ctx$var_id[[node]])
  if (!is.null(target_state)) j$target_state <- abn_json_state_id(ctx, node, target_state)
  j$kind <- kind
  if (!is.null(parent)) {
    j$parent <- ctx$var_id[[parent]]
    if (!is.null(parent_state)) j$parent_state <- abn_json_state_id(ctx, parent, parent_state)
  }
  if (isTRUE(group)) j$group <- ctx$group_id
  if (!is.null(states)) {
    j$states <- lapply(states, function(s) abn_json_state_id(ctx, node, s))
  }
  if (!is.null(scale)) j$scale <- scale
  j["value"] <- list(if (length(value) == 0 || is.na(value)) NULL else value)
  se <- if (is.null(se) || length(se) == 0 || is.na(se)) NULL else unname(as.numeric(se))
  j$uncertainty <- list(standard_error = se)
  list(json = j, node = node, field = field, name = name)
}

# rank mirrors the native order: bayes modes list group.precision before the
# residual precision; grouped MLE order is per-field and rank-independent
abn_json_kind_rank <- c(intercept = 1, coefficient = 2, random_variance = 3,
                        residual_variance = 4, random_covariance = 5)

abn_json_build_parameters <- function(fit, ctx) {
  builder <- switch(abn_json_fit_type(fit),
                    mle = abn_json_params_mle,
                    mle_grouped = abn_json_params_mle_grouped,
                    bayes = abn_json_params_bayes)
  params <- list()
  for (node in ctx$nodes) {
    node_params <- builder(fit, ctx, node)
    rank <- vapply(node_params, function(p) abn_json_kind_rank[[p$json$kind]], numeric(1))
    params <- c(params, node_params[order(rank)])
  }
  for (i in seq_along(params)) {
    params[[i]]$json <- c(list(`_id` = i), params[[i]]$json)
  }
  # reuse stored ids when the parameter sequence is unchanged (stable re-export)
  stored <- ctx$ids$parameters
  if (!is.null(stored) && length(stored) == length(params)) {
    sigs <- vapply(params, function(p) paste(p$node, p$json$kind, p$name, sep = "|"),
                   character(1))
    if (identical(sigs, vapply(stored, function(s) s$sig, character(1)))) {
      ids <- vapply(stored, function(s) s$id, integer(1))
      for (i in seq_along(params)) params[[i]]$json$`_id` <- ids[i]
    }
  }
  params
}

# Design terms of a parent as used in fixed-effect coefficients:
# list of list(parent, state, names = <accepted native names>).
abn_json_parent_terms <- function(ctx, parent, all_levels) {
  dist <- ctx$dists[[parent]]
  lv <- ctx$levels[[parent]]
  if (dist == "multinomial") {
    use <- if (all_levels) lv else lv[-1]
    return(lapply(use, function(l) list(parent = parent, state = l,
                                        names = paste0(parent, l))))
  }
  if (dist == "binomial") {
    return(list(list(parent = parent, state = lv[2],
                     names = c(parent, paste0(parent, lv[2])))))
  }
  list(list(parent = parent, state = NULL, names = parent))
}

abn_json_match_term <- function(terms, name, node) {
  hit <- Filter(function(t) name %in% t$names, terms)
  if (length(hit) != 1) {
    stop("Cannot map coefficient '", name, "' of node '", node,
         "' to a parent term.", call. = FALSE)
  }
  hit[[1]]
}

# Ungrouped MLE: coef / Stderror (+ mse as gaussian residual variance).
abn_json_params_mle <- function(fit, ctx, node) {
  coef <- fit$coef[[node]]
  se <- fit$Stderror[[node]]
  cols <- colnames(coef)
  parents <- ctx$parents[[node]]
  multi_parent <- any(ctx$dists[parents] == "multinomial")
  # design columns in abn order
  design <- list()
  if (!multi_parent) design <- list(list(parent = NULL, state = NULL, intercept = TRUE))
  for (p in parents) {
    for (t in abn_json_parent_terms(ctx, p, all_levels = TRUE)) {
      design[[length(design) + 1L]] <- list(parent = p, state = t$state,
                                            name = if (ctx$dists[[p]] == "multinomial")
                                              t$names else p)
    }
  }
  expected <- list()
  if (ctx$dists[[node]] == "multinomial") {
    child_states <- ctx$levels[[node]][-1]
    for (d in design) {
      for (s in child_states) {
        nm <- if (isTRUE(d$intercept)) paste0(node, "|intercept.", s) else paste0(d$name, s)
        expected[[length(expected) + 1L]] <- c(d, list(target_state = s, native = nm))
      }
    }
  } else {
    for (d in design) {
      nm <- if (isTRUE(d$intercept)) paste0(node, "|intercept") else d$name
      expected[[length(expected) + 1L]] <- c(d, list(target_state = NULL, native = nm))
    }
  }
  natives <- vapply(expected, function(e) e$native, character(1))
  if (ctx$dists[[node]] == "multinomial" && !identical(sort(natives), sort(cols))) {
    # separation fallback (see #268): intercept-only fit despite existing parents
    alt <- lapply(ctx$levels[[node]][-1], function(s) {
      list(intercept = TRUE, parent = NULL, state = NULL, target_state = s,
           native = paste0(node, "|intercept.", s))
    })
    if (identical(sort(vapply(alt, function(e) e$native, character(1))), sort(cols))) {
      expected <- alt
      natives <- vapply(expected, function(e) e$native, character(1))
    }
  }
  if (!identical(sort(natives), sort(cols))) {
    stop("Unexpected coefficient names for node '", node, "': ",
         paste(cols, collapse = ", "), "; expected: ", paste(natives, collapse = ", "),
         call. = FALSE)
  }
  # emit parameters in the order of the original coef columns, which follows
  # the design matrix (e.g. one-hot parents in raw factor order)
  expected <- expected[match(cols, natives)]
  params <- lapply(expected, function(e) {
    kind <- if (isTRUE(e$intercept)) "intercept" else "coefficient"
    abn_json_param(ctx, node, kind, coef[1, e$native], field = "coef", name = e$native,
                   se = se[1, e$native], parent = e$parent, parent_state = e$state,
                   target_state = e$target_state)
  })
  if (ctx$dists[[node]] == "gaussian") {
    params[[length(params) + 1L]] <- abn_json_param(
      ctx, node, "residual_variance", fit$mse[[node]], field = "mse", name = node,
      scale = "variance")
  }
  params
}

# Grouped MLE: mu / betas / sigma / sigma_alpha (lme4, glmmTMB, mblogit).
abn_json_params_mle_grouped <- function(fit, ctx, node) {
  mu <- fit$mu[[node]]
  betas <- fit$betas[[node]]
  sigma <- fit$sigma[[node]]
  sa <- fit$sigma_alpha[[node]]
  terms <- unlist(lapply(ctx$parents[[node]], abn_json_parent_terms, ctx = ctx,
                         all_levels = FALSE), recursive = FALSE)
  params <- list()
  add <- function(p) params[[length(params) + 1L]] <<- p
  present <- function(x) !is.null(x) && length(x) > 0 && !all(is.na(x))

  if (ctx$dists[[node]] == "multinomial") {
    states <- ctx$levels[[node]][-1]
    mu_names <- names(mu) %||% paste0(node, ".", states)
    for (i in seq_along(states)) {
      add(abn_json_param(ctx, node, "intercept", mu[[i]], field = "mu", name = mu_names[i],
                         target_state = states[i]))
    }
    if (present(betas)) {
      betas <- as.matrix(betas)
      for (i in seq_along(states)) {
        for (cn in colnames(betas)) {
          t <- abn_json_match_term(terms, cn, node)
          add(abn_json_param(ctx, node, "coefficient", betas[i, cn], field = "betas",
                             name = cn,
                             parent = t$parent, parent_state = t$state,
                             target_state = states[i]))
        }
      }
    }
    if (present(sa)) {
      sa <- as.matrix(sa)
      for (i in seq_along(states)) {
        for (k in seq_along(states)) {
          if (k < i) next
          nm <- paste(rownames(sa)[i] %||% i, colnames(sa)[k] %||% k, sep = ",")
          if (i == k) {
            add(abn_json_param(ctx, node, "random_variance", sa[i, i], field = "sigma_alpha",
                               name = rownames(sa)[i], scale = "variance", group = TRUE,
                               target_state = states[i]))
          } else {
            add(abn_json_param(ctx, node, "random_covariance", sa[i, k],
                               field = "sigma_alpha", name = nm, scale = "variance",
                               group = TRUE, states = states[c(i, k)]))
          }
        }
      }
    }
    return(params)
  }

  add(abn_json_param(ctx, node, "intercept", mu, field = "mu", name = "(Intercept)"))
  if (present(betas)) {
    for (bn in names(betas)) {
      t <- abn_json_match_term(terms, bn, node)
      add(abn_json_param(ctx, node, "coefficient", betas[[bn]], field = "betas", name = bn,
                         parent = t$parent, parent_state = t$state))
    }
  }
  if (present(sigma) && ctx$dists[[node]] == "gaussian") {
    add(abn_json_param(ctx, node, "residual_variance", sigma[[1]]^2, field = "sigma",
                       name = node, scale = "variance"))
  }
  if (present(sa)) {
    add(abn_json_param(ctx, node, "random_variance", sa[[1]]^2, field = "sigma_alpha",
                       name = node, scale = "variance", group = TRUE))
  }
  params
}

# Bayes: posterior modes (+ quantiles from marginal.quantiles).
abn_json_params_bayes <- function(fit, ctx, node) {
  modes <- fit$modes[[node]]
  lapply(names(modes), function(native) {
    term <- sub(paste0("^", node, "\\|"), "", native, fixed = FALSE)
    args <- list(ctx = ctx, node = node, value = modes[[native]], field = "modes",
                 name = native)
    if (term == "(Intercept)") {
      args$kind <- "intercept"
    } else if (term == "precision") {
      args$kind <- "residual_variance"; args$scale <- "precision"
    } else if (term == "group.precision") {
      args$kind <- "random_variance"; args$scale <- "precision"; args$group <- TRUE
    } else if (term %in% ctx$parents[[node]]) {
      t <- abn_json_parent_terms(ctx, term, all_levels = FALSE)[[1]]
      args$kind <- "coefficient"; args$parent <- term; args$parent_state <- t$state
    } else {
      stop("Cannot map posterior mode '", native, "'.", call. = FALSE)
    }
    p <- do.call(abn_json_param, args)
    q <- fit$marginal.quantiles[[node]][[native]]
    if (!is.null(q)) {
      p$json$uncertainty$posterior_quantiles <- lapply(seq_len(nrow(q)), function(i) {
        list(probability = unname(q[i, "P(X<=x)"]), value = unname(q[i, "x"]))
      })
    }
    p
  })
}

# ---------------------------------------------------------------------------
# Inference
# ---------------------------------------------------------------------------

abn_json_num <- function(x) {
  x <- unname(as.numeric(x))
  if (length(x) == 0 || is.na(x)) NULL else x
}

abn_json_build_inference <- function(fit, ctx, params) {
  bayes <- identical(fit$method, "bayes")
  out <- list(type = if (bayes) "bayesian" else "maximum_likelihood")
  if (bayes && !is.null(fit$priors)) {
    out$priors <- list(
      list(applies_to = "fixed_effects", family = "normal",
           mean = fit$priors$mean, precision = fit$priors$prec),
      list(applies_to = "precisions", family = "gamma",
           shape = fit$priors$loggam.shape, rate = fit$priors$loggam.inv.scale)
    )
  }
  diag <- list(log_marginal_likelihood = abn_json_num(fit$mlik))
  if (!bayes) {
    diag$aic <- abn_json_num(fit$aic)
    diag$bic <- abn_json_num(fit$bic)
  }
  grouped <- !is.null(fit$group.var)
  diag$nodes <- lapply(ctx$nodes, function(node) {
    d <- list(variable = ctx$var_id[[node]],
              log_marginal_likelihood = abn_json_num(fit$mliknode[[node]]))
    if (!bayes) {
      d$aic <- abn_json_num(fit$aicnode[[node]])
      d$bic <- abn_json_num(fit$bicnode[[node]])
      d$mdl <- abn_json_num(fit$mdlnode[[node]])
      d$df <- abn_json_num(fit$df[[node]])
      d$sse <- abn_json_num(fit$sse[[node]])
      # ungrouped gaussian mse is the residual_variance parameter
      if (grouped || ctx$dists[[node]] != "gaussian") d$mse <- abn_json_num(fit$mse[[node]])
    }
    Filter(Negate(is.null), d)
  })
  out$diagnostics <- Filter(Negate(is.null), diag)

  if (bayes && !is.null(fit$marginals)) {
    marginals <- list()
    for (p in params) {
      m <- fit$marginals[[p$node]][[p$name]]
      if (is.null(m)) next
      marginals[[length(marginals) + 1L]] <- list(
        parameter = p$json$`_id`,
        x = as.list(unname(m[, "x"])),
        density = as.list(unname(m[, "f(x)"])))
    }
    out$posterior <- list(marginals = marginals)
  }
  out
}

# ---------------------------------------------------------------------------
# abn extension (abn-internal information only)
# ---------------------------------------------------------------------------

abn_json_build_extension <- function(fit, ctx, params) {
  ext <- list(parameter_names = lapply(params, function(p) {
    list(parameter = p$json$`_id`, field = p$field, name = p$name)
  }))
  if (identical(fit$method, "bayes")) {
    ext$nodes <- lapply(ctx$nodes, function(node) {
      # some of these vectors lose their names (ifelse); order = DAG order
      val <- function(x) {
        x <- if (!is.null(names(x))) x[[node]] else x[[match(node, ctx$nodes)]]
        if (is.null(x) || length(x) == 0 || is.na(x)) NULL else unname(x)
      }
      Filter(Negate(is.null), list(
        variable = ctx$var_id[[node]],
        used_inla = val(fit$used.INLA),
        error_code = val(fit$error.code),
        error_code_desc = val(fit$error.code.desc),
        hessian_accuracy = val(fit$hessian.accuracy)))
    })
  }
  ext
}

# ---------------------------------------------------------------------------
# Structural check of a (parsed or built) network document. Used by export
# before writing and by import before reconstruction.
# ---------------------------------------------------------------------------

abn_json_check_document <- function(doc) {
  fail <- function(...) stop(..., call. = FALSE)
  abn_json_assert_object(doc, "network document")
  abn_json_assert_object(doc[["metadata"]], "metadata", c("schema_version", "issuer"))
  abn_json_assert_string(doc[["metadata"]][["schema_version"]], "metadata.schema_version")
  abn_json_assert_string(doc[["metadata"]][["issuer"]], "metadata.issuer", nonempty = TRUE)
  if (is.null(doc[["metadata"]][["schema_version"]])) {
    fail("Unsupported schema_version ''; expected 'bayesian-network'.")
  }
  if (!identical(doc[["metadata"]][["schema_version"]], "bayesian-network")) {
    fail("Unsupported schema_version '", doc[["metadata"]][["schema_version"]],
         "'; expected 'bayesian-network'.")
  }
  if (!is.character(doc[["metadata"]][["issuer"]]) ||
      length(doc[["metadata"]][["issuer"]]) != 1L ||
      !nzchar(doc[["metadata"]][["issuer"]])) {
    fail("metadata.issuer must be a non-empty string.")
  }
  for (block in c("variables", "groups", "arcs", "parameters", "inference")) {
    if (is.null(doc[[block]])) fail("Missing block '", block, "'.")
    if (block == "inference") {
      abn_json_assert_object(doc[[block]], block, "type")
    } else {
      abn_json_assert_array(doc[[block]], block)
    }
  }

  whole <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= 1 && x <= .Machine$integer.max && x == trunc(x)
  }
  ids <- function(rows, what) {
    if (length(rows) == 0L) return(character(0))
    out <- vapply(rows, function(r) {
      if (!is.list(r)) fail("Every element of ", what, " must be a JSON object.")
      if (!whole(r[["_id"]])) return(NA_character_)
      as.character(round(r[["_id"]]))
    }, character(1))
    if (anyNA(out) || anyDuplicated(out)) {
      fail("Missing, invalid or duplicated _id in ", what, ".")
    }
    out
  }

  var_ids <- ids(doc[["variables"]], "variables")
  param_ids <- ids(doc[["parameters"]], "parameters")
  group_ids <- ids(doc[["groups"]], "groups")
  vars <- stats::setNames(doc[["variables"]], var_ids)
  var_names <- vapply(vars, function(v) {
    if (!is.character(v[["name"]]) || length(v[["name"]]) != 1L ||
        !nzchar(v[["name"]])) {
      fail("Every variable needs a non-empty name.")
    }
    as.character(v[["name"]])
  }, character(1))
  if (anyDuplicated(var_names)) fail("Variable names must be unique.")

  dist_map <- list(
    gaussian = c(type = "continuous", distribution = "gaussian", link = "identity"),
    bernoulli = c(type = "binary", distribution = "bernoulli", link = "logit"),
    poisson = c(type = "count", distribution = "poisson", link = "log"),
    categorical = c(type = "categorical", distribution = "categorical",
                    link = "baseline_logit")
  )
  for (id in var_ids) {
    v <- vars[[id]]
    abn_json_assert_object(v, paste0("variable ", id),
                            c("_id", "name", "type", "distribution", "link"))
    for (field in c("name", "type", "distribution", "link")) {
      abn_json_assert_string(v[[field]], paste0("variable ", id, " ", field), nonempty = TRUE)
    }
    distribution <- v[["distribution"]]
    if (!is.character(distribution) || length(distribution) != 1L ||
        !distribution %in% names(dist_map)) {
      fail("Variable '", v[["name"]], "' has an unsupported distribution.")
    }
    expected <- dist_map[[distribution]]
    observed <- c(type = as.character(v[["type"]]),
                  distribution = distribution,
                  link = as.character(v[["link"]]))
    if (!identical(unname(observed), unname(expected))) {
      fail("Variable '", v[["name"]], "' has incompatible type, distribution and link.")
    }
    states <- v[["states"]]
    if (!is.null(states)) {
      abn_json_assert_array(states, paste0("states of ", v[["name"]]))
      for (state in states) {
        abn_json_assert_object(state, "state", c("_id", "label", "baseline"))
        abn_json_assert_string(state[["label"]], "state label")
        abn_json_assert_boolean(state[["baseline"]], "state baseline")
      }
    }
    if (distribution %in% c("bernoulli", "categorical")) {
      if (!is.list(states) || length(states) < 2L ||
          (distribution == "bernoulli" && length(states) != 2L)) {
        fail("Variable '", v[["name"]], "' needs the correct number of categorical states.")
      }
    } else if (!is.null(states)) {
      fail("Variable '", v[["name"]], "' must not declare states.")
    }
    if (!is.null(v[["transform"]])) {
      if (distribution != "gaussian" || !is.list(v[["transform"]])) {
        fail("Only gaussian variables may declare a transform.")
      }
      abn_json_assert_object(v[["transform"]], "transform", c("center", "scale"))
      center <- v[["transform"]][["center"]]
      scale <- v[["transform"]][["scale"]]
      if (!is.numeric(center) || !is.numeric(scale)) fail("Transform values must be numeric.")
      if (length(center) != 1L || length(scale) != 1L || !is.finite(center) ||
          !is.finite(scale) || scale <= 0) {
        fail("Transform of variable '", v[["name"]], "' needs a finite center and positive scale.")
      }
    }
  }

  state_ids <- lapply(vars, function(v) {
    s <- v[["states"]]
    if (is.null(s)) return(character(0))
    baseline <- vapply(s, function(x) isTRUE(x[["baseline"]]), logical(1))
    if (sum(baseline) != 1L) {
      fail("Variable '", v[["name"]], "' needs exactly one baseline state.")
    }
    labels <- vapply(s, function(x) {
      if (!is.character(x[["label"]]) || length(x[["label"]]) != 1L ||
          is.na(x[["label"]])) {
        fail("Every state of variable '", v[["name"]], "' needs a string label.")
      }
      as.character(x[["label"]])
    }, character(1))
    if (anyDuplicated(labels)) fail("State labels of variable '", v[["name"]], "' must be unique.")
    ids(s, paste0("states of ", v[["name"]]))
  })
  baseline_id <- function(id) {
    s <- vars[[id]][["states"]]
    as.character(round(s[[which(vapply(s, function(x) isTRUE(x[["baseline"]]), logical(1)))]][["_id"]]))
  }

  group_names <- vapply(doc[["groups"]], function(g) {
    if (!is.list(g)) fail("Every group must be a JSON object.")
    if (!is.character(g[["name"]]) || length(g[["name"]]) != 1L ||
        !nzchar(g[["name"]])) fail("Every group needs a non-empty name.")
    as.character(g[["name"]])
  }, character(1))
  if (anyDuplicated(group_names)) fail("Duplicated group names.")
  for (g in doc[["groups"]]) {
    abn_json_assert_object(g, "group", c("_id", "name", "variables"))
    members <- g[["variables"]]
    abn_json_assert_array(members, "group variables")
    if (!is.list(members) || length(members) == 0L) {
      fail("Group '", g[["name"]], "' needs at least one variable.")
    }
    member_ids <- vapply(members, function(id) {
      if (!whole(id)) fail("Group '", g[["name"]], "' has an invalid variable _id.")
      variable <- as.character(round(id))
      if (!variable %in% var_ids) {
        fail("Group '", g[["name"]], "' references unknown variable _id ", id, ".")
      }
      variable
    }, character(1))
    if (anyDuplicated(member_ids)) fail("Group '", g[["name"]], "' has a duplicated variable.")
  }

  ref <- function(id, what) {
    if (!whole(id)) fail("Invalid variable _id referenced by ", what, ".")
    variable <- as.character(round(id))
    if (!variable %in% var_ids) {
      fail("Unknown variable _id ", variable, " referenced by ", what, ".")
    }
    variable
  }
  arc_keys <- character(0)
  for (a in doc[["arcs"]]) {
    abn_json_assert_object(a, "arc", c("source", "target"))
    source <- ref(a[["source"]], "arc")
    target <- ref(a[["target"]], "arc")
    if (identical(source, target)) fail("Arcs must not be self-links.")
    key <- paste(source, target, sep = "->")
    if (key %in% arc_keys) fail("Duplicate arc ", source, " -> ", target, ".")
    arc_keys <- c(arc_keys, key)
  }

  check_scalar <- function(value, what, allow_null = FALSE) {
    if (allow_null && is.null(value)) return(TRUE)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      fail(what, " must be a scalar number", if (allow_null) " or null" else "", ".")
    }
    TRUE
  }
  check_state_ref <- function(id, state_list, what) {
    if (!whole(id)) fail("Invalid ", what, " _id.")
    state <- as.character(round(id))
    if (!state %in% state_list) fail("Unknown ", what, " ", state, ".")
    state
  }

  allowed_kinds <- c("intercept", "coefficient", "residual_variance",
                     "random_variance", "random_covariance")
  variance_kinds <- c("residual_variance", "random_variance", "random_covariance")
  seen_parameters <- character(0)
  for (p in doc[["parameters"]]) {
    abn_json_assert_object(p, "parameter", c("_id", "target", "kind", "value"))
    if (!"value" %in% names(p)) fail("Parameter ", p[["_id"]], " is missing required field 'value'.")
    check_scalar(p[["value"]], "Parameter value", allow_null = TRUE)
    kind <- p[["kind"]]
    if (!is.character(kind) || length(kind) != 1L || !kind %in% allowed_kinds) {
      fail("Unknown parameter kind '", kind %||% "", "' in parameter ", p[["_id"]], ".")
    }
    target <- ref(p[["target"]], "parameter")
    target_dist <- vars[[target]][["distribution"]]
    if (kind == "intercept" && !is.null(p[["parent"]])) {
      fail("Intercept parameter ", p[["_id"]], " must not have a parent.")
    }
    if (kind == "coefficient" && is.null(p[["parent"]])) {
      fail("Coefficient parameter ", p[["_id"]], " needs a parent.")
    }
    if (kind != "coefficient" && kind != "intercept" && !is.null(p[["parent"]])) {
      fail("Only coefficient parameters may have a parent.")
    }
    if (kind %in% variance_kinds &&
        (!is.character(p[["scale"]]) || length(p[["scale"]]) != 1L ||
         !p[["scale"]] %in% c("variance", "precision"))) {
      fail("Parameter ", p[["_id"]], " needs scale 'variance' or 'precision'.")
    }
    if (kind %in% c("residual_variance", "random_variance") && !is.null(p[["value"]])) {
      if (p[["value"]] < 0 || (identical(p[["scale"]], "precision") && p[["value"]] == 0)) {
        fail("Parameter ", p[["_id"]], " needs non-negative variance or strictly positive precision.")
      }
    }
    if (!kind %in% variance_kinds && !is.null(p[["scale"]])) {
      fail("Parameter ", p[["_id"]], " must not have a scale.")
    }
    if (kind == "residual_variance" && !identical(target_dist, "gaussian")) {
      fail("Residual variance parameter ", p[["_id"]], " must target a gaussian variable.")
    }
    if (kind %in% c("random_variance", "random_covariance")) {
      if (!whole(p[["group"]])) fail("Parameter ", p[["_id"]], " needs a valid group _id.")
      if (!as.character(round(p[["group"]])) %in% group_ids) {
        fail("Unknown group ", p[["group"]], " in parameter ", p[["_id"]], ".")
      }
      group <- doc[["groups"]][[match(as.character(p[["group"]]), group_ids)]]
      if (!target %in% vapply(group[["variables"]], as.character, character(1))) {
        fail("Parameter ", p[["_id"]], " targets a variable outside its group.")
      }
    }

    target_state <- NA_character_
    if (!is.null(p[["target_state"]])) {
      if (!target_dist %in% c("bernoulli", "categorical")) {
        fail("Parameter ", p[["_id"]], " has target_state for a non-categorical target.")
      }
      target_state <- check_state_ref(p[["target_state"]], state_ids[[target]],
                                      "target_state")
      if (target_state == baseline_id(target)) {
        fail("Parameter ", p[["_id"]], " uses the baseline target state.")
      }
    }
    if (!is.null(p[["parent_state"]])) {
      if (is.null(p[["parent"]])) fail("Parameter ", p[["_id"]], " needs a parent for parent_state.")
      parent_tmp <- ref(p[["parent"]], "parameter")
      if (!vars[[parent_tmp]][["distribution"]] %in% c("bernoulli", "categorical")) {
        fail("Parameter ", p[["_id"]], " has parent_state for a non-categorical parent.")
      }
    }
    parent <- NA_character_
    parent_state <- NA_character_
    if (!is.null(p[["parent_state"]])) {
      parent_state <- p[["parent_state"]]
    }
    if (!is.null(p[["parent"]])) {
      parent <- ref(p[["parent"]], "parameter")
      if (!paste(parent, target, sep = "->") %in% arc_keys) {
        fail("Parameter ", p[["_id"]], " has parent ", parent,
             " but there is no arc ", parent, " -> ", target, ".")
      }
      if (!is.null(p[["parent_state"]])) {
        parent_state <- check_state_ref(p[["parent_state"]], state_ids[[parent]], "parent_state")
      }
    }
    if (kind == "random_covariance") {
      if (!identical(target_dist, "categorical")) {
        fail("Random covariance parameter ", p[["_id"]], " must target a categorical variable.")
      }
      if (!is.list(p[["states"]]) || length(p[["states"]]) != 2L) {
        fail("Random covariance parameter ", p[["_id"]], " needs exactly two target states.")
      }
      cov_states <- vapply(p[["states"]], check_state_ref, character(1),
                           state_list = state_ids[[target]], what = "random-covariance state")
      if (anyDuplicated(cov_states) || any(cov_states == baseline_id(target))) {
        fail("Random covariance parameter ", p[["_id"]], " needs distinct non-baseline states.")
      }
      if (!is.null(p[["target_state"]])) {
        fail("Random covariance parameter ", p[["_id"]], " must not declare target_state.")
      }
    } else if (!is.null(p[["states"]])) {
      fail("Only random_covariance parameters may declare states.")
    }
    if (!is.null(p[["group"]]) && !kind %in% c("random_variance", "random_covariance")) {
      fail("Only random-effects parameters may declare a group.")
    }
    if (!is.null(p[["uncertainty"]])) {
      abn_json_assert_object(p[["uncertainty"]], "parameter uncertainty")
      se <- p[["uncertainty"]][["standard_error"]]
      if (!is.null(se)) {
        check_scalar(se, "Parameter standard_error")
        if (se < 0) fail("Parameter standard_error must be non-negative.")
      }
      quantiles <- p[["uncertainty"]][["posterior_quantiles"]]
      if (!is.null(quantiles)) {
        abn_json_assert_array(quantiles, "posterior quantiles")
        if (!is.list(quantiles) || length(quantiles) == 0L) {
          fail("Posterior quantiles must be a non-empty JSON array.")
        }
        probabilities <- vapply(quantiles, function(q) {
          abn_json_assert_object(q, "posterior quantile", c("probability", "value"))
          check_scalar(q[["probability"]], "Posterior quantile probability")
          check_scalar(q[["value"]], "Posterior quantile value")
          q[["probability"]]
        }, numeric(1))
        if (any(probabilities < 0 | probabilities > 1) || anyDuplicated(probabilities)) {
          fail("Posterior quantile probabilities must be unique and lie in [0, 1].")
        }
      }
    }
    state_key <- paste(vapply(p[["states"]] %||% list(), as.character, character(1)),
                       collapse = ",")
    key <- paste(target, kind, parent, parent_state, target_state, state_key, sep = "|")
    if (key %in% seen_parameters) fail("Duplicated parameter semantics for target ", target, ".")
    seen_parameters <- c(seen_parameters, key)
  }

  inference <- doc[["inference"]]
  if (!is.list(inference)) fail("inference must be a JSON object.")
  inference_type <- inference[["type"]]
  if (!is.character(inference_type) || length(inference_type) != 1L ||
      !inference_type %in% c("maximum_likelihood", "bayesian")) {
    fail("inference.type must be 'maximum_likelihood' or 'bayesian'.")
  }
  diagnostics <- inference[["diagnostics"]]
  if (!is.null(diagnostics)) {
    abn_json_assert_object(diagnostics, "inference.diagnostics")
    for (field in intersect(names(diagnostics), c("log_marginal_likelihood", "aic", "bic"))) {
      check_scalar(diagnostics[[field]], paste0("Diagnostic field '", field, "'"))
    }
    diag_nodes <- diagnostics[["nodes"]]
    if (!is.null(diag_nodes)) {
      abn_json_assert_array(diag_nodes, "inference.diagnostics.nodes")
      diag_refs <- character(0)
      for (node in diag_nodes) {
        abn_json_assert_object(node, "diagnostics entry", "variable")
        variable <- ref(node[["variable"]], "diagnostics")
        if (variable %in% diag_refs) fail("Duplicated diagnostics variable ", variable, ".")
        diag_refs <- c(diag_refs, variable)
        for (field in intersect(names(node), c("log_marginal_likelihood", "aic", "bic",
                                                "mdl", "df", "sse", "mse"))) {
          check_scalar(node[[field]], paste0("Diagnostic field '", field, "'"))
        }
      }
    }
  }
  posterior <- inference[["posterior"]]
  if (!is.null(posterior)) {
    abn_json_assert_object(posterior, "inference.posterior")
    marginals <- posterior[["marginals"]]
    if (!is.null(marginals)) {
      abn_json_assert_array(marginals, "inference.posterior.marginals")
      marginal_refs <- character(0)
      for (marginal in marginals) {
        abn_json_assert_object(marginal, "posterior marginal", c("parameter", "x", "density"))
        if (!whole(marginal[["parameter"]])) fail("Invalid posterior parameter _id.")
        parameter <- as.character(marginal[["parameter"]])
        if (length(parameter) != 1L || is.na(parameter) || !parameter %in% param_ids) {
          fail("Unknown parameter _id ", marginal[["parameter"]], " in posterior marginals.")
        }
        if (parameter %in% marginal_refs) {
          fail("Duplicated posterior marginal for parameter ", parameter, ".")
        }
        marginal_refs <- c(marginal_refs, parameter)
        x <- marginal[["x"]]; density <- marginal[["density"]]
        abn_json_assert_array(x, "posterior x")
        abn_json_assert_array(density, "posterior density")
        if (!is.list(x) || !is.list(density) || length(x) == 0L ||
            length(x) != length(density)) {
          fail("Posterior marginal ", parameter,
               " needs x and density arrays of equal, non-zero length.")
        }
        for (value in x) check_scalar(value, "Posterior x value")
        density_values <- vapply(density, function(value) {
          check_scalar(value, "Posterior density value")
          value
        }, numeric(1))
        if (any(density_values < 0)) fail("Posterior density values must be non-negative.")
      }
    }
  }

  extensions <- doc[["metadata"]][["extensions"]]
  if (!is.null(extensions)) {
    abn_json_assert_object(extensions, "metadata.extensions")
    ext <- extensions[["abn"]]
    if (!is.null(ext)) {
      abn_json_assert_object(ext, "metadata.extensions.abn")
      parameter_names <- ext[["parameter_names"]]
      if (!is.null(parameter_names)) {
        abn_json_assert_array(parameter_names, "metadata.extensions.abn.parameter_names")
        parameter_refs <- character(0)
        for (entry in parameter_names) {
          abn_json_assert_object(entry, "parameter_names entry", c("parameter", "field", "name"))
          if (!whole(entry[["parameter"]])) fail("Invalid abn extension parameter _id.")
          parameter <- as.character(entry[["parameter"]])
          if (length(parameter) != 1L || is.na(parameter) || !parameter %in% param_ids) {
            fail("Unknown parameter _id ", entry[["parameter"]], " in abn extension.")
          }
          if (parameter %in% parameter_refs) {
            fail("Duplicated parameter_names entry for parameter ", parameter, ".")
          }
          parameter_refs <- c(parameter_refs, parameter)
          for (field in c("field", "name")) {
            if (!is.character(entry[[field]]) || length(entry[[field]]) != 1L ||
                !nzchar(entry[[field]])) {
              fail("metadata.extensions.abn.parameter_names entries need a non-empty ", field, ".")
            }
          }
        }
      }
      nodes <- ext[["nodes"]]
      if (!is.null(nodes)) {
        abn_json_assert_array(nodes, "metadata.extensions.abn.nodes")
        node_refs <- character(0)
        for (node in nodes) {
          abn_json_assert_object(node, "abn nodes extension entry", "variable")
          if (!is.null(node[["used_inla"]])) abn_json_assert_boolean(node[["used_inla"]], "used_inla")
          if (!is.null(node[["error_code_desc"]])) abn_json_assert_string(node[["error_code_desc"]], "error_code_desc")
          for (field in intersect(names(node), c("error_code", "hessian_accuracy"))) {
            check_scalar(node[[field]], field)
          }
          variable <- ref(node[["variable"]], "abn extension")
          if (variable %in% node_refs) {
            fail("Duplicated nodes extension entry for variable ", variable, ".")
          }
          node_refs <- c(node_refs, variable)
        }
      }
    }
  }

  # acyclicity (Kahn)
  edges <- do.call(rbind, lapply(doc[["arcs"]], function(a) {
    c(as.character(round(a[["source"]])), as.character(round(a[["target"]])))
  }))
  if (!is.null(edges)) {
    remaining <- var_ids
    repeat {
      roots <- remaining[!remaining %in% edges[edges[, 1] %in% remaining, 2]]
      if (length(roots) == 0L) break
      remaining <- setdiff(remaining, roots)
    }
    if (length(remaining) > 0L) fail("The arcs contain a cycle.")
  }
  invisible(TRUE)
}
