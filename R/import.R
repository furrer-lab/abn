# Import of bayesian-network JSON documents into abnFit objects.
# See json_plan.md and inst/schemas/bayesian-network.schema.json.

#' Import an abnFit object from a bayesian-network JSON document
#'
#' Reconstructs an \code{abnFit} object from a JSON document as produced by
#' \code{\link{export_abnFit}}. Documents from other issuers are imported from
#' the generic core; abn-specific fields are restored from
#' \code{metadata.extensions.abn} when present.
#'
#' @param file Optional path to a JSON file.
#' @param json Optional JSON string. If both are given, \code{json} wins.
#' @param data Optional observations: a data frame or the path to a
#'   \code{bn-data} data document. The observations are transformed exactly as
#'   \code{\link{fitAbn}} does and attached to \code{abnDag$data.df}; for grouped
#'   fits the grouping column provides \code{group.ids}.
#' @param validate If \code{TRUE}, the document is checked structurally
#'   (schema version, unique ids, resolvable references, acyclic arcs).
#' @return An object of class \code{abnFit}.
#' @seealso \code{\link{export_abnFit}}, \code{\link{fitAbn}},
#'   \code{\link{import_abnData}}
#' @importFrom jsonlite fromJSON
#' @export
import_abnFit <- function(file = NULL, json = NULL, data = NULL, validate = TRUE) {
  if (is.null(file) && is.null(json)) {
    stop("Either 'file' or 'json' must be provided.", call. = FALSE)
  }
  if (!is.null(file) && !is.null(json)) {
    stop("Provide only one of 'file' or 'json'.", call. = FALSE)
  }
  if (!is.null(json)) {
    doc <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  } else {
    if (!file.exists(file)) stop(sprintf("File '%s' does not exist.", file), call. = FALSE)
    doc <- jsonlite::fromJSON(paste(readLines(file, warn = FALSE), collapse = "\n"),
                              simplifyVector = FALSE)
  }
  if (isTRUE(validate)) abn_json_check_document(doc)

  ctx <- abn_json_import_context(doc)
  is_abn <- identical(doc[["metadata"]][["issuer"]], "abn::export_abnFit")
  natives <- abn_json_native_names(doc, ctx, is_abn)

  fit <- switch(ctx$type,
                maximum_likelihood = if (length(doc[["groups"]]) > 0) {
                  abn_json_reconstruct_mle_grouped(doc, ctx, natives)
                } else {
                  abn_json_reconstruct_mle(doc, ctx, natives)
                },
                bayesian = abn_json_reconstruct_bayes(doc, ctx, natives),
                stop("Unsupported inference type: ", ctx$type, call. = FALSE))

  fit <- abn_json_attach_data(fit, ctx, data)
  attr(fit, "abn_json") <- abn_json_import_attr(doc, natives)
  class(fit) <- c("abnFit")
  fit
}

# ---------------------------------------------------------------------------
# Context
# ---------------------------------------------------------------------------

abn_json_import_context <- function(doc) {
  vars <- doc[["variables"]]
  names <- vapply(vars, function(v) v[["name"]], character(1))
  dists_flat <- stats::setNames(vapply(vars, function(v) {
    switch(v[["distribution"]], gaussian = "gaussian", bernoulli = "binomial",
           poisson = "poisson", categorical = "multinomial",
           stop("Unknown distribution: ", v[["distribution"]], call. = FALSE))
  }, character(1)), names)

  levels <- stats::setNames(lapply(vars, function(v) {
    s <- v[["states"]]
    if (is.null(s)) return(NULL)
    base <- which(vapply(s, function(x) isTRUE(x[["baseline"]]), logical(1)))
    labels <- vapply(s, function(x) as.character(x[["label"]]), character(1))
    labels[c(base, seq_along(labels)[-base])]
  }), names)

  centre <- stats::setNames(lapply(vars, function(v) {
    tr <- v[["transform"]]
    if (is.null(tr)) return(NULL)
    c(center = as.numeric(tr[["center"]]), scale = as.numeric(tr[["scale"]]))
  }), names)
  centre <- Filter(Negate(is.null), centre)
  if (length(centre) == 0) centre <- NULL

  dag <- matrix(0, length(names), length(names), dimnames = list(names, names))
  for (a in doc[["arcs"]]) {
    dag[abn_json_var_name(doc, a[["target"]]), abn_json_var_name(doc, a[["source"]])] <- 1
  }

  list(names = names, dists = as.list(dists_flat), dists_flat = dists_flat,
       levels = levels, centre = centre, dag = dag,
       type = doc[["inference"]][["type"]])
}

abn_json_var_name <- function(doc, id) {
  hit <- Filter(function(v) identical(as.character(v[["_id"]]), as.character(id)),
                doc[["variables"]])
  if (length(hit) != 1) stop("Unknown variable _id ", id, ".", call. = FALSE)
  hit[[1]][["name"]]
}

abn_json_state_label <- function(doc, var_id, state_id) {
  v <- Filter(function(x) identical(as.character(x[["_id"]]), as.character(var_id)),
              doc[["variables"]])[[1]]
  hit <- Filter(function(s) identical(as.character(s[["_id"]]), as.character(state_id)),
                v[["states"]] %||% list())
  if (length(hit) != 1) {
    stop("Unknown state ", state_id, " in variable ", v[["name"]], ".", call. = FALSE)
  }
  as.character(hit[[1]][["label"]])
}

# Parameters grouped by target node, in document order.
abn_json_params_by_node <- function(doc, ctx) {
  out <- stats::setNames(vector("list", length(ctx$names)), ctx$names)
  for (p in doc[["parameters"]]) {
    node <- abn_json_var_name(doc, p[["target"]])
    out[[node]][[length(out[[node]]) + 1L]] <- p
  }
  out
}

# ---------------------------------------------------------------------------
# Native names: extension when abn-authored, derived otherwise
# ---------------------------------------------------------------------------

abn_json_native_names <- function(doc, ctx, is_abn) {
  ext <- NULL
  if (is_abn) {
    names_map <- doc[["metadata"]][["extensions"]][["abn"]][["parameter_names"]]
    if (!is.null(names_map)) {
      ext <- stats::setNames(lapply(names_map, function(x) {
        list(field = x[["field"]], name = as.character(x[["name"]]))
      }), vapply(names_map, function(x) as.character(x[["parameter"]]), character(1)))
    }
  }
  out <- list()
  for (p in doc[["parameters"]]) {
    id <- as.character(p[["_id"]])
    out[[id]] <- if (!is.null(ext[[id]])) ext[[id]] else abn_json_derive_name(doc, p)
  }
  out
}

abn_json_derive_name <- function(doc, p) {
  node <- abn_json_var_name(doc, p[["target"]])
  kind <- p[["kind"]]
  bayes <- identical(doc[["inference"]][["type"]], "bayesian")
  grouped <- length(doc[["groups"]]) > 0
  field <- if (kind %in% c("intercept", "coefficient")) {
    if (bayes) "modes" else if (grouped) if (identical(kind, "intercept")) "mu" else "betas"
    else "coef"
  } else if (kind == "residual_variance") {
    if (bayes) "modes" else if (grouped) "sigma" else "mse"
  } else if (kind == "random_variance") {
    if (bayes) "modes" else "sigma_alpha"
  } else {
    "sigma_alpha"
  }
  name <- if (bayes) {
    term <- switch(kind, intercept = "(Intercept)",
                   residual_variance = "precision",
                   random_variance = "group.precision",
                   coefficient = abn_json_var_name(doc, p[["parent"]]),
                   stop("Cannot derive name for kind ", kind, ".", call. = FALSE))
    paste0(node, "|", term)
  } else if (kind %in% c("residual_variance", "random_variance")) {
    node
  } else if (kind == "random_covariance") {
    paste(vapply(p[["states"]], function(s) abn_json_state_label(doc, p[["target"]], s),
                 character(1)), collapse = ",")
  } else if (identical(kind, "intercept")) {
    if (grouped) {
      if (is.null(p[["target_state"]])) "(Intercept)" else
        abn_json_state_label(doc, p[["target"]], p[["target_state"]])
    } else if (!is.null(p[["target_state"]])) {
      paste0(node, "|intercept.", abn_json_state_label(doc, p[["target"]],
                                                       p[["target_state"]]))
    } else {
      paste0(node, "|intercept")
    }
  } else {
    parent <- abn_json_var_name(doc, p[["parent"]])
    grouped <- length(doc[["groups"]]) > 0
    if (!is.null(p[["parent_state"]])) {
      paste0(parent, abn_json_state_label(doc, p[["parent"]], p[["parent_state"]]))
    } else if (!is.null(p[["target_state"]]) && !grouped) {
      # ungrouped multinomial children: abn column names concatenate parent
      # and child state; grouped ones keep the plain term name per row
      paste0(parent, abn_json_state_label(doc, p[["target"]], p[["target_state"]]))
    } else {
      parent
    }
  }
  list(field = field, name = name)
}

# ---------------------------------------------------------------------------
# Diagnostics shared by all fit types
# ---------------------------------------------------------------------------

abn_json_diagnostics <- function(doc, ctx) {
  diag <- doc[["inference"]][["diagnostics"]] %||% list()
  nodes <- diag[["nodes"]] %||% list()
  by_var <- stats::setNames(nodes, vapply(nodes, function(n) {
    abn_json_var_name(doc, n[["variable"]])
  }, character(1)))
  get_num <- function(field, node) {
    v <- by_var[[node]][[field]]
    if (is.null(v)) NA else as.numeric(v)
  }
  # lapply + unlist keeps the native type: all-missing fields stay logical NA
  # (matching e.g. grouped fits where df is a logical NA vector)
  out <- list(
    mliknode = stats::setNames(unlist(lapply(ctx$names, function(n) {
      get_num("log_marginal_likelihood", n)
    })), ctx$names),
    mlik = as.numeric(diag[["log_marginal_likelihood"]])
  )
  if (identical(doc[["inference"]][["type"]], "maximum_likelihood")) {
    # fit field -> diagnostics key; aicnode/bicnode/mdlnode drop the "node" suffix
    key_of <- c(aicnode = "aic", bicnode = "bic", mdlnode = "mdl",
                df = "df", sse = "sse", mse = "mse")
    for (field in names(key_of)) {
      out[[field]] <- stats::setNames(unlist(lapply(ctx$names, function(n) {
        get_num(key_of[[field]], n)
      })), ctx$names)
    }
    out$aic <- as.numeric(diag[["aic"]])
    out$bic <- as.numeric(diag[["bic"]])
    out$mdl <- out$mdlnode  # derived: identical to mdlnode
  }
  out
}

# ---------------------------------------------------------------------------
# Ungrouped MLE
# ---------------------------------------------------------------------------

abn_json_param_value <- function(p) {
  if (is.null(p[["value"]])) NA_real_ else as.numeric(p[["value"]])
}

abn_json_native <- function(natives, p) natives[[as.character(p[["_id"]])]][["name"]]

abn_json_reconstruct_mle <- function(doc, ctx, natives) {
  fit <- list(method = "mle")
  diag <- abn_json_diagnostics(doc, ctx)
  fit[names(diag)] <- diag
  params <- abn_json_params_by_node(doc, ctx)

  fixed_params <- function(node) Filter(function(p) {
    p[["kind"]] %in% c("intercept", "coefficient")
  }, params[[node]])

  fit[["coef"]] <- stats::setNames(lapply(ctx$names, function(node) {
    ps <- fixed_params(node)
    matrix(vapply(ps, abn_json_param_value, numeric(1)), nrow = 1,
           dimnames = list(NULL, vapply(ps, abn_json_native, character(1),
                                        natives = natives)))
  }), ctx$names)
  fit[["Stderror"]] <- stats::setNames(lapply(ctx$names, function(node) {
    ps <- fixed_params(node)
    vals <- vapply(ps, function(p) {
      se <- p[["uncertainty"]][["standard_error"]]
      if (is.null(se)) NA_real_ else as.numeric(se)
    }, numeric(1))
    matrix(vals, nrow = 1, dimnames = list(NULL, vapply(ps, abn_json_native,
                                                         character(1), natives = natives)))
  }), ctx$names)

  # gaussian residual variance from parameters, not from mse diagnostics
  for (node in ctx$names) {
    if (identical(ctx$dists_flat[[node]], "gaussian")) {
      rv <- Filter(function(p) identical(p[["kind"]], "residual_variance"), params[[node]])
      if (length(rv) == 1) fit$mse[[node]] <- as.numeric(rv[[1]][["value"]])
    }
  }

  abn_json_fit_common(fit, ctx)
}

# ---------------------------------------------------------------------------
# Grouped MLE
# ---------------------------------------------------------------------------

abn_json_reconstruct_mle_grouped <- function(doc, ctx, natives) {
  fit <- list(method = "mle")
  diag <- abn_json_diagnostics(doc, ctx)
  fit[names(diag)] <- diag
  params <- abn_json_params_by_node(doc, ctx)
  kind_params <- function(node, kind) Filter(function(p) identical(p[["kind"]], kind),
                                             params[[node]])

  # original mu: unnamed scalars for lme4/glmmTMB nodes, named (coefmat rows)
  # for multinomial nodes
  fit[["mu"]] <- stats::setNames(lapply(ctx$names, function(node) {
    ps <- kind_params(node, "intercept")
    vals <- vapply(ps, abn_json_param_value, numeric(1))
    if (identical(ctx$dists_flat[[node]], "multinomial")) {
      stats::setNames(vals, vapply(ps, abn_json_native, character(1), natives = natives))
    } else {
      unname(vals)
    }
  }), ctx$names)

  fit[["betas"]] <- stats::setNames(lapply(ctx$names, function(node) {
    ps <- kind_params(node, "coefficient")
    if (length(ps) == 0) {
      # original: numeric(0) with empty names attribute for slopeless lme4
      # nodes (fixef[-1] on a length-1 vector), NA for multinomial nodes with
      # a single design column
      return(if (identical(ctx$dists_flat[[node]], "multinomial")) {
        NA
      } else {
        stats::setNames(numeric(0), character(0))
      })
    }
    if (identical(ctx$dists_flat[[node]], "multinomial")) {
      states <- ctx$levels[[node]][-1]
      terms <- unique(vapply(ps, abn_json_native, character(1), natives = natives))
      m <- matrix(NA_real_, length(states), length(terms),
                  dimnames = list(states, terms))
      for (p in ps) {
        i <- match(abn_json_state_label(doc, p[["target"]], p[["target_state"]]), states)
        m[i, abn_json_native(natives, p)] <- abn_json_param_value(p)
      }
      m
    } else {
      stats::setNames(vapply(ps, abn_json_param_value, numeric(1)),
                      vapply(ps, abn_json_native, character(1), natives = natives))
    }
  }), ctx$names)

  # original shapes: multinomial nodes store sigma = NA (logical) and a
  # covariance matrix; glmer nodes store sigma = numeric(0); gaussian nodes
  # store the residual sd
  fit[["sigma"]] <- stats::setNames(lapply(ctx$names, function(node) {
    rv <- kind_params(node, "residual_variance")
    if (length(rv) == 0) {
      if (identical(ctx$dists_flat[[node]], "multinomial")) NA else numeric(0)
    } else {
      unname(sqrt(as.numeric(rv[[1]][["value"]])))
    }
  }), ctx$names)

  fit[["sigma_alpha"]] <- stats::setNames(lapply(ctx$names, function(node) {
    rv <- kind_params(node, "random_variance")
    cv <- kind_params(node, "random_covariance")
    if (length(rv) == 0) return(NA)
    if (length(rv) > 1) {
      # multinomial child: full covariance matrix, values not squared; the
      # original dimnames come from the native names of the rv parameters
      states <- ctx$levels[[node]][-1]
      dn <- vapply(rv, abn_json_native, character(1), natives = natives)
      m <- matrix(NA_real_, length(states), length(states), dimnames = list(dn, dn))
      for (p in rv) {
        i <- match(abn_json_state_label(doc, p[["target"]], p[["target_state"]]), states)
        m[i, i] <- as.numeric(p[["value"]])
      }
      for (p in cv) {
        ss <- vapply(p[["states"]], function(s) {
          abn_json_state_label(doc, p[["target"]], s)
        }, character(1))
        ii <- match(ss[1], states); jj <- match(ss[2], states)
        m[ii, jj] <- as.numeric(p[["value"]]); m[jj, ii] <- as.numeric(p[["value"]])
      }
      m
    } else {
      unname(sqrt(as.numeric(rv[[1]][["value"]])))
    }
  }), ctx$names)

  group <- doc[["groups"]][[1]]
  fit[["group.var"]] <- group[["name"]]
  member_names <- vapply(group[["variables"]], function(id) abn_json_var_name(doc, id),
                         character(1))
  fit[["grouped.vars"]] <- match(member_names, ctx$names)

  abn_json_fit_common(fit, ctx)
}

# ---------------------------------------------------------------------------
# Bayes
# ---------------------------------------------------------------------------

abn_json_reconstruct_bayes <- function(doc, ctx, natives) {
  fit <- list(method = "bayes")
  diag <- abn_json_diagnostics(doc, ctx)
  fit[["mliknode"]] <- diag[["mliknode"]]
  fit[["mlik"]] <- diag[["mlik"]]

  params <- abn_json_params_by_node(doc, ctx)
  fit[["modes"]] <- stats::setNames(lapply(ctx$names, function(node) {
    ps <- params[[node]]
    stats::setNames(vapply(ps, abn_json_param_value, numeric(1)),
                    vapply(ps, abn_json_native, character(1), natives = natives))
  }), ctx$names)

  fit[["coef"]] <- modes2coefs(fit[["modes"]])
  fit[["mse"]] <- getMSEfromModes(fit[["modes"]], ctx$dists)

  priors <- doc[["inference"]][["priors"]]
  if (!is.null(priors)) {
    fixed <- Filter(function(x) identical(x[["applies_to"]], "fixed_effects"), priors)[[1]]
    prec <- Filter(function(x) identical(x[["applies_to"]], "precisions"), priors)[[1]]
    fit[["priors"]] <- list(mean = as.numeric(fixed[["mean"]]),
                            prec = as.numeric(fixed[["precision"]]),
                            loggam.shape = as.numeric(prec[["shape"]]),
                            loggam.inv.scale = as.numeric(prec[["rate"]]))
  }

  marginals <- doc[["inference"]][["posterior"]][["marginals"]]
  if (!is.null(marginals)) {
    by_id <- stats::setNames(marginals, vapply(marginals, function(m) {
      as.character(m[["parameter"]])
    }, character(1)))
    mout <- stats::setNames(vector("list", length(ctx$names)), ctx$names)
    qout <- stats::setNames(vector("list", length(ctx$names)), ctx$names)
    any_q <- FALSE
    for (node in ctx$names) {
      inner_m <- list(); inner_q <- list()
      for (p in params[[node]]) {
        nm <- abn_json_native(natives, p)
        m <- by_id[[as.character(p[["_id"]])]]
        if (!is.null(m)) {
          inner_m[[nm]] <- matrix(c(as.numeric(m[["x"]]), as.numeric(m[["density"]])),
                                  ncol = 2, dimnames = list(NULL, c("x", "f(x)")))
        }
        pq <- p[["uncertainty"]][["posterior_quantiles"]]
        if (!is.null(pq)) {
          any_q <- TRUE
          inner_q[[nm]] <- matrix(
            c(vapply(pq, function(q) as.numeric(q[["probability"]]), numeric(1)),
              vapply(pq, function(q) as.numeric(q[["value"]]), numeric(1))),
            ncol = 2, dimnames = list(NULL, c("P(X<=x)", "x")))
        }
      }
      mout[[node]] <- inner_m
      qout[[node]] <- inner_q
    }
    fit[["marginals"]] <- mout
    if (any_q) fit[["marginal.quantiles"]] <- qout
  }

  if (length(doc[["groups"]]) > 0) {
    group <- doc[["groups"]][[1]]
    fit[["group.var"]] <- group[["name"]]
    member_names <- vapply(group[["variables"]], function(id) abn_json_var_name(doc, id),
                           character(1))
    fit[["grouped.vars"]] <- match(member_names, ctx$names)
  }

  nodes_ext <- doc[["metadata"]][["extensions"]][["abn"]][["nodes"]]
  if (!is.null(nodes_ext)) {
    by_var <- stats::setNames(nodes_ext, vapply(nodes_ext, function(n) {
      abn_json_var_name(doc, n[["variable"]])
    }, character(1)))
    get_flag <- function(field, node) {
      v <- by_var[[node]][[field]]
      if (is.null(v)) NA else v
    }
    # unlist keeps the native type: all-missing vectors stay logical NA
    fit[["used.INLA"]] <- stats::setNames(unlist(lapply(ctx$names, function(n) {
      as.logical(get_flag("used_inla", n))
    })), ctx$names)
    collect <- function(field, coerce) {
      raw <- lapply(ctx$names, function(n) get_flag(field, n))
      if (all(vapply(raw, function(x) length(x) == 1L && isTRUE(is.na(x)), logical(1)))) {
        unlist(raw)
      } else {
        coerce(unlist(raw))
      }
    }
    fit[["error.code"]] <- stats::setNames(collect("error_code", as.numeric), ctx$names)
    fit[["error.code.desc"]] <- stats::setNames(collect("error_code_desc", as.character),
                                                ctx$names)
    fit[["hessian.accuracy"]] <- stats::setNames(collect("hessian_accuracy", as.numeric),
                                                 ctx$names)
  }

  abn_json_fit_common(fit, ctx)
}

# ---------------------------------------------------------------------------
# Common pieces: abnDag, levels, centre
# ---------------------------------------------------------------------------

abn_json_fit_common <- function(fit, ctx) {
  fit[["abnDag"]] <- structure(list(dag = ctx$dag, data.df = abn_json_empty_data(ctx),
                                     data.dists = ctx$dists),
                               class = "abnDag")
  categorical <- ctx$names[ctx$dists_flat[ctx$names] %in% c("binomial", "multinomial")]
  if (length(categorical) > 0) fit[["levels"]] <- ctx$levels[categorical]
  if (!is.null(ctx$centre)) fit[["centre"]] <- ctx$centre
  fit
}

abn_json_empty_data <- function(ctx) {
  columns <- stats::setNames(lapply(ctx$names, function(node) {
    if (ctx$dists_flat[[node]] %in% c("binomial", "multinomial")) {
      factor(character(0), levels = ctx$levels[[node]])
    } else {
      numeric(0)
    }
  }), ctx$names)
  as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
}

# ---------------------------------------------------------------------------
# Data attachment
# ---------------------------------------------------------------------------

abn_json_attach_data <- function(fit, ctx, data) {
  if (is.null(data)) return(fit)
  if (is.character(data) && length(data) == 1L && file.exists(data)) {
    data <- import_abnData(file = data)$data.df
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or a path to a data document.", call. = FALSE)
  }

  missing_cols <- setdiff(ctx$names, colnames(data))
  if (length(missing_cols) > 0) {
    stop("data is missing the columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  for (node in ctx$names) {
    if (ctx$dists_flat[[node]] %in% c("binomial", "multinomial")) {
      if (!setequal(levels(factor(data[[node]])), ctx$levels[[node]])) {
        stop("Levels of '", node, "' in data do not match the document states.",
             call. = FALSE)
      }
    }
  }

  bayes <- identical(ctx$type, "bayesian")
  grouped <- !is.null(fit[["group.var"]])
  if (grouped) {
    if (!(fit[["group.var"]] %in% colnames(data))) {
      stop("data is missing the grouping column '", fit[["group.var"]], "'.",
           call. = FALSE)
    }
    fit[["group.ids"]] <- as.integer(data[[fit[["group.var"]]]])
  }
  data <- data[, ctx$names, drop = FALSE]
  # the data document stores row names as strings; compact default row names
  # are integers in native fits
  if (identical(row.names(data), as.character(seq_len(nrow(data))))) {
    row.names(data) <- NULL
  }

  # mirror the transformations of fitAbn.mle / fitAbn.bayes exactly
  if (!bayes) {
    for (node in names(fit[["centre"]])) {
      tr <- fit[["centre"]][[node]]
      data[[node]] <- (data[[node]] - tr[["center"]]) / tr[["scale"]]
    }
    for (node in ctx$names[ctx$dists_flat[ctx$names] == "binomial"]) {
      if (!inherits(data[[node]], "numeric")) {
        data[[node]] <- as.numeric(factor(data[[node]])) - 1
      }
    }
    for (node in ctx$names[ctx$dists_flat[ctx$names] == "multinomial"]) {
      data[[node]] <- factor(data[[node]])
    }
  } else {
    for (node in ctx$names[ctx$dists_flat[ctx$names] == "binomial"]) {
      data[[node]] <- as.numeric(factor(data[[node]])) - 1
    }
    for (node in names(fit[["centre"]])) {
      tr <- fit[["centre"]][[node]]
      data[[node]] <- (data[[node]] - tr[["center"]]) / tr[["scale"]]
    }
    for (node in colnames(data)) data[[node]] <- as.double(data[[node]])
  }

  fit[["abnDag"]][["data.df"]] <- data
  fit
}

# ---------------------------------------------------------------------------
# Metadata attribute for stable re-export
# ---------------------------------------------------------------------------

abn_json_import_attr <- function(doc, natives) {
  meta <- doc[["metadata"]]
  vars <- doc[["variables"]]
  vnames <- vapply(vars, function(v) v[["name"]], character(1))
  list(
    label = meta[["label"]],
    scenario_id = meta[["scenario_id"]],
    data_reference = meta[["data_reference"]],
    ids = list(
      variables = stats::setNames(vapply(vars, function(v) as.integer(v[["_id"]]),
                                         integer(1)), vnames),
      states = stats::setNames(lapply(vars, function(v) {
        s <- v[["states"]]
        if (is.null(s)) {
          NULL
        } else {
          stats::setNames(vapply(s, function(x) as.integer(x[["_id"]]), integer(1)),
                          vapply(s, function(x) as.character(x[["label"]]), character(1)))
        }
      }), vnames),
      groups = if (length(doc[["groups"]]) > 0) {
        stats::setNames(vapply(doc[["groups"]], function(g) as.integer(g[["_id"]]),
                               integer(1)),
                        vapply(doc[["groups"]], function(g) g[["name"]], character(1)))
      } else {
        NULL
      },
      parameters = lapply(doc[["parameters"]], function(p) {
        list(id = as.integer(p[["_id"]]),
             sig = paste(abn_json_var_name(doc, p[["target"]]), p[["kind"]],
                         natives[[as.character(p[["_id"]])]][["name"]], sep = "|"))
      })
    )
  )
}
