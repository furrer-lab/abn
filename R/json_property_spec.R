# Single source of truth mapping fitAbn fields to the bayesian-network JSON
# document. One row per (field, fit_type):
#   field    : top-level element of an abnFit object
#   fit_type : "mle", "mle_grouped" or "bayes" (see abn_json_fit_type())
#   location : JSON path of the one place the field is stored; NA if not stored
#   reason   : why a field is not stored (NA if stored):
#              "derived"  - recomputed on import from stored values
#              "data"     - per-observation information, restored from the data
#                           document via import_abnFit(data = ...)
#              "excluded" - intentionally not exported (e.g. the call)
# Not fitAbn fields and therefore not listed: pvalue (never set by fitAbn),
# multinomial.states (an artifact of the pre-release importer, not a fitAbn field).

abn_json_spec_rows <- function(fit_type, rows) {
  data.frame(
    field = names(rows),
    fit_type = fit_type,
    location = vapply(rows, function(x) if (x %in% c("derived", "data", "excluded"))
      NA_character_ else x, character(1)),
    reason = vapply(rows, function(x) if (x %in% c("derived", "data", "excluded"))
      x else NA_character_, character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# Fields shared by all fit types.
abn_json_spec_common <- function() {
  c(method = "inference.type",
    abnDag = "arcs",
    centre = "variables[].transform",
    levels = "variables[].states",
    mliknode = "inference.diagnostics.nodes[].log_marginal_likelihood",
    mlik = "inference.diagnostics.log_marginal_likelihood",
    call = "excluded")
}

# MLE score diagnostics (ungrouped and grouped).
abn_json_spec_mle_scores <- function() {
  c(aicnode = "inference.diagnostics.nodes[].aic",
    aic = "inference.diagnostics.aic",
    bicnode = "inference.diagnostics.nodes[].bic",
    bic = "inference.diagnostics.bic",
    mdlnode = "inference.diagnostics.nodes[].mdl",
    mdl = "derived",          # identical to mdlnode
    df = "inference.diagnostics.nodes[].df",
    sse = "inference.diagnostics.nodes[].sse",
    # gaussian nodes of ungrouped fits: residual_variance parameter
    mse = "inference.diagnostics.nodes[].mse")
}

abn_json_spec_grouping <- function() {
  c(group.var = "groups[].name",
    grouped.vars = "groups[].variables",
    group.ids = "data")
}

abn_json_property_spec <- function() {
  rbind(
    abn_json_spec_rows("mle", c(
      abn_json_spec_common(),
      abn_json_spec_mle_scores(),
      coef = "parameters[].value",
      Stderror = "parameters[].uncertainty.standard_error"
    )),
    abn_json_spec_rows("mle_grouped", c(
      abn_json_spec_common(),
      abn_json_spec_mle_scores(),
      abn_json_spec_grouping(),
      mu = "parameters[kind=intercept].value",
      betas = "parameters[kind=coefficient].value",
      sigma = "parameters[kind=residual_variance].value",
      sigma_alpha = "parameters[kind=random_variance|random_covariance].value"
    )),
    abn_json_spec_rows("bayes", c(
      abn_json_spec_common(),
      abn_json_spec_grouping(),
      modes = "parameters[].value",
      coef = "derived",       # modes2coefs(modes)
      mse = "derived",        # 1 / residual precision
      marginal.quantiles = "parameters[].uncertainty.posterior_quantiles",
      marginals = "inference.posterior.marginals",
      priors = "inference.priors",
      used.INLA = "metadata.extensions.abn.nodes[].used_inla",
      error.code = "metadata.extensions.abn.nodes[].error_code",
      error.code.desc = "metadata.extensions.abn.nodes[].error_code_desc",
      hessian.accuracy = "metadata.extensions.abn.nodes[].hessian_accuracy"
    ))
  )
}

# Fit type of an abnFit object as used in abn_json_property_spec().
abn_json_fit_type <- function(fit) {
  if (identical(fit$method, "bayes")) return("bayes")
  if (!is.null(fit$group.var)) return("mle_grouped")
  "mle"
}

# Fields of a fit type; stored = TRUE for fields with a JSON location,
# FALSE for derived/data/excluded fields.
abn_json_fields <- function(fit_type, stored = TRUE) {
  spec <- abn_json_property_spec()
  spec <- spec[spec$fit_type == fit_type, , drop = FALSE]
  spec$field[if (stored) !is.na(spec$location) else is.na(spec$location)]
}

abn_json_check_property_spec <- function(spec = abn_json_property_spec()) {
  if (anyDuplicated(paste(spec$field, spec$fit_type))) {
    stop("JSON property spec: duplicated (field, fit_type).", call. = FALSE)
  }
  stored <- !is.na(spec$location)
  if (any(!is.na(spec$reason[stored])) ||
      any(!spec$reason[!stored] %in% c("derived", "data", "excluded"))) {
    stop("JSON property spec: each row needs either a location or a reason.",
         call. = FALSE)
  }
  for (type in unique(spec$fit_type)) {
    if (anyDuplicated(spec$location[stored & spec$fit_type == type])) {
      stop("JSON property spec: duplicated location for fit type ", type, ".",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

abn_json_check_property_spec()
