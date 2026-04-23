# Structural invariants of the JSON I/O layer.
#
# These tests pin the design contract introduced when the exporter/importer
# were rewritten to encode all node relations purely via structural fields
# (`source.variable_id`, `source.state_id`, `conditions[].parent_variable_id`,
# `conditions[].parent_state_id`). The `name` field on each parameter is
# treated as an opaque, human-readable label and MUST NOT be parsed by the
# importer.
#
# What we verify here:
#   1. No `|` characters leak into any parameter `name` (which would indicate
#      relational info was smuggled into the label string).
#   2. Allowed labels form a small, fixed set.
#   3. Every `parent_variable_id` referenced by a parameter resolves to a real
#      variable in the export.
#   4. Re-exporting an imported model produces byte-identical JSON
#      (idempotence / determinism).
#   5. Mutating every `name` field to a junk string and re-importing yields a
#      model identical to the unmutated import (name-irrelevance).

library(testthat)
library(abn)

# Local copy: `%||%` is internal to the package and not re-exported.
`%||%` <- function(a, b) if (!is.null(a)) a else b

ALLOWED_NAMES <- c(
  "intercept", "beta", "variance", "sigma",
  "sigma_alpha", "random_variance", "random_covariance"
)

# Small helper: parse JSON robustly with simplifyVector = FALSE so that
# parameter / variable lists stay as lists-of-lists.
.parse_json <- function(json_str) {
  jsonlite::fromJSON(json_str, simplifyVector = FALSE)
}

test_that("No '|' character leaks into parameter `name` (standard MLE)", {
  suppressMessages(suppressWarnings({
    fit <- create_test_abnfit_mle()
    parsed <- .parse_json(export_abnFit(fit))
    names_seen <- vapply(parsed$parameters,
                         function(p) as.character(p$name %||% ""),
                         character(1))
    bad <- names_seen[grepl("\\|", names_seen)]
    expect_length(bad, 0,
                  info = paste0("Found relational `|` in name(s): ",
                                paste(unique(bad), collapse = ", ")))
  }))
})

test_that("All parameter `name` values come from the allowed label set", {
  suppressMessages(suppressWarnings({
    fit <- create_test_abnfit_mle()
    parsed <- .parse_json(export_abnFit(fit))
    names_seen <- unique(vapply(parsed$parameters,
                                function(p) as.character(p$name %||% ""),
                                character(1)))
    extras <- setdiff(names_seen, ALLOWED_NAMES)
    expect_length(extras, 0,
                  info = paste0("Unexpected name labels: ",
                                paste(extras, collapse = ", ")))
  }))
})

test_that("Every parent_variable_id resolves to a declared variable", {
  suppressMessages(suppressWarnings({
    fit <- create_test_abnfit_mle()
    parsed <- .parse_json(export_abnFit(fit))
    var_ids <- vapply(parsed$variables,
                      function(v) as.character(v$variable_id),
                      character(1))
    for (p in parsed$parameters) {
      for (coeff in (p$coefficients %||% list())) {
        for (cond in (coeff$conditions %||% list())) {
          pid <- as.character(cond$parent_variable_id %||% NA)
          if (!is.na(pid)) {
            expect_true(pid %in% var_ids,
                        info = paste("Dangling parent_variable_id:", pid))
          }
        }
      }
    }
  }))
})

test_that("Re-export after import is byte-identical (determinism)", {
  suppressMessages(suppressWarnings({
    fit <- create_test_abnfit_mle()
    json1 <- export_abnFit(fit)
    reimported <- import_abnFit(json = json1)
    json2 <- export_abnFit(reimported)
    expect_identical(json1, json2,
                     info = "export(import(export(fit))) must equal export(fit)")
  }))
})

test_that("Mutating parameter `name` does not change the imported model", {
  # Importer must derive everything from structural fields. If we replace
  # every `name` with garbage, the reimported abnFit must still be equal to
  # the original import.
  suppressMessages(suppressWarnings({
    fit <- create_test_abnfit_mle()
    json_str <- export_abnFit(fit)
    base_import <- import_abnFit(json = json_str)

    parsed <- .parse_json(json_str)
    parsed$parameters <- lapply(parsed$parameters, function(p) {
      p$name <- "ZZZ_irrelevant_label"
      p
    })
    mutated_json <- jsonlite::toJSON(parsed, auto_unbox = TRUE,
                                     pretty = FALSE, digits = NA, null = "null")
    mutated_import <- import_abnFit(json = mutated_json)

    expect_true(abnfit_objects_equal(base_import, mutated_import),
                info = "Importer must ignore `name` and rely only on structural fields")
  }))
})
