# Adversarial documents and portable R-adapter-free reconstruction.

test_that("network validation rejects inconsistent variable compatibility", {
  doc <- jdoc_foreign_document()
  doc$variables[[1]]$link <- "logit"
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "incompatible")

  doc <- jdoc_foreign_document()
  doc$variables[[2]]$states[[1]]$baseline <- FALSE
  doc$variables[[2]]$states[[2]]$baseline <- FALSE
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "baseline")
})

test_that("network validation rejects duplicate names, arcs and parameter semantics", {
  doc <- jdoc_foreign_document()
  doc$variables[[2]]$name <- "height"
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "unique")

  doc <- jdoc_foreign_document()
  doc$arcs <- c(doc$arcs, list(list(source = 1, target = 2)))
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "Duplicate arc")

  doc <- jdoc_foreign_document()
  doc$parameters <- c(doc$parameters,
                      list(list(`_id` = 5, target = 2, kind = "intercept", value = 0)))
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "Duplicated parameter")
})

test_that("network validation rejects self-links and misplaced state references", {
  doc <- jdoc_foreign_document()
  doc$arcs <- list(list(source = 1, target = 1))
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "self-link")

  doc <- jdoc_foreign_document()
  baseline <- as.integer(jdoc_state_id(doc$variables[[2]], "no"))
  doc$parameters[[4]]$target_state <- baseline
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "baseline")
})

test_that("network validation rejects unresolved extension and posterior references", {
  doc <- jdoc_foreign_document()
  doc$metadata$extensions <- list(abn = list(parameter_names = list(
    list(parameter = 99, field = "coef", name = "x")
  )))
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "abn extension")

  doc <- jdoc_foreign_document()
  doc$inference$posterior <- list(marginals = list(
    list(parameter = 99, x = list(0, 1), density = list(0.5, 0.5))
  ))
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "posterior")
})

test_that("network validation rejects reference conflicts in diagnostics", {
  doc <- jdoc_foreign_document()
  doc$inference$diagnostics$nodes <- list(
    list(variable = 1, log_marginal_likelihood = 1),
    list(variable = 1, log_marginal_likelihood = 2)
  )
  expect_error(import_abnFit(json = jdoc_serialize(doc)), "Duplicated diagnostics")
})

test_that("special numeric values round-trip through bn-data", {
  dat <- data.frame(g1 = c(1.5, NaN, Inf, -Inf),
                    b1 = factor(c("no", "yes", "no", "yes"), levels = c("no", "yes")))
  json <- export_abnData(dat, data.dists = list(g1 = "gaussian", b1 = "binomial"))
  doc <- jdoc_parse(json)
  expect_equal(doc$data$g1[[2]], "NaN")
  expect_equal(doc$data$g1[[3]], "Infinity")
  expect_equal(doc$data$g1[[4]], "-Infinity")

  back <- import_abnData(json = json)$data.df
  expect_true(is.nan(back$g1[2]))
  expect_identical(back$g1[3:4], c(Inf, -Inf))
})

test_that("ordered and unused categories survive without R adapter metadata", {
  dat <- data.frame(g1 = c("b", NA, "a"),
                    b1 = factor(c("no", "yes", NA), levels = c("no", "maybe", "yes")))
  dat$g1 <- ordered(dat$g1, levels = c("a", "b", "c"))
  json <- export_abnData(dat, data.dists = list(g1 = "multinomial", b1 = "multinomial"))
  back <- import_abnData(json = json)
  generic <- jdoc_parse(json)
  generic$metadata$adapters <- NULL
  generic_json <- jdoc_serialize(generic)

  imported <- import_abnData(json = generic_json)$data.df
  expect_identical(levels(imported$g1), c("a", "b", "c"))
  expect_true(is.ordered(imported$g1))
  expect_identical(levels(imported$b1), c("no", "maybe", "yes"))
  expect_equal(back$data.df, dat)
})

test_that("double-backed count columns retain numeric type", {
  dat <- data.frame(p1 = c(0, 2, NA),
                    b1 = factor(c("no", "yes", "no"), levels = c("no", "yes")))
  json <- export_abnData(dat, data.dists = list(p1 = "poisson", b1 = "binomial"))
  imported <- import_abnData(json = json)$data.df
  expect_identical(class(imported$p1), class(dat$p1))
  expect_identical(imported$p1, dat$p1)
})

test_that("integer count columns retain integer type without adapters", {
  dat <- data.frame(p1 = c(0L, 2L, NA),
                    b1 = factor(c("no", "yes", "no"), levels = c("no", "yes")))
  json <- export_abnData(dat, data.dists = list(p1 = "poisson", b1 = "binomial"))
  generic <- jdoc_parse(json)
  generic$metadata$adapters <- NULL
  imported <- import_abnData(json = jdoc_serialize(generic))$data.df
  expect_type(imported$p1, "integer")
})

test_that("data validation rejects malformed encodings and undeclared labels", {
  dat <- data.frame(p1 = c(1, 2),
                    b1 = factor(c("no", "yes"), levels = c("no", "yes")))
  json <- export_abnData(dat, data.dists = list(p1 = "poisson", b1 = "binomial"))

  bad_count <- jdoc_parse(json)
  bad_count$data$p1[[2]] <- 1.5
  expect_error(import_abnData(json = jdoc_serialize(bad_count)), "non-negative integer")

  bad_label <- jdoc_parse(json)
  bad_label$data$b1[[1]] <- "other"
  expect_error(import_abnData(json = jdoc_serialize(bad_label)), "undeclared category")

  bad_categories <- jdoc_parse(json)
  bad_categories$metadata$columns[[2]]$categories <- c("no", "no")
  expect_error(import_abnData(json = jdoc_serialize(bad_categories)), "unique")

  bad_nullable <- jdoc_parse(json)
  bad_nullable$data$p1[1] <- list(NULL)
  expect_error(import_abnData(json = jdoc_serialize(bad_nullable)), "nullable")
})

test_that("primitive binary columns with one observed category round-trip", {
  dat <- data.frame(b1 = c(0L, 0L, 0L),
                    g1 = c(1.5, 2.5, 3.5))
  json <- export_abnData(dat, data.dists = list(b1 = "binomial", g1 = "gaussian"))
  back <- import_abnData(json = json)$data.df
  expect_identical(as.integer(back$b1), dat$b1)
  expect_type(back$b1, "integer")

  dat_chr <- data.frame(b1 = c("yes", "yes"), g1 = c(1.5, 2.5))
  back_chr <- import_abnData(json = export_abnData(
    dat_chr, data.dists = list(b1 = "binomial", g1 = "gaussian")))$data.df
  expect_identical(back_chr$b1, dat_chr$b1)
})

test_that("declared factor levels survive despite one observed level in raw subsets", {
  dat <- data.frame(b1 = factor(c("yes", "yes"), levels = c("no", "yes")),
                    g1 = c(1.5, 2.5))
  back <- import_abnData(json = export_abnData(
    dat, data.dists = list(b1 = "binomial", g1 = "gaussian")))$data.df
  expect_identical(levels(back$b1), c("no", "yes"))
  expect_identical(class(back$b1), class(dat$b1))
})

test_that("reserved special tokens cannot be used as label categories", {
  dat <- data.frame(g1 = factor(c("a", "NaN"), levels = c("a", "NaN")),
                    b1 = factor(c("no", "yes"), levels = c("no", "yes")))
  expect_error(export_abnData(dat, data.dists = list(g1 = "multinomial", b1 = "binomial")),
               "special numeric tokens")
})

test_that("descriptor inconsistency error identifies the mismatched column", {
  dat <- data.frame(p1 = c(1L, 2L),
                    b1 = factor(c("no", "yes"), levels = c("no", "yes")))
  json <- export_abnData(dat, data.dists = list(p1 = "poisson", b1 = "binomial"))
  bad <- jdoc_parse(json)
  bad$metadata$data_dists$p1 <- "gaussian"
  expect_error(import_abnData(json = jdoc_serialize(bad)), "p1")
})
