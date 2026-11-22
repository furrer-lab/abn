test_that("one-column one-row data exports retain array shapes", {
  dat <- data.frame(C = factor("001", levels = "001"))
  json <- export_abnData(dat, list(C = "multinomial"))
  doc <- jdoc_parse(json)
  expect_type(doc$metadata$column_order, "list")
  expect_length(doc$metadata$column_order, 1L)
  expect_type(doc$metadata$columns, "list")
  expect_type(doc$metadata$columns[[1]]$categories, "list")
  expect_identical(doc$metadata$columns[[1]]$categories[[1]], "001")
  expect_type(doc$metadata$adapters$r$columns$C$class, "list")
  expect_type(doc$metadata$adapters$r$columns$C$levels, "list")
  expect_type(doc$metadata$adapters$r$row_names, "list")
  expect_type(doc$data$C, "list")
  expect_identical(import_abnData(json = json)$data.df, dat)
})

test_that("zero-row datasets preserve declared and unknown category domains", {
  dat <- data.frame(
    known = factor(character(0), levels = c("a", "b", "unused")),
    unknown = character(0),
    numeric = numeric(0)
  )
  dists <- list(known = "multinomial", unknown = "binomial", numeric = "gaussian")
  json <- export_abnData(dat, dists)
  doc <- jdoc_parse(json)
  expect_identical(doc$metadata$columns[[2]]$categories, list())
  expect_identical(doc$metadata$adapters$r$row_names, list())
  expect_equal(doc$metadata$summary$row_count, 0)
  expect_identical(import_abnData(json = json)$data.df, dat)
  doc$metadata$adapters <- NULL
  generic <- import_abnData(json = jdoc_serialize(doc))$data.df
  expect_equal(nrow(generic), 0L)
  expect_identical(levels(generic$known), c("a", "b", "unused"))
  expect_identical(levels(generic$unknown), character(0))
  expect_type(generic$numeric, "double")
})

test_that("all-missing categorical samples need not invent categories", {
  dat <- data.frame(binary = rep(NA_character_, 3),
                    categorical = rep(NA_character_, 3))
  json <- export_abnData(dat, list(binary = "binomial", categorical = "multinomial"))
  doc <- jdoc_parse(json)
  for (column in doc$metadata$columns) expect_identical(column$categories, list())
  expect_identical(import_abnData(json = json)$data.df, dat)
  doc$metadata$adapters <- NULL
  generic <- import_abnData(json = jdoc_serialize(doc))$data.df
  expect_true(all(is.na(generic$binary)))
  expect_true(all(is.na(generic$categorical)))
  expect_identical(levels(generic$binary), character(0))
  expect_identical(levels(generic$categorical), character(0))
})

test_that("single-category raw multinomial samples and unused factor levels round-trip", {
  dat <- data.frame(C = c("only", "only"),
                    known = factor(c("b", "b"), levels = c("a", "b", "c")))
  json <- export_abnData(dat, list(C = "multinomial", known = "multinomial"))
  expect_identical(import_abnData(json = json)$data.df, dat)
  doc <- jdoc_parse(json)
  doc$metadata$adapters <- NULL
  generic <- import_abnData(json = jdoc_serialize(doc))$data.df
  expect_identical(levels(generic$C), "only")
  expect_identical(levels(generic$known), c("a", "b", "c"))
})

test_that("observed labels must still belong to an explicitly declared domain", {
  doc <- jdoc_parse(export_abnData(data.frame(C = c("a", "a")),
                                  list(C = "multinomial")))
  doc$metadata$columns[[1]]$categories <- list()
  expect_error(import_abnData(json = jdoc_serialize(doc)), "undeclared category")
})

test_that("zero-column data are explicitly unsupported", {
  expect_error(export_abnData(data.frame(), list()), "at least one column")
})

test_that("empty and singleton data exports conform to the data schema", {
  samples <- list(data.frame(C = "only"), data.frame(C = character(0)),
                  data.frame(C = rep(NA_character_, 2)))
  for (dat in samples) {
    json <- export_abnData(dat, list(C = "multinomial"))
    jdoc_expect_valid(json, "bn-data.schema.json")
  }
})
