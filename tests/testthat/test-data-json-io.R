data_json_io_specs <- list(
  list(
    dataset = "g2b2c_data",
    rows = 25,
    data.dists = list(
      G1 = "gaussian", B1 = "binomial", B2 = "binomial",
      C = "multinomial", G2 = "gaussian"
    ),
    group.var = NULL,
    expected_semantic_types = c(
      G1 = "continuous", B1 = "binary", B2 = "binary",
      C = "categorical", G2 = "continuous"
    ),
    expected_value_encodings = c(
      G1 = "number", B1 = "labels", B2 = "labels",
      C = "labels", G2 = "number"
    )
  ),
  list(
    dataset = "g2pbcgrp",
    rows = 25,
    data.dists = list(
      G1 = "gaussian", P = "poisson", B = "binomial",
      C = "multinomial", G2 = "gaussian"
    ),
    group.var = "group",
    expected_semantic_types = c(
      G1 = "continuous", P = "count", B = "binary",
      C = "categorical", G2 = "continuous", group = "categorical"
    ),
    expected_value_encodings = c(
      G1 = "number", P = "integer", B = "labels",
      C = "labels", G2 = "number", group = "labels"
    )
  ),
  list(
    dataset = "ex1.dag.data",
    rows = 25,
    data.dists = list(
      b1 = "binomial", p1 = "poisson", g1 = "gaussian",
      b2 = "binomial", p2 = "poisson", b3 = "binomial",
      g2 = "gaussian", b4 = "binomial", b5 = "binomial",
      g3 = "gaussian"
    ),
    group.var = NULL,
    expected_semantic_types = c(
      b1 = "binary", p1 = "count", g1 = "continuous",
      b2 = "binary", p2 = "count", b3 = "binary",
      g2 = "continuous", b4 = "binary", b5 = "binary",
      g3 = "continuous"
    ),
    expected_value_encodings = c(
      b1 = "labels", p1 = "integer", g1 = "number",
      b2 = "labels", p2 = "integer", b3 = "labels",
      g2 = "number", b4 = "labels", b5 = "labels",
      g3 = "number"
    )
  ),
  list(
    dataset = "ex3.dag.data",
    rows = 25,
    data.dists = as.list(stats::setNames(
      rep("binomial", 13),
      paste0("b", 1:13)
    )),
    group.var = "group",
    expected_semantic_types = c(
      stats::setNames(rep("binary", 13), paste0("b", 1:13)),
      group = "categorical"
    ),
    expected_value_encodings = c(
      stats::setNames(rep("labels", 13), paste0("b", 1:13)),
      group = "labels"
    )
  ),
  list(
    dataset = "adg",
    rows = 25,
    data.dists = list(
      AR = "binomial", pneumS = "binomial", female = "binomial",
      livdam = "binomial", eggs = "binomial", wormCount = "poisson",
      age = "gaussian", adg = "gaussian"
    ),
    group.var = "farm",
    expected_semantic_types = c(
      AR = "binary", pneumS = "binary", female = "binary",
      livdam = "binary", eggs = "binary", wormCount = "count",
      age = "continuous", adg = "continuous", farm = "categorical"
    ),
    expected_value_encodings = c(
      AR = "labels", pneumS = "labels", female = "labels",
      livdam = "labels", eggs = "labels", wormCount = "integer",
      age = "number", adg = "number", farm = "labels"
    )
  ),
  list(
    dataset = "FCV",
    rows = 25,
    data.dists = list(
      FCV = "binomial", FHV_1 = "binomial",
      C_felis = "binomial", M_felis = "binomial",
      B_bronchiseptica = "binomial", FeLV = "binomial",
      FIV = "binomial", Gingivostomatitis = "binomial",
      URTD = "binomial", Vaccinated = "binomial",
      Pedigree = "binomial", Outdoor = "binomial",
      Sex = "multinomial", GroupSize = "poisson",
      Age = "gaussian"
    ),
    group.var = NULL,
    expected_semantic_types = c(
      FCV = "binary", FHV_1 = "binary", C_felis = "binary",
      M_felis = "binary", B_bronchiseptica = "binary",
      FeLV = "binary", FIV = "binary", Gingivostomatitis = "binary",
      URTD = "binary", Vaccinated = "binary", Pedigree = "binary",
      Outdoor = "binary", Sex = "categorical", GroupSize = "count",
      Age = "continuous"
    ),
    expected_value_encodings = c(
      FCV = "labels", FHV_1 = "labels", C_felis = "labels",
      M_felis = "labels", B_bronchiseptica = "labels",
      FeLV = "labels", FIV = "labels", Gingivostomatitis = "labels",
      URTD = "labels", Vaccinated = "labels", Pedigree = "labels",
      Outdoor = "labels", Sex = "labels", GroupSize = "integer",
      Age = "number"
    )
  )
)

load_data_json_fixture <- function(spec) {
  fixture_env <- new.env(parent = emptyenv())
  suppressWarnings(data(list = spec$dataset, package = "abn", envir = fixture_env))
  if (!exists(spec$dataset, envir = fixture_env, inherits = FALSE)) {
    data_file <- data_json_fixture_file(spec$dataset)
    if (file.exists(data_file)) load(data_file, envir = fixture_env)
  }
  expect_true(exists(spec$dataset, envir = fixture_env, inherits = FALSE),
              info = paste("attached dataset not found:", spec$dataset))
  data.df <- get(spec$dataset, envir = fixture_env, inherits = FALSE)
  data.df[seq_len(min(spec$rows, nrow(data.df))), , drop = FALSE]
}

data_json_fixture_file <- function(dataset) {
  dataset_file <- switch(
    dataset,
    "ex1.dag.data" = "ex1data.RData",
    "ex3.dag.data" = "ex3data.RData",
    paste0(dataset, ".RData")
  )
  file.path(testthat::test_path("..", "..", "data"), dataset_file)
}

data_json_column_metadata <- function(parsed, column_name) {
  columns <- parsed$metadata$columns
  matches <- vapply(columns, function(column) {
    identical(as.character(column$name), column_name)
  }, logical(1))
  expect_equal(sum(matches), 1L, info = paste("column metadata for", column_name))
  columns[[which(matches)]]
}

expect_data_json_metadata <- function(parsed, data.df, spec) {
  expect_true("metadata" %in% names(parsed))
  expect_true("data" %in% names(parsed))
  expect_equal(parsed$metadata$schema_version, "bn-data-v1")
  expect_equal(parsed$metadata$issuer, "abn::export_abnData")
  expect_equal(parsed$metadata$orientation, "columnar")
  expect_equal(unlist(parsed$metadata$column_order), colnames(data.df))
  expect_equal(names(parsed$data), colnames(data.df))
  expect_equal(length(parsed$metadata$columns), ncol(data.df))

  expected_node_names <- names(spec$data.dists)
  expected_group_name <- spec$group.var
  for (column_name in colnames(data.df)) {
    column <- data_json_column_metadata(parsed, column_name)
    expect_equal(column$name, column_name)
    expect_equal(
      column$semantic_type,
      unname(spec$expected_semantic_types[[column_name]])
    )
    expect_equal(
      column$value_encoding,
      unname(spec$expected_value_encodings[[column_name]])
    )
    expect_true("nullable" %in% names(column))

    if (identical(column_name, expected_group_name)) {
      expect_equal(column$role, "group")
      expect_null(column$model_distribution)
    } else {
      expect_true(column_name %in% expected_node_names)
      expect_equal(column$role, "node")
      expect_equal(column$model_distribution, unname(spec$data.dists[[column_name]]))
    }

    if (is.factor(data.df[[column_name]])) {
      expect_equal(unlist(column$categories), levels(data.df[[column_name]]))
      expect_equal(column$ordered, is.ordered(data.df[[column_name]]))
    }
  }
}

expect_data_json_values <- function(parsed, data.df) {
  for (column_name in colnames(data.df)) {
    values <- parsed$data[[column_name]]
    expect_equal(length(values), nrow(data.df), info = column_name)
    expected_values <- data.df[[column_name]]
    if (is.factor(expected_values)) expected_values <- as.character(expected_values)
    expected_values <- as.list(expected_values)
    expect_equal(values, expected_values, info = column_name)
  }
}

expect_imported_abn_data <- function(imported, data.df, spec) {
  expect_true(is.list(imported))
  expect_true("data.df" %in% names(imported))
  expect_true("data.dists" %in% names(imported))
  expect_true("group.var" %in% names(imported))
  expect_s3_class(imported$data.df, "data.frame")
  expect_equal(colnames(imported$data.df), colnames(data.df))
  expect_equal(nrow(imported$data.df), nrow(data.df))
  expect_equal(imported$data.dists, spec$data.dists)
  expect_equal(imported$group.var, spec$group.var)

  for (column_name in colnames(data.df)) {
    original <- data.df[[column_name]]
    roundtrip <- imported$data.df[[column_name]]
    expect_equal(class(roundtrip), class(original), info = column_name)
    if (is.factor(original)) {
      expect_equal(levels(roundtrip), levels(original), info = column_name)
      expect_equal(is.ordered(roundtrip), is.ordered(original), info = column_name)
    }
    expect_equal(roundtrip, original, info = column_name)
  }
}

test_that("export_abnData serializes exemplary attached datasets", {
  for (spec in data_json_io_specs) {
    data.df <- load_data_json_fixture(spec)

    json <- export_abnData(
      data.df = data.df,
      data.dists = spec$data.dists,
      group.var = spec$group.var,
      include_summary = TRUE
    )

    expect_type(json, "character")
    expect_length(json, 1)
    expect_true(jsonlite::validate(json), info = spec$dataset)

    parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
    expect_data_json_metadata(parsed, data.df, spec)
    expect_data_json_values(parsed, data.df)

    summary_result <- validate_data_json_summary(parsed, on_error = "return")
    expect_true(summary_result$valid, info = spec$dataset)
    expect_equal(parsed$metadata$summary, summary_result$summary_expected,
                 info = spec$dataset)
  }
})

test_that("import_abnData round-trips exemplary attached datasets", {
  for (spec in data_json_io_specs) {
    data.df <- load_data_json_fixture(spec)
    json <- export_abnData(
      data.df = data.df,
      data.dists = spec$data.dists,
      group.var = spec$group.var,
      include_summary = TRUE
    )

    imported <- import_abnData(json = json)

    expect_imported_abn_data(imported, data.df, spec)
  }
})

test_that("export_abnData writes JSON files that import_abnData can read", {
  spec <- data_json_io_specs[[1]]
  data.df <- load_data_json_fixture(spec)
  temp_file <- tempfile(fileext = ".json")
  on.exit(unlink(temp_file), add = TRUE)

  result <- export_abnData(
    data.df = data.df,
    data.dists = spec$data.dists,
    file = temp_file,
    include_summary = TRUE
  )

  expect_equal(result, temp_file)
  expect_true(file.exists(temp_file))
  imported <- import_abnData(file = temp_file)
  expect_imported_abn_data(imported, data.df, spec)
})
