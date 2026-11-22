row_label_data <- function() {
  data.frame(x = c(0.5, 1.5, 2.5, 3.5), y = c(4, 3, 2, 1))
}

row_label_export <- function(dat) {
  export_abnData(dat, data.dists = list(x = "gaussian", y = "gaussian"))
}

test_that("R adapter preserves string row labels without guessing numeric types", {
  for (labels in list(c("001", "1", "999999999999999999999999", "sample-a"),
                      c("001", "002", "003", "004"),
                      c("1", "2", "3", "4"))) {
    dat <- row_label_data()
    row.names(dat) <- labels
    json <- row_label_export(dat)
    doc <- jsonlite::fromJSON(json, simplifyVector = FALSE)
    expect_identical(doc$metadata$adapters$r$row_names_type, "character")
    expect_true(all(vapply(doc$metadata$adapters$r$row_names, is.character, logical(1))))
    back <- import_abnData(json = json)$data.df
    expect_identical(row.names(back), labels)
    expect_identical(typeof(attr(back, "row.names")), "character")
    expect_identical(back, dat)
  }
})

test_that("explicit R adapter metadata restores automatic and integer row names", {
  automatic <- row_label_data()
  integer <- row_label_data()
  row.names(integer) <- c(11L, 12L, 13L, 14L)
  for (type in c("automatic", "integer")) {
    dat <- if (type == "automatic") automatic else integer
    json <- row_label_export(dat)
    doc <- jsonlite::fromJSON(json, simplifyVector = FALSE)
    expect_identical(doc$metadata$adapters$r$row_names_type, type)
    back <- import_abnData(json = json)$data.df
    expect_identical(back, dat)
    expect_identical(typeof(attr(back, "row.names")), "integer")
    expect_identical(.row_names_info(back, type = 1L), .row_names_info(dat, type = 1L))
  }
})

test_that("legacy row-label metadata without a type marker defaults to strings", {
  dat <- row_label_data()
  # Include a canonical sequence as well as leading zeros: neither is evidence
  # of an integer type when the adapter has not explicitly supplied that type.
  for (labels in list(c("001", "002", "003", "004"), c("1", "2", "3", "4"))) {
    row.names(dat) <- labels
    doc <- jsonlite::fromJSON(row_label_export(dat), simplifyVector = FALSE)
    doc$metadata$adapters$r$row_names_type <- NULL
    json <- jsonlite::toJSON(doc, auto_unbox = TRUE, null = "null", digits = NA)
    back <- import_abnData(json = json)$data.df
    expect_identical(row.names(back), labels)
    expect_identical(typeof(attr(back, "row.names")), "character")
  }
})

test_that("without the R adapter custom row labels are absent and defaults are generated", {
  dat <- row_label_data()
  row.names(dat) <- c("001", "1", "999999999999999999999999", "sample-a")
  for (remove_all_adapters in c(FALSE, TRUE)) {
    doc <- jsonlite::fromJSON(row_label_export(dat), simplifyVector = FALSE)
    if (remove_all_adapters) {
      doc$metadata$adapters <- NULL
    } else {
      doc$metadata$adapters$r <- NULL
    }
    json <- jsonlite::toJSON(doc, auto_unbox = TRUE, null = "null", digits = NA)
    back <- import_abnData(json = json)$data.df
    expect_identical(back$x, dat$x)
    expect_identical(back$y, dat$y)
    expect_identical(row.names(back), as.character(seq_len(nrow(dat))))
    expect_true(.row_names_info(back, type = 1L) < 0L)
  }
})

test_that("inconsistent R row-name types cannot silently change string labels", {
  dat <- row_label_data()
  row.names(dat) <- c("001", "002", "003", "004")
  doc <- jsonlite::fromJSON(row_label_export(dat), simplifyVector = FALSE)
  doc$metadata$adapters$r$row_names_type <- "integer"
  json <- jsonlite::toJSON(doc, auto_unbox = TRUE, null = "null", digits = NA)
  expect_error(import_abnData(json = json), "exact, representable integer labels")

  doc$metadata$adapters$r$row_names_type <- "automatic"
  json <- jsonlite::toJSON(doc, auto_unbox = TRUE, null = "null", digits = NA)
  expect_error(import_abnData(json = json), "default row sequence")
})
