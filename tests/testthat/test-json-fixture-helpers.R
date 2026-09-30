test_that("JSON fixture helper compares unnamed grouped scalars positionally", {
  expect_no_error(
    json_fixture_compare_named_list(list(a = 1), list(a = 1))
  )
})

test_that("JSON fixture helper compares named vectors by name", {
  expect_no_error(
    json_fixture_compare_named_list(
      list(a = c(x = 1, y = 2)),
      list(a = c(y = 2, x = 1))
    )
  )
})

test_that("JSON fixture helper rejects name loss for multi-value vectors", {
  expect_error(
    json_fixture_compare_named_list(
      list(a = c(x = 1, y = 2)),
      list(a = c(1, 2))
    ),
    "lost names"
  )
})

test_that("JSON fixture helper accepts Bayesian empty SE imported as NA", {
  original <- matrix(numeric(0), nrow = 0, ncol = 0)
  imported <- matrix(NA_real_, nrow = 1, ncol = 1,
                     dimnames = list(NULL, "b1|intercept"))

  expect_no_error(
    json_fixture_compare_matrix(original, imported, allow_imported_all_na = TRUE)
  )
})

test_that("JSON fixture helper treats empty and NA grouped fields as absent", {
  expect_no_error(
    json_fixture_compare_named_list(list(a = numeric(0)), list(a = NA_real_))
  )
})
