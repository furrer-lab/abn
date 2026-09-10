test_that("export_abnFit works for MLE method without group.var", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      result <- export_abnFit(abn_fit)

      # ASSERT
      expect_type(result, "character")
      expect_true(jsonlite::validate(result))
    })
  })
})

test_that("export_abnFit writes to file when 'file' argument is provided", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()
      temp_file <- tempfile(fileext = ".json")
      on.exit(unlink(temp_file), add = TRUE)

      # ACT
      result <- export_abnFit(abn_fit, file = temp_file)

      # ASSERT
      expect_equal(result, temp_file) # Function should return file path invisibly
      expect_true(file.exists(temp_file))
      expect_true(file.info(temp_file)$size > 0)

      # Verify file contains valid JSON
      file_content <- readLines(temp_file, warn = FALSE)
      file_json <- paste(file_content, collapse = "\n")
      expect_true(jsonlite::validate(file_json))
    })
  })
})

test_that("exported JSON has correct top-level structure", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      json_str <- export_abnFit(abn_fit)
      parsed <- jsonlite::fromJSON(json_str)

      # ASSERT
      expect_true("variables" %in% names(parsed))
      expect_true("parameters" %in% names(parsed))
      expect_true("arcs" %in% names(parsed))
      expect_null(parsed$scenario_id) # Default should be NULL
      expect_null(parsed$label) # Default should be NULL
    })
  })
})

test_that("scenario_id and label can be set", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      json_str <- export_abnFit(abn_fit, scenario_id = "test_scenario_1", label = "Test Model")
      parsed <- jsonlite::fromJSON(json_str)

      # ASSERT
      expect_equal(parsed$scenario_id, "test_scenario_1")
      expect_equal(parsed$label, "Test Model")
    })
  })
})

test_that("variables array contains all expected nodes with correct structure", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      expect_true(is.data.frame(parsed$variables) || is.list(parsed$variables))
      expect_equal(nrow(parsed$variables), 8)  # 8 nodes

      # Check required fields in variables
      expect_true("variable_id" %in% names(parsed$variables))
      expect_true("attribute_name" %in% names(parsed$variables))
      expect_true("model_type" %in% names(parsed$variables))

      # Check all node IDs are present
      expected_nodes <- c("b1", "b2", "g1", "g2", "m1", "m2", "p1", "p2")
      expect_true(all(expected_nodes %in% parsed$variables$variable_id))

      # Check model types are valid
      valid_types <- c("gaussian", "binomial", "poisson", "multinomial")
      expect_true(all(parsed$variables$model_type %in% valid_types))
    })
  })
})

test_that("multinomial variables have states array", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      # Find multinomial variables
      multinomial_vars <- parsed$variables[parsed$variables$model_type == "multinomial", ]

      # Each multinomial variable should have states
      for (i in seq_len(nrow(multinomial_vars))) {
        states <- multinomial_vars$states[[i]]
        expect_true(!is.null(states))
        expect_true(is.data.frame(states) || is.list(states))
        expect_true("state_id" %in% names(states))
        expect_true("value_name" %in% names(states))
        expect_true("is_baseline" %in% names(states))
      }
    })
  })
})

test_that("parameters array has correct structure", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      expect_true(is.data.frame(parsed$parameters) || is.list(parsed$parameters))
      expect_true(nrow(parsed$parameters) > 0 || length(parsed$parameters) > 0)

      # Check required fields in parameters
      expect_true("parameter_id" %in% names(parsed$parameters))
      expect_true("name" %in% names(parsed$parameters))
      expect_true("link_function_name" %in% names(parsed$parameters))
      expect_true("source" %in% names(parsed$parameters))
      expect_true("coefficients" %in% names(parsed$parameters))
    })
  })
})

test_that("parameter source has correct structure", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit), simplifyVector = FALSE)

      # ASSERT
      # Check first parameter's source
      first_param <- parsed$parameters[[1]]
      source <- first_param$source[1]

      expect_true("variable_id" %in% names(source))
      # state_id is optional, so we don't require it
    })
  })
})

test_that("parameter coefficients have correct structure", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      # Check first parameter's coefficients
      first_param <- parsed$parameters[1, ]
      coefficients <- first_param$coefficients[[1]]

      expect_true(is.data.frame(coefficients) || is.list(coefficients))
      expect_true("value" %in% names(coefficients))
      expect_true("stderr" %in% names(coefficients))
      expect_true("condition_type" %in% names(coefficients))
      expect_true("conditions" %in% names(coefficients))

      # Check condition_type is valid
      valid_condition_types <- c("intercept", "linear_term", "CPT_combination",
                                 "variance", "random_variance", "random_covariance")
      expect_true(coefficients$condition_type[1] %in% valid_condition_types)
    })
  })
})

test_that("arcs array has correct structure", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      expect_true("arcs" %in% names(parsed))

      if (is.data.frame(parsed$arcs) && nrow(parsed$arcs) > 0) {
        expect_true("source_variable_id" %in% names(parsed$arcs))
        expect_true("target_variable_id" %in% names(parsed$arcs))
      } else if (is.list(parsed$arcs) && length(parsed$arcs) > 0) {
        # Check first arc
        first_arc <- parsed$arcs[[1]]
        expect_true("source_variable_id" %in% names(first_arc))
        expect_true("target_variable_id" %in% names(first_arc))
      }
    })
  })
})

test_that("export_abnFit works for MLE method with group.var", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      # create_test_abnfit_mle_groups() # Uncomment to regenerate the test object
      # load(file = "tests/testthat/testdata/abnfit_mle_groups.Rdata") # Load pre-saved fitted object
      load(file = "testdata/abnfit_mle_groups.Rdata") # Load pre-saved fitted object

      # ACT
      result <- export_abnFit(abn_fit)

      # ASSERT
      expect_type(result, "character")
      expect_true(jsonlite::validate(result))
    })
  })
})

test_that("grouped model export has variance parameters", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      # load(file = "tests/testthat/testdata/abnfit_mle_groups.Rdata")
      load(file = "testdata/abnfit_mle_groups.Rdata")

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      # Check that we have variance-related parameters
      variance_params <- parsed$parameters[
        parsed$parameters$coefficients[[1]]$condition_type %in%
          c("variance", "random_variance", "random_covariance"),
      ]

      # At least some nodes should have variance parameters
      expect_true(nrow(variance_params) > 0 || length(variance_params) > 0)
    })
  })
})

test_that("intercept parameters have empty conditions", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      # Find intercept parameters
      intercept_params <- parsed$parameters[parsed$parameters$name == "intercept", ]

      if (nrow(intercept_params) > 0) {
        # Check first intercept parameter
        first_intercept <- intercept_params[1, ]
        coef <- first_intercept$coefficients[[1]]
        expect_equal(coef$condition_type[1], "intercept")

        # Conditions should be empty list for intercepts
        conditions <- coef$conditions[[1]]
        expect_true(is.list(conditions))
        expect_equal(length(conditions), 0)
      }
    })
  })
})

test_that("linear_term parameters have parent in conditions", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      # Find parameters with condition_type = "linear_term"
      for (i in seq_len(nrow(parsed$parameters))) {
        param <- parsed$parameters[i, ]
        coef <- param$coefficients[[1]]

        if (coef$condition_type[1] == "linear_term") {
          conditions <- coef$conditions[[1]]
          expect_true(is.data.frame(conditions) || is.list(conditions))
          expect_true("parent_variable_id" %in% names(conditions))
        }
      }
    })
  })
})

test_that("link functions are correctly assigned", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      # Check that link functions match distribution types
      valid_links <- c("identity", "logit", "log")
      expect_true(all(parsed$parameters$link_function_name %in% valid_links))
    })
  })
})
