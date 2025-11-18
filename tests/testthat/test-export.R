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

test_that("extract_parameters_by_distribution_grouping extracts parameters correctly for grouped nodes", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE - Gaussian node with grouping
      mu_gauss <- 1.5
      betas_gauss <- c(parent1 = 0.5, parent2 = -0.3)
      sigma_gauss <- 0.2
      sigma_alpha_gauss <- 0.1
      distribution_gauss <- "gaussian"
      node_id_gauss <- "g1"

      # ACT
      result_gauss <- extract_parameters_by_distribution_grouping(
        mu_gauss, betas_gauss, sigma_gauss, sigma_alpha_gauss, distribution_gauss, node_id_gauss
      )

      # ASSERT - Basic structure
      expect_type(result_gauss, "list")
      expect_true("fixed_effects" %in% names(result_gauss))
      expect_true("random_effects" %in% names(result_gauss))

      # ASSERT - Fixed effects structure
      expect_equal(result_gauss$fixed_effects$intercept, mu_gauss)
      expect_equal(result_gauss$fixed_effects$coefficients, as.list(betas_gauss))

      # ASSERT - Random effects structure
      expect_equal(result_gauss$random_effects$sigma, sigma_gauss)
      expect_equal(result_gauss$random_effects$sigma_alpha, sigma_alpha_gauss)
    })
  })
})

test_that("extract_parameters_by_distribution_grouping handles multinomial nodes correctly", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE - Multinomial node with grouping
      mu_mult <- c("m1.2" = 0.5, "m1.3" = 1.0)
      betas_mult <- matrix(c(0.1, 0.3, -0.2, 0.4), nrow = 2,
                           dimnames = list(c("2", "3"), c("parent1", "parent2")))
      sigma_mult <- matrix(c(0.02, 0.01, 0.01, 0.03), nrow = 2,
                           dimnames = list(c("m1.2~1", "m1.3~1"), c("m1.2~1", "m1.3~1")))
      sigma_alpha_mult <- matrix(c(0.02, 0.01, 0.01, 0.03), nrow = 2,
                                 dimnames = list(c("m1.2~1", "m1.3~1"), c("m1.2~1", "m1.3~1")))
      distribution_mult <- "multinomial"
      node_id_mult <- "m1"

      # ACT
      result_mult <- extract_parameters_by_distribution_grouping(
        mu_mult, betas_mult, sigma_mult, sigma_alpha_mult, distribution_mult, node_id_mult
      )

      # ASSERT - Structure for multinomial
      expect_type(result_mult, "list")
      expect_true("fixed_effects" %in% names(result_mult))
      expect_true("random_effects" %in% names(result_mult))
      expect_true("categories" %in% names(result_mult$fixed_effects))

      # ASSERT - Category structure
      expect_equal(length(result_mult$fixed_effects$categories), length(mu_mult))
      expect_true("category_2" %in% names(result_mult$fixed_effects$categories))
      expect_true("category_3" %in% names(result_mult$fixed_effects$categories))

      # ASSERT - Random effects for multinomial
      expect_equal(result_mult$random_effects$sigma, sigma_mult)
      expect_type(result_mult$random_effects$sigma_alpha, "list")
    })
  })
})

test_that("extract_parameters_by_distribution_grouping handles nodes with no parent coefficients", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE - Node with no parents (empty betas)
      mu_empty <- 2.0
      betas_empty <- numeric(0)  # No parents
      sigma_empty <- 0.5
      sigma_alpha_empty <- 0.2
      distribution_empty <- "binomial"
      node_id_empty <- "b1"

      # ACT
      result_empty <- extract_parameters_by_distribution_grouping(
        mu_empty, betas_empty, sigma_empty, sigma_alpha_empty, distribution_empty, node_id_empty
      )

      # ASSERT
      expect_type(result_empty, "list")
      expect_equal(result_empty$fixed_effects$intercept, mu_empty)
      expect_equal(result_empty$fixed_effects$coefficients, list())  # Empty list for no coefficients
      expect_equal(result_empty$random_effects$sigma, sigma_empty)
      expect_equal(result_empty$random_effects$sigma_alpha, sigma_alpha_empty)
    })
  })
})

test_that("extract_parameters_by_distribution_grouping handles matrix betas correctly", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE - Node with matrix coefficients
      mu_matrix <- 1.2
      betas_matrix <- matrix(c(0.3, -0.1), nrow = 1, ncol = 2,
                             dimnames = list(NULL, c("parent1", "parent2")))
      sigma_matrix <- 0.1
      sigma_alpha_matrix <- 0.05
      distribution_matrix <- "poisson"
      node_id_matrix <- "p1"

      # ACT
      result_matrix <- extract_parameters_by_distribution_grouping(
        mu_matrix, betas_matrix, sigma_matrix, sigma_alpha_matrix, distribution_matrix, node_id_matrix
      )

      # ASSERT
      expect_type(result_matrix, "list")
      expect_equal(result_matrix$fixed_effects$intercept, mu_matrix)
      expect_type(result_matrix$fixed_effects$coefficients, "list")
      expect_equal(length(result_matrix$fixed_effects$coefficients), ncol(betas_matrix))
      expect_true("parent1" %in% names(result_matrix$fixed_effects$coefficients))
      expect_true("parent2" %in% names(result_matrix$fixed_effects$coefficients))
    })
  })
})

