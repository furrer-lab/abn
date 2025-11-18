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

test_that("nodes section contains all expected nodes with correct structure", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      expected_nodes <- c("b1", "b2", "g1", "g2", "m1", "m2", "p1", "p2")
      expect_equal(sort(names(parsed$nodes)), sort(expected_nodes))

      # Check each node has required components
      for (node_name in expected_nodes) {
        node <- parsed$nodes[[node_name]]
        expect_true("label" %in% names(node))
        expect_true("distribution" %in% names(node))
        expect_true("parameterisation" %in% names(node))
        expect_equal(node$label, node_name)
        expect_true(node$distribution %in% c("gaussian", "binomial", "poisson", "multinomial"))
      }
    })
  })
})

test_that("parameterisation structure differs by distribution type", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      # Gaussian node should have intercept and coefficients
      g1_params <- parsed$nodes$g1$parameterisation
      expect_true("intercept" %in% names(g1_params))
      expect_true("coefficients" %in% names(g1_params))

      # Binomial node (b1) should have intercept only
      b1_params <- parsed$nodes$b1$parameterisation
      expect_true("intercept" %in% names(b1_params))

      # Multinomial node should have categories
      m1_params <- parsed$nodes$m1$parameterisation
      expect_true("categories" %in% names(m1_params))
    })
  })
})

test_that("arcs section has correct structure", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      expect_true("arcs" %in% names(parsed))
      arcs <- parsed$arcs
      expect_true("source" %in% names(arcs))
      expect_true("target" %in% names(arcs))
      expect_true("frequency" %in% names(arcs))
      expect_true("significance" %in% names(arcs))
      expect_true("constraint" %in% names(arcs))

      # Should be empty list for these fields as per the current implementation
      expect_equal(length(arcs$frequency), 0)
      expect_equal(length(arcs$significance), 0)
      expect_equal(length(arcs$constraint), 0)
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

