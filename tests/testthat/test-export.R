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
      expect_true("graph" %in% names(parsed))
      expect_true("nodes" %in% names(parsed))
      expect_true("arcs" %in% names(parsed))
    })
  })
})

test_that("graph metadata is correctly exported for MLE", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      abn_fit <- create_test_abnfit_mle()

      # ACT
      parsed <- jsonlite::fromJSON(export_abnFit(abn_fit))

      # ASSERT
      expect_equal(parsed$graph$method, "mle")
      expect_equal(parsed$graph$n_nodes, 8)
      expect_equal(parsed$graph$n_observations, 250)

      # Check scores are present
      expect_true("scores" %in% names(parsed$graph))
      expect_true("aic" %in% names(parsed$graph$scores))
      expect_true("bic" %in% names(parsed$graph$scores))
      expect_true("mlik" %in% names(parsed$graph$scores))
      expect_true("mdl" %in% names(parsed$graph$scores))
      expect_true("groupVar" %in% names(parsed$graph))
      expect_true("groupedVariables" %in% names(parsed$graph))
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
      create_test_abnfit_mle_groups() # Uncomment to regenerate the test object
      load(file = "tests/testthat/testdata/abnfit_mle_groups.Rdata") # Load pre-saved fitted object
      # load(file = "/testdata/abnfit_mle_groups.Rdata") # Load pre-saved fitted object

      # ACT
      result <- export_abnFit(abn_fit)

      # ASSERT
      expect_type(result, "character")
      expect_true(jsonlite::validate(result))
    })
  })
})

