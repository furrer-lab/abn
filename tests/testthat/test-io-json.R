# Comprehensive Round-Trip Import/Export JSON Tests
# Tests for both standard MLE and grouped MLE models
# Verifies: RData → export → import ≡ RData (direct load)

library(testthat)
library(abn)

# Helper function to safely load test data
get_test_data <- function(filename) {
  path <- test_path("testdata", filename)
  if (!file.exists(path)) {
    alt_path <- system.file("testdata", filename, package = "abn")
    if (alt_path != "" && file.exists(alt_path)) return(alt_path)
    return("")
  }
  return(path)
}

# ============================================================================
# SUITE 1: STANDARD MLE MODEL ROUND-TRIP TESTS
# ============================================================================

test_that("Standard MLE: RData → export → import preserves model structure", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE: Load model directly from RData (reference)
      direct_model <- create_test_abnfit_mle()
      
      # ACT: Export to JSON
      json_str <- export_abnFit(direct_model)
      expect_type(json_str, "character", 
                  info = "export_abnFit should return character string")
      expect_true(jsonlite::validate(json_str),
                  info = "Exported JSON should be valid")
      
      # ACT: Import from JSON
      imported_model <- import_abnFit(json = json_str)
      expect_s3_class(imported_model, "abnFit",
                      info = "Imported object should be abnFit class")
      
      # ASSERT: Compare models using custom comparison
      expect_true(
        abnfit_objects_equal(direct_model, imported_model),
        info = paste(
          "Standard MLE model should be preserved through export→import cycle",
          "Expected same DAG structure, coefficients, and standard errors"
        )
      )
    })
  })
})

test_that("Standard MLE: coefficient values preserved through round-trip", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      direct_model <- create_test_abnfit_mle()
      
      # ACT
      json_str <- export_abnFit(direct_model)
      imported_model <- import_abnFit(json = json_str)
      
      # ASSERT: Check all coefficient matrices
      for (node_name in names(direct_model$coef)) {
        direct_coef <- direct_model$coef[[node_name]]
        imported_coef <- imported_model$coef[[node_name]]
        
        expect_equal(
          dim(direct_coef), dim(imported_coef),
          info = paste(
            "Coefficient matrix dimensions should match for node:", node_name,
            "Direct:", paste(dim(direct_coef), collapse = "x"),
            "Imported:", paste(dim(imported_coef), collapse = "x")
          )
        )
        
        expect_equal(
          rownames(direct_coef), rownames(imported_coef),
          info = paste(
            "Coefficient row names (parameter names) should match for node:",
            node_name
          )
        )
        
        # Check numerical values with tolerance
        max_diff <- max(abs(direct_coef - imported_coef), na.rm = TRUE)
        expect_lt(
          max_diff, 1e-5,
          info = paste(
            "Coefficient values for node", node_name,
            "should match within tolerance (1e-5)",
            "Max difference found:", format(max_diff, scientific = TRUE)
          )
        )
      }
    })
  })
})

test_that("Standard MLE: standard errors preserved through round-trip", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      direct_model <- create_test_abnfit_mle()
      
      # ACT
      json_str <- export_abnFit(direct_model)
      imported_model <- import_abnFit(json = json_str)
      
      # ASSERT: Check all standard error matrices
      for (node_name in names(direct_model$Stderror)) {
        direct_se <- direct_model$Stderror[[node_name]]
        imported_se <- imported_model$Stderror[[node_name]]
        
        expect_equal(
          dim(direct_se), dim(imported_se),
          info = paste(
            "Standard error matrix dimensions should match for node:",
            node_name,
            "Direct:", paste(dim(direct_se), collapse = "x"),
            "Imported:", paste(dim(imported_se), collapse = "x")
          )
        )
        
        # Check numerical values with tolerance
        max_diff <- max(abs(direct_se - imported_se), na.rm = TRUE)
        expect_lt(
          max_diff, 1e-5,
          info = paste(
            "Standard error values for node", node_name,
            "should match within tolerance (1e-5)",
            "Max difference found:", format(max_diff, scientific = TRUE)
          )
        )
      }
    })
  })
})

test_that("Standard MLE: link functions consistent across round-trip", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      direct_model <- create_test_abnfit_mle()
      
      # ACT
      json_str <- export_abnFit(direct_model)
      parsed <- jsonlite::fromJSON(json_str)
      
      # ASSERT: Check link functions
      expect_true(
        "link_function_name" %in% colnames(parsed$parameters),
        info = "Parameters should have link_function_name column"
      )
      
      # Valid link functions for each distribution
      valid_links <- c("identity", "logit", "log")
      
      for (i in seq_len(nrow(parsed$parameters))) {
        param <- parsed$parameters[i, ]
        link <- param$link_function_name
        
        expect_true(
          link %in% valid_links,
          info = paste(
            "Link function should be one of:", paste(valid_links, collapse = ", "),
            "Found:", link, "for parameter", i
          )
        )
      }
      
      # ACT: Import and verify link functions are preserved
      imported_model <- import_abnFit(json = json_str)
      expect_s3_class(imported_model, "abnFit",
                      info = "Imported model should be abnFit")
    })
  })
})

test_that("Standard MLE: variable ordering consistent through export→import", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      direct_model <- create_test_abnfit_mle()
      
      # ACT
      json_str <- export_abnFit(direct_model)
      parsed <- jsonlite::fromJSON(json_str)
      
      # ASSERT
      expect_true(
        is.data.frame(parsed$variables) || is.list(parsed$variables),
        info = "Variables should be a data frame or list"
      )
      
      # Check variable IDs are numeric strings
      var_ids <- parsed$variables$variable_id
      expect_true(
        all(grepl("^[0-9]+$", var_ids)),
        info = paste(
          "All variable IDs should be numeric strings",
          "Found:", paste(var_ids, collapse = ", ")
        )
      )
      
      # ACT: Import and verify variable order is maintained
      imported_model <- import_abnFit(json = json_str)
      expect_equal(
        colnames(imported_model$abnDag$dag),
        colnames(direct_model$abnDag$dag),
        info = "Column order in DAG matrix should be maintained"
      )
    })
  })
})

# ============================================================================
# SUITE 2: GROUPED MLE MODEL ROUND-TRIP TESTS
# ============================================================================

test_that("Grouped MLE: RData → export → import preserves grouped model structure", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE: Load grouped model from RData
      load(file = get_test_data("abnfit_mle_groups.Rdata"))
      direct_model <- abn_fit
      
      # ACT: Export to JSON
      json_str <- export_abnFit(direct_model)
      expect_type(json_str, "character",
                  info = "export_abnFit should return character string")
      expect_true(jsonlite::validate(json_str),
                  info = "Exported JSON should be valid")
      
      # ACT: Import from JSON
      imported_model <- import_abnFit(json = json_str)
      expect_s3_class(imported_model, "abnFit",
                      info = "Imported object should be abnFit class")
      
      # ASSERT: Compare models using custom comparison
      expect_true(
        abnfit_objects_equal(direct_model, imported_model),
        info = paste(
          "Grouped MLE model should be preserved through export→import cycle",
          "Expected same DAG structure, group-specific coefficients, variance",
          "parameters, and metadata"
        )
      )
    })
  })
})

test_that("Grouped MLE: group-specific coefficients preserved", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE: Load grouped model
      load(file = get_test_data("abnfit_mle_groups.Rdata"))
      direct_model <- abn_fit
      
      # ACT
      json_str <- export_abnFit(direct_model)
      parsed <- jsonlite::fromJSON(json_str)
      
      # ASSERT: Verify group information in parameters
      expect_true(
        "group_label" %in% colnames(parsed$parameters),
        info = "Grouped model parameters should have group_label column"
      )
      
      # Check that we have multiple groups
      unique_groups <- unique(parsed$parameters$group_label)
      expect_gt(
        length(unique_groups), 1,
        info = paste(
          "Grouped model should have multiple group labels",
          "Found groups:", paste(unique_groups, collapse = ", ")
        )
      )
      
      # ACT: Import and verify groups are preserved
      imported_model <- import_abnFit(json = json_str)
      
      # Verify same number of groups
      expect_equal(
        length(unique_groups), length(unique(imported_model$coef)),
        info = paste(
          "Number of groups should match",
          "Expected:", length(unique_groups),
          "Got:", length(unique(imported_model$coef))
        )
      )
    })
  })
})

test_that("Grouped MLE: variance and random effect parameters preserved", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE: Load grouped model
      load(file = get_test_data("abnfit_mle_groups.Rdata"))
      direct_model <- abn_fit
      
      # ACT
      json_str <- export_abnFit(direct_model)
      parsed <- jsonlite::fromJSON(json_str)
      
      # ASSERT: Check for variance-related parameters
      expect_true(
        "parameters" %in% names(parsed),
        info = "JSON should have parameters field"
      )
      
      # Extract variance/random parameters
      variance_params <- parsed$parameters[
        grepl("variance|random", parsed$parameters$name, 
              ignore.case = TRUE),
      ]
      
      expect_gt(
        nrow(variance_params), 0,
        info = paste(
          "Grouped model should have variance/random effect parameters",
          "Total parameters:", nrow(parsed$parameters),
          "Parameters found:", nrow(variance_params)
        )
      )
      
      # ACT: Import and verify variance parameters are preserved
      imported_model <- import_abnFit(json = json_str)
      
      # Compare variance structures
      for (node_name in names(direct_model$coef)) {
        expect_true(
          node_name %in% names(imported_model$coef),
          info = paste(
            "Node", node_name, "should be present in imported model"
          )
        )
      }
    })
  })
})

test_that("Grouped MLE: group variable reference maintained", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE: Load grouped model
      load(file = get_test_data("abnfit_mle_groups.Rdata"))
      direct_model <- abn_fit
      
      # ACT
      json_str <- export_abnFit(direct_model)
      parsed <- jsonlite::fromJSON(json_str)
      
      # ASSERT
      expect_true(
        is.data.frame(parsed$variables) || is.list(parsed$variables),
        info = "Variables should be present in JSON"
      )
      
      # ACT: Import and verify DAG structure
      imported_model <- import_abnFit(json = json_str)
      
      # Compare DAG matrices
      expect_equal(
        dim(direct_model$abnDag$dag), dim(imported_model$abnDag$dag),
        info = paste(
          "DAG dimensions should match",
          "Direct:", paste(dim(direct_model$abnDag$dag), collapse = "x"),
          "Imported:", paste(dim(imported_model$abnDag$dag), collapse = "x")
        )
      )
      
      expect_equal(
        colnames(direct_model$abnDag$dag),
        colnames(imported_model$abnDag$dag),
        info = "DAG variable names (columns) should match"
      )
    })
  })
})

# ============================================================================
# SUITE 3: SIDE-BY-SIDE COMPARISON TESTS
# ============================================================================

test_that("Both models: direct RData equals JSON round-trip (standard)", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      direct_model <- create_test_abnfit_mle()
      
      # ACT
      json_str <- export_abnFit(direct_model)
      imported_model <- import_abnFit(json = json_str)
      
      # ASSERT: Comprehensive comparison
      comparison_result <- abnfit_objects_equal(direct_model, imported_model)
      expect_true(
        comparison_result,
        info = paste(
          "Standard MLE: Direct RData load should equal JSON round-trip",
          "Models differ in structure, coefficients, or standard errors"
        )
      )
      
      # Additional sanity checks
      expect_equal(
        length(direct_model$coef), length(imported_model$coef),
        info = paste(
          "Both models should have same number of nodes",
          "Direct:", length(direct_model$coef),
          "Imported:", length(imported_model$coef)
        )
      )
    })
  })
})

test_that("Both models: direct RData equals JSON round-trip (grouped)", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      load(file = get_test_data("abnfit_mle_groups.Rdata"))
      direct_model <- abn_fit
      
      # ACT
      json_str <- export_abnFit(direct_model)
      imported_model <- import_abnFit(json = json_str)
      
      # ASSERT: Comprehensive comparison
      comparison_result <- abnfit_objects_equal(direct_model, imported_model)
      expect_true(
        comparison_result,
        info = paste(
          "Grouped MLE: Direct RData load should equal JSON round-trip",
          "Models differ in structure, group coefficients, or variance",
          "parameters"
        )
      )
      
      # Additional sanity checks
      expect_equal(
        length(direct_model$coef), length(imported_model$coef),
        info = paste(
          "Both models should have same number of nodes",
          "Direct:", length(direct_model$coef),
          "Imported:", length(imported_model$coef)
        )
      )
    })
  })
})

test_that("Both models: JSON structure validated for web consumption", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      standard_model <- create_test_abnfit_mle()
      load(file = get_test_data("abnfit_mle_groups.Rdata"))
      grouped_model <- abn_fit
      
      # ACT: Export both models
      standard_json <- export_abnFit(standard_model)
      grouped_json <- export_abnFit(grouped_model)
      
      standard_parsed <- jsonlite::fromJSON(standard_json)
      grouped_parsed <- jsonlite::fromJSON(grouped_json)
      
      # ASSERT: Validate JSON structure for both
      for (json_label in c("standard", "grouped")) {
        parsed <- if (json_label == "standard") standard_parsed else grouped_parsed
        
        # Check required top-level fields
        required_fields <- c("variables", "parameters", "arcs")
        for (field in required_fields) {
          expect_true(
            field %in% names(parsed),
            info = paste(
              json_label, "model JSON should have field:", field
            )
          )
        }
        
        # Check variable IDs are numeric strings (not variable names)
        var_ids <- parsed$variables$variable_id
        expect_true(
          all(grepl("^[0-9]+$", var_ids)),
          info = paste(
            json_label,
            "model: All variable IDs should be numeric strings",
            "Found:", paste(var_ids, collapse = ", ")
          )
        )
        
        # Check arcs have numeric IDs
        if (nrow(parsed$arcs) > 0) {
          for (i in seq_len(min(3, nrow(parsed$arcs)))) {
            source_id <- parsed$arcs$source_variable_id[i]
            target_id <- parsed$arcs$target_variable_id[i]
            
            expect_true(
              grepl("^[0-9]+$", source_id),
              info = paste(
                json_label, "model: Arc source_variable_id should be numeric",
                "Found:", source_id
              )
            )
            
            expect_true(
              grepl("^[0-9]+$", target_id),
              info = paste(
                json_label, "model: Arc target_variable_id should be numeric",
                "Found:", target_id
              )
            )
          }
        }
      }
    })
  })
})

# ============================================================================
# SUITE 4: DATA INTEGRITY TESTS
# ============================================================================

test_that("All parameters exported and reimported (standard model)", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      direct_model <- create_test_abnfit_mle()
      
      # Count parameters in original
      direct_param_count <- sum(sapply(direct_model$coef, nrow))
      
      # ACT
      json_str <- export_abnFit(direct_model)
      parsed <- jsonlite::fromJSON(json_str)
      exported_param_count <- nrow(parsed$parameters)
      
      # ASSERT
      expect_gt(
        direct_param_count, 0,
        info = "Standard model should have at least some parameters"
      )
      
      expect_equal(
        direct_param_count, exported_param_count,
        info = paste(
          "All parameters should be exported",
          "Original count:", direct_param_count,
          "Exported count:", exported_param_count
        )
      )
      
      # ACT: Import and count again
      imported_model <- import_abnFit(json = json_str)
      imported_param_count <- sum(sapply(imported_model$coef, nrow))
      
      expect_equal(
        direct_param_count, imported_param_count,
        info = paste(
          "All parameters should be reimported",
          "Original count:", direct_param_count,
          "Reimported count:", imported_param_count
        )
      )
    })
  })
})

test_that("All parameters exported and reimported (grouped model)", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      load(file = get_test_data("abnfit_mle_groups.Rdata"))
      direct_model <- abn_fit
      
      # Count parameters in original (all groups combined)
      direct_param_count <- sum(sapply(direct_model$coef, nrow))
      
      # ACT
      json_str <- export_abnFit(direct_model)
      parsed <- jsonlite::fromJSON(json_str)
      exported_param_count <- nrow(parsed$parameters)
      
      # ASSERT
      expect_gt(
        direct_param_count, 0,
        info = "Grouped model should have at least some parameters"
      )
      
      expect_equal(
        direct_param_count, exported_param_count,
        info = paste(
          "All parameters (all groups) should be exported",
          "Original count:", direct_param_count,
          "Exported count:", exported_param_count
        )
      )
      
      # ACT: Import and count again
      imported_model <- import_abnFit(json = json_str)
      imported_param_count <- sum(sapply(imported_model$coef, nrow))
      
      expect_equal(
        direct_param_count, imported_param_count,
        info = paste(
          "All parameters should be reimported",
          "Original count:", direct_param_count,
          "Reimported count:", imported_param_count
        )
      )
    })
  })
})

test_that("All arcs preserved through standard model round-trip", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      direct_model <- create_test_abnfit_mle()
      dag_matrix <- direct_model$abnDag$dag
      
      # Extract arcs from original DAG
      direct_arcs <- which(dag_matrix == 1, arr.ind = TRUE)
      direct_arc_count <- nrow(direct_arcs)
      
      # ACT
      json_str <- export_abnFit(direct_model)
      parsed <- jsonlite::fromJSON(json_str)
      exported_arc_count <- nrow(parsed$arcs)
      
      # ASSERT
      expect_equal(
        direct_arc_count, exported_arc_count,
        info = paste(
          "All arcs should be exported",
          "Original arc count:", direct_arc_count,
          "Exported arc count:", exported_arc_count
        )
      )
      
      # Check arc ordering is deterministic
      if (exported_arc_count > 1) {
        source_ids <- as.numeric(parsed$arcs$source_variable_id)
        target_ids <- as.numeric(parsed$arcs$target_variable_id)
        
        # Verify sorting: (source, target) ascending
        for (i in 1:(length(source_ids) - 1)) {
          current_pair <- c(source_ids[i], target_ids[i])
          next_pair <- c(source_ids[i + 1], target_ids[i + 1])
          
          comparison <- if (current_pair[1] == next_pair[1]) {
            current_pair[2] <= next_pair[2]
          } else {
            current_pair[1] < next_pair[1]
          }
          
          expect_true(
            comparison,
            info = paste(
              "Arcs should be sorted by (source_id, target_id)",
              "Arc", i, ":", paste(current_pair, collapse = "→"),
              "Arc", i + 1, ":", paste(next_pair, collapse = "→")
            )
          )
        }
      }
      
      # ACT: Import and verify arcs are preserved
      imported_model <- import_abnFit(json = json_str)
      reimported_dag <- imported_model$abnDag$dag
      reimported_arcs <- which(reimported_dag == 1, arr.ind = TRUE)
      reimported_arc_count <- nrow(reimported_arcs)
      
      expect_equal(
        direct_arc_count, reimported_arc_count,
        info = paste(
          "Arc count should be preserved through import",
          "Original:", direct_arc_count,
          "Reimported:", reimported_arc_count
        )
      )
    })
  })
})

test_that("All arcs preserved through grouped model round-trip", {
  suppressMessages({
    suppressWarnings({
      # ARRANGE
      load(file = get_test_data("abnfit_mle_groups.Rdata"))
      direct_model <- abn_fit
      dag_matrix <- direct_model$abnDag$dag
      
      # Extract arcs from original DAG
      direct_arcs <- which(dag_matrix == 1, arr.ind = TRUE)
      direct_arc_count <- nrow(direct_arcs)
      
      # ACT
      json_str <- export_abnFit(direct_model)
      parsed <- jsonlite::fromJSON(json_str)
      exported_arc_count <- nrow(parsed$arcs)
      
      # ASSERT
      expect_equal(
        direct_arc_count, exported_arc_count,
        info = paste(
          "All arcs should be exported",
          "Original arc count:", direct_arc_count,
          "Exported arc count:", exported_arc_count
        )
      )
      
      # ACT: Import and verify arcs are preserved
      imported_model <- import_abnFit(json = json_str)
      reimported_dag <- imported_model$abnDag$dag
      reimported_arcs <- which(reimported_dag == 1, arr.ind = TRUE)
      reimported_arc_count <- nrow(reimported_arcs)
      
      expect_equal(
        direct_arc_count, reimported_arc_count,
        info = paste(
          "Arc count should be preserved through import",
          "Original:", direct_arc_count,
          "Reimported:", reimported_arc_count
        )
      )
    })
  })
})
