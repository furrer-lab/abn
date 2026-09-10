# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(abn)

# Set the number of threads for BLAS and OpenMP if possible
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(2)
} else {
  # E.g., on CRAN noSuggests environment, RhpcBLASctl is not available
  message("RhpcBLASctl is not available, BLAS thread control will not be set.")
}

test_check("abn")
