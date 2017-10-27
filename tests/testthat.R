Sys.setenv("R_TESTS" = "")
library(testthat)
library(pcr)

test_check("pcr")
