context("test analyses functions")

test_that("pcr_ddct in average_ct mode", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_ddct(ct1,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain')

  expect_identical(res, rel_express1)
})

test_that("pcr_ddct in average_dct mode", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_ddct(ct2,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain',
                  mode = 'average_dct')

  expect_identical(res, rel_express2)
})

test_that("pcr_dct calculates and returns the correct values", {
   # make a data.frame of two identical columns
   pcr_hk <- data.frame(
     GAPDH1 = ct1$GAPDH,
     GAPDH2 = ct1$GAPDH
     )

   ## add grouping variable
   group_var <- rep(c('brain', 'kidney'), each = 6)

   # calculate caliberation
   res <- pcr_dct(pcr_hk,
                  group_var = group_var,
                  reference_group = 'brain')
   expect_identical(res, hk_express)
})

test_that("pcr_curve in average_amounts mode", {
  # make a vector of RNA amounts
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate standard curve
  standard_curve <- pcr_assess(ct3,
                               amount = amount,
                               method = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate standard amounts and error
  res <- pcr_curve(ct1,
                   group_var = group_var,
                   reference_gene = 'GAPDH',
                   reference_group = 'brain',
                   intercept = intercept,
                   slope = slope)

  expect_identical(res, standarized1)
})

test_that("pcr_curve in average_normalized mode", {
  # make a vector of RNA amounts
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate standard curve
  standard_curve <- pcr_assess(ct3,
                               amount = amount,
                               method = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate standard amounts and error
  res <- pcr_curve(ct2,
                   group_var = group_var,
                   reference_gene = 'GAPDH',
                   reference_group = 'brain',
                   intercept = intercept,
                   slope = slope,
                   mode = 'average_normalized')
  expect_identical(res, standarized2)
})

test_that("pcr_analyze calls the right methods", {
  # default method ddct
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_analyze(ct1,
                     group_var = group_var,
                     reference_gene = 'GAPDH',
                     reference_group = 'brain')

  expect_identical(res, rel_express1)

  # make a data.frame of two identical columns
  pcr_hk <- data.frame(
    GAPDH1 = ct1$GAPDH,
    GAPDH2 = ct1$GAPDH
  )

  # method: delta_ct
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate caliberation
  res <- pcr_analyze(pcr_hk,
                      group_var = group_var,
                      reference_group = 'brain',
                      method = 'delta_ct')
  expect_identical(res, hk_express)

  # method: relative_curve
  # make a vector of RNA amounts
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate standard curve
  standard_curve <- pcr_assess(ct3,
                               amount = amount,
                               method = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate standard amounts and error
  res <- pcr_analyze(ct1,
                    group_var = group_var,
                    reference_gene = 'GAPDH',
                    reference_group = 'brain',
                    intercept = intercept,
                    slope = slope,
                    method = 'relative_curve')
  expect_identical(res, standarized1)

})


