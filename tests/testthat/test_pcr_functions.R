context("test normalization functions")

test_that("pcr_ave returns the proper averages in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_ave(ct1, group_var = group_var)

  expect_identical(ave[, 2:3], rel_express1[, 2:3])
})

test_that("pcr_norm returns the proper dct values in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_ave(ct1, group_var = group_var)
  dct = pcr_norm(ave, 'GAPDH')

  expect_equal(rel_express1$norm, dct$c_myc)
})

test_that("pcr_calib retruns the proper ddct values in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_ave(ct1, group_var = group_var)
  dct = pcr_norm(ave, 'GAPDH')
  ddct = pcr_calib(dct, 'brain')

  expect_equal(rel_express1$norm_calib, ddct$c_myc)
})

test_that("pcr_error returns the proper errors in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  sds = pcr_sd(ct1, group_var = group_var)
  errors = pcr_error(sds, reference_gene = 'GAPDH')

  expect_equal(rel_express1$error, errors$c_myc)
})

test_that("pcr_cv calculates the correct cv", {
  # input dilution variable
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate curve
  standard_curve <- pcr_assess(ct3,
                               amount = amount,
                               mode = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  # calculate amounts
  amounts2 <- pcr_amount(ct1, intercept, slope)

  # calculate cv
  group = rep(c('brain', 'kidney'), each = 6)
  cv <- pcr_cv(amounts2,
               group_var = group,
               reference_gene = 'GAPDH')

  expect_equal(round(pcr_amounts$cv, 2), round(cv$c_myc, 2))
})

test_that("pcr_amount calculates amount and returns the right formate", {
  # input dilution variable
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate curve
  standard_curve <- pcr_assess(ct3, amount = amount, mode = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  # calculate manually
  amounts1 <- ct1 %>%
    mutate(c_myc = 10 ^ ((c_myc - intercept[1]) / slope[1]),
           GAPDH = 10 ^ ((GAPDH - intercept[2]) / slope[2]))

  # test
  amounts2 <- pcr_amount(ct1, intercept, slope)
  expect_equal(amounts1, amounts2)
})

test_that("pcr_assess retruns the proper values in the right formate", {
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # default mode "effeciency"
  eff <- pcr_assess(ct3,
                    amount = amount)
  expect_s3_class(eff, 'data.frame')

  # testing problem!!
  # standard curve mode
  co_eff1 <- coefficients(lm(ct3$c_myc ~ log10(amount)))
  co_eff2 <- coefficients(lm(ct3$GAPDH ~ log10(amount)))

  sc <- pcr_assess(ct3,
                   amount = amount,
                   mode = 'standard_curve')
  expect_s3_class(sc, 'data.frame')
  # testing problem!!
})

test_that("pcr_analyze returns the proper values in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  norm_rel <- pcr_analyze(ct1,
                            group_var = rep(c('brain', 'kidney'), each = 6),
                            reference_gene = 'GAPDH',
                            reference_group = 'brain')

  expect_equal(rel_express1$norm_rel, norm_rel$norm_rel)
  expect_equal(rel_express1$int_upper, norm_rel$int_upper)
  expect_equal(rel_express1$int_lower, norm_rel$int_lower)

  norm_rel <- pcr_analyze(ct2,
                            group_var = rep(c('brain', 'kidney'), each = 6),
                            reference_gene = 'GAPDH',
                            reference_group = 'brain',
                            mode = 'average_dct')

  expect_equal(rel_express2$norm_rel, norm_rel$norm_rel)
  expect_equal(rel_express2$int_upper, norm_rel$int_upper)
  expect_equal(rel_express2$int_lower, norm_rel$int_lower)

  group_var <- rep(c('brain', 'kidney'), each = 6)
  rel <- pcr_analyze(housekeeping,
                     group_var = group_var,
                     reference_group = 'brain',
                     method = 'delta_ct')

  expect_identical(hk_express, rel)

  ## calculate standard curve
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  standard_curve <- pcr_assess(ct3,
                             amount = amount,
                             mode = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  ## calculate standard amounts and error
  standarized <- pcr_analyze(ct1,
                             group_var = group_var,
                             reference_gene = 'GAPDH',
                             reference_group = 'brain',
                             intercept = intercept,
                             slope = slope,
                             method = 'relative_curve')
  # testing problem!!
})




