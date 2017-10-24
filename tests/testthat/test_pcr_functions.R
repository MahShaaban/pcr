context("test normalization functions")

test_that("pcr_ave returns the proper averages in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_ave(pcr1_ct, group_var = group_var)

  expect_identical(ave[, 2:3], pcr1_norm[, 2:3])
})

test_that("pcr_norm returns the proper dct values in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_ave(pcr1_ct, group_var = group_var)
  dct = pcr_norm(ave, 'GAPDH')

  expect_equal(pcr1_norm$norm, dct$c_myc)
})

test_that("pcr_calib retruns the proper ddct values in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_ave(pcr1_ct, group_var = group_var)
  dct = pcr_norm(ave, 'GAPDH')
  ddct = pcr_calib(dct, 'brain')

  expect_equal(pcr1_norm$norm_calib, ddct$c_myc)
})

test_that("pcr_error returns the proper errors in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  sds = pcr_sd(pcr1_ct, group_var = group_var)
  errors = pcr_error(sds, reference_gene = 'GAPDH')

  expect_equal(pcr1_norm$error, errors$c_myc)
})

test_that("pcr_cv calculates the correct cv", {
  # input dilution variable
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate curve
  standard_curve <- pcr_assess(pcr_dilute,
                               amount = amount,
                               mode = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  # calculate amounts
  amounts2 <- pcr_amount(pcr1_ct, intercept, slope)

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
  standard_curve <- pcr_assess(pcr_dilute, amount = amount, mode = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  # calculate manually
  amounts1 <- pcr1_ct %>%
    mutate(c_myc = 10 ^ ((c_myc - intercept[1]) / slope[1]),
           GAPDH = 10 ^ ((GAPDH - intercept[2]) / slope[2]))

  # test
  amounts2 <- pcr_amount(pcr1_ct, intercept, slope)
  expect_equal(amounts1, amounts2)
})

test_that("pcr_assess retruns the proper values in the right formate", {
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # default mode "effeciency"
  eff <- pcr_assess(pcr_dilute,
                    amount = amount)
  expect_s3_class(eff, 'data.frame')

  # testing problem!!
  # standard curve mode
  co_eff1 <- coefficients(lm(pcr_dilute$c_myc ~ log10(amount)))
  co_eff2 <- coefficients(lm(pcr_dilute$GAPDH ~ log10(amount)))

  sc <- pcr_assess(pcr_dilute,
                   amount = amount,
                   mode = 'standard_curve')
  expect_s3_class(sc, 'data.frame')
  # testing problem!!
})

test_that("pcr_analyze returns the proper values in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  norm_rel <- pcr_analyze(pcr1_ct,
                            group_var = rep(c('brain', 'kidney'), each = 6),
                            reference_gene = 'GAPDH',
                            reference_group = 'brain')

  expect_equal(pcr1_norm$norm_rel, norm_rel$norm_rel)
  expect_equal(pcr1_norm$int_upper, norm_rel$int_upper)
  expect_equal(pcr1_norm$int_lower, norm_rel$int_lower)

  norm_rel <- pcr_analyze(pcr2_ct,
                            group_var = rep(c('brain', 'kidney'), each = 6),
                            reference_gene = 'GAPDH',
                            reference_group = 'brain',
                            mode = 'average_dct')

  expect_equal(pcr2_norm$norm_rel, norm_rel$norm_rel)
  expect_equal(pcr2_norm$int_upper, norm_rel$int_upper)
  expect_equal(pcr2_norm$int_lower, norm_rel$int_lower)

  group_var <- rep(c('brain', 'kidney'), each = 6)
  rel <- pcr_analyze(pcr_hk,
                     group_var = group_var,
                     reference_group = 'brain',
                     method = 'delta_ct')

  expect_identical(pcr_hk_calib, rel)

  ## calculate standard curve
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  standard_curve <- pcr_assess(pcr_dilute,
                             amount = amount,
                             mode = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  ## calculate standard amounts and error
  standarized <- pcr_analyze(pcr1_ct,
                             group_var = group_var,
                             reference_gene = 'GAPDH',
                             reference_group = 'brain',
                             intercept = intercept,
                             slope = slope,
                             method = 'relative_curve')
  # testing problem!!
})




