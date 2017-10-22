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
})
