context("test normalization functions")

test_that("pcr_ave returns the proper averages in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_ave(pcr1_ct, group_var = group_var)

  expect_identical(ave[, 2:3], pcr1_norm[, 2:3])
})

test_that("pcr_dct returns the proper dct values in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_ave(pcr1_ct, group_var = group_var)
  dct = pcr_dct(ave, 'GAPDH')

  expect_equal(pcr1_norm$dct, dct$c_myc)
})

test_that("pcr_ddct retruns the proper ddct values in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_ave(pcr1_ct, group_var = group_var)
  dct = pcr_dct(ave, 'GAPDH')
  ddct = pcr_ddct(dct, 'brain')

  expect_equal(pcr1_norm$ddct, ddct$c_myc)
})

test_that("pcr_error returns the proper errors in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  sds = pcr_sd(pcr1_ct, group_var = group_var)
  errors = pcr_error(sds, reference_gene = 'GAPDH')

  expect_equal(pcr1_norm$error, errors$c_myc)
})

test_that("pcr_normalize returns the proper values in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  norm_rel <- pcr_normalize(pcr1_ct,
                            group_var = rep(c('brain', 'kidney'), each = 6),
                            reference_gene = 'GAPDH',
                            reference_group = 'brain')

  expect_equal(pcr1_norm$norm_rel, norm_rel$norm_rel)
  expect_equal(pcr1_norm$int_upper, norm_rel$int_upper)
  expect_equal(pcr1_norm$int_lower, norm_rel$int_lower)

  norm_rel <- pcr_normalize(pcr2_ct,
                            group_var = rep(c('brain', 'kidney'), each = 6),
                            reference_gene = 'GAPDH',
                            reference_group = 'brain',
                            mode = 'average_dct')

  expect_equal(pcr2_norm$norm_rel, norm_rel$norm_rel)
  expect_equal(pcr2_norm$int_upper, norm_rel$int_upper)
  expect_equal(pcr2_norm$int_lower, norm_rel$int_lower)
})
