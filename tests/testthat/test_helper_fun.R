context("test helper functions")

test_that("pcr_average averages by group_var", {
  # average by group_var
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave <- pcr_average(ct1, group_var = group_var)

  ave2 <- mutate(ct1, group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) mean(x))

  expect_identical(ave, ave2)
})

test_that("pcr_average averages and returns a tidy data.frame", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave <- pcr_average(ct1, group_var = group_var, tidy = TRUE)

  ave2 <- mutate(ct1, group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) mean(x)) %>%
    gather(gene, average, -group)

  expect_identical(ave, ave2)
})

test_that("pcr_average averages by amount", {
  # average by amount
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
  ave <- pcr_average(ct3, amount = amount)

  ave2 <- mutate(ct3, amount = amount) %>%
    group_by(amount) %>%
    summarise_all(function(x) mean(x))

  expect_identical(ave, ave2)
})

test_that("pcr_normalize normalizes by subtraction", {

  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_average(ct1, group_var = group_var)
  dct = pcr_normalize(ave, 'GAPDH')

  expect_equal(rel_express1$normalized, dct$c_myc)
})

test_that("pcr_normalize normalizes by subtraction", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_average(ct1, group_var = group_var)
  norm = pcr_normalize(ave, 'GAPDH', mode = 'divide')

  norm2 <- ave$c_myc / ave$GAPDH

  expect_equal(norm$c_myc, norm2)
})

test_that("pcr_caliberate calculates the caliberated expression by subtraction", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_average(ct1, group_var = group_var)
  dct = pcr_normalize(ave, 'GAPDH')
  ddct = pcr_caliberate(dct, 'brain')

  expect_equal(rel_express1$caliberated, ddct$c_myc)
})

test_that("pcr_caliberate calculates the caliberated expression by division", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave = pcr_average(ct1, group_var = group_var)
  norm = pcr_normalize(ave, 'GAPDH')
  calib = pcr_caliberate(norm, 'brain', mode = 'divide')

  calib2 <- norm$c_myc / norm$c_myc[1]
  expect_equal(calib$c_myc, calib2)
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
                               method = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  # calculate amounts
  amounts2 <- pcr_amount(ct1, intercept, slope)

  # calculate cv
  group = rep(c('brain', 'kidney'), each = 6)
  cv <- pcr_cv(amounts2,
               group_var = group,
               reference_gene = 'GAPDH')

})


