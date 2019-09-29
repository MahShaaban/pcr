context("test helper functions")

test_that(".pcr_average averages by group_var", {
  # average by group_var
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave <- .pcr_average(ct1, var = group_var)

  ave2 <- aggregate(ct1, by = list(var = group_var), mean)

  expect_identical(ave, ave2)
})

test_that(".pcr_average averages and returns a tidy data.frame", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave <- .pcr_average(ct1, var = group_var, tidy = TRUE)

  ave2 <- aggregate(ct1, by = list(var = group_var), mean)

  ave2 <- data.frame(
    var = rep(ave2$var, ncol(ct1)),
    gene = rep(names(ct1), each = ncol(ct1)),
    average = unlist(ave2[, names(ct1)], use.names = FALSE)
  )

  expect_identical(ave, ave2)
})

test_that(".pcr_average averages by amount", {
  # average by amount
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
  ave <- .pcr_average(ct3, var = amount)

  ave2 <- aggregate(ct3, by = list(var = amount), mean)

  data.frame(
    var = rep(ave2$var, ncol(ct3)),
    gene = rep(names(ct3), each = nrow(ct3)),
    average = unlist(ave2[, names(ct3)], use.names = FALSE)
  )
})

test_that(".pcr_normalize normalizes by subtraction", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave <- .pcr_average(ct1, var = group_var)
  norm <- .pcr_normalize(ave[, -1], 'GAPDH', mode = 'subtract')

  norm2 <- ave$c_myc - ave$GAPDH

  expect_equal(norm$c_myc, norm2)
})

test_that(".pcr_normalize normalizes by subtraction", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave <- .pcr_average(ct1, var = group_var)
  norm <- .pcr_normalize(ave[, -1], 'GAPDH', mode = 'divide')

  norm2 <- ave$c_myc / ave$GAPDH

  expect_equal(norm$c_myc, norm2)
})

test_that(".pcr_calibrate calculates calibrated expression by subtraction", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave <- .pcr_average(ct1, var = group_var)
  norm <- .pcr_normalize(ave[, -1], 'GAPDH')
  calib <- .pcr_calibrate(norm,
                          group_var = unique(group_var),
                          reference_group = 'brain',
                          mode = 'subtract')

  calib2 <- norm$c_myc - norm$c_myc[1]
  expect_equal(calib$c_myc, calib2)
})

test_that(".pcr_calibrate calculates calibrated expression of two genes", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ct <- ct1
  ct$c_myc2 <- ct$c_myc

  ave <- .pcr_average(ct, var = group_var)
  norm <- .pcr_normalize(ave[, -1], 'GAPDH')
  calib <- .pcr_calibrate(norm,
                          group_var = unique(group_var),
                          reference_group = 'brain',
                          mode = 'subtract')

  calib2 <- norm$c_myc - norm$c_myc[1]
  expect_equal(calib$c_myc, calib2)
})

test_that(".pcr_calibrate calculates the calibrated expression by division", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave <- .pcr_average(ct1, var = group_var)
  norm <- .pcr_normalize(ave[, -1], 'GAPDH')
  calib <- .pcr_calibrate(norm,
                          group_var = unique(group_var),
                          reference_group = 'brain',
                          mode = 'divide')

  calib2 <- norm$c_myc / norm$c_myc[1]
  expect_equal(calib$c_myc, calib2)
})

test_that(".pcr_error returns the proper errors in the right format", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  sds <- .pcr_sd(ct1, group_var = group_var)
  errors <- .pcr_error(sds[, -1], reference_gene = 'GAPDH')

  errors2 <- sqrt((sds$c_myc^2) + (sds$GAPDH^2))
  expect_equal(errors$c_myc, errors2)
})

test_that(".pcr_cv calculates the correct cv", {
  # input dilution variable
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate curve
  standard_curve <- pcr_assess(ct3,
                               amount = amount,
                               method = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  # calculate amounts
  amounts2 <- .pcr_amount(ct1, intercept, slope)

  # calculate cv
  group <- rep(c('brain', 'kidney'), each = 6)
  cv <- .pcr_cv(amounts2,
               group_var = group)

})
