context("test helper functions")

test_that(".pcr_average works", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  ave <- apply(ct1, 2, .pcr_average, var = group_var)

  ave2 <- aggregate(ct1, by = list(var = group_var), mean)

  expect_equal(ave[,1], ave2$c_myc)
  expect_equal(ave[,2], ave2$GAPDH)
})

test_that(".pcr_sd works", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  sds <- apply(ct1, 2, .pcr_sd, var = group_var)

  sds2 <- aggregate(ct1, by = list(var = group_var), sd)

  expect_equal(sds[,1], sds2$c_myc)
  expect_equal(sds[,2], sds2$GAPDH)
})

test_that(".pcr_cv works", {
  group_var <- rep(c('brain', 'kidney'), each = 6)
  cv <- apply(ct1, 2, .pcr_cv, var = group_var)

  cv2 <- aggregate(ct1, by = list(var = group_var), function(x) sd(x)/mean(x))

  expect_equal(cv[,1], cv2$c_myc)
  expect_equal(cv[,2], cv2$GAPDH)
})

test_that(".pcr_normalize works", {
  norm1 <- .pcr_normalize(ct1$c_myc, ct1$GAPDH)
  norm2 <- .pcr_normalize(ct1$c_myc, ct1$GAPDH, mode = 'divide')

  expect_equal(norm1, ct1$c_myc - ct1$GAPDH)
  expect_equal(norm2, ct1$c_myc / ct1$GAPDH)
})

test_that(".pcr_error works", {
  error <- .pcr_error(sd(ct1$c_myc), sd(ct1$GAPDH))
  error2 <- sqrt(sd(ct1$c_myc)^2 + sd(ct1$GAPDH)^2)
  expect_equal(error, error2)
})

test_that(".pcr_amount works", {
  line <- .pcr_amount(ct1$c_myc, 1, .1)
  line2 <- 10^((ct1$c_myc - 1) / .1)
  expect_equal(line, line2)
})

test_that(".pcr_relative works", {
  vec <- .pcr_relative(ct1$c_myc)
  vec2 <- 2^-ct1$c_myc
  expect_equal(vec, vec2)
})

test_that(".pcr_rsquared works", {
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
  rr <- .pcr_rsquared(ct3$c_myc, amount)
  rr2 <- cor(ct3$c_myc, log10(amount))^2

  expect_equal(rr, rr2)
})
