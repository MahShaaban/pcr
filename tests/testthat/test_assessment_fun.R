context("test assessment functions")

test_that("pcr_efficiency calculates the correct intercept and slope", {
  # make amount/dilution variable
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate the standard curve
  res <- pcr_efficiency(ct3,
                        amount = amount,
                        reference_gene = 'GAPDH')

  log_amount <- log10(amount)
  x <- with(ct3, c_myc - GAPDH)
  c <- coef(lm(x ~ log_amount))

  expect_equal(unname(c), unlist(res[, 2:3], use.names = FALSE))
})

test_that("pcr_efficiency retruns a plot", {
  # make amount/dilution variable
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate the standard curve
  gg <- pcr_efficiency(ct3,
                       amount = amount,
                       reference_gene = 'GAPDH',
                       plot = TRUE)
  expect_identical(class(gg), c("gg", "ggplot"))
})

test_that("pcr_standard calculates the correct intercept and slope", {
  # make amount/dilution variable
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate the standard curve
  res <- pcr_standard(ct3,
                      amount = amount)

  log_amount <- log10(amount)

  c <- coef(lm(ct3$c_myc ~ log_amount))

  #expect_equal(unlist(res[1, 2:3], use.names = FALSE), unname(c))

  c <- coef(lm(ct3$GAPDH ~ log_amount))

  #expect_equal(unlist(res[2, 2:3], use.names = FALSE), unname(c))
})

test_that("pcr_standard retruns a plot", {
  # make amount/dilution variable
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate the standard curve
  gg <- pcr_standard(ct3,
                     amount = amount,
                     plot = TRUE)

  expect_identical(class(gg), c("gg", "ggplot"))
})

test_that("pcr_assess calls the correct methods", {
  # default: standard_curve
  # make amount/dilution variable
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate the standard curve
  res <- pcr_standard(ct3,
                      amount = amount)

  log_amount <- log10(amount)

  c <- coef(lm(ct3$c_myc ~ log_amount))

  #expect_equal(unlist(res[1, 2:3], use.names = FALSE), unname(c))


  # method: efficiency
  # make amount/dilution variable
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate the standard curve
  res <- pcr_assess(ct3,
                    amount = amount,
                    reference_gene = 'GAPDH',
                    method = 'efficiency')

  log_amount <- log10(amount)
  x <- with(ct3, c_myc - GAPDH)
  c <- coef(lm(x ~ log_amount))

  expect_equal(unname(c), unlist(res[, 2:3], use.names = FALSE))
})
