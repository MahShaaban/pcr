context("test assessment functions")

# loading required libraries
library(readr)
library(dplyr)

# delta_ct method for efficiency calculations
# locate and read file
fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
ct3 <- read_csv(fl)

## create amounts/diluion variable
amount <- c(1, .5, .2, .1, .05, .02, .01)

# apply delta ct method
# for amplification efficiency calculations
## add a log10 amount/dilution variable
## normalize by subtracting a reference gene
## calculate error using sd of gene and control
## calculate intervals by adding/subtracting errors from normalized values
efficiency <- ct3 %>%
  mutate(log_amount = rep(log10(amount), each = 3)) %>%
  group_by(log_amount) %>%
  summarise(normalized = mean(c_myc) - mean(GAPDH),
            error = sqrt((sd(c_myc)^2) + (sd(GAPDH)^2)),
            lower = normalized - error,
            upper = normalized + error)

# start testing

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

  # testing problem

  c <- coef(lm(ct3$GAPDH ~ log_amount))

  # testing problem
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

  # testing problem

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
