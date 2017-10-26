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


