context("test analyses functions")

test_that("pcr_ddct in separate_tube mode", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_ddct(ct1,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain')

  res2 <- ct1 %>%
    mutate(group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) mean(x))

  norm <- res2$c_myc - res2$GAPDH
  expect_equal(res$normalized, norm)

  calib <- norm[2] - norm[1]
  expect_equal(res$calibrated[2], calib)
})

test_that("pcr_ddct in same_tube mode", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_ddct(ct2,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain',
                  mode = 'same_tube')

  res2 <- data_frame(c_myc = ct2$c_myc - ct2$GAPDH) %>%
    mutate(group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) mean(x))

  expect_equal(res$normalized, res2$c_myc)

  calib <- res2$c_myc[2] - res2$c_myc[1]
  expect_equal(res$calibrated[2], calib)
})

test_that("pcr_dct calculates and returns the correct values", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # make a data.frame of two identical columns
   pcr_hk <- data.frame(
     GAPDH = ct1$GAPDH
     )

   ## add grouping variable
   group_var <- rep(c('brain', 'kidney'), each = 6)

   # calculate caliberation
   res <- pcr_dct(pcr_hk,
                  group_var = group_var,
                  reference_group = 'brain')
  res2 <- pcr_hk %>%
    mutate(group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) mean(x))

  calib <- res2$GAPDH[2] - res2$GAPDH[1]
  expect_equal(res$calibrated[2], calib)
})

test_that("pcr_curve in separate_tube mode", {
  # make a vector of RNA amounts
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate standard curve
  standard_curve <- pcr_assess(ct3,
                               amount = amount,
                               method = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate standard amounts and error
  res <- pcr_curve(ct1,
                   group_var = group_var,
                   reference_gene = 'GAPDH',
                   reference_group = 'brain',
                   intercept = intercept,
                   slope = slope)

  res2 <- ct1 %>%
    mutate(c_myc = 10 ^ ((c_myc - intercept[1]) / slope[1]),
           GAPDH = 10 ^ ((GAPDH - intercept[2]) / slope[2])) %>%
    mutate(group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) mean(x))

  norm <- res2$c_myc / res2$GAPDH
  expect_equal(res$normalized, norm)
})

test_that("pcr_curve in same_tube mode", {
  # make a vector of RNA amounts
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate standard curve
  standard_curve <- pcr_assess(ct3,
                               amount = amount,
                               method = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate standard amounts and error
  res <- pcr_curve(ct2,
                   group_var = group_var,
                   reference_gene = 'GAPDH',
                   reference_group = 'brain',
                   intercept = intercept,
                   slope = slope,
                   mode = 'same_tube')

  res2 <- ct2 %>%
    mutate(c_myc = 10 ^ ((c_myc - intercept[1]) / slope[1]),
           GAPDH = 10 ^ ((GAPDH - intercept[2]) / slope[2]))
  norm <- data_frame(group = group_var,
                     c_myc = res2$c_myc / res2$GAPDH) %>%
    group_by(group) %>%
    summarise(c_myc = mean(c_myc))
  expect_equal(res$normalized, norm$c_myc)
})

test_that("pcr_analyze calls the right methods", {
  # default method ddct
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_analyze(ct1,
                     group_var = group_var,
                     reference_gene = 'GAPDH',
                     reference_group = 'brain')

  res2 <- ct1 %>%
    mutate(group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) mean(x))

  norm <- res2$c_myc - res2$GAPDH
  expect_equal(res$normalized, norm)

  calib <- norm[2] - norm[1]
  expect_equal(res$calibrated[2], calib)

  # make a data.frame of two identical columns
  pcr_hk <- data.frame(
    GAPDH = ct1$GAPDH
  )

  # method: delta_ct
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate caliberation
  res <- pcr_analyze(pcr_hk,
                      group_var = group_var,
                      reference_group = 'brain',
                      method = 'delta_ct')
  res2 <- pcr_hk %>%
    mutate(group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) mean(x))

  calib <- res2$GAPDH[2] - res2$GAPDH[1]
  expect_equal(res$calibrated[2], calib)

  # method: relative_curve
  # make a vector of RNA amounts
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate standard curve
  standard_curve <- pcr_assess(ct3,
                               amount = amount,
                               method = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate standard amounts and error
  res <- pcr_analyze(ct1,
                    group_var = group_var,
                    reference_gene = 'GAPDH',
                    reference_group = 'brain',
                    intercept = intercept,
                    slope = slope,
                    method = 'relative_curve')

  res2 <- ct1 %>%
    mutate(c_myc = 10 ^ ((c_myc - intercept[1]) / slope[1]),
           GAPDH = 10 ^ ((GAPDH - intercept[2]) / slope[2])) %>%
    mutate(group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) mean(x))

  norm <- res2$c_myc / res2$GAPDH
  expect_equal(res$normalized, norm)
  })

test_that("pcr_ddct returns a plot", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_ddct(ct1,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain',
                  plot = TRUE)

  expect_identical(class(res), c("gg", "ggplot"))
})

test_that("pcr_dct returns a plot", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # make a data.frame of two identical columns
  pcr_hk <- data.frame(
    GAPDH = ct1$GAPDH
  )

  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate caliberation
  res <- pcr_dct(pcr_hk,
                 group_var = group_var,
                 reference_group = 'brain',
                 plot = TRUE)

  expect_identical(class(res), c("gg", "ggplot"))
})


test_that("pcr_curve returns a plot", {
  # make a vector of RNA amounts
  amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

  # calculate standard curve
  standard_curve <- pcr_assess(ct3,
                               amount = amount,
                               method = 'standard_curve')
  intercept <- standard_curve$intercept
  slope <- standard_curve$slope

  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate standard amounts and error
  res <- pcr_curve(ct2,
                   group_var = group_var,
                   reference_gene = 'GAPDH',
                   reference_group = 'brain',
                   intercept = intercept,
                   slope = slope,
                   mode = 'same_tube',
                   plot = TRUE)

  expect_identical(class(res), c("gg", "ggplot"))
})
