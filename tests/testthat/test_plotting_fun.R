context("testing plotting functions")

test_that("pcr_plot returns a ggplot", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_ddct(ct1,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain')
  gg <- pcr_plot(res, method = 'delta_delta_ct')

  expect_identical(class(gg), c("gg", "ggplot"))
})

test_that("pcr_plot returns a dodge columns plot", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_ddct(ct1,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain')

  res <- bind_rows(res,
                   mutate(res, gene = 'c_myc2'))

  gg <- pcr_plot(res, method = 'delta_delta_ct')

  expect_identical(class(gg), c("gg", "ggplot"))
})

test_that("pcr_plot returns a faceted plot", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_ddct(ct1,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain')

  res <- bind_rows(res,
                   mutate(res, gene = 'c_myc2'))

  gg <- pcr_plot(res, method = 'delta_delta_ct', facets = TRUE)

  expect_identical(class(gg), c("gg", "ggplot"))
})
