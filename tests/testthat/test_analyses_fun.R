context("test analyses functions")

test_that("pcr_ddct works", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_ddct(ct1,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain')

  expect_s3_class(res, 'data.frame')
  expect_equal(ncol(res), 8)
  expect_equal(nrow(res), length(unique(group_var)))
  expect_equal(unique(res$gene), setdiff(names(ct1), 'GAPDH'))


  res2 <- pcr_ddct(ct1,
                   group_var = group_var,
                   reference_gene = 'GAPDH',
                   reference_group = 'brain',
                   mode = 'same_tube')

  expect_s3_class(res2, 'data.frame')
  expect_equal(ncol(res2), 8)
  expect_equal(nrow(res2), length(unique(group_var)))
  expect_equal(unique(res2$gene), setdiff(names(ct1), 'GAPDH'))

  expect_error(
    pcr_ddct(ct1,
             group_var = group_var,
             reference_gene = 'GAPDH',
             reference_group = 'brain',
             mode = 'same')
  )

  gg <- pcr_ddct(ct1,
                 group_var = group_var,
                 reference_gene = 'GAPDH',
                 reference_group = 'brain',
                 plot = TRUE)

  expect_identical(class(gg), c("gg", "ggplot"))
})

test_that("pcr_dct works", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_dct(ct1,
                 group_var = group_var,
                 reference_group = 'brain')

  expect_s3_class(res, 'data.frame')
  expect_equal(ncol(res), 7)
  expect_equal(nrow(res), length(unique(group_var))*length(unique(names(ct1))))
  expect_equal(unique(res$gene), names(ct1))

  res2 <- pcr_dct(ct1,
                  group_var = group_var,
                  reference_group = 'brain',
                  mode = 'same_tube')

  expect_s3_class(res2, 'data.frame')
  expect_equal(ncol(res2), 7)
  expect_equal(nrow(res2), length(unique(group_var))*length(unique(names(ct1))))
  expect_equal(unique(res2$gene), names(ct1))

  expect_error(
    pcr_dct(ct1,
            group_var = group_var,
            reference_group = 'brain',
            mode = 'same')
  )

  gg <- pcr_dct(ct1,
                group_var = group_var,
                reference_group = 'brain',
                plot = TRUE)

  expect_identical(class(gg), c("gg", "ggplot"))
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

  expect_s3_class(res, 'data.frame')
  expect_equal(ncol(res), 7)
  expect_equal(nrow(res), length(unique(group_var)))
  expect_equal(unique(res$gene), setdiff(names(ct1), 'GAPDH'))

  res2 <- pcr_curve(ct1,
                    group_var = group_var,
                    reference_gene = 'GAPDH',
                    reference_group = 'brain',
                    intercept = intercept,
                    slope = slope,
                    mode = 'same_tube')

  expect_s3_class(res2, 'data.frame')
  expect_equal(ncol(res2), 7)
  expect_equal(nrow(res2), length(unique(group_var)))
  expect_equal(unique(res2$gene), setdiff(names(ct1), 'GAPDH'))

  expect_error(
    pcr_curve(ct1,
              group_var = group_var,
              reference_gene = 'GAPDH',
              reference_group = 'brain',
              intercept = intercept,
              slope = slope,
              mode = 'same')
  )

  gg <- pcr_curve(ct1,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain',
                  intercept = intercept,
                  slope = slope,
                  plot = TRUE)

  expect_identical(class(gg), c("gg", "ggplot"))
})
