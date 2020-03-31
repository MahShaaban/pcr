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

test_that('pcr_analyze work for more than two genes', {
  set.seed(123)
  ct <- data.frame(
    geneA = rnorm(9, 30, 1),
    geneC = rnorm(9, 25, 1),
    geneB = rnorm(9, 20, 1)
  )
  var <- rep(c('group1', 'group2', 'group3'), 3)

  res <- pcr_analyze(ct,
                     group_var = var,
                     reference_gene = 'geneA',
                     reference_group = 'group1')
  a <- pcr_analyze(ct[, c('geneA', 'geneC')],
                   group_var = var,
                   reference_gene = 'geneA',
                   reference_group = 'group1')
  b <- pcr_analyze(ct[, c('geneA', 'geneB')],
                   group_var = var,
                   reference_gene = 'geneA',
                   reference_group = 'group1')

  expect_equal(rbind(a, b), res)
})

test_that('pcr_analyze respects order', {
  # create a data.frame
  set.seed(123)
  ct <- data.frame(
    geneA = rnorm(9, 30, 1),
    geneB = rnorm(9, 25, 1),
    geneC = rnorm(9, 20, 1)
  )

  # create a grouping variable
  var <- rep(c('group1', 'group2', 'group3'), each = 3)

  # run pcr_analyze
  res1 <- pcr_analyze(ct,
                      group_var = var,
                      reference_gene = 'geneA',
                      reference_group = 'group1')

  # change the columns order and retun pcr_analyze
  ind <- sample(1:ncol(ct))
  ct <- ct[, ind]

  res2 <- pcr_analyze(ct,
                      group_var = var,
                      reference_gene = 'geneA',
                      reference_group = 'group1')

  expect_true(all(res1[order(res1$group, res1$gene),]==
                    res2[order(res2$group, res2$gene),]))

  # change the rows order and rerun pcr_analyze
  ind <- sample(1:nrow(ct))
  ct <- ct[ind,]
  var <- var[ind]

  res3 <- pcr_analyze(ct,
                      group_var = var,
                      reference_gene = 'geneA',
                      reference_group = 'group1')

  expect_true(all(res1[order(res1$group, res1$gene),]==
                    res3[order(res3$group, res3$gene),]))
})

