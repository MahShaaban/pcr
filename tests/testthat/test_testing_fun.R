context("test testing functions")

test_that("pcr_test runs the t.test correctly", {
  fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
  ct4 <- readr::read_csv(fl)

  # make group variable
  group <- rep(c('control', 'treatment'), each = 12)

  # test using pcr_test
  res <- pcr_test(ct4,
                  group_var = group,
                  reference_gene = 'ref',
                  reference_group = 'control',
                  test = 't.test')

  # test using pcr_ttest
  res2 <- pcr_ttest(ct4,
                    group_var = group,
                    reference_gene = 'ref',
                    reference_group = 'control')

  expect_identical(res, res2)

  norm <- ct4$target - ct4$ref

  group <- relevel(factor(group), ref = 'treatment')

  tt <- t.test(norm ~ group)

  expect_equal(res$estimate, unname(tt$estimate[1] - tt$estimate[2]))
  expect_equal(res$p_value, tt$p.value)
  expect_equal(res$lower, tt$conf.int[1])
  expect_equal(res$upper, tt$conf.int[2])
})


test_that("pcr_test runs the wilcox.test correctly", {
  fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
  ct4 <- readr::read_csv(fl)

  # make group variable
  group <- rep(c('control', 'treatment'), each = 12)

  # test using pcr_test
  res <- pcr_test(ct4,
                  group_var = group,
                  reference_gene = 'ref',
                  reference_group = 'control',
                  test = 'wilcox.test')

  # test using pcr_wilcox
  res2 <- pcr_wilcox(ct4,
                     group_var = group,
                     reference_gene = 'ref',
                     reference_group = 'control' )

  expect_identical(res, res2)

  norm <- ct4$target - ct4$ref

  group <- relevel(factor(group), ref = 'treatment')

  wt <- wilcox.test(norm ~ group, conf.int = TRUE)

  expect_equal(res$estimate, unname( wt$estimate[1] - wt$estimate[2]))
  expect_equal(res$p_value, wt$p.value)
  expect_equal(res$lower, wt$conf.int[1])
  expect_equal(res$upper, wt$conf.int[2])
})

test_that("pcr_test runs the lm correctly", {
  fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
  ct4 <- readr::read_csv(fl)

  # make group variable
  group <- rep(c('control', 'treatment'), each = 12)


  # test using pcr_test
  res <- pcr_test(ct4,
                  group_var = group,
                  reference_gene = 'ref',
                  reference_group = 'control',
                  test = 'lm')

  # test using pcr_test
  res2 <- pcr_lm(ct4,
                 group_var = group,
                 reference_gene = 'ref',
                 reference_group = 'control')

  expect_identical(res, res2)

  norm <- ct4$target - ct4$ref

  group <- relevel(factor(group), ref = 'control')

  ll <- lm(norm ~ group)

  conf_int <- confint(lm(norm ~ group))

  expect_equal(res$estimate, unname(ll$coefficients[2]))
  expect_equal(res$lower, conf_int[-1, 1])
  expect_equal(res$upper, conf_int[-1, 2])
})

test_that("pcr_test runs the lm correctly with multiple groups", {
  fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
  ct4 <- readr::read_csv(fl)

  # make group variable
  group <- rep(c('control', 'treatment'), each = 12)

  # subset treatment data
  ct_treated <- ct4[ group == 'treatment',]

  # make dose variable
  dose <- rep(c(100, 80, 60, 40), each = 3)

  # test using lm
  res <- pcr_test(ct_treated,
                  group_var = as.character(dose),
                  reference_gene = 'ref',
                  reference_group = '40',
                  test = 'lm')

  norm <- ct_treated$target - ct_treated$ref

  group <- relevel(factor(dose), ref = '40')

  ll <- lm(norm ~ group)

  expect_equal(res$estimate, unname(ll$coefficients[c(4, 2, 3)]))
})

test_that("pcr_test runs the lm correctly with a model matrix", {
  fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
  ct4 <- readr::read_csv(fl)

  # make a model matrix
  group <- rep(c('control', 'treatment'), each = 12)
  group <- relevel(factor(group), ref = 'control')
  dose <- rep(c(100, 80, 60, 40), each = 3, times = 2)
  mm <- model.matrix(~group+dose+group:dose, data = data.frame(group, dose))

  # test using lm
  res <- pcr_test(ct4,
                  reference_gene = 'ref',
                  model_matrix = mm,
                  test = 'lm')

  norm <- ct4$target - ct4$ref
  ll <- lm(norm ~ mm + 0)

  expect_equal(res$estimate, unname(ll$coefficients)[-1])
})

test_that("pcr_test runs the lm to adjust for separate runs", {
  fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
  ct4 <- readr::read_csv(fl)

  # make a model matrix
  group <- rep(c('control', 'treatment'), each = 12)
  group <- relevel(factor(group), ref = 'control')
  set.seed(1234)
  run <- factor(rep(c(1:3), 8))
  mm <- model.matrix(~group + group:run, data = data.frame(group, run))

  # test using lm
  res <- pcr_test(ct4,
                  reference_gene = 'ref',
                  model_matrix = mm,
                  test = 'lm')

  norm <- ct4$target - ct4$ref
  ll <- lm(norm ~ mm + 0)

  expect_equal(res$estimate, unname(ll$coefficients)[-1])
})

test_that("pcr_test runs the lm to adjust for rna quality", {
  fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
  ct4 <- readr::read_csv(fl)

  # make a model matrix
  group <- rep(c('control', 'treatment'), each = 12)
  group <- relevel(factor(group), ref = 'control')
  set.seed(1234)
  quality <- scale(rnorm(n = 24, mean = 1.9, sd = .1))
  mm <- model.matrix(~group + group:quality, data = data.frame(group, quality))

  # test using lm
  res <- pcr_test(ct4,
                  reference_gene = 'ref',
                  model_matrix = mm,
                  test = 'lm')

  norm <- ct4$target - ct4$ref
  ll <- lm(norm ~ mm + 0)

  expect_equal(res$estimate, unname(ll$coefficients)[-1])
})
