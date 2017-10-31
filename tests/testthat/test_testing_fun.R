context("test testing functions")

test_that("pcr_test runs the t.test correctly", {
  fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
  ct4 <- readr::read_csv(fl)

  # make group variable
  group <- rep(c('control', 'treatment'), each = 12)

  # test using t-test
  res <- pcr_test(ct4,
                  group_var = group,
                  reference_gene = 'ref',
                  reference_group = 'control',
                  test = 't.test')

  norm <- ct4$target - ct4$ref

  group <- relevel(factor(group), ref = 'treatment')

  tt <- tidy(t.test(norm ~ group))

  expect_equal(res$estimate, tt$estimate)
  expect_equal(res$p_value, tt$p.value)
  expect_equal(res$lower, tt$conf.low)
  expect_equal(res$upper, tt$conf.high)
})


test_that("pcr_test runs the wilcox.test correctly", {
  fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
  ct4 <- readr::read_csv(fl)

  # make group variable
  group <- rep(c('control', 'treatment'), each = 12)

  # test using wilcox.test
  res <- pcr_test(ct4,
                  group_var = group,
                  reference_gene = 'ref',
                  reference_group = 'control',
                  test = 'wilcox.test')

  norm <- ct4$target - ct4$ref

  group <- relevel(factor(group), ref = 'treatment')

  wt <- tidy(wilcox.test(norm ~ group, conf.int = TRUE))

  expect_equal(res$estimate, wt$estimate)
  expect_equal(res$p_value, wt$p.value)
  expect_equal(res$lower, wt$conf.low)
  expect_equal(res$upper, wt$conf.high)
})

test_that("pcr_test runs the lm correctly", {
  fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
  ct4 <- readr::read_csv(fl)

  # make group variable
  group <- rep(c('control', 'treatment'), each = 12)

  # test using lm
  res <- pcr_test(ct4,
                  group_var = group,
                  reference_gene = 'ref',
                  reference_group = 'control',
                  test = 'lm')

  norm <- ct4$target - ct4$ref

  group <- relevel(factor(group), ref = 'control')

  ll <- tidy(lm(norm ~ group))[-1,]

  conf_int <- tidy(confint(lm(norm ~ group)))[-1,]

  expect_equal(res$estimate, ll$estimate)
  expect_equal(res$p_value, ll$p.value)
  expect_equal(res$lower, conf_int$X2.5..)
  expect_equal(res$upper, conf_int$X97.5..)
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

  ll <- tidy(lm(norm ~ group))[-1,][c(3,1,2),]
  conf_int <- tidy(confint(lm(norm ~ group)))[-1,][c(3,1,2),]


  expect_equal(res$estimate, ll$estimate)
  expect_equal(res$p_value, ll$p.value)
  expect_equal(res$lower, conf_int$X2.5..)
  expect_equal(res$upper, conf_int$X97.5..)
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
  ll <- tidy(lm(norm ~ mm + 0))[-1,]
  conf_int <- confint(lm(norm ~ mm + 0))[-1,]

  expect_equal(res$estimate, ll$estimate)
  expect_equal(res$p_value, ll$p.value)
  expect_equal(res$lower, unname(conf_int[,1]))
  expect_equal(res$upper, unname(conf_int[,2]))
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
  ll <- tidy(lm(norm ~ mm + 0))[-1,]
  conf_int <- confint(lm(norm ~ mm + 0))[-1,]

  expect_equal(res$estimate, ll$estimate)
  expect_equal(res$p_value, ll$p.value)
  expect_equal(res$lower, unname(conf_int[,1]))
  expect_equal(res$upper, unname(conf_int[,2]))
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
  ll <- tidy(lm(norm ~ mm + 0))[-1,]
  conf_int <- confint(lm(norm ~ mm + 0))[-1,]

  expect_equal(res$estimate, ll$estimate)
  expect_equal(res$p_value, ll$p.value)
  expect_equal(res$lower, unname(conf_int[,1]))
  expect_equal(res$upper, unname(conf_int[,2]))
})
