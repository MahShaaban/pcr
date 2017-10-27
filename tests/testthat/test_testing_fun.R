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

  # test using t-test
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
