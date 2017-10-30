context("test analyses functions")

# load required libraries
library(readr)
library(dplyr)
library(tidyr)

# double delta method / separate_tube mode
## locate and read file
fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
ct1 <- read_csv(fl)

## add grouping variable
ct1$group <- rep(c('brain', 'kidney'), each = 6)

## average the ct values
## normalize by substracting the referece gene
## calibrate by subtracting the reference group
## calculate relative expression
df1 <- ct1 %>%
  group_by(group) %>%
  summarise_all(function(x) mean(x)) %>%
  group_by(group) %>%
  summarise(c_myc = c_myc - GAPDH) %>%
  gather(gene, normalized, -group) %>%
  mutate(calibrated = normalized - .$normalized[1]) %>%
  mutate(relative_expression = 2 ^ -calibrated)

## calculate standard deviations
## calculate the error with the reference gene
df2 <- ct1 %>%
  group_by(group) %>%
  summarise_all(function(x) sd(x)) %>%
  group_by(group) %>%
  summarise(c_myc = sqrt((c_myc^2) + (GAPDH^2))) %>%
  gather(gene, error, -group)

## join the two data.frames
## calculate the intervals
rel_express1 <- left_join(df1, df2) %>%
  mutate(lower = 2 ^ -(calibrated + error),
         upper = 2 ^ -(calibrated - error))


# # double delta method / same_tube mode
# locate and read file
fl <- system.file('extdata', 'ct2.csv', package = 'pcr')
ct2 <- read_csv(fl)

# add grouping variable
group <- rep(c('brain', 'kidney'), each = 6)

# apply the comparative double delta ct method
## normalize by substracting the referece gene
## average the normalized values
## calibrate by subtracting the reference group
## calculate relative expression
df1 <- data_frame(c_myc = with(ct2, c_myc - GAPDH)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise(c_myc = mean(c_myc)) %>%
  gather(gene, normalized, -group) %>%
  mutate(calibrated = normalized - .$normalized[1]) %>%
  mutate(relative_expression = 2 ^ -calibrated)

## calculate standard deviations
df2 <- data_frame(c_myc = with(ct2, c_myc - GAPDH)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise(c_myc = sd(c_myc)) %>%
  gather(gene, error, -group)

## join the two data.frames
## calculate the intervals
rel_express2 <- left_join(df1, df2) %>%
  mutate(lower = 2 ^ -(calibrated + error),
         upper = 2 ^ -(calibrated - error))

# make a data frame of identical columns
# locate and read file
fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
ct1 <- read_csv(fl)

houskeeping <- data_frame(
  GAPDH1 = ct1$GAPDH,
  GAPDH2 = ct1$GAPDH
)

# add grouping variable
houskeeping$group <- rep(c('brain', 'kidney'), each = 6)

# apply the comparative delta ct method
## average the ct values
## calibrate by subtracting the reference group
df1 <- group_by(houskeeping, group) %>%
  summarise_all(function(x) mean(x)) %>%
  mutate_if(is.numeric, function(x) x - x[1]) %>%
  gather(gene, calibrated, -group) %>%
  mutate(fold_change = 2 ^ -calibrated)

# calculate the standard deviations
df2 <- group_by(houskeeping, group) %>%
  summarise_all(function(x) sd(x)) %>%
  gather(gene, error, -group)

## join the two data.frames
## calculate the intervals
hk_express <- full_join(df1, df2) %>%
  mutate(lower = 2 ^ -(calibrated + error),
         upper = 2 ^ -(calibrated - error))

## calculate intercept and slopes
## locate and read data
fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
ct3 <- read_csv(fl)

### make a log amount variable
amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
log_amount <- log10(amount)

### calculated the linear trend and coeffecients
### these are calculated for each gene separately
ll <- lm(ct3$c_myc ~ log_amount)
intercept1 <- ll$coefficients[1]
slope1 <- ll$coefficients[2]

ll <- lm(ct3$GAPDH ~ log_amount)
intercept2 <- ll$coefficients[1]
slope2 <- ll$coefficients[2]

# applying the standard curve method/ average amounts mode
## locate and read file
fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
ct1 <- read_csv(fl)

# make a grouping variable
group <- rep(c('brain', 'kidney'), each = 6)

# calculate the amounts of rna in samples
## calculate the input amounts 10 ^ ((ct - b) / m)
## average the amounts
## normalize by division on referece gene
## calibrate by division on reference group
df1 <- ct1 %>%
  mutate(c_myc = 10 ^ ((c_myc - intercept1) / slope1),
         GAPDH = 10 ^ ((GAPDH - intercept2) / slope2)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise_all(function(x) mean(x)) %>%
  group_by(group) %>%
  summarise(c_myc = c_myc / GAPDH) %>%
  gather(gene, normalized, -group) %>%
  mutate(calibrated = normalized / .$normalized[1])

## calculate standard deviations
## calculate cv
df2 <- ct1 %>%
  mutate(c_myc = 10 ^ ((c_myc - intercept1) / slope1),
         GAPDH = 10 ^ ((GAPDH - intercept2) / slope2)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise_all(function(x) sd(x)/mean(x)) %>%
  group_by(group) %>%
  summarise(c_myc = sqrt((c_myc^2) + (GAPDH^2))) %>%
  gather(gene, error, -group)


## join data.frames
## calculate intervals
## modify error = error(cv) * normalized
standarized1 <- full_join(df1, df2) %>%
  mutate(lower = calibrated - error,
         upper = calibrated + error,
         error = error * normalized)

# applying the standard curve method/ average normalized mode
# locate and read files
fl <- system.file('extdata', 'ct2.csv', package = 'pcr')
ct2 <- read_csv(fl)

# calculate the amounts of rna in samples
## calculate the input amounts 10 ^ ((ct - b) / m)
## average the amounts
## normalize by division on referece gene
## average the normalized amounts
## calibrate by division on reference group
df1 <- ct2 %>%
  mutate(c_myc = 10 ^ ((c_myc - intercept1) / slope1),
         GAPDH = 10 ^ ((GAPDH - intercept2) / slope2)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise(c_myc = mean(c_myc / GAPDH)) %>%
  gather(gene, normalized, -group) %>%
  mutate(calibrated = normalized / .$normalized[1])

# calculate cv as mean/sd
df2 <- ct2 %>%
  mutate(c_myc = 10 ^ ((c_myc - intercept1) / slope1),
         GAPDH = 10 ^ ((GAPDH - intercept2) / slope2)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise(c_myc = sd(c_myc / GAPDH) / mean(c_myc / GAPDH)) %>%
  gather(gene, error, -group)

## join data.frames
## calculate intervals
## modify error = error(cv) * normalized
standarized2 <- full_join(df1, df2) %>%
  mutate(lower = calibrated - error,
         upper = calibrated + error,
         error = error * normalized)

# start testing

test_that("pcr_ddct in separate_tube mode", {
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate all values and errors in one step
  res <- pcr_ddct(ct1,
                  group_var = group_var,
                  reference_gene = 'GAPDH',
                  reference_group = 'brain')

  expect_identical(res, rel_express1)
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

  expect_identical(res, rel_express2)
})

test_that("pcr_dct calculates and returns the correct values", {
  # make a data.frame of two identical columns
   pcr_hk <- data.frame(
     GAPDH1 = ct1$GAPDH,
     GAPDH2 = ct1$GAPDH
     )

   ## add grouping variable
   group_var <- rep(c('brain', 'kidney'), each = 6)

   # calculate caliberation
   res <- pcr_dct(pcr_hk,
                  group_var = group_var,
                  reference_group = 'brain')
   expect_identical(res, hk_express)
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

  # testing problem
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
  # testing problem!!
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

  expect_identical(res, rel_express1)

  # make a data.frame of two identical columns
  pcr_hk <- data.frame(
    GAPDH1 = ct1$GAPDH,
    GAPDH2 = ct1$GAPDH
  )

  # method: delta_ct
  ## add grouping variable
  group_var <- rep(c('brain', 'kidney'), each = 6)

  # calculate caliberation
  res <- pcr_analyze(pcr_hk,
                      group_var = group_var,
                      reference_group = 'brain',
                      method = 'delta_ct')
  expect_identical(res, hk_express)

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
  # testing problem
  })


