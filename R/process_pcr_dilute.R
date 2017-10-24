# load libraries
library(readr)
library(dplyr)
library(tidyr)
library(devtools)

# locate averages and sd files of dilution experiments
fl1 <- system.file('extdata', 'pcr_dilute_ave.csv', package = 'pcr')
fl2 <- system.file('extdata', 'pcr_dilute_error.csv', package = 'pcr')

# create amounts variable
amount <- c(1, .5, .2, .1, .05, .02, .01)

# read averages and sds and regenerate ct values
set.seed(1234)
generated_dat <- list(ave = read_csv(fl1),
                   sd = read_csv(fl2)) %>%
  bind_rows(.id = 'stat') %>%
  mutate(amount = rep(amount, 2)) %>%
  gather(gene, ct, -stat, -amount) %>%
  spread(stat, ct) %>%
  group_by(amount, gene) %>%
  mutate(ct = paste(rnorm(3, mean = ave, sd = sd), collapse = ',')) %>%
  transform(ct = strsplit(ct, ',')) %>%
  unnest()

# create a data.frame of generated ct values
pcr_dilute <- with(generated_dat, split(ct, gene)) %>%
  bind_cols() %>%
  mutate_all(function(x) as.numeric(x)) %>%
  arrange(-row_number())

write_csv(pcr_dilute, path = 'inst/extdata/pcr_dilute.csv')

# locate and read file
fl <- system.file('extdata', 'pcr_dilute.csv', package = 'pcr')
pcr_dilute <- read_csv(fl)
use_data(pcr_dilute, overwrite = TRUE)

# calculate dct and errors
pcr_effeciency <- pcr_dilute %>%
  mutate(log_amount = rep(log10(amount), each = 3)) %>%
  group_by(log_amount) %>%
  summarise(dct = mean(c_myc) - mean(GAPDH),
            error = sqrt((sd(c_myc)^2) + (sd(GAPDH)^2)),
            int_lower = dct - error,
            int_upper = dct + error)

use_data(pcr_effeciency, overwrite = TRUE)


# calculate standard curve
## calculate intercept and slopes
ll <- lm(pcr_dilute$c_myc ~ rep(log10(amount), each = 3))
intercept1 <- ll$coefficients[1]
slope1 <- ll$coefficients[2]

ll <- lm(pcr_dilute$GAPDH ~ rep(log10(amount), each = 3))
intercept2 <- ll$coefficients[1]
slope2 <- ll$coefficients[2]

# make a grouping variable
group <- rep(c('brain', 'kidney'), each = 6)

# calculate input amounts, averages and normalization
input_amounts <- pcr1_ct %>%
  mutate(c_myc = 10 ^ ((c_myc - intercept1) / slope1),
         GAPDH = 10 ^ ((GAPDH - intercept2) / slope2)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise_all(function(x) mean(x)) %>%
  mutate(norm = c_myc / GAPDH) %>%
  mutate(calib = norm / .$norm[1])

# calculate cv and errors
ave <- input_amounts[,1:3] %>%
  gather(gene, ave, -group)

error <- pcr1_ct %>%
  mutate(c_myc = 10 ^ ((c_myc - intercept1) / slope1),
         GAPDH = 10 ^ ((GAPDH - intercept2) / slope2)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise_all(function(x) sd(x)) %>%
  gather(gene, sd, -group) %>%
  full_join(ave) %>%
  group_by(group, gene) %>%
  summarise(cv = sd/ave) %>%
  spread(gene, cv) %>%
  mutate(cv = sqrt((c_myc^2) + (GAPDH^2))) %>%
  select(group, cv)

pcr_amounts <- full_join(input_amounts, error) %>%
  mutate(sd = norm * cv)

use_data(pcr_amounts, overwrite = TRUE)
