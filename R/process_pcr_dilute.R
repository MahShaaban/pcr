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
  mutate_all(function(x) as.numeric(x))

write_csv(pcr_dilute, path = 'inst/extdata/pcr_dilute.csv')
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
