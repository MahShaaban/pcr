# load required libraries
library(readr)
library(dplyr)
library(devtools)

# locate and read file
fl <- system.file('extdata', 'pcr2_ct.csv', package = 'pcr')
pcr2_ct <- read_csv(fl)

devtools::use_data(pcr2_ct, overwrite = TRUE)

# normalize as same wells
ave <- pcr2_ct %>%
  mutate(dct = c_myc - GAPDH) %>%
  select(-GAPDH) %>%
  mutate(group = rep(c('brain', 'kidney'), each = 6)) %>%
  group_by(group) %>%
  summarise(ave = mean(dct)) %>%
  mutate(ddct = ave - .$ave[1]) %>%
  mutate(norm_rel = 2 ^ -ddct)

## calculate errors
errors <- pcr2_ct %>%
  mutate(dct = c_myc - GAPDH) %>%
  select(-GAPDH) %>%
  mutate(group = rep(c('brain', 'kidney'), each = 6)) %>%
  group_by(group) %>%
  summarise(error = sd(dct))

## propagate errors
pcr2_norm <- left_join(ave, errors) %>%
  mutate(int_lower = 2 ^ -(ddct + error),
         int_upper = 2 ^ -(ddct - error))
use_data(pcr2_norm, overwrite = TRUE)
