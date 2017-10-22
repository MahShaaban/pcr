# load required libraries
library(readr)
library(dplyr)
library(devtools)

# locate and read file
fl <- system.file('extdata', 'pcr2_ct.csv', package = 'pcr')
pcr2_ct <- read_csv(fl)

devtools::use_data(pcr2_ct, overwrite = TRUE)

# normalize as same wells
norm_rel <- pcr2_ct %>%
  mutate(norm = c_myc - GAPDH) %>%
  select(-GAPDH) %>%
  mutate(group = rep(c('brain', 'kidney'), each = 6)) %>%
  group_by(group) %>%
  summarise(ave = mean(norm)) %>%
  mutate(norm_calib = ave - .$ave[1]) %>%
  mutate(norm_rel = 2 ^ -norm_calib)

## calculate errors
errors <- pcr2_ct %>%
  mutate(norm = c_myc - GAPDH) %>%
  select(-GAPDH) %>%
  mutate(group = rep(c('brain', 'kidney'), each = 6)) %>%
  group_by(group) %>%
  summarise(error = sd(norm))

## propagate errors
pcr2_norm <- left_join(norm_rel, errors) %>%
  mutate(int_lower = 2 ^ -(norm_calib + error),
         int_upper = 2 ^ -(norm_calib - error))
use_data(pcr2_norm, overwrite = TRUE)
