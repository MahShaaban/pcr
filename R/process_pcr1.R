# load required libraries
library(readr)
library(dplyr)
library(devtools)

# locate and read file
fl <- system.file('extdata', 'pcr1_ct.csv', package = 'pcr')
pcr1_ct <- read_csv(fl)

use_data(pcr1_ct, overwrite = TRUE)

# add grouping variable
pcr1_ct$group <- rep(c('brain', 'kidney'), each = 6)

# normalize as separate wells
## calculate averages, ddct and normalized relative values
norm_rel <- pcr1_ct %>%
  group_by(group) %>%
  summarise_all(function(x) mean(x)) %>%
  mutate(norm = c_myc - GAPDH) %>%
  mutate(norm_calib = norm - .$norm[1]) %>%
  mutate(norm_rel = 2 ^ -norm_calib)

## calculate errors
errors <- pcr1_ct %>%
  group_by(group) %>%
  summarise_all(function(x) sd(x)) %>%
  group_by(group) %>%
  summarise(error = sqrt((c_myc^2) + (GAPDH^2)))

## propagate errors
pcr1_norm <- left_join(norm_rel, errors) %>%
  mutate(int_lower = 2 ^ -(norm_calib + error),
         int_upper = 2 ^ -(norm_calib - error))
use_data(pcr1_norm, overwrite = TRUE)
