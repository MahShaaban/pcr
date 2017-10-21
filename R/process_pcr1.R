# load required libraries
library(readr)
library(dplyr)

# locate and read file
fl <- system.file('extdata', 'pcr1_ct.csv', package = 'pcr')
pcr1_ct <- read_csv(fl)

devtools::use_data(pcr1_ct, overwrite = TRUE)

# add grouping variable
pcr1_ct$group <- rep(c('brain', 'kidney'), each = 6)

# normalize as separate wells
## calculate averages, ddct and normalized relative values
ave <- pcr1_ct %>%
  group_by(group) %>%
  summarise_all(function(x) mean(x)) %>%
  mutate(dct = c_myc - GAPDH) %>%
  mutate(ddct = dct - .$dct[1]) %>%
  mutate(norm_rel = 2 ^ -ddct)

## calculate errors
errors <- pcr1_ct %>%
  group_by(group) %>%
  summarise_all(function(x) sd(x)) %>%
  group_by(group) %>%
  summarise(error = sqrt((c_myc^2) + (GAPDH^2)))

## propagate errors
pcr1_norm <- left_join(ave, errors) %>%
  mutate(int_lower = 2 ^ -(ddct + error),
         int_upper = 2 ^ -(ddct - error))
devtools::use_data(pcr1_norm, overwrite = TRUE)
