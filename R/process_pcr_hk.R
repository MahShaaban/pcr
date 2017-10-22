# load required libraries
library(dplyr)
library(tidyr)
library(readr)
library(devtools)

# make a data frame of identical columns
# locate and read file
fl <- system.file('extdata', 'pcr1_ct.csv', package = 'pcr')
pcr1_ct <- read_csv(fl)

pcr_hk <- data_frame(
  GAPDH1 = pcr1_ct$GAPDH,
  GAPDH2 = pcr1_ct$GAPDH
)
use_data(pcr_hk, overwrite = TRUE)

# add grouping variable
pcr_hk$group <- rep(c('brain', 'kidney'), each = 6)

# caliberate the two genes
calib <- group_by(pcr_hk, group) %>%
  summarise_all(function(x) mean(x)) %>%
  mutate_if(is.numeric, function(x) x - x[1]) %>%
  gather(gene, calib, -group) %>%
  mutate(rel = 2 ^ -calib)

# calculate errors
errors <- group_by(pcr_hk, group) %>%
  summarise_all(function(x) sd(x)) %>%
  gather(gene, sd, -group)

# calculate intervals
pcr_hk_calib <- full_join(calib, errors) %>%
  mutate(int_lower = 2 ^ -(calib + sd),
         int_upper = 2 ^ -(calib - sd))
use_data(pcr_hk_calib, overwrite = TRUE)
