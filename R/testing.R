# load libraries
library(readr)
library(dplyr)
library(devtools)

# locate and read file
fl <- system.file('extdata', 'ct4_raw.tsv', package = 'pcr')
ct4 <- read_tsv(fl, col_names = c('concentraion', 'ct', 'group'), n_max = 48)

# make groups
groups <- data_frame(
  group = 1:4,
  gene_type = rep(c('target', 'ref'), each = 2),
  sample_type = rep(c('control', 'treatment'), 2)
)

# join ct values and groups
ct4 <- full_join(ct4, groups)

# reshape the ct values
ct4 <- with(ct4, split(ct, gene_type)) %>% bind_rows

# write csv of only ct values
write_csv(ct4, path = 'inst/extdata/ct4.csv')

use_data(ct4, overwrite = TRUE)
