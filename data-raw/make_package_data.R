# load libraries
library(readr)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)

# ct3
# regenerate the dilution experiment ct values
## locate averages and sd files of dilution experiments
fl1 <- system.file('extdata', 'ct3_ave.csv', package = 'pcr')
fl2 <- system.file('extdata', 'ct3_sd.csv', package = 'pcr')

## create amounts/diluion variable
amount <- c(1, .5, .2, .1, .05, .02, .01)

## generate random numbers around ave and sd
### read averages and sds and regenerate ct values in a list
### bind list items in a tidy data.frame
### generate ct values using rnorm on ave and sd
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

### create a data.frame of generated numeric ct values
ct3 <- with(generated_dat, split(ct, gene)) %>%
  bind_cols() %>%
  mutate_all(function(x) as.numeric(x)) %>%
  arrange(-row_number())

### write csve and use data
write_csv(ct3, path = 'inst/extdata/ct3.csv')

# ct4
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

# write csv of only ct values and use data
write_csv(ct4, path = 'inst/extdata/ct4.csv')

# ct5
# locate and read gnorm_example.xlsx
fl <- system.file('extdata', 'genorm_example.xlsx', package = 'pcr')

raw_data <- read_xlsx(fl, sheet = 1)
ct5 <- raw_data %>%
  fill(Target) %>%
  select(Target, Sample, `Cq value`) %>%
  setNames(c('gene', 'group', 'ct')) %>%
  group_by(gene, group) %>%
  mutate(row_id = row_number()) %>%
  spread(gene, ct) %>%
  ungroup() %>%
  select(-row_id, -group)

write_csv(ct5, path = 'inst/extdata/ct5.csv')

#
# https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2002-3-7-research0034/MediaObjects/13059_2001_453_MOESM1_ESM.txt
fl <- system.file('extdata', 'hk_tissues.txt', package = 'pcr')
raw_data <- read_tsv(fl, skip = 1, n_max = 85)
hk_tissue <- select(raw_data, -1)
hk_tissue$tissue <- str_replace_all(raw_data$X1, '\\d+', '')
write_csv(hk_tissue, path = 'inst/extdata/hk_tissue.csv')
