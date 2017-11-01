# load libraries
library(readr)
library(dplyr)
library(tidyr)

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
