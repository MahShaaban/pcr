###############################################################################
# File name: standard_curve.R
# Description: contains the code to load and process ct raw values using the
#   standard curve method
# Details: loads the raw data and make it available within the package for the
# puroposes of testing and documentations. When this file excuted withing the
# package it should do the following:
#   + regenerate the dilution experiment ct values 'ct3.csv' and an object with
#     from averages and the same name from the averages and standard deviations
#     provided in ct3_ave.csv' and 'ct3_sd.csv'
#   + calculate the standard curve (slope and intercept) using 'ct3.csv'
#   + apply the standard curve method/ average amounts mode on 'ct1.csv' and
#     produce the 'standarized1' object
#   + apply the standard curve method/ average normalized mode on 'ct2.csv' and
#     produce the 'standarized2' object
# Author: Mahmoud Ahmed <mahmoud.s.fahmy@students.kasralainy.edu.eg>
###############################################################################

# load libraries
library(readr)
library(dplyr)
library(tidyr)
library(devtools)

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
generated_dat <- with(generated_dat, split(ct, gene)) %>%
  bind_cols() %>%
  mutate_all(function(x) as.numeric(x)) %>%
  arrange(-row_number())

write_csv(generated_dat, path = 'inst/extdata/ct3.csv')


## calculate intercept and slopes
## locate and read data
fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
ct3 <- read_csv(fl)

use_data(ct3, overwrite = TRUE)

### make a log amount variable
log_amount <- rep(log10(amount), each = 3)

### calculated the linear trend and coeffecients
### these are calculated for each gene separately
ll <- lm(ct3$c_myc ~ log_amount)
intercept1 <- ll$coefficients[1]
slope1 <- ll$coefficients[2]

ll <- lm(ct3$GAPDH ~ log_amount)
intercept2 <- ll$coefficients[1]
slope2 <- ll$coefficients[2]

# applying the standard curve method/ average amounts mode
## locate and read file
fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
ct1 <- read_csv(fl)

# make a grouping variable
group <- rep(c('brain', 'kidney'), each = 6)

# calculate the amounts of rna in samples
## calculate the input amounts 10 ^ ((ct - b) / m)
## average the amounts
## normalize by division on referece gene
## caliberate by division on reference group
df1 <- ct1 %>%
  mutate(c_myc = 10 ^ ((c_myc - intercept1) / slope1),
         GAPDH = 10 ^ ((GAPDH - intercept2) / slope2)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise_all(function(x) mean(x)) %>%
  group_by(group) %>%
  summarise(c_myc = c_myc / GAPDH) %>%
  gather(gene, normalized, -group) %>%
  mutate(caliberated = normalized / .$normalized[1])

## calculate standard deviations
## calculate cv
df2 <- ct1 %>%
  mutate(c_myc = 10 ^ ((c_myc - intercept1) / slope1),
         GAPDH = 10 ^ ((GAPDH - intercept2) / slope2)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise_all(function(x) sd(x)/mean(x)) %>%
  group_by(group) %>%
  summarise(c_myc = sqrt((c_myc^2) + (GAPDH^2))) %>%
  gather(gene, error, -group)


## join data.frames
## calculate intervals
## modify error = error(cv) * normalized
standarized1 <- full_join(df1, df2) %>%
  mutate(lower = caliberated - error,
         upper = caliberated + error,
         error = error * normalized)

use_data(standarized1, overwrite = TRUE)


# applying the standard curve method/ average normalized mode
# locate and read files
fl <- system.file('extdata', 'ct2.csv', package = 'pcr')
ct2 <- read_csv(fl)

# calculate the amounts of rna in samples
## calculate the input amounts 10 ^ ((ct - b) / m)
## average the amounts
## normalize by division on referece gene
## average the normalized amounts
## caliberate by division on reference group
df1 <- ct2 %>%
  mutate(c_myc = 10 ^ ((c_myc - intercept1) / slope1),
         GAPDH = 10 ^ ((GAPDH - intercept2) / slope2)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise(c_myc = mean(c_myc / GAPDH)) %>%
  gather(gene, normalized, -group) %>%
  mutate(caliberated = normalized / .$normalized[1])

# calculate cv as mean/sd
df2 <- ct2 %>%
  mutate(c_myc = 10 ^ ((c_myc - intercept1) / slope1),
         GAPDH = 10 ^ ((GAPDH - intercept2) / slope2)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise(c_myc = sd(c_myc / GAPDH) / mean(c_myc / GAPDH)) %>%
  gather(gene, error, -group)

## join data.frames
## calculate intervals
## modify error = error(cv) * normalized
standarized2 <- full_join(df1, df2) %>%
  mutate(lower = caliberated - error,
         upper = caliberated + error,
         error = error * normalized)
use_data(standarized2, overwrite = TRUE)
