###############################################################################
# File name: comparative.R
# Description: contains the code to load and process ct raw values using
#   different comparative methods.
# Details: loads the raw data and make it available within the package for the
# puroposes of testing and documentations. When this file excuted withing the
# package it should do the following:
#   + apply comparative ct method in separate_tube mode on 'ct1.csv' and produce
#     'rel_express1' object; a data.frame containing calculated values
#   + apply comparative ct method in same_tube mode on 'ct2.csv' and
#     produce 'rel_express2' object; a data.frame containing calculated values
#   + apply the delat_ct method to calculate efficiency in a dilution experiment
#     'ct3.csv' and produces 'efficiency'
#   + use the housekeeping gene data from 'ct1.csv' to apply the delta_ct method
#     to calculate the fold change of the house keeping genes 'hk_express'
# Author: Mahmoud Ahmed <mahmoud.s.fahmy@students.kasralainy.edu.eg>
###############################################################################

# load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(devtools)

# double delta method / separate_tube mode
## locate and read file
fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
ct1 <- read_csv(fl)

use_data(ct1, overwrite = TRUE)

## add grouping variable
ct1$group <- rep(c('brain', 'kidney'), each = 6)

## average the ct values
## normalize by substracting the referece gene
## calibrate by subtracting the reference group
## calculate relative expression
df1 <- ct1 %>%
  group_by(group) %>%
  summarise_all(function(x) mean(x)) %>%
  group_by(group) %>%
  summarise(c_myc = c_myc - GAPDH) %>%
  gather(gene, normalized, -group) %>%
  mutate(calibrated = normalized - .$normalized[1]) %>%
  mutate(relative_expression = 2 ^ -calibrated)

## calculate standard deviations
## calculate the error with the reference gene
df2 <- ct1 %>%
  group_by(group) %>%
  summarise_all(function(x) sd(x)) %>%
  group_by(group) %>%
  summarise(c_myc = sqrt((c_myc^2) + (GAPDH^2))) %>%
  gather(gene, error, -group)

## join the two data.frames
## calculate the intervals
rel_express1 <- left_join(df1, df2) %>%
  mutate(lower = 2 ^ -(calibrated + error),
         upper = 2 ^ -(calibrated - error))

use_data(rel_express1, overwrite = TRUE)



# # double delta method / same_tube mode
# locate and read file
fl <- system.file('extdata', 'ct2.csv', package = 'pcr')
ct2 <- read_csv(fl)

devtools::use_data(ct2, overwrite = TRUE)

# add grouping variable
group <- rep(c('brain', 'kidney'), each = 6)

# apply the comparative double delta ct method
## normalize by substracting the referece gene
## average the normalized values
## calibrate by subtracting the reference group
## calculate relative expression
df1 <- data_frame(c_myc = with(ct2, c_myc - GAPDH)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise(c_myc = mean(c_myc)) %>%
  gather(gene, normalized, -group) %>%
  mutate(calibrated = normalized - .$normalized[1]) %>%
  mutate(relative_expression = 2 ^ -calibrated)

## calculate standard deviations
df2 <- data_frame(c_myc = with(ct2, c_myc - GAPDH)) %>%
  mutate(group = group) %>%
  group_by(group) %>%
  summarise(c_myc = sd(c_myc)) %>%
  gather(gene, error, -group)

## join the two data.frames
## calculate the intervals
rel_express2 <- left_join(df1, df2) %>%
  mutate(lower = 2 ^ -(calibrated + error),
         upper = 2 ^ -(calibrated - error))

use_data(rel_express2, overwrite = TRUE)

# make a data frame of identical columns
# locate and read file
fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
ct1 <- read_csv(fl)

houskeeping <- data_frame(
  GAPDH1 = ct1$GAPDH,
  GAPDH2 = ct1$GAPDH
)

# add grouping variable
houskeeping$group <- rep(c('brain', 'kidney'), each = 6)

# apply the comparative delta ct method
## average the ct values
## calibrate by subtracting the reference group
df1 <- group_by(houskeeping, group) %>%
  summarise_all(function(x) mean(x)) %>%
  mutate_if(is.numeric, function(x) x - x[1]) %>%
  gather(gene, calibrated, -group) %>%
  mutate(fold_change = 2 ^ -calibrated)

# calculate the standard deviations
df2 <- group_by(houskeeping, group) %>%
  summarise_all(function(x) sd(x)) %>%
  gather(gene, error, -group)

## join the two data.frames
## calculate the intervals
hk_express <- full_join(df1, df2) %>%
  mutate(lower = 2 ^ -(calibrated + error),
         upper = 2 ^ -(calibrated - error))

use_data(hk_express, overwrite = TRUE)


# delta_ct method for efficiency calculations
# locate and read file
fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
ct3 <- read_csv(fl)

## create amounts/diluion variable
amount <- c(1, .5, .2, .1, .05, .02, .01)

# apply delta ct method
# for amplification efficiency calculations
## add a log10 amount/dilution variable
## normalize by subtracting a reference gene
## calculate error using sd of gene and control
## calculate intervals by adding/subtracting errors from normalized values
efficiency <- ct3 %>%
  mutate(log_amount = rep(log10(amount), each = 3)) %>%
  group_by(log_amount) %>%
  summarise(normalized = mean(c_myc) - mean(GAPDH),
            error = sqrt((sd(c_myc)^2) + (sd(GAPDH)^2)),
            lower = normalized - error,
            upper = normalized + error)

use_data(efficiency, overwrite = TRUE)
