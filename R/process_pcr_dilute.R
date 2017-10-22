# load libraries
library(readr)
library(dplyr)
library(devtools)

# locate and read file
fl <- system.file('extdata', 'pcr_dilute_ave.csv', package = 'pcr')
pcr_dilute_ave <- read_csv(fl)

fl <- system.file('extdata', 'pcr_dilute_error.csv', package = 'pcr')
pcr_dilute_error <- read_csv(fl)

use_data(pcr_dilute_ave, overwrite = TRUE)
use_data(pcr_dilute_error, overwrite = TRUE)

# diluations
amount <- c(1, .5, .2, .1, .05, .02, .01)

# calculate dct and errors
pcr_effeciency <- data_frame(log_amount = log10(amount),
                             dct = pcr_dilute_ave$c_myc - pcr_dilute_ave$GAPDH,
                             error = sqrt((pcr_dilute_error$c_myc^2) + (pcr_dilute_error$GAPDH^2)),
                             int_lower = dct - error,
                             int_upper = dct + error)
use_data(pcr_effeciency, overwrite = TRUE)
