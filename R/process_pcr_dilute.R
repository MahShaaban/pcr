# load libraries
library(readr)
library(dplyr)
library(devtools)

# locate and read file
fl <- system.file('extdata', 'pcr_dilute.csv', package = 'pcr')
pcr_dilute <- read_csv(fl)

use_data(pcr_dilute, overwrite = TRUE)

# calculate delta ct and errors
pcr_effeciency <- pcr_dilute %>%
  mutate(log_amount = log10(amount),
         dct = ave_cmyc - ave_GAPDH,
         error = sqrt((sd_cmyc^2) + (sd_GAPDH^2)),
         int_lower = dct - error,
         int_upper = dct + error)

use_data(pcr_effeciency, overwrite = TRUE)
