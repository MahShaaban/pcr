# load required libraries
library(readr)
library(devtools)

# ct1
fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
if(file.exists(fl)) {
  ct1 <- read_csv(fl)
  use_data(ct1, overwrite = TRUE)
} else {
  message(paste("File", fl, "doesn't exist"))
}

# ct2
fl <- system.file('extdata', 'ct2.csv', package = 'pcr')
if(file.exists(fl)) {
  ct2 <- read_csv(fl)
  use_data(ct2, overwrite = TRUE)
} else {
  message(paste("File", fl, "doesn't exist"))
}

# ct3
fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
if(file.exists(fl)) {
  ct3 <- read_csv(fl)
  use_data(ct3, overwrite = TRUE)
} else {
  message(paste("File", fl, "doesn't exist"))
}

# ct4
fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
if(file.exists(fl)) {
  ct4 <- read_csv(fl)
  use_data(ct4, overwrite = TRUE)
} else {
  message(paste("File", fl, "doesn't exist"))
}

# ct5
fl <- system.file('extdata', 'ct5.csv', package = 'pcr')
if(file.exists(fl)) {
  ct5 <- read_csv(fl)
  use_data(ct5, overwrite = TRUE)
} else {
  message(paste("File", fl, "doesn't exist"))
}
