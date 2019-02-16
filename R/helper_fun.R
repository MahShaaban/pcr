#' Average values by a variable
#'
#' Uses a group_by statement to average values in a data.fram by a variable
#'
#' @param df A data.frame
#' @param group_var A vector or lenght equals the number or rows of df
#' @param amount A vector or lenght equals the number or rows of df
#' @param tidy A logical, default FALSE. When TRUE returns a tidy data.frame
#'
#' @details Used to average ct or input amounts (or their averages) by a
#' grouping variable; group_var for experimental groups or amount for serial
#' dilutions.
#'
#' @return A data.frame with a column for the grouping variable; group_var or
#' amount and one for each of the original columns. When tidy is TRUE the
#' returns a tidy data.frame with 3 columns; group/amount, gene and average.
#'
#' @examples
#' # using a group_var variabale
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate averages
#' .pcr_average(ct1, group_var = group_var)
#'
#' # calculate averages and return a tidy data.frame
#' .pcr_average(ct1, group_var = group_var, tidy = TRUE)
#'
#' # using a amount variable
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
#'
#' # add amount variable
#' amount <- amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # calculate averages
#' .pcr_average(ct3, amount = amount)
#'
#' # calculate averages and return a tidy data.frame
#' .pcr_average(ct3, amount = amount, tidy = TRUE)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate group_by summarise_all
#' @importFrom tidyr gather
#'
#' @keywords internal
#'
#' @export
.pcr_average <- function(df, group_var, amount, tidy = FALSE) {
  # when group_var is provided
  # calculate the averages using group_var in group_by
  if(!missing(group_var)) {
    ave <- mutate(df, group = group_var) %>%
      group_by(group) %>%
      summarise_all(function(x) mean(x))

    # when tidy is TRUE return a tidy data.frame
    if(tidy == TRUE) {
      ave <- ave %>%
        gather(gene, average, -group)
    }
  }

  # when group_var is not provided and amount is
  # calculate the averages using amount in group_by
  if(missing(group_var)) {
    ave <- mutate(df, amount = amount) %>%
      group_by(amount) %>%
      summarise_all(function(x) mean(x))

    # when tidy is TRUE return a tidy data.frame
    if(tidy == TRUE) {
      ave <- ave %>%
        gather(gene, average, -amount)
    }
  }

  return(ave)
}

#' Normalize values by a column
#'
#' Uses subtraction or division to normalize values in all columns to a certain
#' specified column
#'
#' @inheritParams .pcr_average
#' @param reference_gene A character string of the name of the column
#' corresponding to the reference gene
#' @param mode A character string of the normalization mode to be used. Default
#' is 'subtract'. Other possible modes include 'divide'
#'
#' @details Used to normalize ct or input amounts (or their averages) by a
#' a reference_gene/column
#'
#' @return A data.frame with a column for each of the original columns after
#' subtraction or division to a reference_gene/column which is dropped.
#' The function returns ignores non numeric columns. When tidy is TRUE the
#' returns a tidy data.frame with the columns: gene and average as
#' well as any non numeric columns such as a grouping variable group/amount.
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
#'
#' # normalize the ct values
#' .pcr_normalize(ct1, reference_gene = 'GAPDH')
#'
#' # normalize by division
#' .pcr_normalize(ct1, reference_gene = 'GAPDH', mode = 'divide')
#'
#' # add grouping variable and average first
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#' ave <- .pcr_average(ct1, group_var = group_var)
#'
#' # normalize by subtraction
#' .pcr_normalize(ave, 'GAPDH')
#'
#' # normalize by division
#' .pcr_normalize(ave, 'GAPDH', mode = 'divide')
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate_if starts_with
#'
#' @keywords internal
#'
#' @export
.pcr_normalize <- function(df, reference_gene, mode = 'subtract',
                           tidy = FALSE) {
  # get the reference_gene column and unlist
  ref <- select(df, reference_gene) %>% unlist(use.names = FALSE)

  if(mode == 'subtract') {
    # drop the reference_gene columns
    # ignore non numeric columns
    # subtract reference_gene from all other columns
    norm <- select(df, -starts_with(reference_gene)) %>%
      mutate_if(is.numeric, function(x) x - ref)
  } else if(mode == 'divide') {
    # drop the reference_gene columns
    # ignore non numeric columns
    # divide all other columns by reference_gene
    norm <- select(df, -starts_with(reference_gene)) %>%
      mutate_if(is.numeric, function(x) x / ref)
  }

  # retrun a tidy data.frame when tidy == TRUE
  if(tidy == TRUE) {
    norm <- gather(norm, gene, normalized, -group)
  }
  return(norm)
}

#' Calibrate values by a row
#'
#' Uses subtraction or division to caliberate values in all rows to a sepcified
#' row
#'
#' @inheritParams .pcr_average
#' @inheritParams .pcr_normalize
#' @param reference_group A character string of the the entery in the rows of a
#' grouping variable
#'
#' @details Used to calibrate average ct or input amounts by a reference_group/row
#'
#' @return A data.frame of the same dimensions after subtracting or dividing by
#' a reference_group/row. The function returns ignores non numeric columns.
#' When tidy is TRUE returns a tidy data.frame with the columns: gene and
#' calibrated as well as any non numeric columns such as a grouping variable
#' group/amount.
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate averages
#' ave <- .pcr_average(ct1, group_var = group_var)
#'
#' # calculate delta ct
#' dct <- .pcr_normalize(ave, 'GAPDH')
#'
#' # calculate delta delta ct
#' .pcr_calibrate(dct, 'brain', tidy = TRUE)
#'
#' # calculate delta delta ct and return a tidy data.frame
#' .pcr_calibrate(dct, 'brain', tidy = TRUE)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate_if
#' @importFrom tidyr gather
#'
#' @keywords internal
#'
#' @export
.pcr_calibrate <- function(df, reference_group, mode = 'subtract',
                           tidy = FALSE) {
  # get the row index of the reference group
  ind <- which(df$group == reference_group)

  if(mode == 'subtract') {
    # ignore non numeric columns
    # subtract reference_group from all other rows
    calib <- mutate_if(df, is.numeric, function(x) x - x[ind])
  } else if(mode == 'divide') {
    # ignore non numeric columns
    # divide all other rows by the reference_group
    calib <- mutate_if(df, is.numeric, function(x) x / x[ind])
  }

  # return a tidy data.frame when tidy == TRUE
  if(tidy == TRUE) {
    calib <- calib %>%
      gather(gene, calibrated, -group)
  }
  return(calib)
}

#' Calculate standard deviation
#'
#' Uses a group_by statement to calculate the standard deviations of values in
#' a data.fram grouped by a variable
#'
#' @inheritParams .pcr_average
#'
#' @details Used to calculate the standard devaitions of ct after grouping by a
#' group_var for the experimental groups
#'
#' @return A data.frame with a column for the grouping variable; group_var and
#' one for each of the original columns. When tidy is TRUE the
#' returns a tidy data.frame with 3 columns; group/amount, gene and error.
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate standard deviations
#' .pcr_sd(ct1, group_var = group_var)
#'
#' # calculate standard deviations and return a tidy data.frame
#' .pcr_sd(ct1, group_var = group_var, tidy = TRUE)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate group_by summarise_all
#' @importFrom tidyr gather
#' @importFrom stats sd
#'
#' @keywords internal
#'
#' @export
.pcr_sd <- function(df, group_var, tidy = FALSE) {
  # group_by the group_var
  # calculate standard deviation
  sd <- mutate(df, group = group_var) %>%
   group_by(group) %>%
   summarise_all(function(x) sd(x))

  # return a tidy data.frame when tidy == TRUE
  if(tidy == TRUE) {
   sd <- sd %>%
     gather(gene, error, -group)
  }

  return(sd)
}

#' Calculate error terms
#'
#' Uses a specified column as a reference to calculate the error between it
#' and another column.
#'
#' @inheritParams .pcr_normalize
#' @inheritParams .pcr_average
#'
#' @details Used to sum the error of a gene and a reference_gene/column
#'
#' @return A data.frame with a column for each of the original columns after
#' taking the square root of the squared standard deviations of a target gene
#' and a reference_gene/column which is dropped. The function ignores
#' non numeric columns. When tidy is TRUE the returns a tidy data.frame with
#' the columns: gene and error as well as any non numeric columns such as a
#' grouping variable group/amount.
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate standard deviations
#' sds <- .pcr_sd(ct1, group_var = group_var)
#'
#' # calculate errors
#' .pcr_error(sds, reference_gene = 'GAPDH')
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select starts_with mutate_if
#' @importFrom tidyr gather
#'
#' @keywords internal
#'
#' @export
.pcr_error <- function(df, reference_gene, tidy = FALSE) {
  # get the reference_gene column and unlist
  ref <- select(df, reference_gene) %>%
    unlist(use.names = FALSE)

  # drop the reference gene column
  # ignore non numeric columns
  # calculate the error term by takine the squar root of the sum squares
  error <- select(df, -starts_with(reference_gene)) %>%
    mutate_if(is.numeric, function(x) sqrt((x^2) + (ref)^2))

  # return a tidy data.frame when tidy == TRUE
  if(tidy == TRUE) {
    error <- error %>%
      gather(gene, error, -group)
  }

  return(error)
}

#' Calculates the coefficient of variation
#'
#' Calculates the coefficient of variation of a gene by deviding standard
#' deviations of each group by their averages
#'
#' @param amounts A data.frame of the calculated input amounts returned by
#' \code{\link{.pcr_amount}}
#' @inheritParams .pcr_average
#'
#' @details Used to calculate the coefficient of variation of the input amounts
#'  after grouping by a grouping variable; group_var for experimental groups.
#'
#' @return A data.frame with a column for the grouping variable; group_var and
#' one for each of the original columns. When tidy is TRUE the
#' returns a tidy data.frame with 3 columns; group/amount, gene and error
#'
#' @examples
#' # locate and read data
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
#'
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
#'
#' # make a vector of RNA amounts
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # calculate curve
#' standard_curve <- pcr_assess(ct3, amount = amount, method = 'standard_curve')
#' intercept <- standard_curve$intercept
#' slope <- standard_curve$slope
#'
#' # calculate amounts
#' input_amounts <- .pcr_amount(ct1,
#'                             intercept = intercept,
#'                             slope = slope)
#'
#' # make grouping variable
#' group <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate cv errors
#' .pcr_cv(input_amounts,
#'        group_var = group)
#'
#' @importFrom tidyr gather spread
#' @importFrom dplyr group_by summarise ungroup
#'
#' @keywords internal
#'
#' @export
.pcr_cv <- function(amounts, group_var, tidy = FALSE) {
  # group_by group_var and calculate cv
  cv <- mutate(amounts, group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) sd(x)/mean(x))

  # return a tidy data.frame when tidy == TRUE
  if(tidy == TRUE) {
    cv <- gather(cv, gene, error, -group)
  }

  return(cv)
}

#' Calculate PCR RNA amounts
#'
#'
#' @inheritParams .pcr_average
#' @param intercept A numeric vector of length equals the number of columns of df,
#'  one for each gene such as the output of \link{pcr_assess}
#' @param slope A numeric vector of length equals the number of columns of df,
#' one for each gene such as the output of \link{pcr_assess}
#'
#' @details Used to alculate the amount of RNA in a PCR experimental sample using the
#' information provided by the standard curve, namely the slope and the
#' intercept calculated in advance for each gene in a similar expermiment.
#'
#' @return A data.frame with the same dimensions as df containing the
#' amount of RNA in each sample in each gene
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr data_frame full_join group_by mutate bind_cols
#' @importFrom tidyr gather
#'
#' @examples
#' # locate and read data
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
#'
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
#'
#' # make a vector of RNA amounts
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # calculate curve
#' standard_curve <- pcr_assess(ct3, amount = amount, method = 'standard_curve')
#' intercept <- standard_curve$intercept
#' slope <- standard_curve$slope
#'
#' # calculate amounts
#' .pcr_amount(ct1,
#'            intercept = intercept,
#'            slope = slope)
#'
#' @keywords internal
#'
#' @export
.pcr_amount <- function(df, intercept, slope) {
  # make a data.frame of intercept, slope ana gene names
  curve <- data_frame(intercept,
                      slope,
                      gene = names(df))
  # tidy the data.frame
  ct <- df %>%
    gather(gene, ct)

  # calculate input amounts
  amounts <- full_join(ct, curve) %>%
    group_by(gene) %>%
    mutate(amount = 10 ^ ((ct - intercept)/slope))

  # reshap input amounts
  with(amounts, split(amount, gene)) %>%
    bind_cols()
}

#' Calculate the linear trend
#'
#' Calculates the linear trend; intercept and slope between two variables
#'
#' @param df A data.frame of raw ct values or the delta ct values calculated
#' by \link{.pcr_normalize}
#' @param amount A numeric vector input amounts/dilutions of legnth equals the
#' number of the rows of df.
#'
#' @details Used to calculate the linear trend; intercept and slope for a line
#' between each column of ct or delta ct values and the log10 input amount
#'
#' @return A data.frame of 4 columns
#' \itemize{
#'   \item gene The column names of df
#'   \item intercept The intercept of the line
#'   \item slope The slope of teh line
#'   \item r_squared The squared correlation
#' }
#'
#'
#' @importFrom purrr map
#' @importFrom stats cor lm coefficients
#' @importFrom dplyr data_frame bind_rows
#'
#' @examples
#' # locate and read file
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
#'
#' # make amount/dilution variable
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # calculate trend
#' .pcr_trend(ct3, amount = amount)
#'
#' @keywords internal
#'
#' @export
.pcr_trend <- function(df, amount) {
  # make a trend line using linear regression
  trend_line <- map(df, function(x) {
    # calculate the r squared
    r_squared <- cor(x, log10(amount))^2

    # calculate the model
    ll <- lm(x ~ log10(amount))

    # get coeffecients
    coeff <- coefficients(ll)

    # make a data.frame of intercep, slope and rsquared
    data_frame(
      intercept = coeff[1],
      slope = coeff[2],
      r_squared = r_squared
    )
  })

  trend_line <- bind_rows(trend_line, .id = 'gene')

  return(trend_line)
}
