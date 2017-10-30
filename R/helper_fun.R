#' Calculate average ct value
#'
#' Takes a data.frame of raw ct values and returns the averages in different
#' groups or dilutions
#'
#' @param df A data.frame of numeric ct value of n genes in columns and m
#' samples in raws
#' @param group_var A vector of length m as a grouping variables of samples
#' @param amount A numeric vector of the input dilutions or amounts
#' @param tidy A logical, default FALSE. When TRUE returns a tidy data.frame
#'
#' @return A data.frame with a column for each and an addition colum for the,
#' grouping variable that was provided, and a raw for each unique the  entry
#' in the group_var or amount. When tidy is TRUE the function returns a tidy
#' data.frame with three columns; group/amount, gene and average.
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
#' pcr_average(ct1, group_var = group_var)
#'
#' # calculate averages and return a tidy data.frame
#' pcr_average(ct1, group_var = group_var, tidy = TRUE)
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
#' pcr_average(ct3, amount = amount)
#'
#' # calculate averages and return a tidy data.frame
#' pcr_average(ct3, amount = amount, tidy = TRUE)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate group_by summarise_all
#' @importFrom tidyr gather
#'
#' @export
pcr_average <- function(df, group_var, amount, tidy = FALSE) {
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

#' Normalize ct values to a reference gene
#'
#' Takes a data.frame of ct values or average ct values for each gene in each
#' condition and returns the ct values normalized to a reference gene.
#'
#' @param df A data.frame of such as that returned by \link{pcr_average}
#' @param reference_gene A character string of the name of the column
#' corresponding to the reference gene
#' @param mode A character string of the normalization mode to be used. Default
#' is 'subtract'. Other possible modes include 'divide'
#' @inheritParams  pcr_average
#'
#' @return A data.frame of normalized ct values for each gene and a
#' grouping variable. When mode is 'subtract' (default) each the reference gene
#'  is subtracted from each column, and division is used instead when mode is
#' 'divide'. The function returns ignores non numeric columns and drops the one
#' that corresponds to the reference_gene
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
#'
#' # normalize the ct values
#' pcr_normalize(ct1, reference_gene = 'GAPDH')
#'
#' # normalize by division
#' pcr_normalize(ct1, reference_gene = 'GAPDH', mode = 'divide')
#'
#' # add grouping variable and average first
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#' ave <- pcr_average(ct1, group_var = group_var)
#'
#' # normalize by subtraction
#' pcr_normalize(ave, 'GAPDH')
#'
#' # normalize by division
#' pcr_normalize(ave, 'GAPDH', mode = 'divide')
#'
#' @importFrom dplyr select mutate_if starts_with
#'
#' @export
pcr_normalize <- function(df, reference_gene, mode = 'subtract', tidy = FALSE) {
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

#' calibrate ct values to a reference group
#'
#' Takes a data.frame of the ct or normalized ct values for each gene in each
#' condition and returns a data.frame of the calibrated ct value of each gene
#' e.x. (delta ct of gene of interest - delta ct value of same gene in
#' reference group)
#'
#' @param df A data.frame of ct or normalized ct values for each gene and a
#' grouping variable such as the output of the \link{pcr_normalize}
#' @param reference_group A character string of the reference group as it is
#' recorded in the grouping variable
#' @inheritParams pcr_normalize
#' @inheritParams pcr_average
#'
#' @return A data.frame of the calibrated ct values for each gene in a grouping
#' variable by a reference group
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
#' ave <- pcr_average(ct1, group_var = group_var)
#'
#' # calculate delta ct
#' dct <- pcr_normalize(ave, 'GAPDH')
#'
#' # calculate delta delta ct
#' pcr_calibrate(dct, 'brain', tidy = TRUE)
#'
#' # calculate delta delta ct and return a tidy data.frame
#' pcr_calibrate(dct, 'brain', tidy = TRUE)
#'
#' @importFrom dplyr filter select mutate_if
#' @importFrom tidyr gather
#'
#' @export
pcr_calibrate <- function(df, reference_group, mode = 'subtract', tidy = FALSE) {
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

#' Calculate standard error value
#'
#' @inheritParams pcr_average
#'
#' @return A data.frame of ncol n and nrow equals the number of unique
#' grouping variables containing the standard error values of each gene in each
#' group
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
#' pcr_sd(ct1, group_var = group_var)
#'
#' # calculate standard deviations and return a tidy data.frame
#' pcr_sd(ct1, group_var = group_var, tidy = TRUE)
#'
#' @importFrom dplyr mutate group_by summarise_all
#' @importFrom tidyr gather
#' @importFrom stats sd
#'
#' @export
pcr_sd <- function(df, group_var, tidy = FALSE) {
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

#' Calculate error value
#'
#' @param df A data.frame of ncol n and nrow equals the number of unique
#' grouping variables containing the standard error values of each gene in each
#' group
#' @inheritParams pcr_normalize
#' @inheritParams pcr_average
#'
#' @return A data.frame of error values for each gene and a grouping
#' variable. The column corresponding to the reference gene is dropped.
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
#' sds <- pcr_sd(ct1, group_var = group_var)
#'
#' # calculate errors
#' pcr_error(sds, reference_gene = 'GAPDH')
#'
#' @importFrom dplyr select starts_with mutate_if
#' @importFrom tidyr gather
#'
#' @export
pcr_error <- function(df, reference_gene, tidy = FALSE) {
  # get the reference_gene column and unlist
  ref <- select(df, reference_gene) %>% unlist(use.names = FALSE)

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

#' Calculate PCR RNA amounts
#'
#' Calculate the amount of RNA in a PCR experimental sample using the
#' information provided by the standard curve, namely the slope and the
#' intercept calculated in advance for each gene in a similar expermiment.
#'
#' @inheritParams pcr_average
#' @param intercept A numeric vector of length equals ncol(df), one for each
#' gene such as the output of \link{pcr_assess}
#' @param slope A numeric vector of length ncol(df), one for each gene such as
#' the output of \link{pcr_assess}
#'
#' @return A data.frame of dimensions equal that of the input df containing the
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
#' pcr_amount(ct1,
#'            intercept = intercept,
#'            slope = slope)
#'
#' @export
pcr_amount <- function(df, intercept, slope) {
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

#' Calculates cv
#'
#' Error terms for the standard curve relative quantification
#'
#' @param amounts A data.frame of calculated input RNA amounts
#' @param mode A character string; 'separate_tube' or 'same_tube'
#' @inheritParams pcr_average
#' @inheritParams pcr_normalize
#'
#' @return A data.frame of the cv term for each gene in each sample
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
#' input_amounts <- pcr_amount(ct1,
#'                             intercept = intercept,
#'                             slope = slope)
#'
#' # make grouping variable
#' group <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate cv errors
#' pcr_cv(input_amounts,
#'        group_var = group,
#'        reference_gene = 'GAPDH')
#'
#' @importFrom tidyr gather spread
#' @importFrom dplyr group_by summarise ungroup
#'
#' @export
pcr_cv <- function(amounts, group_var, reference_gene, mode = 'average_amount',
                   tidy = FALSE) {
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

#' Calculate the linear trend for the standard curve
#'
#' @inheritParams pcr_average
#'
#' @return A data.frame of 4 columns
#' \itemize{
#'   \item gene
#'   \item intercept
#'   \item slope
#'   \item r_squared
#' }
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
#' pcr_trend(ct3, amount = amount)
#'
#' @importFrom purrr map
#' @importFrom stats cor lm coefficients
#' @importFrom dplyr data_frame bind_rows
#'
#' @export
pcr_trend <- function(df, amount) {
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
