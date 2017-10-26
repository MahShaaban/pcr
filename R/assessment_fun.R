#' Assess PCR data quality
#'
#' Assess the quality of a PCR experiment by examining the amplification
#' effeciency when using the comparative ct method or calculating the standard
#' curve
#'
#' @param df A data.frame of one or more columns containing ct values of a
#' reference gene and another from an experiment run with differen dilutions
#' (RNA amounts)
#' @param amount A numeric vector of length equals nrow df with RNA amounts
#' @param mode A character string of assessment mode. Default "effeciency"
#' @param plot A logical default FALSE of whether plot or return the data
#'
#' @return When mode == 'effeciency' return a data.frame of 5 columns or a plot
#' . when mode == 'standard_curve' retruns a data.frame of intercept and slope
#' or a plot for each of the input columns.
#'
#' @examples
#' # locate and read data
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
#'
#' # make a vector of RNA amounts
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # calculate effeciencey
#' pcr_assess(ct3,
#'            amount = amount)
#'
#' # calculate standard curve
#' pcr_assess(ct3,
#'            amount = amount,
#'            mode = 'standard_curve')
#'
#' @importFrom dplyr mutate summarise_all select data_frame bind_rows
#' @importFrom purrr map
#' @importFrom stats lm coefficients
#'
#' @export
pcr_assess <- function(df, amount, mode = 'effeciency', plot = FALSE) {
  if(mode == 'effeciency') {
    log_amount <- unique(log10(amount))

    ave <- df %>%
      mutate(amount = amount) %>%
      group_by(amount) %>%
      summarise_all(function(x) mean(x)) %>%
      select(-amount)

    dct <- unlist(ave[, 1] - ave[, 2], use.names = FALSE)

    error <- df %>%
      mutate(amount = amount) %>%
      group_by(amount) %>%
      summarise_all(function(x) sd(x)) %>%
      select(-amount)

    error <- as.numeric(sqrt((error[, 1] ^ 2) + error[, 2] ^ 2))

    int_lower = unlist(dct - error)
    int_upper = unlist(dct + error)

    effeciency <- data_frame(
      log_amount = log_amount,
      dct = dct,
      error = error,
      int_lower = int_lower,
      int_upper = int_upper
    )

    return(effeciency)
  } else if(mode == 'standard_curve') {
    standard_cruve <- map(df, function(x) {
      ll <- lm(x ~ log10(amount))
      coeff <- coefficients(ll)
      data_frame(intercept = coeff[1],
                 slope = coeff[2])
    }) %>%
      bind_rows(.id = 'gene')

    return(standard_cruve)
  }
}

