#' Calculate amplification efficiency
#'
#' @inheritParams pcr_average
#' @inheritParams pcr_normalize
#' @param plot A logical (default FALSE)
#'
#' @return A data.frame of 4 columns when plot == FALSE. Returned by calling
#' \link{pcr_trend}
#'
#' @examples
#' # locate and read file
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
#'
#' # make amount/dilution variable
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # calculate amplification efficiency
#' pcr_efficiency(ct3,
#'                amount = amount,
#'                reference_gene = 'GAPDH')
#'
#' # plot amplification efficiency
#' pcr_efficiency(ct3,
#'                amount = amount,
#'                reference_gene = 'GAPDH',
#'                plot = TRUE)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr full_join mutate
#' @import ggplot2
#'
#' @export
pcr_efficiency <- function(df, amount, reference_gene, plot = FALSE) {
  # calculate delta_ct
  dct <- pcr_normalize(df, reference_gene = reference_gene)
  # calculate trend; intercep, slop and r_squared
  trend <- pcr_trend(dct, amount = amount)

  # return data when plot is false
  if(plot == TRUE) {
    # calculate average of normalizex values
    ave <- pcr_normalize(ct3, reference_gene = reference_gene) %>%
      pcr_average(group_var = amount, tidy = TRUE)

    # calculate the standard deviation
    sd <- pcr_normalize(ct3, reference_gene = reference_gene) %>%
      pcr_sd(group_var = amount, tidy = TRUE)

    # merge data.frames and calculate intervals
    dat <- full_join(ave, sd) %>%
      full_join(trend) %>%
      mutate(group = log10(group),
             lower = average - error,
             upper = average + error)

    # make efficiency plot
    gg <- ggplot(dat, aes(x = group, y = average)) +
      geom_point() +
      geom_errorbar(aes(ymin = lower, ymax = upper)) +
      facet_wrap(~gene) +
      geom_smooth(method = 'lm')

    return(gg)
  } else if(plot == FALSE) {
    return(trend)
  }
}

#' Calculate the standard curve
#'
#' @inheritParams pcr_average
#' @inheritParams pcr_normalize
#' @param plot A logical (default FALSE)
#'
#' @return A data.frame of 4 columns when plot == FALSE. Returned by calling
#' \link{pcr_trend}
#'
#' @examples
#' # locate and read file
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
#'
#' # make amount/dilution variable
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # calculate the standard curve
#' pcr_efficiency(ct3,
#'                amount = amount,
#'                reference_gene = 'GAPDH')
#'
#' # plot the standard curve
#' pcr_standard(ct3,
#'              amount = amount,
#'              plot = TRUE)
#'
#' @importFrom dplyr mutate full_join
#' @importFrom tidyr gather
#' @import ggplot2
#'
#' @export
pcr_standard <- function(df, amount, plot = FALSE) {
  # return data when plot is false
  # calculate trend; intercep, slop and r_squared
  trend <- pcr_trend(df, amount)

  # when plot == TRUE
  # plot a standard curve for each gene
  if(plot == TRUE) {
    # make a data.frame of data
    dat <- mutate(df, log_amount = log10(amount)) %>%
      gather(gene, ct, -log_amount) %>%
      full_join(trend)
    gg <- ggplot(dat, aes(x = log_amount, y = ct)) +
      geom_point() +
      facet_wrap(~gene)
    return(gg)
    } else if(plot == FALSE) {
      return(trend)
    }
}

#' Assess qPCR data quality
#'
#' @inheritParams pcr_average
#' @param method A character string; 'standard_curve' defalult or 'efficiency'
#' @param ... Arguments passed to other methods
#'
#' @return A data.frame or a plot. For details; \link{pcr_standard} and
#' \link{pcr_efficiency}
#'
#' @export
pcr_assess <- function(df, method = 'standard_curve', ...) {
  switch(method,
         'standard_curve' = pcr_standard(df, ...),
         'efficiency' = pcr_efficiency(df, ...))
}
