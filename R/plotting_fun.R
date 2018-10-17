#' Plotting function
#'
#' @param df A data.frame such as this returned by \link{pcr_analyze}
#' @param method A character string. Possible input includes 'delta_delta_ct',
#' 'delta_ct' or 'relative_curve'
#' @param facets A logical of whether or not to use facets when applicable
#'
#' @return A ggplot object. A bar graph of caculated ge
#'
#' @importFrom ggplot2 ggplot geom_col geom_errorbar aes_string facet_wrap
#'
#' @examples
#' ## locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate all delta_delta_ct model
#' df <- pcr_ddct(ct1,
#'                group_var = group_var,
#'                reference_gene = 'GAPDH',
#'                reference_group = 'brain')
#'
#' # make a plot
#' .pcr_plot_analyze(df, method = 'delta_delta_ct')
#'
#' # make a data.frame of two identical columns
#' pcr_hk <- data.frame(
#'   GAPDH1 = ct1$GAPDH,
#'   GAPDH2 = ct1$GAPDH
#'   )
#'
#' # calculate delta_ct model
#' df <- pcr_dct(pcr_hk,
#'               group_var = group_var,
#'               reference_group = 'brain')
#'
#' # make a plot
#' .pcr_plot_analyze(df, method = 'delta_ct')
#' .pcr_plot_analyze(df, method = 'delta_ct', facet = TRUE)
#'
#' # calculate curve
#' # locate and read data
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
#'
#' # make a vector of RNA amounts
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' standard_curve <- pcr_assess(ct3,
#'                              amount = amount,
#'                              method = 'standard_curve')
#' intercept <- standard_curve$intercept
#' slope <- standard_curve$slope
#'
#' # calculate the rellative_curve model
#' df <- pcr_curve(ct1,
#'                 group_var = group_var,
#'                 reference_gene = 'GAPDH',
#'                 reference_group = 'brain',
#'                 intercept = intercept,
#'                 slope = slope)
#'
#' # make a plot
#' .pcr_plot_analyze(df, method = 'relative_curve')
#'
.pcr_plot_analyze <- function(df, method, facets = FALSE) {
  y <- switch (method,
               'delta_delta_ct' = 'relative_expression',
               'delta_ct' = 'fold_change',
               'relative_curve' = 'calibrated'
  )

  if(length(unique(df$gene)) == 1) {
    gg <- ggplot(df, aes_string(x = 'group', y = y)) +
      geom_col() +
      geom_errorbar(aes_string(ymin = 'lower', ymax = 'upper'))
  } else if(facets == TRUE) {
    gg <- ggplot(df, aes_string(x = 'group', y = y)) +
      geom_col() +
      geom_errorbar(aes_string(ymin = 'lower', ymax = 'upper')) +
      facet_wrap(~ gene)
  } else {
    gg <- ggplot(df, aes_string(x = 'group', y = y, fill = 'gene')) +
      geom_col(position = 'dodge') +
      geom_errorbar(aes_string(ymin = 'lower', ymax = 'upper'),
                    position = 'dodge')
  }

  return(gg)
}

#' Plot quality assessment graphs
#'
#' @param df A data.frame
#' @param amount A numeric vector
#' @param reference_gene A character string
#' @param method A character string; 'efficiency' or 'standard_curve'
#'
#' @return A plot
#'
#' @examples
#' # locate and read file
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
#'
#' # make amount/dilution variable
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # plot the standard curves
#' .pcr_plot_assess(ct3,
#'                 amount = amount,
#'                 reference_gene = 'GAPDH',
#'                 method = 'standard_curve')
#'
#' # plot amplification efficiency
#' .pcr_plot_assess(ct3,
#'                 amount = amount,
#'                 reference_gene = 'GAPDH',
#'                 method = 'efficiency')
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr full_join mutate
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar facet_wrap geom_smooth
#'
.pcr_plot_assess <- function(df, amount, reference_gene, method) {
  switch (method,
    'efficiency' = {
      # calculate delta_ct
      dct <- .pcr_normalize(df, reference_gene = reference_gene)

      # calculate trend; intercep, slop and r_squared
      trend <- .pcr_trend(dct, amount = amount)

      # return data when plot is false
      # calculate trend; intercep, slop and r_squared
      trend <- .pcr_trend(dct, amount = amount)

      # calculate average of normalizex values
      ave <- .pcr_normalize(df, reference_gene = reference_gene) %>%
        .pcr_average(group_var = amount, tidy = TRUE)

      # calculate the standard deviation
      sd <- .pcr_normalize(df, reference_gene = reference_gene) %>%
        .pcr_sd(group_var = amount, tidy = TRUE)

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
    },
    'standard_curve' = {
      # calculate trend; intercep, slop and r_squared
      trend <- .pcr_trend(df, amount)

      # make a data.frame of data
      dat <- mutate(df, log_amount = log10(amount)) %>%
        gather(gene, ct, -log_amount) %>%
        full_join(trend)
      gg <- ggplot(dat, aes(x = log_amount, y = ct)) +
        geom_point() +
        facet_wrap(~gene)

      return(gg)
    }
  )
}
