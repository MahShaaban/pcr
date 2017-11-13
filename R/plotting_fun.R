#' Plotting function
#'
#' @param df A data.frame
#' @param method A character string
#' @param facets A character string
#'
#' @return A plot
#'
#' @importFrom ggplot2 ggplot geom_col geom_errorbar aes_string facet_wrap
#'
#' @export
pcr_plot_analyze <- function(df, method, facets = FALSE) {
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
      geom_errorbar(aes_string(ymin = 'lower', ymax = 'upper'))
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
#' @importFrom magrittr %>%
#' @importFrom dplyr full_join mutate
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar facet_wrap geom_smooth
#'
#' @export
pcr_plot_assess <- function(df, amount, reference_gene, method) {
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
