#' Plotting function
#'
#' @param df A data.frame such as this returned by \link{pcr_analyze}
#' @param method A character string. Possible input includes 'delta_delta_ct',
#' 'delta_ct' or 'relative_curve'
#' @param facets A logical of whether or not to use facets when applicable
#'
#' @return A ggplot object. A bar graph of caculated ge
#'
#' @keywords internal
#'
#' @examples
#' ## locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- read.csv(fl)
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
#' pcr:::.pcr_plot_analyze(df, method = 'delta_delta_ct')
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
#' pcr:::.pcr_plot_analyze(df, method = 'delta_ct')
#' pcr:::.pcr_plot_analyze(df, method = 'delta_ct', facet = TRUE)
#'
#' # calculate curve
#' # locate and read data
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- read.csv(fl)
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
#' pcr:::.pcr_plot_analyze(df, method = 'relative_curve')
#'
#' @importFrom ggplot2 ggplot geom_col geom_errorbar aes_string facet_wrap

.pcr_plot_analyze <- function(df, method, facets = FALSE) {
  # switch to a value/column to plot
  y <- switch (method,
               'delta_delta_ct' = 'relative_expression',
               'delta_ct' = 'fold_change',
               'relative_curve' = 'calibrated'
  )

  # make plot
  if (length(unique(df$gene)) == 1) {
    gg <- ggplot(df, aes_string(x = 'group', y = y)) +
      geom_col() +
      geom_errorbar(aes_string(ymin = 'lower', ymax = 'upper'))
  } else if (facets == TRUE) {
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
#' @keywords internal
#'
#' @examples
#' # locate and read file
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- read.csv(fl)
#'
#' # make amount/dilution variable
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # plot the standard curves
#' pcr:::.pcr_plot_assess(ct3,
#'                 amount = amount,
#'                 reference_gene = 'GAPDH',
#'                 method = 'standard_curve')
#'
#' # plot amplification efficiency
#' pcr:::.pcr_plot_assess(ct3,
#'                        amount = amount,
#'                        reference_gene = 'GAPDH',
#'                        method = 'efficiency')
#'
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap geom_smooth

.pcr_plot_assess <- function(df, amount, reference_gene, method) {
  switch (method,
    'efficiency' = {
      # extract the reference gene and genes of interest
      ref <- subset(df, select = reference_gene, drop = TRUE)
      goi <- subset(df, select = names(df) != reference_gene)

      # normalize the genes of interest by the reference gene
      dct <- apply(goi,
                   MARGIN = 2,
                   FUN = function(x) {
                     .pcr_normalize(x, ref)
                   })

      # make a data.frame from amounts, dct values and gene names
      df <- data.frame(x = rep(log10(amount), times = ncol(dct)),
                       y = unlist(as.numeric(dct), use.names = FALSE),
                       gene = rep(colnames(dct), each = nrow(dct)))

      # make plot
      gg <- ggplot(df, aes(x = df$x, y = df$y)) +
        geom_point() +
        facet_wrap(~df$gene) +
        geom_smooth(method = 'lm')

      return(gg)
    },
    'standard_curve' = {
      # make a data.frame of amounts, ct values and gene names
      df <- data.frame(x = rep(log10(amount), times = ncol(df)),
                       y = unlist(df, use.names = FALSE),
                       gene = rep(colnames(df), each = nrow(df)))

      # make plot
      gg <- ggplot(df, aes(x = df$x, y = df$y)) +
        geom_point() +
        facet_wrap(~df$gene)

      return(gg)
    }
  )
}
