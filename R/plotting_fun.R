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
pcr_plot <- function(df, method, facets = FALSE) {
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
