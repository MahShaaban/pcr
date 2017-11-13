#' Calculate amplification efficiency
#'
#' Uses the \eqn{C_T} values from a serial dilution experiment to calculate the
#' amplification efficiency of a PCR reaction.
#'
#' @param df A data.frame of \eqn{C_T} values with genes in the columns and samples
#' in rows rows. Each sample are replicates of a known input/dilution given by amount
#' @param amount A numeric vector of the input amounts or dilutions. The length
#' of this vector should equal the row number of df
#' @param reference_gene A character string of the column name of a control gene
#' @param plot A logical (default FALSE) to indicate whether to return a data.frame
#' or a plot
#'
#' @return When plot is FALSE returns a data.frame of 4 columns describing the line
#' between the \eqn{\Delta C_T} of target genes and the log of amount
#' \itemize{
#'   \item gene The column names of df. reference_gene is dropped
#'   \item intercept The intercept of the line
#'   \item slope The slope of the line
#'   \item r_squared The squared correlation
#' }
#'
#' When plot is TRUE returns a graph instead shows the average and
#' standard deviation of of the \eqn{\Delta C_T} at different input amounts.
#' In addition, a linear trend line is drawn.
#'
#' @details Fortunately, regardless of the method used in the analysis of qPCR
#' data, The quality assessment are done in a similar way. It requires an
#' experiment similar to that of calculating the standard curve. Serial
#' dilutions of the genes of interest and controls are used as input to the
#' reaction and different calculations are made. The amplification efficiency is
#' approximated be the linear trend between the difference between the \eqn{C_T}
#' value of a gene of interest and a control/reference (\eqn{\Delta C_T}) and
#' the log input amount. This piece of information is required when using the
#' \eqn{\Delta \Delta C_T} model. Typically, the slope of the curve should be very
#' small and the \eqn{R^2} value should be very close to one. Other analysis
#' methods are recommended when this is not the case.
#'
#' @references Livak, Kenneth J, and Thomas D Schmittgen. 2001. “Analysis of
#' Relative Gene Expression Data Using Real-Time Quantitative PCR and the
#' Double Delta CT Method.” Methods 25 (4). ELSEVIER.
#' doi:10.1006/meth.2001.1262.
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
#' @export
pcr_efficiency <- function(df, amount, reference_gene, plot = FALSE) {

  # return data when plot is false
  if(plot == TRUE) {
    gg <- pcr_plot_assess(df, amount, reference_gene, method = 'efficiency')

    return(gg)
  } else if(plot == FALSE) {
    # calculate delta_ct
    dct <- .pcr_normalize(df, reference_gene = reference_gene)

    # calculate trend; intercep, slop and r_squared
    trend <- .pcr_trend(dct, amount = amount)

    return(trend)
  }
}

#' Calculate the standard curve
#'
#' Uses the \eqn{C_T} values from a serial dilution experiment to calculate the
#' a curve for each gene and the log of the input amount
#'
#' @inheritParams pcr_efficiency
#'
#' @return When plot is FALSE returns a data.frame of 4 columns describing the line
#' between the \eqn{C_T} of each gene and the log of amount
#' \itemize{
#'   \item gene The column names of df
#'   \item intercept The intercept of the line
#'   \item slope The slope of the line
#'   \item r_squared The squared correlation
#' }
#'
#' When plot is TRUE returns a graph instead shows the average and
#' standard deviation of of the \eqn{C_T} at different input amounts.
#'
#' @details Fortunately, regardless of the method used in the analysis of qPCR
#' data, The quality assessment are done in a similar way. It requires an
#' experiment similar to that of calculating the standard curve. Serial
#' dilutions of the genes of interest and controls are used as input to
#' the reaction and different calculations are made.
#' Curves are required for each gene using the $C_T$ value and the log of the
#' input amount. In this case, a separate slope and intercept are required for
#'  the calculation of the relative expression when applying the standard curve
#'  model.
#'
#' @references Livak, Kenneth J, and Thomas D Schmittgen. 2001. “Analysis of
#' Relative Gene Expression Data Using Real-Time Quantitative PCR and the
#' Double Delta CT Method.” Methods 25 (4). ELSEVIER.
#' doi:10.1006/meth.2001.1262.
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
#' pcr_standard(ct3,
#'              amount = amount)
#'
#' # plot the standard curve
#' pcr_standard(ct3,
#'              amount = amount,
#'              plot = TRUE)
#'
#'
#' @export
pcr_standard <- function(df, amount, plot = FALSE) {
  # return data when plot is false
  # when plot == TRUE
  # plot a standard curve for each gene
  if(plot == TRUE) {
    pcr_plot_assess(df, amount, method = 'standard_curve')
    } else if(plot == FALSE) {
      # calculate trend; intercep, slop and r_squared
      trend <- .pcr_trend(df, amount)

      return(trend)
    }
}

#' Assess qPCR data quality
#'
#' A unified interface to invoke different quality assessment methods of qPCR data.
#'
#' @inheritParams pcr_efficiency
#' @param method A character string; 'standard_curve' (default) or 'efficiency'
#' for invoking a certain quality assessment model
#' @param ... Arguments passed to the methods
#'
#' @return A data.frame or a plot. For details; \link{pcr_standard} and
#' \link{pcr_efficiency}
#'
#' @details The different quality assessment methods can be invoked using the
#' argument method with 'standard_curve' or 'efficiency'. Alternatively, the
#' same methods can be applied by using the corresponding functions:
#' \link{pcr_standard} or \link{pcr_efficiency} for calculating the
#' amplification efficiency of a PCR reaction or the individual standard
#' curves respectively. Unlike the amplification efficiency calculation when,
#' using the double delta ct model, the standard curves are required in
#' calculating the standard curve analysis model.
#'
#' @examples
#' #' # locate and read file
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
#'
#' # make amount/dilution variable
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # calculate the standard curve
#' pcr_assess(ct3,
#'            amount = amount,
#'            method = 'standard_curve')
#'
#' # calculate amplification efficiency
#' pcr_assess(ct3,
#'            amount = amount,
#'            reference_gene = 'GAPDH',
#'            method = 'efficiency')
#'
#' @export
pcr_assess <- function(df, method = 'standard_curve', ...) {
  switch(method,
         'standard_curve' = pcr_standard(df, ...),
         'efficiency' = pcr_efficiency(df, ...))
}
