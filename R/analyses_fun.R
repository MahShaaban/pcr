#' Calculate the delta_delta_ct model
#'
#' Uses the \eqn{C_T} values and a reference gene and a group to calculate the
#' delta delta \eqn{C_T} model to estimate the normalized relative expression
#' of target genes.
#'
#' @param df A data.frame of \eqn{C_T} values with genes in the columns and
#' samples in rows rows
#' @param group_var A character vector of a grouping variable. The length of
#' this variable should equal the number of rows of df
#' @param reference_gene A character string of the column name of a control
#' gene
#' @param reference_group A character string of the control group in group_var
#' @param mode A character string of; 'separate_tube' (default) or 'same_tube'.
#' This is to indicate whether the different genes were run in separate or the
#' same PCR tube
#' @param plot A logical (default is FALSE)
#' @param ... Arguments passed to customize plot
#'
#' @return A data.frame of 8 columns:
#' \itemize{
#'   \item group The unique entries in group_var
#'   \item gene The column names of df. reference_gene is dropped
#'   \item normalized The \eqn{C_T} value (or the average \eqn{C_T} value) of
#'   target genes after subtracting that of the reference_gene
#'   \item calibrated The normalized average \eqn{C_T} value of target genes
#'   after subtracting that of the reference_group
#'   \item relative_expression The expression of target genes normalized by
#'   a reference_gene and calibrated by a reference_group
#'   \item error The standard deviation of the relative_expression
#'   \item lower The lower interval of the relative_expression
#'   \item upper The upper interval of the relative_expression
#' }
#' When \code{plot} is TRUE, returns a bar graph of the relative expression of
#' the genes in the column and the groups in the column group. Error bars are
#' drawn using the columns lower and upper. When more one gene are plotted the
#' default in dodge bars. When the argument facet is TRUE a separate panel is
#'  drawn for each gene.
#'
#' @details The comparative \eqn{C_T} methods assume that the cDNA templates of
#' the gene/s of interest as well as the control/reference gene have similar
#' amplification efficiency. And that this amplification efficiency is near
#' perfect. Meaning, at a certain threshold during the linear portion of the
#' PCR reaction, the amount of the gene of the interest and the control double
#' each cycle. Another assumptions is that, the expression difference between
#' two genes or two samples can be captured by subtracting one (gene or sample
#' of interest) from another (reference).  This final assumption requires also
#' that these references don't change with the treatment or the course in
#' question.
#'
#' @examples
#' ## locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- read.csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate all values and errors in one step
#' pcr_ddct(ct1,
#'          group_var = group_var,
#'          reference_gene = 'GAPDH',
#'          reference_group = 'brain')
#'
#' # return a plot
#' pcr_ddct(ct1,
#'          group_var = group_var,
#'          reference_gene = 'GAPDH',
#'          reference_group = 'brain',
#'          plot = TRUE)
#'
#' @export
pcr_ddct <- function(df, group_var, reference_gene, reference_group,
                     mode = 'separate_tube', plot = FALSE, ...) {
  # extract the reference gene and genes of interest
  ref <- subset(df, select = reference_gene, drop = TRUE)
  goi <- subset(df, select = names(df) != reference_gene)

  # apply the calculations
  res <- apply(goi,
               MARGIN = 2,
               FUN = function(x) {
                 if (mode == 'separate_tube') {
                   # calculate the averages
                   x_ave <- .pcr_average(x, group_var)
                   ref_ave <- .pcr_average(ref, group_var)

                   # calculate the dct
                   dct <- .pcr_normalize(x_ave, ref_ave)

                   # calculate the standard deviations
                   x_sd <- .pcr_sd(x, group_var)
                   ref_sd <- .pcr_sd(ref, group_var)

                   # calculate the error terms
                   error <- .pcr_error(x_sd, ref_sd)
                 } else if (mode == 'same_tube') {
                   # normalize first
                   norm <- .pcr_normalize(x, ref)

                   # calculate the averages
                   dct <- .pcr_average(norm, group_var)

                   # calculate the error
                   error <- .pcr_sd(norm, group_var)
                 } else {
                   stop("mode should be one of 'separate_tube' or 'same_tube'.")
                 }

                 # calculate the ddct
                 group_ref <- subset(dct, unique(group_var) == reference_group)
                 ddct <- .pcr_normalize(dct, group_ref)

                 # calculate the relative expression
                 rel_expr <- .pcr_relative(ddct)

                 # calculate the error bars
                 upper <- .pcr_relative(ddct - error)
                 lower <- .pcr_relative(ddct + error)

                 # make a data.frame
                 data.frame(
                   group = unique(group_var),
                   gene = '',
                   normalized = dct,
                   calibrated = ddct,
                   relative_expression = rel_expr,
                   error = error,
                   lower = lower,
                   upper = upper
                 )
               })

  # format results into a data.frame
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  res$gene <- rep(names(goi), each = length(unique(group_var)))

  # return
  # return plot when plot == TRUE
  if (plot) {
    gg <- .pcr_plot_analyze(res, method = 'delta_delta_ct', ...)
    return(gg)
  } else {
    return(res)
  }
}

#' Calculate the delta_ct model
#'
#' Uses the \eqn{C_T} values and a reference group to calculate the delta \eqn{C_T}
#' model to estimate the relative fold change of a gene between groups
#'
#' @inheritParams pcr_ddct
#'
#' @return A data.frame of 7 columns
#' \itemize{
#'   \item group The unique entries in group_var
#'   \item gene The column names of df
#'   \item calibrated The average \eqn{C_T} value of target genes after
#'   subtracting that of the reference_group
#'   \item fold_change The fold change of genes relative to a reference_group
#'   \item error The standard deviation of the fold_change
#'   \item lower The lower interval of the fold_change
#'   \item upper The upper interval of the fold_change
#' }
#' When \code{plot} is TRUE, returns a bar graph of the fold change of
#' the genes in the column and the groups in the column group. Error bars are
#' drawn using the columns lower and upper. When more one gene are plotted the
#' default in dodge bars. When the argument facet is TRUE a separate panel is
#'  drawn for each gene.
#'
#' @details This method is a variation of the double delta \eqn{C_T} model,
#' \code{\link{pcr_ddct}}. It can be used to calculate the fold change
#' of in one sample relative to the others. For example, it can be used to
#' compare and choosing a control/reference genes.
#'
#' @references Livak, Kenneth J, and Thomas D Schmittgen. 2001. “Analysis of
#' Relative Gene Expression Data Using Real-Time Quantitative PCR and the
#' Double Delta CT Method.” Methods 25 (4). ELSEVIER.
#' doi:10.1006/meth.2001.1262.
#'
#' @examples
#' # locate and read file
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- read.csv(fl)
#'
#' # make a data.frame of two identical columns
#' pcr_hk <- data.frame(
#'   GAPDH1 = ct1$GAPDH,
#'   GAPDH2 = ct1$GAPDH
#'   )
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate caliberation
#' pcr_dct(pcr_hk,
#'         group_var = group_var,
#'         reference_group = 'brain')
#'
#' # returns a plot
#' pcr_dct(pcr_hk,
#'         group_var = group_var,
#'         reference_group = 'brain',
#'         plot = TRUE)
#'
#' # returns a plot with facets
#' pcr_dct(pcr_hk,
#'         group_var = group_var,
#'         reference_group = 'brain',
#'         plot = TRUE,
#'         facet = TRUE)
#'
#' @export
pcr_dct <- function(df, group_var, reference_group,
                    mode = 'separate_tube', plot = FALSE, ...) {
  # apply the calculations
  res <- apply(df,
               MARGIN = 2,
               FUN = function(x) {
                 if (mode == 'separate_tube') {
                   # calculate the averages
                   x_ave <- .pcr_average(x, group_var)
                   group_ref <- subset(x_ave, unique(group_var) == reference_group)
                   dct <- .pcr_normalize(x_ave, group_ref)

                   # calculate the standard deviations
                   error <- .pcr_sd(x, group_var)

                 } else if (mode == 'same_tube') {
                   # normalize first
                   group_ref <- subset(x, group_var == reference_group)
                   norm <- .pcr_normalize(x, group_ref)

                   # calculate the averages
                   dct <- .pcr_average(norm, group_var)

                   # calculate the error
                   error <- .pcr_sd(norm, group_var)
                 } else {
                   stop("mode should be one of 'separate_tube' or 'same_tube'.")
                 }

                 # calculate the relative expression
                 rel_expr <- .pcr_relative(dct)

                 # calculate the error bars
                 upper <- .pcr_relative(dct - error)
                 lower <- .pcr_relative(dct + error)

                 # make a data.frame
                 data.frame(
                   group = unique(group_var),
                   gene = '',
                   calibrated = dct,
                   fold_change = rel_expr,
                   error = error,
                   lower = lower,
                   upper = upper
                 )
               })

  # format results into a data.frame
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  res$gene <- rep(names(df), each = length(unique(group_var)))

  # return
  # return plot when plot == TRUE
  if (plot) {
    gg <- .pcr_plot_analyze(res, method = 'delta_ct', ...)
    return(gg)
  } else {
    return(res)
  }
}

#' Calculate the standard curve model
#'
#' Uses the \eqn{C_T} values and a reference gene and a group, in addition to the
#' intercept and slope of each gene form a serial dilution experiment, to calculate
#' the standard curve model and estimate the normalized relative expression of the
#' target genes.
#'
#' @inheritParams pcr_ddct
#' @param intercept A numeric vector of intercept and length equals the number of genes
#' @param slope A numeric vector of slopes length equals the number of genes
#'
#' @return A data.frame of 7 columns
#' \itemize{
#'   \item group The unique entries in group_var
#'   \item gene The column names of df
#'   \item normalized The normalized expression of target genes relative to a reference_gene
#'   \item calibrated The calibrated expression of target genes relative to a reference_group
#'   \item error The standard deviation of normalized relative expression
#'   \item lower The lower interval of the normalized relative expression
#'   \item upper The upper interval of the normalized relative expression
#' }
#' When \code{plot} is TRUE, returns a bar graph of the calibrated expression
#' of the genes in the column and the groups in the column group. Error bars
#' are drawn using the columns lower and upper. When more one gene are plotted
#' the default in dodge bars. When the argument facet is TRUE a separate
#' panel is drawn for each gene.
#'
#' @details this model doesn't assume perfect amplification but rather actively
#' use the amplification in calculating the relative expression. So when the
#' amplification efficiency of all genes are 100\% both methods should give
#' similar results. The standard curve method is applied using two steps.
#' First, serial dilutions of the mRNAs from the samples of interest are used
#' as input to the PCR reaction. The linear trend of the log input amount and
#' the resulting \eqn{C_T} values for each gene are used to calculate an intercept
#' and a slope. Secondly, these intercepts and slopes are used to calculate the
#'  amounts of mRNA of the genes of interest and the control/reference in the
#'  samples of interest and the control sample/reference. These amounts are
#'  finally used to calculate the relative expression.
#'
#' @references Livak, Kenneth J, and Thomas D Schmittgen. 2001. “Analysis of
#' Relative Gene Expression Data Using Real-Time Quantitative PCR and the
#' Double Delta CT Method.” Methods 25 (4). ELSEVIER.
#' doi:10.1006/meth.2001.1262.
#'
#' @examples
#' # locate and read file
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- read.csv(fl)
#'
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- read.csv(fl)
#'
#' # make a vector of RNA amounts
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # calculate curve
#' standard_curve <- pcr_assess(ct3, amount = amount, method = 'standard_curve')
#' intercept <- standard_curve$intercept
#' slope <- standard_curve$slope
#'
#' # make grouping variable
#' group <- rep(c('brain', 'kidney'), each = 6)
#'
#' # apply the standard curve method
#' pcr_curve(ct1,
#'           group_var = group,
#'           reference_gene = 'GAPDH',
#'           reference_group = 'brain',
#'           intercept = intercept,
#'           slope = slope)
#'
#' # returns a plot
#' pcr_curve(ct1,
#'           group_var = group,
#'           reference_gene = 'GAPDH',
#'           reference_group = 'brain',
#'           intercept = intercept,
#'           slope = slope,
#'           plot = TRUE)
#'
#' @export
pcr_curve <- function(df, group_var, reference_gene, reference_group,
                      mode = 'separate_tube', intercept, slope,
                      plot = FALSE, ...) {
  amounts <- mapply(function(d, a, b) .pcr_amount(d, a, b),
                    d = df, a = intercept, b = slope)
  amounts <- as.data.frame(amounts)

  # extract the reference gene and genes of interest
  ref <- subset(amounts, select = reference_gene, drop = TRUE)
  goi <- subset(amounts, select = names(amounts) != reference_gene)

  # apply the calculations
  res <- apply(goi,
               MARGIN = 2,
               FUN = function(x) {
                 if (mode == 'separate_tube') {
                   # calculate the averages
                   x_ave <- .pcr_average(x, group_var)
                   ref_ave <- .pcr_average(ref, group_var)

                   # calculate the normalized values
                   norm <- .pcr_normalize(x_ave, ref_ave, mode = 'divide')

                   # calculate the standard deviations
                   x_cv <- .pcr_cv(x, group_var)
                   ref_cv <- .pcr_cv(ref, group_var)

                   # calculate the error terms
                   error <- .pcr_error(x_cv, ref_cv)
                 } else if (mode == 'same_tube') {
                   # normalize first
                   norm <- .pcr_normalize(x, ref, mode = 'divide')

                   # calculate the error
                   error <- .pcr_cv(norm, group_var)
                 } else {
                   stop("mode should be one of 'separate_tube' or 'same_tube'.")
                 }

                 # calculate the calibrated values
                 group_ref <- subset(norm, unique(group_var) == reference_group)
                 calib <- .pcr_normalize(norm, group_ref, mode = 'divide')

                 # calculate the error bars
                 upper <- calib + error
                 lower <- calib - error
                 error <- error * norm

                 # make a data.frame
                 data.frame(
                   group = unique(group_var),
                   gene = '',
                   normalized = norm,
                   calibrated = calib,
                   error = error,
                   lower = lower,
                   upper = upper
                 )
               })

  # format results into a data.frame
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  res$gene <- rep(names(goi), each = length(unique(group_var)))

  # return
  # return plot when plot == TRUE
  if (plot) {
    gg <- .pcr_plot_analyze(res, method = 'relative_curve', ...)
    return(gg)
  } else {
    return(res)
  }
}

#' Apply qPCR analysis methods
#'
#' A unified interface to invoke different analysis methods of qPCR data.
#'
#' @inheritParams pcr_ddct
#' @inheritParams pcr_curve
#' @param method A character string; 'delta_delta_ct' default, 'delta_ct' or
#' 'relative_curve' for invoking a certain analysis model
#' @param ... Arguments passed to the methods
#'
#' @return A data.frame by default, when \code{plot} is TRUE returns a plot.
#' For details; \link{pcr_ddct}, \link{pcr_dct} and \link{pcr_curve}.
#'
#' @details The different analysis methods can be invoked using the
#' argument method with 'delta_delta_ct' default, 'delta_ct' or
#' 'relative_curve' for the double delta \eqn{C_T}, delta ct or the standard curve
#' model respectively. Alternatively, the same methods can be applied by using
#' the corresponding functions directly: \link{pcr_ddct}, \link{pcr_dct} or
#' \link{pcr_curve}
#'
#' @references Livak, Kenneth J, and Thomas D Schmittgen. 2001. “Analysis of
#' Relative Gene Expression Data Using Real-Time Quantitative PCR and the
#' Double Delta CT Method.” Methods 25 (4). ELSEVIER.
#' doi:10.1006/meth.2001.1262.
#'
#' @examples
#' # applying the delta delta ct method
#' ## locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- read.csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate all values and errors in one step
#' pcr_analyze(ct1,
#'             group_var = group_var,
#'             reference_gene = 'GAPDH',
#'             reference_group = 'brain',
#'             method = 'delta_delta_ct')
#'
#' # return a plot
#' pcr_analyze(ct1,
#'             group_var = group_var,
#'             reference_gene = 'GAPDH',
#'             reference_group = 'brain',
#'             method = 'delta_delta_ct',
#'             plot = TRUE)
#'
#' # applying the delta ct method
#' # make a data.frame of two identical columns
#' pcr_hk <- data.frame(
#'   GAPDH1 = ct1$GAPDH,
#'   GAPDH2 = ct1$GAPDH
#'   )
#'
#' # calculate fold change
#' pcr_analyze(pcr_hk,
#'             group_var = group_var,
#'             reference_group = 'brain',
#'             method = 'delta_ct')
#'
#' # return a plot
#' pcr_analyze(pcr_hk,
#'             group_var = group_var,
#'             reference_group = 'brain',
#'             method = 'delta_ct',
#'             plot = TRUE)
#'
#' # applying the standard curve method
#' # locate and read file
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- read.csv(fl)
#'
#' # make a vector of RNA amounts
#' amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)
#'
#' # calculate curve
#' standard_curve <- pcr_assess(ct3, amount = amount, method = 'standard_curve')
#' intercept <- standard_curve$intercept
#' slope <- standard_curve$slope
#'
#' # apply the standard curve method
#' pcr_analyze(ct1,
#'            group_var = group_var,
#'            reference_gene = 'GAPDH',
#'            reference_group = 'brain',
#'            intercept = intercept,
#'            slope = slope,
#'            method = 'relative_curve')
#'
#' # return a plot
#' pcr_analyze(ct1,
#'            group_var = group_var,
#'            reference_gene = 'GAPDH',
#'            reference_group = 'brain',
#'            intercept = intercept,
#'            slope = slope,
#'            method = 'relative_curve',
#'            plot = TRUE)
#'
#' @export
pcr_analyze <- function(df, method = 'delta_delta_ct', ...) {
  switch(method,
         'delta_delta_ct' = pcr_ddct(df, ...),
         'delta_ct' = pcr_dct(df, ...),
         'relative_curve' = pcr_curve(df, ...))
}
