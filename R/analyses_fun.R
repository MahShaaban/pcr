#' Calculate the delta_delta_ct model
#'
#' Uses the \eqn{C_T} values and a reference gene and a group to calculate the delta
#' delta \eqn{C_T} model to estimate the normalized relative expresssion of target
#' genes.
#'
#' @param df A data.frame of \eqn{C_T} values with genes in the columns and samples
#' in rows rows
#' @param group_var A character vector of a grouping variable. The length of
#' this variable should equal the number of rows of df
#' @param reference_gene A character string of the column name of a control gene
#' @param reference_group A character string of the control group in group_var
#' @param mode A character string of; 'separate_tube' (default) or 'same_tube'.
#' This is to indicate whether the different genes were run in separate or the
#' same PCR tube
#'
#' @return A data.frame of 8 columns:
#' \itemize{
#'   \item group The unique entries in group_var
#'   \item gene The column names of df. reference_gene is dropped
#'   \item normalized The \eqn{C_T} value (or the average \eqn{C_T} value) of target genes
#'   after subtracting that of the reference_gene
#'   \item calibrated The normalized average \eqn{C_T} value of target genes after
#'   subtracting that of the reference_group
#'   \item relative_expression The expression of target genes normalized by
#'   a reference_gene and caliberated by a reference_group
#'   \item error The standard deviation of the relative_xpression
#'   \item lower The lower interval of the relative_expression
#'   \item upper The upper interval of the relative_expression
#' }
#'
#' @details The comparitive \eqn{C_T} methods assume that the cDNA templates of the
#'  gene/s of interest as well as the control/reference gene have similar
#'  amplification efficiency. And that this amplification efficiency is near
#'  perfect. Meaning, at a certain threshold during the linear portion of the
#'  PCR reaction, the amount of the gene of the interest and the control double
#'   each cycle. Another assumptios is that, the expression difference between
#'   two genes or two samples can be captured by subtracting one (gene or
#'   sample of interest) from another (reference).  This final assumption
#'   requires also that these references doesnt change with the treatment or
#'   the course in question.
#'
#' @examples
#' ## locate and read raw ct data
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
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
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr mutate full_join
#'
#' @export
pcr_ddct <- function(df, group_var, reference_gene, reference_group, mode = 'separate_tube') {
  # calculate the delta_ct
  if(mode == 'separate_tube') {
    # calculate average ct and normalize
    ave <- .pcr_average(df, group_var = group_var)
    dct <- .pcr_normalize(ave, reference_gene = reference_gene)
    } else if(mode == 'same_tube') {
      # normalize and average normalized ct values
      dct <- .pcr_normalize(df, reference_gene = reference_gene)
      dct <- .pcr_average(dct, group_var = group_var)
    }
  # retain the normalized ct
  delta_ct <- gather(dct, gene, normalized, -group)

  # calculate the delta_delta_ct
  ddct <- .pcr_calibrate(dct, reference_group = reference_group, tidy = TRUE)

  # calculate the relative expression
  norm_rel <- mutate(ddct, relative_expression = 2 ^ -calibrated)

  if(mode == 'separate_tube') {
    # calculate the error from ct values
    sds <- .pcr_sd(df, group_var = group_var)
    error <- .pcr_error(sds, reference_gene = reference_gene, tidy = TRUE)
    } else if(mode == 'same_tube') {
      # calculate error from normalized ct values
      dct <- .pcr_normalize(df, reference_gene = reference_gene)
      error <- .pcr_sd(dct, group_var = group_var, tidy = TRUE)
    }

  # merge data.frames and calculate intervals
  res <- full_join(delta_ct, ddct) %>%
    full_join(norm_rel) %>%
    full_join(error) %>%
    mutate(lower = 2 ^ -(calibrated + error),
           upper = 2 ^ -(calibrated - error))

  return(res)
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
#'
#' @details This method is a variation of the double delta \eqn{C_T} model,
#' \code{\link{pcr_ddct}}. It can be used to caculate the fold change
#' of in one sample relative to the others. For example, it can be used to
#' compare and choosing a control/referece genes.
#'
#' @references Livak, Kenneth J, and Thomas D Schmittgen. 2001. “Analysis of
#' Relative Gene Expression Data Using Real-Time Quantitative PCR and the
#' Double Delta CT Method.” Methods 25 (4). ELSEVIER.
#' doi:10.1006/meth.2001.1262.
#'
#' @examples
#' # locate and read file
#' fl <- system.file('extdata', 'ct1.csv', package = 'pcr')
#' ct1 <- readr::read_csv(fl)
#'
#' # make a data.frame of two identical columns
#' pcr_hk <- data.frame(
#'   GAPDH1 = ct1$GAPDH,
#'   GAPDH2 = ct1$GAPDH
#'   )
#'
## add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate caliberation
#' pcr_dct(pcr_hk,
#'         group_var = group_var,
#'         reference_group = 'brain')
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr mutate full_join
#'
#' @export
pcr_dct <- function(df, group_var, reference_gene, reference_group, mode = 'separate_tube') {
  if(mode == 'separate_tube') {
    # average ct and calibrate to a reference group
    ave <- .pcr_average(df, group_var = group_var)
    dct <- .pcr_calibrate(ave, reference_group = reference_group)
  } else if(mode == 'same_tube') {
    # calibrate ct and average
    dct <- .pcr_calibrate(df, reference_group = reference_group)
    dct <- .pcr_average(dct, group_var = group_var)
  }

  # retain calibrated values
  # calculate the fold change
  calib <- gather(dct, gene, calibrated, -group) %>%
    mutate(fold_change = 2 ^ -calibrated)

  if(mode == 'separate_tube') {
    # calculate the standard deviation from ct values
    sds <- .pcr_sd(df, group_var = group_var, tidy = TRUE)
  } else if(mode == 'same_tube') {
    # calibrate ct values to a reference group
    # calculated sd from calibrated values
    dct <- .pcr_calibrate(df, reference_group = reference_group)
    sds <- .pcr_sd(dct, group_var = group_var, tidy = TRUE)
  }

  # join data frame and calculate intervals
  res <- full_join(calib, sds) %>%
    mutate(lower = 2 ^ -(calibrated + error),
           upper = 2 ^ -(calibrated - error))

  return(res)
}

#' Calculate the standard curve model
#'
#' Uses the \eqn{C_T} values and a reference gene and a group, in addition to the
#' intercept and slope of each gene form a serial dilution experiment, to calculate
#' the standard curve model and estimate the normalized relative expresssion of the
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
#'   \item normalized The normalized expresison of target genes relative to a reference_gene
#'   \item calibrated The calibrated expression of target genes relative to a reference_group
#'   \item error The standard deviation of normalized relative expression
#'   \item lower The lower interval of the normalized relative expression
#'   \item upper The upper interval of the normalized relative expression
#' }
#'
#' @details this model doesn't assume perfect amplificaion but rather actively
#' use the amplificaion in calculating the relative expresion. So when the
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
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr full_join mutate
#' @export
pcr_curve <- function(df, group_var, reference_gene, reference_group, mode = 'separate_tube',
                      intercept, slope) {
  # calculate the amount of rna in samples
  amounts <- .pcr_amount(df,
                        intercept = intercept,
                        slope = slope)
  if(mode == 'separate_tube') {
    # average amounts and normalize by a reference_gene
    ave <- .pcr_average(amounts, group_var = group_var)
    norm <- .pcr_normalize(ave, reference_gene = reference_gene, mode = 'divide')
  } else if(mode == 'same_tube') {
    # normalize amounts and average
    norm <- .pcr_normalize(amounts, reference_gene = reference_gene, mode = 'divide')
    norm <- .pcr_average(norm, group_var = group_var)
  }

  # retain normalized amounts
  normalized <- gather(norm, gene, normalized, -group)

  # calibrate to a reference_group
  calib <- .pcr_calibrate(norm, reference_group = reference_group,
                          mode = 'divide', tidy = TRUE)

  if(mode == 'separate_tube') {
    # calculate cv from amounts
    cv <- .pcr_cv(amounts, group_var = group_var)
    error <- .pcr_error(cv, reference_gene = reference_gene, tidy = TRUE)
  } else if(mode == 'same_tube') {
    # calculate cv from normalized amounts
    norm <- .pcr_normalize(amounts, reference_gene = reference_gene, mode = 'divide')
    error <- .pcr_cv(norm, group_var = group_var, tidy = TRUE)
  }
  # join data.frames and calculate intervals
  res <- full_join(normalized, calib) %>%
    full_join(error) %>%
    mutate(lower = calibrated - error,
           upper = calibrated + error,
           error = error * normalized)

  return(res)
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
#' @return A data.frame. For details; \link{pcr_ddct}, \link{pcr_dct} and
#' \link{pcr_curve}
#'
#' @details The different analysis methods can be invoked using the
#' argument method with 'delta_delta_ct' default, 'delta_ct' or
#' 'relative_curve' for the double delta \eqn{C_T}, delta ct or the standard curve
#' model respectively. Alternativel, the same methods can be applied by using
#' the corresponding functions direcly: \link{pcr_ddct}, \link{pcr_dct} or
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
#' ct1 <- readr::read_csv(fl)
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
#' # applying the delta ct method
#' # make a data.frame of two identical columns
#' pcr_hk <- data.frame(
#'   GAPDH1 = ct1$GAPDH,
#'   GAPDH2 = ct1$GAPDH
#'   )
#'
#' # calculate caliberation
#' pcr_analyze(pcr_hk,
#'             group_var = group_var,
#'             reference_group = 'brain',
#'             method = 'delta_ct')
#'
#' # applying the standard curve method
#' # locate and read file
#' fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
#' ct3 <- readr::read_csv(fl)
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
#' @export
pcr_analyze <- function(df, method = 'delta_delta_ct', ...) {
  switch(method,
         'delta_delta_ct' = pcr_ddct(df, ...),
         'delta_ct' = pcr_dct(df, ...),
         'relative_curve' = pcr_curve(df, ...))
}

