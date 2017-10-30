#' Calculate the delta_delta_ct
#'
#' @param df A data.frame of ct values
#' @param group_var A character vector of a grouping variable
#' @param reference_gene A character string of the column name of a control gene
#' @param reference_group A character string of the control group
#' @param mode A character string of; 'separate_tube' or 'same_tube'
#'
#' @return A data.frame of 8 columns
#' \itemize{
#'   \item group
#'   \item gene
#'   \item normalized
#'   \item calibrated
#'   \item relative_expression
#'   \item error
#'   \item lower
#'   \item upper
#' }
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
#' @importFrom tidyr gather
#' @importFrom dplyr mutate full_join
#'
#' @export
pcr_ddct <- function(df, group_var, reference_gene, reference_group, mode = 'separate_tube') {
  # calculate the delta_ct
  if(mode == 'separate_tube') {
    # calculate average ct and normalize
    ave <- pcr_average(df, group_var = group_var)
    dct <- pcr_normalize(ave, reference_gene = reference_gene)
    } else if(mode == 'same_tube') {
      # normalize and average normalized ct values
      dct <- pcr_normalize(df, reference_gene = reference_gene)
      dct <- pcr_average(dct, group_var = group_var)
    }
  # retain the normalized ct
  delta_ct <- gather(dct, gene, normalized, -group)

  # calculate the delta_delta_ct
  ddct <- pcr_calibrate(dct, reference_group = reference_group, tidy = TRUE)

  # calculate the relative expression
  norm_rel <- mutate(ddct, relative_expression = 2 ^ -calibrated)

  if(mode == 'separate_tube') {
    # calculate the error from ct values
    sds <- pcr_sd(df, group_var = group_var)
    error <- pcr_error(sds, reference_gene = reference_gene, tidy = TRUE)
    } else if(mode == 'same_tube') {
      # calculate error from normalized ct values
      dct <- pcr_normalize(df, reference_gene = reference_gene)
      error <- pcr_sd(dct, group_var = group_var, tidy = TRUE)
    }

  # merge data.frames and calculate intervals
  res <- full_join(delta_ct, ddct) %>%
    full_join(norm_rel) %>%
    full_join(error) %>%
    mutate(lower = 2 ^ -(calibrated + error),
           upper = 2 ^ -(calibrated - error))

  return(res)
}

#' Calculate the delta_ct
#'
#' @inheritParams pcr_ddct
#'
#' @return A data.frame of 7 columns
#' \itemize{
#'   \item group
#'   \item gene
#'   \item calibrated
#'   \item fold_change
#'   \item error
#'   \item lower
#'   \item upper
#' }
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
#' @importFrom tidyr gather
#' @importFrom dplyr mutate full_join
#'
#' @export
pcr_dct <- function(df, group_var, reference_gene, reference_group, mode = 'separate_tube') {
  if(mode == 'separate_tube') {
    # average ct and calibrate to a reference group
    ave <- pcr_average(df, group_var = group_var)
    dct <- pcr_calibrate(ave, reference_group = reference_group)
  } else if(mode == 'same_tube') {
    # calibrate ct and average
    dct <- pcr_calibrate(df, reference_group = reference_group)
    dct <- pcr_average(dct, group_var = group_var)
  }

  # retain calibrated values
  # calculate the fold change
  calib <- gather(dct, gene, calibrated, -group) %>%
    mutate(fold_change = 2 ^ -calibrated)

  if(mode == 'separate_tube') {
    # calculate the standard deviation from ct values
    sds <- pcr_sd(df, group_var = group_var, tidy = TRUE)
  } else if(mode == 'same_tube') {
    # calibrate ct values to a reference group
    # calculated sd from calibrated values
    dct <- pcr_calibrate(df, reference_group = reference_group)
    sds <- pcr_sd(dct, group_var = group_var, tidy = TRUE)
  }

  # join data frame and calculate intervals
  res <- full_join(calib, sds) %>%
    mutate(lower = 2 ^ -(calibrated + error),
           upper = 2 ^ -(calibrated - error))

  return(res)
}

#' Calculate the standard_curve
#'
#' @inheritParams pcr_ddct
#' @param mode A character string; 'separate_tube' or 'same_tube'
#' @param intercept A numeric vector of intercept and length equals the number of genes
#' @param slope A numeric vector of slopes length equals the number of genes
#'
#' @return A data.frame of 7 columns
#' \itemize{
#'   \item group
#'   \item gene
#'   \item normalized
#'   \item calibrated
#'   \item error
#'   \item lower
#'   \item upper
#' }
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
#' pcr_ddct(ct1,
#'          group_var = group,
#'          reference_gene = 'GAPDH',
#'          reference_group = 'brain')
#'
#' @importFrom tidyr gather
#' @importFrom dplyr full_join mutate
#' @export
pcr_curve <- function(df, group_var, reference_gene, reference_group, mode = 'separate_tube',
                      intercept, slope) {
  # calculate the amount of rna in samples
  amounts <- pcr_amount(df,
                        intercept = intercept,
                        slope = slope)
  if(mode == 'separate_tube') {
    # average amounts and normalize by a reference_gene
    ave <- pcr_average(amounts, group_var = group_var)
    norm <- pcr_normalize(ave, reference_gene = reference_gene, mode = 'divide')
  } else if(mode == 'same_tube') {
    # normalize amounts and average
    norm <- pcr_normalize(amounts, reference_gene = reference_gene, mode = 'divide')
    norm <- pcr_average(norm, group_var = group_var)
  }

  # retain normalized amounts
  normalized <- gather(norm, gene, normalized, -group)

  # calibrate to a reference_group
  calib <- pcr_calibrate(norm, reference_group = reference_group,
                          mode = 'divide', tidy = TRUE)

  if(mode == 'separate_tube') {
    # calculate cv from amounts
    cv <- pcr_cv(amounts, group_var = group_var, mode = mode)
    error <- pcr_error(cv, reference_gene = reference_gene, tidy = TRUE)
  } else if(mode == 'same_tube') {
    # calculate cv from normalized amounts
    norm <- pcr_normalize(amounts, reference_gene = reference_gene, mode = 'divide')
    error <- pcr_cv(norm, group_var = group_var, tidy = TRUE)
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
#' @inheritParams pcr_average
#' @param method A character string; 'delta_delta_ct' default, 'delta_ct' or
#' 'relative_curve'
#' @param ... Arguments passed to other methods
#'
#' @return A data.frame or a plot. For details; \link{pcr_ddct}, \link{pcr_dct} and
#' \link{pcr_curve}
#'
#' @export
pcr_analyze <- function(df, method = 'delta_delta_ct', ...) {
  switch(method,
         'delta_delta_ct' = pcr_ddct(df, ...),
         'delta_ct' = pcr_dct(df, ...),
         'relative_curve' = pcr_curve(df, ...))
}

