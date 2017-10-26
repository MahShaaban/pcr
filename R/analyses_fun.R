#' Calculate the delta_delta_ct
#'
#' @param df A data.frame of ct values
#' @param group_var A character vector of a grouping variable
#' @param reference_gene A character string of the column name of a control gene
#' @param reference_group A character string of the control group
#' @param mode A character string of; 'average_ct' or 'average_dct'
#'
#' @return A data.frame of 8 columns
#' \itemize{
#'   \item group
#'   \item gene
#'   \item normalized
#'   \item caliberated
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
pcr_ddct <- function(df, group_var, reference_gene, reference_group, mode = 'average_ct') {
  # calculate the delta_ct
  if(mode == 'average_ct') {
    # calculate average ct and normalize
    ave <- pcr_average(df, group_var = group_var)
    dct <- pcr_normalize(ave, reference_gene = reference_gene)
    } else if(mode == 'average_dct') {
      # normalize and average normalized ct values
      dct <- pcr_normalize(df, reference_gene = reference_gene)
      dct <- pcr_average(dct, group_var = group_var)
    }
  # retain the normalized ct
  delta_ct <- gather(dct, gene, normalized, -group)

  # calculate the delta_delta_ct
  ddct <- pcr_caliberate(dct, reference_group = reference_group, tidy = TRUE)

  # calculate the relative expression
  norm_rel <- mutate(ddct, relative_expression = 2 ^ -caliberated)

  if(mode == 'average_ct') {
    # calculate the error from ct values
    sds <- pcr_sd(df, group_var = group_var)
    error <- pcr_error(sds, reference_gene = reference_gene, tidy = TRUE)
    } else if(mode == 'average_dct') {
      # calculate error from normalized ct values
      dct <- pcr_normalize(df, reference_gene = reference_gene)
      error <- pcr_sd(dct, group_var = group_var, tidy = TRUE)
    }

  # merge data.frames and calculate intervals
  res <- full_join(delta_ct, ddct) %>%
    full_join(norm_rel) %>%
    full_join(error) %>%
    mutate(lower = 2 ^ -(caliberated + error),
           upper = 2 ^ -(caliberated - error))

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
#'   \item caliberated
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
pcr_dct <- function(df, group_var, reference_gene, reference_group, mode = 'average_ct') {
  if(mode == 'average_ct') {
    # average ct and calibrate to a reference group
    ave <- pcr_average(df, group_var = group_var)
    dct <- pcr_caliberate(ave, reference_group = reference_group)
  } else if(mode == 'average_dct') {
    # caliberate ct and average
    dct <- pcr_caliberate(df, reference_group = reference_group)
    dct <- pcr_average(dct, group_var = group_var)
  }

  # retain caliberated values
  # calculate the fold change
  calib <- gather(dct, gene, caliberated, -group) %>%
    mutate(fold_change = 2 ^ -caliberated)

  if(mode == 'average_ct') {
    # calculate the standard deviation from ct values
    sds <- pcr_sd(df, group_var = group_var, tidy = TRUE)
  } else if(mode == 'average_dct') {
    # caliberate ct values to a reference group
    # calculated sd from caliberated values
    dct <- pcr_caliberate(df, reference_group = reference_group)
    sds <- pcr_sd(dct, group_var = group_var, tidy = TRUE)
  }

  # join data frame and calculate intervals
  res <- full_join(calib, sds) %>%
    mutate(lower = 2 ^ -(caliberated + error),
           upper = 2 ^ -(caliberated - error))

  return(res)
}

#' Calculate the standard_curve
#'
#' @inheritParams pcr_ddct
#' @param mode A character string; 'average_amounts' or 'average_normalized'
#' @param intercept A numerice vector of intercepst and length equals the number of genes
#' @param slope A numerice vector of slopes length equals the number of genes
#'
#' @return A data.frame of 7 columns
#' \itemize{
#'   \item group
#'   \item gene
#'   \item normalized
#'   \item caliberated
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
#' standard_curve <- pcr_assess(ct3, amount = amount, mode = 'standard_curve')
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
pcr_curve <- function(df, group_var, reference_gene, reference_group, mode = 'average_amounts',
                      intercept, slope) {
  # calculate the amount of rna in samples
  amounts <- pcr_amount(df,
                        intercept = intercept,
                        slope = slope)
  if(mode == 'average_amounts') {
    # average amounts and normalize by a reference_gene
    ave <- pcr_average(amounts, group_var = group_var)
    norm <- pcr_normalize(ave, reference_gene = reference_gene, mode = 'divide')
  } else if(mode == 'average_normalized') {
    # normalize amounts and average
    norm <- pcr_normalize(amounts, reference_gene = reference_gene, mode = 'divide')
    norm <- pcr_average(norm, group_var = group_var)
  }

  # retain normalized amounts
  normalized <- gather(norm, gene, normalized, -group)

  # caliberate to a reference_group
  calib <- pcr_caliberate(norm, reference_group = reference_group,
                          mode = 'divide', tidy = TRUE)

  if(mode == 'average_amounts') {
    # calculate cv from amounts
    cv <- pcr_cv(amounts, group_var = group_var, mode = mode)
    error <- pcr_error(cv, reference_gene = reference_gene, tidy = TRUE)
  } else if(mode == 'average_normalized') {
    # calculate cv from normalized amounts
    norm <- pcr_normalize(amounts, reference_gene = reference_gene, mode = 'divide')
    error <- pcr_cv(norm, group_var = group_var, tidy = TRUE)
  }
  # join data.frames and calculate intervals
  res <- full_join(normalized, calib) %>%
    full_join(error) %>%
    mutate(lower = caliberated - error,
           upper = caliberated + error,
           error = error * normalized)

  return(res)
}
