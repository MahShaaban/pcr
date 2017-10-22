#' Calculate average ct value
#'
#' Takes a data.frame of raw ct values and returns the averages in different
#' groups
#'
#' @param df A data.frame of numeric ct value of n genes in columns and m
#' samples in raws.
#' @param group_var A vector of length m as a grouping variables of samples
#'
#' @return A data.frame of ncol n and nrow equalse the number of unique
#' grouping variables containing the averages of ct values of each gene in each
#' group
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'pcr1_ct.csv', package = 'pcr')
#' pcr1_ct <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate averages
#' pcr_ave(pcr1_ct, group_var = group_var)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate group_by summarise_all
#'
#' @export
pcr_ave <- function(df, group_var) {
  if(length(group_var) != nrow(df)) {
    stop("group_var should be a vector of length equals nrow of df")
  }

  mutate(df, group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) mean(x))
}

#' Normalize ct values to a reference gene
#'
#' Takes a data.frame of ct values or average ct values for each gene in each
#' condition and returns the ct values normalized to a reference gene.
#'
#' @param df A data.frame of ct or average values for each gene in each
#' conditon in addition to a grouping variable such as the output of
#' \link{pcr_ave}
#' @param reference_gene A character string of the name of the column
#' correspoinding to the reference gene
#'
#' @return A data.frame of normalized ct values for each gene and a
#' grouping variable. The column corresponding to the reference gene is
#' dropped.
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'pcr1_ct.csv', package = 'pcr')
#' pcr1_ct <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate averages
#' ave <- pcr_ave(pcr1_ct, group_var = group_var)
#'
#' # calculate delta ct
#' pcr_norm(ave, 'GAPDH')
#'
#' @importFrom dplyr select mutate_if starts_with
#'
#' @export
pcr_norm <- function(df, reference_gene) {
  if(!reference_gene %in% names(df)) {
    stop("reference gene should be a character string and one of the column names")
  }

  ref <- select(df, reference_gene) %>% unlist(use.names = FALSE)

  select(df, -starts_with(reference_gene)) %>%
    mutate_if(is.numeric, function(x) x - ref)
}

#' Caliberate ct values to a reference group
#'
#' Takes a data.frame of the ct or normalized ct values for each gene in each
#' condition and returns a data.frame of the caliberated ct value of each gene
#' e.x. (delat ct of gene of interest - delta ct value of same gene in
#' reference group)
#'
#' @param df A data.frame of ct or normalized ct values for each gene and a
#' grouping variable such as the output of the \link{pcr_norm}
#' @param reference_group A character string of the reference group as it is
#' recorded in the grouping variable
#'
#' @return A data.fram of the caliberated ct values for each gene in a grouping
#' variable by a reference group
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'pcr1_ct.csv', package = 'pcr')
#' pcr1_ct <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate averages
#' ave <- pcr_ave(pcr1_ct, group_var = group_var)
#'
#' # calculate delta ct
#' dct <- pcr_norm(ave, 'GAPDH')
#'
#' # calculate delta delta ct
#' pcr_calib(dct, 'brain')
#'
#' @importFrom dplyr filter select mutate_if
#' @export
pcr_calib <- function(df, reference_group) {
  if(!reference_group %in% unlist(df$group)) {
    stop("reference_group should be a character string and one of the group variable entries")
  }

  ref <- filter(df, group == reference_group) %>%
    select(-group) %>%
    unlist(use.names = FALSE)

  mutate_if(df, is.numeric, function(x) x - ref)
}

#' Calculate standard error value
#'
#' @inheritParams pcr_ave
#'
#' @return A data.frame of ncol n and nrow equalse the number of unique
#' grouping variables containing the standard error values of each gene in each
#' group
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'pcr1_ct.csv', package = 'pcr')
#' pcr1_ct <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate standard deviations
#' pcr_sd(pcr1_ct, group_var = group_var)
#'
#' @importFrom dplyr mutate group_by summarise_all
#' @importFrom stats sd
#'
#' @export
pcr_sd <- function(df, group_var) {
  if(length(group_var) != nrow(df)) {
    stop("group_var should be a vector of length equals nrow of df")
  }

  mutate(df, group = group_var) %>%
    group_by(group) %>%
    summarise_all(function(x) sd(x))
}

#' Calculate error value
#'
#' @param df A data.frame of ncol n and nrow equalse the number of unique
#' grouping variables containing the standard error values of each gene in each
#' group
#' @inheritParams pcr_norm
#'
#' @return A data.frame of error values for each gene and a grouping
#' variable. The column corresponding to the reference gene is dropped.
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'pcr1_ct.csv', package = 'pcr')
#' pcr1_ct <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate standard deviations
#' sds <- pcr_sd(pcr1_ct, group_var = group_var)
#'
#' # calculate errors
#' pcr_error(sds, reference_gene = 'GAPDH)
#'
#' @importFrom dplyr select starts_with mutate_if
#'
#' @export
pcr_error <- function(df, reference_gene) {
  if(!reference_gene %in% names(df)) {
    stop("reference gene should be a character string and one of the column names")
  }

  ref <- select(df, reference_gene) %>% unlist(use.names = FALSE)

  select(df, -starts_with(reference_gene)) %>%
    mutate_if(is.numeric, function(x) sqrt((x^2) + (ref)^2))
}

#' pcr_analyze
#'
#' A wrapper function to perform different analysis methods and modes
#'
#' @inheritParams pcr_ave
#' @inheritParams pcr_norm
#' @inheritParams pcr_calib
#' @param intervals A logical (default TRUE) to whether or not to calculate the
#' error intervals
#' @param mode A character string. Possible inputs are "average_ct" or
#' "average_dct"
#' @param method A character string. Default is "delta_delta_ct"
#'
#' @return A tidy data.frame of calculated expression values and errors
#'
#' @examples
#' # locate and read raw ct data
#' fl <- system.file('extdata', 'pcr1_ct.csv', package = 'pcr')
#' pcr1_ct <- readr::read_csv(fl)
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate all values and errors in one step
#' pcr_analyze(pcr1_ct,
#' group_var = rep(c('brain', 'kidney'), each = 6),
#' reference_gene = 'GAPDH',
#' reference_group = 'brain')
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr full_join mutate
#'
#' @export
pcr_analyze <- function(df, group_var, reference_gene, reference_group,
                        intervals = TRUE, method = 'delta_delta_ct',
                        mode = 'average_ct') {
  if(method == 'delta_delta_ct') {
    if(mode == 'average_ct') {
      ave <- pcr_ave(df, group_var = group_var)
      dct <- pcr_norm(ave, reference_gene = reference_gene)
    } else if(mode == 'average_dct') {
      dct <- pcr_norm(df, reference_gene = reference_gene)
      dct <- pcr_ave(dct, group_var = group_var)
    }

    ddct <- pcr_calib(dct, reference_group = reference_group)

    norm_rel <- gather(ddct, gene, ddct, -group) %>%
      mutate(norm_rel = 2 ^ -ddct)

    if(intervals == TRUE) {
      if(mode == 'average_ct') {
        sds <- pcr_sd(df, group_var = group_var)
        errors <- pcr_error(sds, reference_gene = reference_gene)
      } else if(mode == 'average_dct') {
        dct <- pcr_norm(df, reference_gene = reference_gene)
        errors <- pcr_sd(dct, group_var = group_var)
      }
      errors <- gather(errors, gene, error, -group)

      norm_rel <- norm_rel %>%
        full_join(errors) %>%
        mutate(int_lower = 2 ^ - (ddct + error),
               int_upper = 2 ^ - (ddct - error))
    }
    return(norm_rel)
  }
}

#' Caliberate genes to a reference group
#'
#' Caliberate ct values of genes to a reference group and obtain the fold-
#' change usually for reference genes
#'
#' @inheritParams pcr_ave
#' @inheritParams pcr_norm
#' @inheritParams pcr_calib
#' @inheritParams pcr_analyze
#' @param method A character string of calculation method. Default "delta_ct"
#'
#' @return A data.frame of caliberated, relative values and errors
#'
#' @examples
#' # locate and read file
#' fl <- system.file('extdata', 'pcr1_ct.csv', package = 'pcr')
#' pcr1_ct <- readr::read_csv(fl)
#'
#' # make a data.frame of two identical columns
#' pcr_hk <- data.frame(
#'   GAPDH1 = pcr1_ct$GAPDH,
#'   GAPDH2 = pcr1_ct$GAPDH
#'   )
#'
#' # add grouping variable
#' group_var <- rep(c('brain', 'kidney'), each = 6)
#'
#' # calculate caliberation
#' pcr_caliberate(pcr_hk,
#'                group_var = group_var,
#'                reference_group = 'brain')
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr mutate full_join
#'
#' @export
pcr_caliberate <- function(df, group_var, reference_group, intervals = TRUE,
                           method = 'delta_ct', mode = 'average_ct') {
  if(method == 'delta_ct') {
    if(mode == 'average_ct') {
      ave <- pcr_ave(df, group_var = group_var)
      dct <- pcr_calib(ave, reference_group = reference_group)
    } else if(mode == 'average_dct') {
      dct <- pcr_calib(df, reference_group = reference_group)
      dct <- pcr_ave(dct, group_var = group_var)
    }

    rel <- gather(dct, gene, calib, -group) %>%
      mutate(rel = 2 ^ -calib)

    if(intervals == TRUE) {
      if(mode == 'average_ct') {
        sds <- pcr_sd(df, group_var = group_var)
      } else if(mode == 'average_dct') {
        dct <- pcr_calib(df, reference_group = reference_group)
        sds <- pcr_sd(dct, group_var = group_var)
      }
      sds <- gather(sds, gene, sd, -group)

      rel <- rel %>%
        full_join(sds) %>%
        mutate(int_lower = 2 ^ - (calib + sd),
               int_upper = 2 ^ - (calib - sd))
    }
    return(rel)
  }
}

#' Assess PCR data quality
#'
#' @param df A data.frame of exactly two columns containing average ct values
#' of a reference gene and another from an experiment run with differen
#' dilutions (RNA amounts)
#' @param error An optional data.frame of dimentions equals that of df and
#' contains an error measure from the same experiments
#' @param amount A numeric vector of length equals nrow df with RNA amounts
#' @param mode A character string of assessment mode. Default "effeciency"
#' @param plot A logical default FALSE of whether plot or return the data
#'
#' @return A data.frame of 5 columns
#'
#' @examples
#' # locate and read data
#' fl <- system.file('extdata', 'pcr_dilute_ave.csv', package = 'pcr')
#' pcr_dilute_ave <- readr::read_csv(fl)
#'
#' fl <- system.file('extdata', 'pcr_dilute_error.csv', package = 'pcr')
#' pcr_dilute_error <- readr::read_csv(fl)
#'
#' # make a vector of RNA amounts
#' amount <- c(1, .5, .2, .1, .05, .02, .01)
#'
#' # calculate effeciencey
#' pcr_assess(pcr_dilute_ave,
#'            error = pcr_dilute_error,
#'            amount = amount)
#'
#' @importFrom dplyr data_frame
#'
#' @export
pcr_assess <- function(df, error = NULL, amount, mode = 'effeciency',
                      plot = FALSE) {

  log_amount <- log10(amount)

  dct <- unlist(df[, 1] - df[, 2], use.names = FALSE)

  if(!is.null(error)) {
    error <- as.numeric(sqrt((error[, 1] ^ 2) + error[, 2] ^ 2))
    int_lower = unlist(dct - error)
    int_upper = unlist(dct + error)

    effeciency <- data_frame(
      log_amount = log_amount,
      dct = dct,
      error = error,
      int_lower = int_lower,
      int_upper = int_upper
    )
  } else {
    effeciency <- data_frame(
      log_amount = log_amount,
      dct = dct
    )
  }
  return(effeciency)
}
