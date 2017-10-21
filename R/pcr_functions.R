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

#' Calculate delta ct value
#'
#' Takes a data.frame of average ct values for each gene in each condition and
#' returns the delta ct value (average gene of interest - average reference
#' gene)
#'
#' @param df A data.frame of average ct values for each gene in each conditon
#' in addition to a grouping variable such as the output of \link{pcr_ave}
#' @param reference_gene A character string of the name of the column
#' correspoinding to the reference gene
#'
#' @return A data.frame of delta ct values for each gene and a grouping
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
#' # calculate averages
#' ave <- pcr_ave(pcr1_ct, group_var = group_var)
#'
#' # calculate delta ct
#' pcr_dct(ave, 'GAPDH')
#'
#' @importFrom dplyr select mutate_if starts_with
#'
#' @export
pcr_dct <- function(df, reference_gene) {
  if(!reference_gene %in% names(df)) {
    stop("reference gene should be a character string and one of the column names")
  }

  ref <- select(df, reference_gene) %>% unlist(use.names = FALSE)

  select(df, -starts_with(reference_gene)) %>%
    mutate_if(is.numeric, function(x) x - ref)
}

#' Calculate delta delta ct value
#'
#' Takes a data.frame of the delta ct values for each gene in each condition
#' and returns a data.frame of the delta delta value of each gene (delat ct of
#' gene of interest - delta ct value of same gene in reference group)
#'
#' @param df A data.frame of delta ct values for each gene and a grouping
#' variable such as the output of the \link{pcr_dct}
#' @param reference_group A character string of the reference group as it is
#' recorded in the grouping variable
#'
#' @return A data.fram of the delta delta values for each gene in a grouping
#' variable
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
#' dct <- pcr_dct(ave, 'GAPDH')
#'
#' # calculate delta delta ct
#' pcr_ddct(dct, 'brain')
#'
#' @importFrom dplyr filter select mutate_if
#' @export
pcr_ddct <- function(df, reference_group) {
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
#' @inheritParams pcr_dct
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

#' Normalize PCR experment data
#'
#' A wrapper function to calculate the normalized relative values of genes of
#' interest to a reference gene and a reference group.
#'
#' @inheritParams pcr_ave
#' @inheritParams pcr_dct
#' @inheritParams pcr_ddct
#' @param intervals A logical (default TRUE) to whether or not to calculate the
#' error intervals
#' @param mode A character string. Possible inputs are "average_ct" or
#' "average_dct"
#'
#' @return A tidy data.frame of calculated normalized relative values and
#' errors
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
#' pcr_normalize(pcr1_ct,
#' group_var = rep(c('brain', 'kidney'), each = 6),
#' reference_gene = 'GAPDH',
#' reference_group = 'brain')
#'
#' @importFrom dplyr full_join mutate
#' @importFrom tidyr gather
#'
#' @export
pcr_normalize <- function(df, group_var, reference_gene, reference_group,
                          intervals = TRUE, mode = 'average_ct') {
  if(mode == 'average_ct') {
    ave <- pcr_ave(df, group_var = group_var)
    dct <- pcr_dct(ave, reference_gene = reference_gene)
  } else if(mode == 'average_dct') {
    dct <- pcr_dct(df, reference_gene = reference_gene)
    dct <- pcr_ave(dct, group_var = group_var)
  }

  ddct <- pcr_ddct(dct, reference_group = reference_group)

  norm_rel <- gather(ddct, gene, ddct, -group) %>%
    mutate(norm_rel = 2 ^ -ddct)

  if(intervals == TRUE) {
    if(mode == 'average_ct') {
      sds <- pcr_sd(df, group_var = group_var)
      errors <- pcr_error(sds, reference_gene = reference_gene)
    } else if(mode == 'average_dct') {
      dct <- pcr_dct(df, reference_gene = reference_gene)
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
