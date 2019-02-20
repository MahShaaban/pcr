#' Calculate m values for reference genes
#'
#' Calculates the stability measure for a number of reference genes
#'
#' @param df A \code{data.frame} of the expression values.
#' @param group_var A \code{character} vector of the groups of samples of
#' length equals the number of rows in \code{df}.
#' @param reference_gene A \code{character} vector of the names of reference
#' genes in \code{df} column names.
#'
#' @return A \code{data.frame} of a number of colums equals the length of
#' \code{reference_gene} plus a \code{group} column and a number of rows equals
#' the number of unique values in \code{group_var}.
#'
#' @examples
#' ## locate and read relative expression of reference genes
#' fl <- system.file('extdata', 'hk_tissue.csv', package = 'pcr')
#' hk_tissue <- readr::read_csv(fl)
#'
#' # make a grouping variable and names of reference genes
#' group_var <- ifelse(grepl('POOL', hk_tissue$tissue),
#'     'Normal',
#'      hk_tissue$tissue)
#' reference_gene <- names(hk_tissue)[-11]
#'
#' # calculate m values for all reference genes in all groups
#' pcr_m(hk_tissue, group_var, reference_gene)
#'
#' @importFrom stats aggregate
#'
#' @export
pcr_m <- function(df, group_var, reference_gene) {
  # define values for the loop
  ## l is an empty list to collect m values
  l <- list()

  # loop over all reference genes
  ## calculate the m value for each reference gene
  for(i in seq_along(reference_gene)) {
    # get name of a reference gene
    cl <- reference_gene[i]

    # get names of other reference genes
    ot <- setdiff(reference_gene, cl)

    # calculate the log2 of the expression ratio
    r <- log2(df[, ot]/unlist(df[, cl]))

    # l is an empty list to collect values
    ll <- list()

    # loop over ratios and calculate sd
    for(j in seq_along(names(r))) {
      # get the names of the numerator
      c <- names(r)[j]

      # calculate sd for the ratio in each group
      s <- aggregate(r[, c], list(group_var), sd)[, 2]

      # add sd to the empty list
      ll[[c]] <- s
    }
    l[[cl]] <- rowMeans(as.data.frame(ll))
  }

  # reshape results as data.frame
  res <- as.data.frame(l)

  # add a group column with the group names
  res$group <- unique(group_var)[order(unique(group_var))]

  # return results
  return(res)
}

#' Calculate stability of a combination of reference genes
#'
#' Calculate the stability of a decreasing number of combinations of
#' reference genes by removing the least stable reference gene across all
#' groups.
#'
#' @inheritParams pcr_m
#'
#' @return A \code{data.frame} of a number of columns equals the number of
#' uniqure entries in \code{group_var} plus two other columns; \code{gene}
#' and \code{n} with the names and numbers of the combinations of
#' reference genes used.
#'
#' @examples
#' ## locate and read relative expression of reference genes
#' fl <- system.file('extdata', 'hk_tissue.csv', package = 'pcr')
#' hk_tissue <- readr::read_csv(fl)
#'
#' # make a grouping variable and names of reference genes
#' group_var <- ifelse(grepl('POOL', hk_tissue$tissue),
#'     'Normal',
#'      hk_tissue$tissue)
#' reference_gene <- names(hk_tissue)[-11]
#'
#' # calculate m values for all reference genes in all groups
#' pcr_stability(hk_tissue, group_var, reference_gene)
#'
#' @export
pcr_stability <- function(df, group_var, reference_gene) {
  # stop if the length of reference_gene less than two
  if(length(reference_gene) < 2) {
    stop("pcr_stability require two or more reference genes.")
  }

  # define values for the loop
  ## i as the number of reference genes
  i <- length(reference_gene)

  ## left as the names of reference genes
  left <- reference_gene

  ## ll an empty list to collect stability measures
  ll <- list()

  # loop over reference genes
  ## calculate m values for each group
  ### remove least stable
  #### return average stability for decreasing number of combinations
  while (i > 1) {
    # calculate stability for a given number of reference genes
    m <- pcr_m(df,
               group_var = group_var,
               reference_gene = left)

    # decrease number of reference genes by one
    i <- i-1

    # get average stability for each reference gene in all groups
    ## remove least stable
    ### reformate as data.frame and,
    #### add group as column names and,
    #### add remaining genes in a column
    v <- rowMeans(m[, -i-2])
    v <- as.data.frame(t(v))
    names(v) <- m$group
    v$gene <- paste(left, collapse = '/')

    # return average stability for the number of reference genes
    ll[[i]] <- v

    # remove the reference gene with least stability in all groups
    left <- names(colMeans(m[, -i-2]))[-1]
  }

  # bind lists and reshape as data.frame
  ## add the number of reference genes to a column n
  res <- do.call('rbind', ll)
  res <- as.data.frame(res)
  res$n <- 2:(nrow(res)+1)

  # return data.frame
  return(res)
}

pcr_sne <- function(df, reference_gene) {
  l1 <- list()
  for(j in 1:ncol(df)) {
    v1 <- df[, j]
    v2 <- df[, -j]

    r1 <- v2/unlist(v1)
    l1[[j]] <- r1
  }

  l2 <- list()
  for(p in 1:nrow(df)) {
    v1 <- df[p,]
    v2 <- df[-p,]

    r2 <- apply(v2, 1, function(x) x/unlist(v1))
    l2[[p]] <- r2
  }

  return(l2)
}
