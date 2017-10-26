#' Statistical testing of PCR data
#'
#' @inheritParams pcr_average
#' @param test A character string; 't.test' default or 'wilcox.test'
#'
#' @return A data.frame of 5 columns
#' \itemize{
#'   \item gene
#'   \item estimate
#'   \item p_value
#'   \item lower
#'   \item upper
#' }
#'
#' @examples
#' # locate and read data
#' fl <- system.file('extdata', 'ct4_raw.tsv', package = 'pcr')
#' ct4 <- readr::read_csv(fl)(fl)
#'
#' # make group variable
#' group <- rep(c('control', 'treatment'), each = 12)
#' group <- relevel(factor(group), ref = 'treatment')
#'
#' # test using t-test
#' pcr_test(ct4,
#'          group_var = group,
#'          reference_gene = 'ref',
#'          test = 't.test')
#'
#' # test using wilcox.test
#' pcr_test(ct4,
#'          group_var = group,
#'          reference_gene = 'ref',
#'          test = 'wilcox.test')
#'
#' @importFrom purrr map
#' @importFrom dplyr data_frame
#' @importFrom broom tidy
#' @importFrom stats t.test wilcox.test
#'
#' @export
pcr_test <- function(df, group_var, reference_gene, test = 't.test') {
  # calculate the delta_ct values
  norm <- pcr_normalize(df, reference_gene = reference_gene)

  dat <- switch(test,
                't.test' = {
                  map(norm, function(x) {
                    t_test <- t.test(x ~ group_var)
                    t_test <- tidy(t_test)
                    t_test <- with(t_test,
                                   data_frame(
                                     estimate = estimate,
                                     p_value = p.value,
                                     lower = conf.low,
                                     upper = conf.high
                                   ))
                    })
                  },
                'wilcox.test' = {
                  map(norm, function(x) {
                    wilcox_test <- wilcox.test(x ~ group_var, conf.int = TRUE)
                    wilcox_test <- tidy(wilcox_test)
                    wilcox_test <- with(wilcox_test,
                                   data_frame(
                                     estimate = estimate,
                                     p_value = p.value,
                                     lower = conf.low,
                                     upper = conf.high
                                   ))
                  })
                })
  bind_rows(dat, .id = 'gene')
}
