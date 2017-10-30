#' Statistical testing of PCR data
#'
#' @inheritParams pcr_average
#' @inheritParams pcr_normalize
#' @inheritParams pcr_calibrate
#' @param test A character string; 't.test' default, 'wilcox.test' or lm
#'
#' @return A data.frame of 5 columns in addition to term when test == 'lm'
#' \itemize{
#'   \item term
#'   \item gene
#'   \item estimate
#'   \item p_value
#'   \item lower
#'   \item upper
#' }
#'
#' @examples
#' # locate and read data
#' fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
#' ct4 <- readr::read_csv(fl)
#'
#' # make group variable
#' group <- rep(c('control', 'treatment'), each = 12)
#'
#' # test using t-test
#' pcr_test(ct4,
#'          group_var = group,
#'          reference_gene = 'ref',
#'          reference_group = 'control',
#'          test = 't.test')
#'
#' # test using wilcox.test
#' pcr_test(ct4,
#'          group_var = group,
#'          reference_gene = 'ref',
#'          reference_group = 'control',
#'          test = 'wilcox.test')
#' # testing using lm
#' pcr_test(ct4,
#'          group_var = group,
#'          reference_gene = 'ref',
#'          reference_group = 'control',
#'          test = 'lm')
#'
#' @importFrom purrr map
#' @importFrom dplyr data_frame
#' @importFrom broom tidy
#' @importFrom stats t.test wilcox.test lm confint relevel
#'
#' @export
pcr_test <- function(df, group_var, reference_gene, reference_group,
                     test = 't.test') {
  # calculate the delta_ct values
  norm <- pcr_normalize(df, reference_gene = reference_gene)

  # adjust group_var for formuls
  if(test == 'lm') {
    group_var <- relevel(factor(group_var), ref = reference_group)
  } else {
    group_levels <- unique(group_var)
    if(length(group_levels) != 2) {
      stop('t.test and wilcox are only applied to two group comparisons')
    }
    ref <- group_levels[group_levels != reference_group]
    group_var <- relevel(factor(group_var), ref = ref)
  }
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
                },
                'lm' = {
                  map(norm, function(x) {
                    linear_model <- lm(x ~ group_var)
                    conf_int <- confint(linear_model)

                    linear_model <- tidy(linear_model)[-1,]
                    conf_int <- tidy(conf_int)[-1,]

                    data_frame(
                      term = linear_model$term,
                      estimate = linear_model$estimate,
                      p_value = linear_model$p.value,
                      lower = conf_int$`X2.5..`,
                      upper = conf_int$`X97.5..`
                    )
                  })
                })
  bind_rows(dat, .id = 'gene')
}
