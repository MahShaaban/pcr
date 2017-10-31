#' Statistical testing of PCR data
#'
#' A unified interface to different statistical significance tests for qPCR data
#'
#' @inheritParams pcr_ddct
#' @param test A character string; 't.test' default, 'wilcox.test' or 'lm'
#' @param model_matrix A model matrix for advanced experimental design. for
#' constructing such a matrix with different variables check
#' \code{\link[stats]{model.matrix}}
#' @param ... Other arguments for the testing methods
#'
#' @return A data.frame of 5 columns in addition to term when test == 'lm'
#' \itemize{
#'   \item term The linear regression comparison terms
#'   \item gene The column names of df. reference_gene is dropped
#'   \item estimate The estimate for each term
#'   \item p_value The p-value fot each term
#'   \item lower The low 95\% confidence interval
#'   \item upper The high 95\% confidence interval
#' }
#' For details about the test methods themselves and different parameters,
#' consult \code{\link[stats]{t.test}}, \code{\link[stats]{wilcox.test}}
#' and \code{\link[stats]{lm}}
#'
#' @details The simple t-test can be used to test the significance of the
#' difference between two condtions \eqn{\Delta C_T}. t-test assumes in addition,
#'  that the input \eqn{C_T} values are normally distributed and the variance
#'  between condtions are comprable.
#' Wilcoxon test can be used when sample size is samll and those two last a
#' ssumpiton are hard to achieve.
#'
#' Two use the linear regression here. A null hypothesis is formulated as following,
#' \deqn{
#'   C_{T, target, treatment} - C_{T, control, treatment} =
#'   C_{T, target, control} - C_{T, control, control}
#'   \quad \textrm{or} \quad  \Delta\Delta C_T
#' }
#' This is exactly the \eqn{\Delta\Delta C_T} as explained earlier. So the
#' \eqn{\Delta\Delta C_T} is estimated and the null is rejected when
#' \eqn{\Delta\Delta C_T \ne 0}.
#'
#' @references Yuan, Joshua S, Ann Reed, Feng Chen, and Neal Stewart. 2006.
#' “Statistical Analysis of Real-Time PCR Data.” BMC Bioinformatics 7 (85).
#' BioMed Central. doi:10.1186/1471-2105-7-85.
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
#'
#' # testing using lm
#' pcr_test(ct4,
#'          group_var = group,
#'          reference_gene = 'ref',
#'          reference_group = 'control',
#'          test = 'lm')
#'
#' # testing advanced designs using a model matrix
#' # make a model matrix
#' group <- relevel(factor(group), ref = 'control')
#' dose <- rep(c(10, 2, 0.4, .08), each = 3, times = 2)
#' mm <- model.matrix(~group:dose, data = data.frame(group, dose))
#'
#' # test using t-test
#' pcr_test(ct4,
#'          reference_gene = 'ref',
#'          model_matrix = mm,
#'          test = 'lm')
#'
#' # using linear models to check the effect of RNA quality
#' # make a model matrix
#' group <- relevel(factor(group), ref = 'control')
#' set.seed(1234)
#' quality <- scale(rnorm(n = 24, mean = 1.9, sd = .1))
#' mm <- model.matrix(~group + group:quality, data = data.frame(group, quality))
#'
#' # testing using lm
#' pcr_test(ct4,
#'          reference_gene = 'ref',
#'          model_matrix = mm,
#'          test = 'lm')
#' # using linear model to check the effects of mixing separate runs
#' # make a model matrix
#' group <- relevel(factor(group), ref = 'control')
#' run <- factor(rep(c(1:3), 8))
#' mm <- model.matrix(~group + group:run, data = data.frame(group, run))
#'
#' # test using lm
#' pcr_test(ct4,
#'          reference_gene = 'ref',
#'          model_matrix = mm,
#'          test = 'lm')
#'
#' @importFrom purrr map
#' @importFrom dplyr data_frame
#' @importFrom broom tidy
#' @importFrom stats t.test wilcox.test lm confint relevel model.matrix
#'
#' @export
pcr_test <- function(df, group_var, reference_gene, reference_group,
                     test = 't.test', model_matrix = NULL, ...) {
  # calculate the delta_ct values
  norm <- .pcr_normalize(df, reference_gene = reference_gene)

  # adjust group_var for formuls
  if(test == 'lm' & is.null(model_matrix)) {
    group_var <- relevel(factor(group_var), ref = reference_group)
  } else if (test != 'lm') {
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
                    t_test <- t.test(x ~ group_var, ...)
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
                    wilcox_test <- wilcox.test(x ~ group_var, conf.int = TRUE, ...)
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
                    if(is.null(model_matrix)) {
                      linear_model <- lm(x ~ group_var, ...)
                      conf_int <- confint(linear_model, ...)
                    } else {
                      linear_model <- lm(x ~ model_matrix + 0, ...)
                      conf_int <- confint(linear_model, ...)
                    }

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
