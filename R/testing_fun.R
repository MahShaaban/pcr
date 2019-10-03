#' Statistical testing of PCR data
#'
#' A unified interface to different statistical significance tests for qPCR
#' data
#'
#' @inheritParams pcr_ddct
#' @param test A character string; 't.test' default, 'wilcox.test' or 'lm'
#' @param ... Other arguments for the testing methods
#'
#' @return A data.frame of 5 columns in addition to term when test == 'lm'
#' \itemize{
#'   \item term The linear regression comparison terms
#'   \item gene The column names of df. reference_gene is dropped
#'   \item estimate The estimate for each term
#'   \item p_value The p-value for each term
#'   \item lower The low 95\% confidence interval
#'   \item upper The high 95\% confidence interval
#' }
#' For details about the test methods themselves and different parameters,
#' consult \code{\link[stats]{t.test}}, \code{\link[stats]{wilcox.test}}
#' and \code{\link[stats]{lm}}
#'
#' @details The simple t-test can be used to test the significance of the
#' difference between two conditions \eqn{\Delta C_T}. t-test assumes in
#' addition, that the input \eqn{C_T} values are normally distributed and the
#' variance between conditions are comparable. Wilcoxon test can be used when
#' sample size is small and those two last assumptions are hard to achieve.
#'
#' Two use the linear regression here. A null hypothesis is formulated as
#' following,
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
#' ct4 <- read.csv(fl)
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
#' dose <- rep(c(100, 80, 60, 40), each = 3, times = 2)
#' mm <- model.matrix(~group:dose, data = data.frame(group, dose))
#'
#' # test using lm
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
#'
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
#' @export
pcr_test <- function(df, test = 't.test', ...) {
  switch (test,
    't.test' = pcr_ttest(df, ...),
    'wilcox.test' = pcr_wilcox(df, ...),
    'lm' = pcr_lm(df, ...)
  )
}

#' t-test qPCR data
#'
#' @inheritParams pcr_ddct
#' @param tidy A \code{logical} whether to return a \code{list} of \code{htest}
#' or a tidy \code{data.frame}. Default TRUE.
#' @param ... Other arguments to \code{\link[stats]{t.test}}
#'
#' @return A data.frame of 5 columns
#' \itemize{
#'   \item gene The column names of df. reference_gene is dropped
#'   \item estimate The estimate for each term
#'   \item p_value The p-value for each term
#'   \item lower The low 95\% confidence interval
#'   \item upper The high 95\% confidence interval
#' }
#' When \code{tidy} is FALSE, returns a \code{list} of \code{htest} objects.
#'
#' @examples
#' # locate and read data
#' fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
#' ct4 <- read.csv(fl)
#'
#' # make group variable
#' group <- rep(c('control', 'treatment'), each = 12)
#'
#' # test
#' pcr_ttest(ct4,
#'           group_var = group,
#'           reference_gene = 'ref',
#'           reference_group = 'control')
#'
#' # test using t.test method
#' pcr_test(ct4,
#'          group_var = group,
#'          reference_gene = 'ref',
#'          reference_group = 'control',
#'          test = 't.test')
#'
#' @importFrom stats t.test relevel
#'
#' @export
pcr_ttest <- function(df, group_var, reference_gene, reference_group,
                      tidy = TRUE, ...) {
  # only if group_var is not a factor
  if (!is.factor(group_var)) {
    # adjust the reference group
    group_levels <- unique(group_var)

    if(length(group_levels) != 2) {
      stop('t.test is only applied to two group comparisons.')
    }

    ref_group <- group_levels[group_levels != reference_group]
    group_var <- relevel(factor(group_var), ref = ref_group)
  }

  # extract the reference gene and genes of interest
  ref <- subset(df, select = reference_gene, drop = TRUE)
  goi <- subset(df, select = names(df) != reference_gene)

  # apply the calculations
  tst <- apply(goi,
               MARGIN = 2,
               FUN = function(x) {
                 norm <- .pcr_normalize(x, ref)
                 t.test(norm ~ group_var, ...)
               })

  # make a tidy data.frame or return htest object
  if (tidy) {
    tst <- lapply(tst,
                  FUN = function(x) {
                    data.frame(
                      gene = '',
                      estimate = unname(x$estimate[1] - x$estimate[2]),
                      p_value = x$p.value,
                      lower = x$conf.int[1],
                      upper = x$conf.int[2]
                    )
                  })
    tst <- do.call(rbind, tst)
    tst$gene <- names(goi)
    rownames(tst) <- NULL

  }

  # return
  return(tst)
}

#' Wilcoxon test qPCR data
#'
#' @inheritParams pcr_ddct
#' @param tidy A \code{logical} whether to return a \code{list} of \code{htest}
#' or a tidy \code{data.frame}. Default TRUE.
#' @param ... Other arguments to \code{\link[stats]{wilcox.test}}
#'
#' @return A data.frame of 5 columns
#' \itemize{
#'   \item gene The column names of df. reference_gene is dropped
#'   \item estimate The estimate for each term
#'   \item p_value The p-value for each term
#'   \item lower The low 95\% confidence interval
#'   \item upper The high 95\% confidence interval
#' }
#'
#' When \code{tidy} is FALSE, returns a \code{list} of \code{htest} objects.
#'
#' @examples
#' # locate and read data
#' fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
#' ct4 <- read.csv(fl)
#'
#' # make group variable
#' group <- rep(c('control', 'treatment'), each = 12)
#'
#' # test
#' pcr_wilcox(ct4,
#'            group_var = group,
#'            reference_gene = 'ref',
#'            reference_group = 'control')
#'
#' # test using wilcox.test method
#' pcr_test(ct4,
#'          group_var = group,
#'          reference_gene = 'ref',
#'          reference_group = 'control',
#'          test = 'wilcox.test')
#'
#' @importFrom stats wilcox.test relevel
#'
#' @export
pcr_wilcox <- function(df, group_var, reference_gene, reference_group,
                       tidy = TRUE, ...) {
  # only if group_var is not a factor
  if (!is.factor(group_var)) {
    # adjust the reference group
    group_levels <- unique(group_var)

    if(length(group_levels) != 2) {
      stop('wilcox.test is only applied to two group comparisons.')
    }

    ref_group <- group_levels[group_levels != reference_group]
    group_var <- relevel(factor(group_var), ref = ref_group)
  }

  # extract the reference gene and genes of interest
  ref <- subset(df, select = reference_gene, drop = TRUE)
  goi <- subset(df, select = names(df) != reference_gene)

  # apply the calculations
  tst <- apply(goi,
               MARGIN = 2,
               FUN = function(x) {
                 norm <- .pcr_normalize(x, ref)
                 wilcox.test(norm ~ group_var, conf.int = TRUE, ...)
               })

  # make a tidy data.frame or return htest object
  if (tidy) {
    tst <- lapply(tst,
                  FUN = function(x) {
                    data.frame(
                      gene = '',
                      estimate = unname(x$estimate),
                      p_value = x$p.value,
                      lower = x$conf.int[1],
                      upper = x$conf.int[2]
                    )
                  })
    tst <- do.call(rbind, tst)
    tst$gene <- names(goi)
    rownames(tst) <- NULL
  }

  # return
  return(tst)
}

#' Linear regression qPCR data
#'
#' @inheritParams pcr_ddct
#' @param model_matrix A model matrix for advanced experimental design. for
#' constructing such a matrix with different variables check
#' \code{\link[stats]{model.matrix}}
#' @param mode A character string for the normalization mode. Possible values
#' are "subtract" (default) or "divide".
#' @param tidy A \code{logical} whether to return a \code{list} of
#' \code{\link[stats]{lm}} or a tidy \code{data.frame}. Default TRUE.
#' @param ... Other arguments to \code{\link[stats]{lm}}
#'
#' @return A data.frame of 6 columns
#' \itemize{
#'   \item term The term being tested
#'   \item gene The column names of df. reference_gene is dropped
#'   \item estimate The estimate for each term
#'   \item p_value The p-value for each term
#'   \item lower The low 95\% confidence interval
#'   \item upper The high 95\% confidence interval
#' }
#' When \code{tidy} is FALSE, returns a \code{list} of \code{\link[stats]{lm}}
#' objects.
#'
#' @examples
#' # locate and read data
#' fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
#' ct4 <- read.csv(fl)
#'
#' # make group variable
#' group <- rep(c('control', 'treatment'), each = 12)
#'
#' # test
#' pcr_lm(ct4,
#'        group_var = group,
#'        reference_gene = 'ref',
#'        reference_group = 'control')
#'
#' # testing using lm method
#' pcr_test(ct4,
#'          group_var = group,
#'          reference_gene = 'ref',
#'          reference_group = 'control',
#'          test = 'lm')
#'
#' @importFrom stats lm confint relevel
#'
#' @export
pcr_lm <- function(df, group_var, reference_gene, reference_group,
                   model_matrix = NULL, mode = 'subtract', tidy = TRUE, ...) {
  # adjust group_var for formuls
  if (missing(model_matrix)) {
    if (!is.factor(group_var)) {
      group_var <- relevel(factor(group_var), ref = reference_group)
    }
  }

  # extract the reference gene and genes of interest
  ref <- subset(df, select = reference_gene, drop = TRUE)
  goi <- subset(df, select = names(df) != reference_gene)

  # apply the calculations
  tst <- apply(goi,
               MARGIN = 2,
               FUN = function(x) {
                 norm <- .pcr_normalize(x, ref)
                 if (is.null(model_matrix)) {
                   lm(norm ~ group_var, ...)
                 } else {
                   lm(norm ~ model_matrix + 0, ...)
                 }
               })
  if (tidy) {
    tst <- lapply(tst,
                  FUN = function(x) {
                    mod <- x
                    conf_int <- confint(mod)
                    data.frame(
                      gene = '',
                      term = names(mod$coefficients)[-1],
                      estimate = unname(mod$coefficients)[-1],
                      p_value = summary(mod)$coefficients[-1, 4],
                      lower = conf_int[-1, 1],
                      upper = conf_int[-1, 2]
                    )
                  })
    tst <- do.call(rbind, tst)
    tst$gene <- names(goi)
    rownames(tst) <- NULL
  }

  # return
  return(tst)
}
