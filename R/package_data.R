#' \eqn{C_T} values from qPCR (separate tubes)
#'
#' A dataset containing the \eqn{C_T} values of two genes from a qPCR experiment.
#' Samples were prepared from human tissues; Brain and kidney (n = 6) each.
#' Primers for each genes were run in separate reaction tubes.
#'
#' @format A data.frame with 12 rows and 2 variables:
#' \describe{
#'   \item{c_myc}{\eqn{C_T} values of the target gene \strong{c-myc}}
#'   \item{GAPDH}{\eqn{C_T} values of the control gene \strong{GAPDH}}
#' }
#'
#' @seealso \code{\link{ct2}}
#' @seealso \code{\link{ct3}}
#'
#' @source \url{http://www3.appliedbiosystems.com/cms/groups/mcb_support/documents/generaldocuments/cms_040980.pdf}
"ct1"

#' \eqn{C_T} values from qPCR (same tubes)
#'
#' A dataset containing the \eqn{C_T} values of two genes from a qPCR experiment.
#' Samples were prepared from human tissues; Brain and kidney (n = 6) each.
#' Primers for both genes were run in the same tubes with different reporting dyes.
#'
#' @format A data.frame with 12 rows and 2 variables:
#' \describe{
#'   \item{c_myc}{\eqn{C_T} values of the target gene \strong{c-myc}}
#'   \item{GAPDH}{\eqn{C_T} values of the control gene \strong{GAPDH}}
#' }
#'
#' @seealso \code{\link{ct1}}
#' @seealso \code{\link{ct3}}
#'
#' @source \url{http://www3.appliedbiosystems.com/cms/groups/mcb_support/documents/generaldocuments/cms_040980.pdf}
"ct2"

#' \eqn{C_T} values from qPCR (Serial dilutions)
#'
#' A dataset containing the \eqn{C_T} values of two genes from a serial dilution
#' qPCR experiment. The original dataset shows only the averages and standard
#' deviations of each of the 7 different diluions (1, .5, .2, .1, .05, .02 and .01).
#' These summarise were used to regenerate 3 replicates for each of the dilutions
#' to be used in testing and examples of the different functions.
#'
#'
#' @format A data.frame with 21 rows and 2 variables:
#' \describe{
#'   \item{c_myc}{\eqn{C_T} values of the target gene \strong{c-myc}}
#'   \item{GAPDH}{\eqn{C_T} values of the control gene \strong{GAPDH}}
#' }
#'
#' @seealso \code{\link{ct1}}
#' @seealso \code{\link{ct2}}
#'
#' @source \url{http://www3.appliedbiosystems.com/cms/groups/mcb_support/documents/generaldocuments/cms_040980.pdf}
"ct3"

#' \eqn{C_T} values from qPCR (Serial dilutions)
#'
#' A dataset containing the \eqn{C_T} values of two genes from a controlled serial
#'  dilution qPCR experiment. The data were prepaired from four different dilutions
#'  (10, 2, 0.4 and 0.08) and two control groups; control and treatment (n = 12) each.
#'
#'
#' @format A data.frame with 24 rows and 2 variables:
#' \describe{
#'   \item{ref}{\eqn{C_T} values of the reference gene}
#'   \item{target}{\eqn{C_T} values of the target gene}
#' }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1395339/}
"ct4"
