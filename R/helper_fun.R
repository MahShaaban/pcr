#' Get vector average by a variable
#'
#' @param vec A vector of numerics
#' @param var A grouping variable
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c('group1', 'group2'), 3)
#' pcr:::.pcr_average(vec, var)
#'
#' @importFrom stats aggregate

.pcr_average <- function(vec, var) {
  res <- aggregate(vec,
                   by = list(var),
                   FUN = mean)
  res <- res[order(unique(var)),]
  return(res$x)
}

#' Get vector standard deviation by a variable
#'
#' @inheritParams .pcr_average
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c('group1', 'group2'), 3)
#' pcr:::.pcr_sd(vec, var)
#'
#' @importFrom stats aggregate sd

.pcr_sd <- function(vec, var) {
  res <- aggregate(vec,
                   by = list(var),
                   FUN = sd)
  res <- res[order(unique(var)),]
  return(res$x)
}

#' Get vector coefficient of variance by a variable
#'
#' @inheritParams .pcr_average
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c('group1', 'group2'), 3)
#' pcr:::.pcr_cv(vec, var)
#'
#' @importFrom stats aggregate sd

.pcr_cv <- function(vec, var) {
  res <- aggregate(vec,
                   by = list(var),
                   FUN = function(x) sd(x)/mean(x))
  res <- res[order(unique(var)),]               
  return(res$x)
}

#' Normalize vector by another
#'
#' @inheritParams .pcr_average
#' @param ref A numeric vector
#' @param mode Either 'subtract' or 'divide'
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' ref <- rnorm(6, 30, .1)
#' pcr:::.pcr_normalize(vec, ref)
#' pcr:::.pcr_normalize(vec, ref, mode = 'divide')
#'
#' @importFrom stats aggregate sd

.pcr_normalize <- function(vec, ref, mode = 'subtract') {
  if (mode == 'subtract') {
    res <- vec - ref
  } else if (mode == 'divide') {
    res <- vec / ref
  } else {
    stop("mode can be one of 'subtract' or 'divide'.")
  }
  return(res)
}

#' Propage two vectors
#'
#' @inheritParams .pcr_average
#' @inheritParams .pcr_normalize
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' ref <- rnorm(6, 30, .1)
#' pcr:::.pcr_error(vec, ref)

.pcr_error <- function(vec, ref) {
  res <- sqrt(vec^2 + ref^2)
  return(res)
}

#' Calculate the amounts
#'
#' @inheritParams .pcr_average
#' @param a A numeric
#' @param b A numeric
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' pcr:::.pcr_amount(vec, 1, 1)

.pcr_amount <- function(vec, a, b) {
  res <- 10 ^ ((vec - a)/b)
  return(res)
}

#' Raise two to a vector power
#'
#' @inheritParams .pcr_average
#'
#' @return A numeric
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' pcr:::.pcr_relative(vec)

.pcr_relative <- function(vec, amp_eff, col_name) {
  if (amp_eff) {
      a_e <- amp_eff %>% select(contains(col_name))
      res <- a_e[1,1] ^ (-vec)
      return(res)
  } else {
      res <- 2 ^ (-vec)
      return(res)
  }
}                  
                   
#' Calculate R squared
#'
#' @inheritParams .pcr_average
#'
#' @return A numeric
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c(.1, .5), 3)
#' pcr:::.pcr_rsquared(vec, var)
#'
#' @importFrom stats cor

.pcr_rsquared <- function(vec, var) {
  if(anyNA(vec)){
    warning(paste0(sum(is.na(vec)),
                   " NAs detected. ",
                   "Ensure samples are still in the dynamic range"))
  }
  res <- cor(vec,
             log10(var),
             use = "complete.obs")^2
  return(res)
}

#' Calculate the intercept of a line
#'
#' @inheritParams .pcr_average
#'
#' @return A numeric
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c(.1, .5), 3)
#' pcr:::.pcr_intercept(vec, var)
#'
#' @importFrom stats lm coefficients

.pcr_intercept <- function(vec, var) {
  ll <- lm(vec ~ log10(var))
  res <- coefficients(ll)
  return(res[[1]])
}

#' Calculate the slope of a line
#'
#' @inheritParams .pcr_average
#'
#' @return A numeric
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c(.1, .5), 3)
#' pcr:::.pcr_slope(vec, var)
#'
#' @importFrom stats lm coefficients

.pcr_slope <- function(vec, var) {
  ll <- lm(vec ~ log10(var))
  res <- coefficients(ll)
  return(res[[2]])
}
