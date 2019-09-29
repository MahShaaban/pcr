.pcr_average <- function(vec, var) {
  res <- aggregate(vec,
                   by = list(var),
                   FUN = mean)
  return(res$x)
}

vec <- rnorm(6, 30, 1)
var <- rep(c('group1', 'group2'), 3)
.pcr_average(vec, var)

.pcr_sd <- function(vec, var) {
  res <- aggregate(vec,
                   by = list(var),
                   FUN = sd)
  return(res$x)
}

vec <- rnorm(6, 30, 1)
var <- rep(c('group1', 'group2'), 3)
.pcr_sd(vec, var)

.pcr_cv <- function(vec, var) {
  res <- aggregate(vec,
                   by = list(var),
                   FUN = function(x) sd(x)/mean(x))
  return(res$x)
}

vec <- rnorm(6, 30, 1)
var <- rep(c('group1', 'group2'), 3)
.pcr_cv(vec, var)

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

vec <- rnorm(6, 30, 1)
ref <- rnorm(6, 30, .1)
.pcr_normalize(vec, ref)
.pcr_normalize(vec, ref, mode = 'divide')

.pcr_error <- function(vec, ref) {
  res <- sqrt(vec^2 + ref^2)
  return(res)
}

vec <- rnorm(2, 30, 1)
ref <- rnorm(2, 30, .1)
.pcr_error(vec, ref)

.pcr_amount <- function(vec, a, b) {
  res <- 10 ^ ((vec - a)/b)
  return(res)
}

vec <- rnorm(6, 30, 1)
.pcr_amount(vec, 1, 1)

.pcr_relative <- function(vec) {
  res <- 2 ^ (-vec)
  return(res)
}

vec <- rnorm(6, 30, 1)
.pcr_amount(vec, 1, 1)

.pcr_rsquared <- function(vec, var) {
  res <- cor(vec, log10(var))^2
  return(res)
}

vec <- rnorm(6, 30, 1)
var <- rep(c(.1, .5), 3)
.pcr_rsquared(vec, var)

.pcr_intercept <- function(vec, var) {
  ll <- lm(vec ~ log10(var))
  res <- coefficients(ll)
  return(res[[1]])
}

vec <- rnorm(6, 30, 1)
var <- rep(c(.1, .5), 3)
.pcr_intercept(vec, var)

.pcr_slop <- function(vec, var) {
  ll <- lm(vec ~ log10(var))
  res <- coefficients(ll)
  return(res[[2]])
}

vec <- rnorm(6, 30, 1)
var <- rep(c(.1, .5), 3)
.pcr_slop(vec, var)
