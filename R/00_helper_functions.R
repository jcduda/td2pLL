


#' @title Accessing EC50 values of a td2pLL model
#'
#' @description For parameters of a td2pLL model
#' and given `time` values, the EC50 value is retrieved.
#' @param coefs Named numeric vector containing the parameters of a td2pLL
#' model containing the coefficients h, delta, gamma and c0.
#' For an object of class td2pLL, these can be accessed via `coef()`.
#' @param times Numeric vector of the times for which the EC50 shall be retrieved.
#' Entries must be greater than 0.
#' @return A data.frame with columns time and EC50 containing the EC50s at
#'  `times`.


get_EC50s <- function(coefs, times) {
  stopifnot(length(setdiff(names(coefs), c("h", "delta", "gamma", "c0"))) == 0)
  stopifnot(is.numeric(times) & length(times) >= 1)
  stopifnot(all(times >= 0))

  EC50s <- sapply(times, function(time) coefs["delta"] * time^(-coefs["gamma"]) + coefs["c0"])

  res <- data.frame(time = times, EC50 = EC50s)

  return(res)
}
