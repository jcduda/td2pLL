
#' @export
#' @title Accessing ED50 values of a td2pLL model
#'
#' @description For parameters of a td2pLL model
#'  and given `time` values, the ED50 value is retrieved, i.e.
#'  the dose at which 50 percent of the maximal effect is reached.
#' @param coefs (Named `numeric(4)`)\cr
#'  Contains the parameters of a td2pLL, namely `h`, `delta`, `gamma` and `c0`.
#'  For an object of class `td2pLL`, these can be accessed via `coef()`.
#' @param times (`numeric()`)\cr
#'  Times at which the ED50 shall be retrieved.
#'  Entries must be greater than 0.
#' @return (`data.frame()`)\cr
#'  A data.frame with columns `time` and `ED50`.


get_ED50s <- function(coefs, times) {
  stopifnot(length(setdiff(names(coefs), c("h", "delta", "gamma", "c0"))) == 0)
  stopifnot(is.numeric(times) & length(times) >= 1)
  stopifnot(all(times >= 0))

  ED50s <- sapply(times, function(time) coefs["delta"] * time^(-coefs["gamma"]) + coefs["c0"])

  res <- data.frame(time = times, ED50 = ED50s)

  return(res)
}
