
# function: td2pLL
#' @export
#' @title Calculate response value for given td2pLL model
#' @description `td2pLL` returns the response value of a fully specified
#' time-dose two-parameter log-logistic model at a certain time and dose
#' value:
#' \deqn{f(d,t)=100-100\frac{d^h}{EC_{50}(t)^h + d^h}}
#'  with
#'  \deqn{EC_{50}(t) = \Delta \cdot t^{-\gamma} + C_0}
#' @param time numeric for the time value where the response shall be calculated
#' @param dose numeric for the dose (or concentration) value where the response
#' shall be calculated
#' @param h The `h` parameter of the model
#' @param gamma The `gamma` parameter of the model
#' @param c0 The `c0` parameter of the model
#' @param delta The `delta` parameter of the model
#' @return The response value of the model _in percent_ at the given `time` and
#' `dose` value.
#' @examples
#' td2pLL(time = 4, dose = 0.1, h = 2, gamma = 2.5, c0 = 0.1, delta = 0.3)
td2pLL <- function(time, dose, h, gamma, c0, delta) {
  ED50 <- delta * time^(-gamma) + c0
  resp <- 100 - 100 * (dose^h) / (ED50^h + dose^h)
  return(resp = resp)
}
