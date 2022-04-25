#' @export
#' @title Gradient function of the td2pLL model
#' @description For details on the td2pLL model, see [fit_td2pLL]
#' @param time (`numeric(1)`)
#' @param dose (`numeric(1)`)
#' @param h (`numeric(1)`)
#' @param delta (`numeric(1)`)
#' @param gamma (`numeric(1)`)
#' @param c0 (`numeric(1)`)
#' @return `numeric(5)` Gradient vector of the td2pLL model at given
#'  parameters. \cr
#'  Note that the first entry is a 1 for the intercept. In the the [td2pLL]
#'  model, it is assumed to be 100. However, for initial normalization,
#'  it is important to consider this parameter as unknown for optimal
#'  design calculations.

grad_td2pLL <- function(time, dose, h, delta, gamma, c0){
  # To avoid log(0) throwing an error
  lg2 <- function(x) ifelse(x == 0, 0, log(x))

  B <- delta * time^(-gamma) + c0
  const1 <- -100*dose^h / (B^h + dose^h)^2

  # h
  g1 <- const1 * B^h * lg2(dose / B)
  # delta
  g2 <- const1 * (-1) * h * B^(h-1) * time^(-gamma)
  # gamma
  g3 <- const1 * B^(h-1) * h * delta * lg2(time) * time^(-gamma)
  # c0
  g4 <- const1 * (-1) * B^(h-1) * h

  cbind(e0 = 1, h = g1, delta = g2, gamma = g3, c0 = g4)
}

#' @export
#' @title Gradient function of the td2pLL model when Emax=-100 is NOT known
#' @description For details on the td2pLL model, see [fit_td2pLL]
#' @param time (`numeric(1)`)
#' @param dose (`numeric(1)`)
#' @param h (`numeric(1)`)
#' @param delta (`numeric(1)`)
#' @param gamma (`numeric(1)`)
#' @param c0 (`numeric(1)`)
#' @param Emax (`numeric(1)`)
#' @return `numeric(6)` Gradient vector of the td2pLL model at given
#'  parameters. \cr
#'  Note that the first entry is a 1 for the intercept. In the the [td2pLL]
#'  model, it is assumed to be 100. However, for initial normalization,
#'  it is important to consider this parameter as unknown for optimal
#'  design calculations.
#'  The last entry is not the derivate by `Emax` which is assumed fixed at -100
#'  in the [td2pLL] model.
#'  It is know assumed unknown so that in the optimal design, measurements
#'  for larger doses are also recommended.

grad_td2pLL_noEmax <- function(time, dose, h, delta, gamma, c0, Emax){
  # To avoid log(0) throwing an error
  lg2 <- function(x) ifelse(x == 0, 0, log(x))

  B <- delta * time^(-gamma) + c0
  const1 <- Emax*dose^h / (B^h + dose^h)^2

  # h
  g1 <- const1 * B^h * lg2(dose / B)
  # delta
  g2 <- const1 * (-1) * h * B^(h-1) * time^(-gamma)
  # gamma
  g3 <- const1 * B^(h-1) * h * delta * lg2(time) * time^(-gamma)
  # c0
  g4 <- const1 * (-1) * B^(h-1) * h
  # Emax
  g5 <- dose^h / (B^h + dose^h)

  cbind(e0 = 1, h = g1, delta = g2, gamma = g3, c0 = g4, Emax = g5)
}



