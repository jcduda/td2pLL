#' @title Gradient function of the td2pLL model
#' @description For details on the td2pLL model, see [fit_td2pLL]
#' @param time (`numeric(1)`)
#' @param dose (`numeric(1)`)
#' @param h (`numeric(1)`)
#' @param delta (`numeric(1)`)
#' @param gamma (`numeric(1)`)
#' @param c0 (`numeric(1)`)
#' @return `numeric(4)` Gradient vector of the td2pLL model at given
#'  parameters.

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

  cbind(h = g1, delta = g2, gamma = g3, c0 = g4)
}


#' @title Gradient of td2pLL for PSO
#' @description Wrapper function of [grad_td2pLL] to pass
#' it to the [psoOptDes] function as argument gradient.

grad_td2pLL_pso <- function(x, theta){
  t(grad_td2pLL(time = x[1], dose = x[2],
                h = theta[1], delta = theta[2], gamma = theta[3], c0 = theta[4])
  )
}
