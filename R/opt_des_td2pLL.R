#' @title Calculate D-optimal design for td2pLL model
#' @description Calculate a (local) D-optimal design for a td2pLL model using
#'  particle swarm optimization.
#' @inheritParams psoOptDesign
#' @param theta (`numeric(4)`) \cr
#'  Parameter vector of the (non-linear) td2pLL model parameters in the order
#'   h, delta, gamma and c0.
#' @examples
#' res <- opt_des_td2pLL(theta =  c(h=2, delta=0.2, gamma=1.3, c0=0.2))
#' res
#' @export
opt_des_td2pLL <- function(nPoints = 9, theta, control = NULL, Lb = c(1, 0),
                           Ub = c(10, 1), wFixed = NULL) {

  if(is.null(control)) control <- list(numIt = 200, numPart = 300, setProgressBar= TRUE)

  res <- psoOptDesign(crit = Dcrit, nPoints = nPoints, dimension = 2, control = control,
                      wFixed = wFixed, Lb = Lb, Ub = Ub, gradient = grad_td2pLL_pso,
                      theta = theta)

  return(res)

}
