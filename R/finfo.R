#' @export
#' @title Fisher Information Matrix for multidimensional designs
#' @details Calculates the fisher information matrix for a multidimensional
#' optimal design in dependency on weights, the design, the gradient and
#' theta (as we assume a non-linear model).
#' @param w (`numeric`) \cr
#'  Weight vector with length same as rows of `M`. Weights must be within 0 and 1
#'  and must sum up to 1.
#' @param M (`matrix()`) \cr
#'  Design with support points as rows.
#' @param gradient (`function(x, theta)`)
#'   Gradient function of the model with input parameters `x` and `theta`.
#'   `x` is a `numeric()` of length same as number of columns of `M`
#'    representing one support point and \cr
#'   `theta` is a `numeric()` containing the parameters of the model.
#' @param theta `numeric()` \cr
#'  Parameter vector of the model.
#' @return Fisher information matrix as a symmetric matrix with `length(theta)`
#'  rows.
#' @examples
#' finfo(w = c(1/3, 1/3, 1/3),
#'      M = rbind(c(1.5, 0.4),
#'                c(2, 0.3),
#'                c(5, 0.9)),
#'      gradient = grad_td2pLL_pso,
#'      theta = c(h=2, delta=0.2, gamma=1.3, c0=0.2))

finfo <- function(w, M, gradient, theta){
  stopifnot(length(w) == nrow(M))
  # Add linear parameter e0
  dimTheta <- length(theta) + 1
  tmp <- matrix(0, dimTheta, dimTheta)
  for(i in 1:nrow(M)){
    tmp = tmp + w[i] * (gradient(M[i, ], theta) %*% t(gradient(M[i, ], theta)))
  }

  return(tmp)
}

#' @export
#' @title Determinant of Fisher Information Matrix
#' @description Can be used as `crit` in [psoOptDesign].
#' @inheritParams finfo
#' @return The determinant of the fisher information matrix.
#' @examples
#' Dcrit(w = c(1/4, 1/4, 1/4, 1/4),
#' M = rbind(c(1.5, 0.4),
#'          c(2,   0.3),
#'          c(5, 0.9),
#'          c(7, 1)),
#' gradient = grad_td2pLL_pso,
#' theta = c(h=2, delta=0.2, gamma=1.3, c0=0.2))

Dcrit <- function(w, M, gradient, theta){
  return((det(finfo(w, M, gradient, theta))))
}
