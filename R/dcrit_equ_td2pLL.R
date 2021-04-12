#' @title D-optimality equivalence theorem for multidimensional design
#' @param M_x (`numeric()`) \cr
#' A single support point
#' @param des (`list()`) \cr
#' A calculated design structured as the [opt_des_td2pLL] output
#' @param finfo_inv (optimal `matrix()`)
#'  Inverse of Fisher Information matrix.
#' @inheritParams opt_des_td2pLL

dcrit_equ_td2pLL = function(M_x, des, theta, finfo_inv = NULL)
{
  w = des$weights
  M = des$supPoints
  gradient <- grad_td2pLL_pso

  if(is.null(dim(M_x)) & length(M_x) == 2) M_x <- matrix(M_x, nrow = 1)

  if(!is.null(finfo_inv)) {
    minfo0 <- finfo_inv
  } else {
  info0 = finfo(w, M, gradient, theta)
  minfo0 = solve(info0)
  }

  result = t(gradient(M_x,theta)) %*% minfo0 %*% gradient(M_x,theta) -nrow(minfo0)
  return(result)
}

