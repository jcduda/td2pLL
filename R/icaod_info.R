#' @title Information matrix for [locally()] assuming unknown e0 and known eMax
#' @description If e0 is assumed unknown and eMax is assumed known in the td2pLL model,
#'  the optimal design will propose a support point in dose=0, but not in the maximal dose.
#' To assume additionally that eMax is unknown, use [icaod_info_noEmax].
#' @param x (`numeric(2*7)`) \cr
#'  Contains the current support points under consideration where first all
#'  time coordinates are stored in x and then all dose coordinates.
#' @param w (`numeric(7)`) \cr
#'  Weights of the support points. `sum(w)` must be 1.
#' @param param (`numeric(5)`) \cr
#'  Contains assumed true parameters for `e0`, `h`, `delta`, `gamma`, `c0` where
#'  the value of `e0` has no effect on the optimal design calculation as it
#'  is a linear parameter.
#' @return Information matrix of dimension 5 times 5.
#'
icaod_info <- function(x, w, param){

  times <- x[1:(length(x)/2)]
  doses <- x[(length(x)/2 + 1):length(x)]
  if(is.null(names(param)))
    names(param) <- c("e0", "h", "delta", "gamma", "c0")
  M <- array(0, c(5, 5))
  for(i in 1:length(doses)){
    a <- grad_td2pLL(time = times[i], dose = doses[i], h = param["h"],
                     delta = param["delta"], gamma = param["gamma"], c0 = param["c0"])
    M <- M + t(a) %*% a * w[i]
  }

  # Unfortunately much slower:
  # M <- matrix(rowSums(apply(data.frame(time = times, dose = doses, weight = w), 1, function(x){
  #   a <- grad_td2pLL(time = x[1], dose = x[2], h = param["h"],
  #               delta = param["delta"], gamma = param["gamma"], c0 = param["c0"])
  #   t(a) %*% a * x[3]
  # })), nrow = 5, ncol = 5)

  return(M)
}


#' @title Information matrix for [locally()] assuming unknown e0 and unknown eMax
#' @description If e0 is assumed unknown and eMax is assumed unknown in the td2pLL model
#' (where usually we assume e0=100 and eMax=100 for dytotoxicity assays)
#'  the optimal design will propose a support point in dose=0 and in the maximal dose.
#' To assume that eMax is known, use [icaod_info].
#' @param x (`numeric(2*8)`) \cr
#'  Contains the current support points under consideration where first all
#'  time coordinates are stored in x and then all dose coordinates.
#' @param w (`numeric(8)`) \cr
#'  Weights of the support points. `sum(w)` must be 1.
#' @param param (`numeric(6)`) \cr
#'  Contains assumed true parameters for `e0`, `h`, `delta`, `gamma`, `c0`, `Emax` where
#'  the values of `e0` and `Emax` have no effect on the optimal design calculation
#'  because they are linear parameter.
#' @return Information matrix of dimension 6 times 6.
#'
icaod_info_noEmax <- function(x, w, param){

  times <- x[1:(length(x)/2)]
  doses <- x[(length(x)/2 + 1):length(x)]
  if(is.null(names(param)))
    names(param) <- c("e0", "h", "delta", "gamma", "c0", "Emax")
  M <- array(0, c(6, 6))
  for(i in 1:length(doses)){
    a <- grad_td2pLL_noEmax(time = times[i], dose = doses[i],
                            h = param["h"], delta = param["delta"], gamma = param["gamma"],
                            c0 = param["c0"],
                            Emax = param["Emax"])
    M <- M + t(a) %*% a * w[i]
  }
  return(M)
}

