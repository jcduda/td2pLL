

#' @export
#' @title Calculate optimal locally D-optimal design for td2pLL model
#' @description If prior knowledge on the (rough) model shape is available,
#'  an optimal experimental design for a td2pLL model can be calculated.
#'  The numerical back bone is the Imperialist Competitive Algorithm
#'  implemented by Masoudi et al (2017).
#' @param param (`numeric(4)`) \cr
#'  Assumed parameters for h, delta, gamma and c0 in this order.
#' @param Emax_known (`logical(1)`) \cr
#'  Whether or not the maximal effect is known and hence, if it is required
#'  to test for larger doses.
#'  If `FALSE` (default) then the maximal effect (which is -100 in a td2pLL model
#'  for typical cytotoxicity data) is assumed unknown and hence the optimal design
#'  will propose to make measurement at the largest possible dose in order to
#'  gather information on the unknown, maximal effect `Emax`. \cr
#'  Note that the mean control level `e0` which is assumed `e0=100` in the
#'  typical td2pLL model is, for optima design planning, always assumed to be
#'  unknown as otherwise no control measurements in dose=0 would be porposed
#'  by the optimal design.
#' @param lx (`numeric(2)`) \cr
#'   Minimal value for time and dose. Note that for dose, you cannot give
#'   a lower bound of zero due to technical reasons.
#'   Instead, provide a very small lower bound (default is exp(-8) = 0.00034)
#'   that represents zero.
#' @param ux (`numeric(2)`) \cr
#'   Maximal value for time and dose.
#' @param iter (`integer(1)`) \cr
#'   Number of iterations in the Imperialist Competitive Algorithm. Default is 1000.
#' @param ICA.control (`list()`) \cr
#'   Additional control parameters passed to the [locally()] function.
#'   See [ICA.control()] for details.
#' @return An object of class `minimax` generated with [locally()].
#'   For `Emax_known` = `TRUE`, the design has 7 support points. Otherwise,
#'   it has 8.
#' @examples
#'  td_opt_1 <- td_opt(param = c(h = 2, delta = 0.2, gamma = 1.3, c0 = 0.2),
#'    Emax_known = TRUE,
#'     ICA.control = list(ncount = 300, rseed = 1905, trace = FALSE),
#'     iter = 600)
#'  td_opt_1
#'  plot_td_des(td_opt_1)
#'  plot_td_dcrit_equ(td_opt_1)
#'  td_opt_2 <- td_opt(param = c(h = 2, delta = 0.2, gamma = 1.3, c0 = 0.2),
#'    Emax_known = FALSE,
#'     ICA.control = list(ncount = 300, rseed = 1905, revol_rate = 0.5, trace = FALSE),
#'     iter = 600)
#'  td_opt_2
#'  plot_td_des(td_opt_2)
#'  plot_td_dcrit_equ(td_opt_2)
#'  td_opt_3 <- td_opt(param = c(h = 1, delta = 0.5, gamma = 2, c0 = 0.01),
#'    Emax_known = TRUE, ICA.control = list(ncount = 300, rseed = 1905, trace = FALSE),
#'     iter = 1000)
#'  plot_td_des(td_opt_3)
#'  plot_td_dcrit_equ(td_opt_3, plot_theta = 100, n_grid = 100, dose_lim = c(0, 0.2))


td_opt <- function(param, Emax_known = FALSE, lx = c(1, exp(-8)), ux = c(10, 1),
                   iter = 1000, ICA.control = NULL){

  stopifnot(all(names(param) == c("h", "delta", "gamma", "c0")))
  stopifnot(is.numeric(param))

  if(Emax_known == TRUE){
    fimfunc <- icaod_info
    parvars <- c("e0", "h", "delta", "gamma", "c0")
    inipars <- c(e0 = 100, param )
  } else {
    fimfunc <- icaod_info_noEmax
    parvars <- c("e0", "h", "delta", "gamma", "c0", "Emax")
    inipars <- c(e0 = 100, param, Emax = -100 )
  }

  if(is.null(ICA.control)) ICA.control <- list(ncount = 100)

  icaod_opt <- locally(predvars = c("time", "dose"),
                              parvars = parvars,
                              fimfunc= fimfunc,
                              npar = ifelse(Emax_known, 5, 6),
                              ICA.control = ICA.control,
                              inipars = inipars,
                              lx = lx,
                              ux = ux,
                              k = ifelse(Emax_known, 7, 8),
                              iter = iter)

  return(icaod_opt)
}


