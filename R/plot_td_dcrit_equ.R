
#' @export
#' @title Plot sensitivity function of calculated optimal design to verify optimality
#' @description The values at each point (time, dose) is the gradient at the point
#' times the inverse of the FIM of the calculated optimal design (at the assumes
#' parameters) times the gradient.
#' If the resulting surface is below zero, the design is optimal and in the support
#' points the value is zero. \cr
#' Due to computational inaccuracies, the maximal value is usually not exactly zero.
#' @inheritParams plot_td_des
#' @param time_lim (`numeric(2)`) \cr
#'  Lower and upper bound of the time component.
#' @param dose_lim (`numeric(2)`) \cr
#'  Lower and upper bound of the dose component.
#' @param n_grid (`numeric(1)`) \cr
#'  The values will be calculates on a `n_grid` time `n_grid` grid.
#' @param title (`character(1)`) \cr
#'  Optional plot title. By default, the maximal value and its position are
#'  displayed in the title.
#' @param return_values (`logical(1)`) \cr
#'  Defines if the values of the points are returned.
#' @param plot_theta,plot_phi,plot_r
#'  Angles theta, phi and r passed to [persp]

plot_td_dcrit_equ <- function(td_opt_des,
                              time_lim = NULL,
                              dose_lim = NULL, n_grid = 61, title = NULL,
                              return_values = FALSE,
                              plot_theta = 120, plot_phi = 30, plot_r = sqrt(3)){

  des <- get_td_opt_des(td_opt_des)
  param <- td_opt_des$arg$inipars

  if(is.null(time_lim)){
    time_lim <- c(td_opt_des$arg$lx[1], td_opt_des$arg$ux[1])
  }

  if(is.null(dose_lim)){
    dose_lim <- c(td_opt_des$arg$lx[2], td_opt_des$arg$ux[2])
  }

  #if(log10_dose){
  #  if(dose_lim[1] <=0) stop("Lower dose limit cannot be <=0 if log10_dose = TRUE")
  #  dose_values <- 10^(seq(log10(dose_lim[1]), log10(dose_lim[2]), length.out = n_grid))
  #} else {
  dose_values <- seq(dose_lim[1], dose_lim[2], length = n_grid)
  #}
  time_values <- seq(time_lim[1], time_lim[2], length = n_grid)

  df_values <- expand.grid(time = time_values, dose = dose_values)

  stopifnot(length(param) %in% c(5, 6))
  if(length(param) == 5){
    tmp_info <- icaod_info
    tmp_grad <- grad_td2pLL
  } else {
    tmp_info <- icaod_info_noEmax
    tmp_grad <- grad_td2pLL_noEmax
  }

  finfo_inv <- solve(tmp_info(x = c(des$time, des$dose), w = des$weight, param = param))

  dcrit_equ <- function(x, finfo_inv, param){
    if(length(param) == 5){
      grad_val <- unlist(tmp_grad(time = x[1], dose = x[2],
                                  h = param["h"], delta = param["delta"],
                                  gamma = param["gamma"],
                                  c0 = param["c0"]))
    } else {
      grad_val <- unlist(tmp_grad(time = x[1], dose = x[2],
                                  h = param["h"], delta = param["delta"],
                                  gamma = param["gamma"],
                                  c0 = param["c0"],
                                  Emax = param["Emax"]))
    }

    return(grad_val %*% finfo_inv %*% t(grad_val) - nrow(finfo_inv))
  }


  df_values$eq <- apply(df_values, 1, function(x){dcrit_equ(x = x,
                                                            finfo_inv = finfo_inv,
                                                            param = param)
  })

  eq_values_matrix <- matrix(df_values$eq, nrow = n_grid, ncol = n_grid)

  if(is.null(title)) title <- paste0("Max: ", round(max(df_values$eq), 3)
                                     ," at t=", round(df_values$time[which.max(df_values$eq)], 3)
                                     ,", dose=", round(df_values$dose[which.max(df_values$eq)], 3)
  )

  graphics::persp(time_values, dose_values, eq_values_matrix,
                  main = title,
                  theta = plot_theta, phi = plot_phi, r = plot_r,
                  ticktype = "detailed", xlab="time", ylab="dose",
                  zlab= "Equ. Thm. Value",
                  box = T,
                  col = "white",
                  border = "black",
                  shade = NA)

  if(return_values) return(list(values = df_values))

}
