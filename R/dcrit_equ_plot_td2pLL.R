

#' @title D-optimality equivalence theorem plot for td2pLL model
#' @description Using the [dcrit_equ_td2pLL] function to calculate the
#' equivalence theorem inequality at one design points, the corresponding
#' equivalence theorem plot can be visualized with this function.
#' Note that, if the entire surface is equal to or below zero, the design
#' is D-optimal.
#' @note The plot object itself cannot be saved as it is generated with [persp].
#' @inheritParams dcrit_equ_td2pLL
#' @inheritParams opt_des_td2pLL
#' @param time_lim (`numeric(2)`) \cr
#' Lower and upper bound for time
#' @param dose_lim (`numeric(2)`) \cr
#' Lower and upper bound of dose
#' @param n_grid (`integer(1)`) \cr
#' `n_grid` times `n_grid` spans the grid for evaluating the equivalence theorem
#' within the time and dose boundaries
#' @param return_values (`logical(1)`) \cr
#'  If the equivalence theorem values at each grid point shall be returned,
#'  default is `TRUE`.
#' @param title (`character(1)`) \cr
#'  Plot title passed to `persp`
#' @param plot_theta (`numeric(1)`) \cr
#'  Theta angle passed to `persp`. Default is 45.
#' @param plot_phi (`numeric(1)`) \cr
#'  Phi angle passed to `persp`. Default is 15.
#' @param plot_r (`numeric(1)`) \cr
#'  Rotation passed to persp`. Default is sqrt(3).
#' @examples
#'  theta_0 = c(h=2, delta=0.2, gamma=1.3, c0=0.2)
#'  set.seed(1905)
#'  des <- opt_des_td2pLL(theta =  theta_0)
#'  equ_values <- dcrit_equ_plot_td2pLL(des = des, theta = theta_0)$values
#'  # From other perspectives:
#'  dcrit_equ_plot_td2pLL(des = des, theta = theta_0, return_values = FALSE,
#'   plot_theta = 100, plot_phi = 30)
#'  dcrit_equ_plot_td2pLL(des = des, theta = theta_0, return_values = FALSE,
#'   plot_theta = 250, plot_phi = 60)
#'  dcrit_equ_plot_td2pLL(des = des, theta = theta_0, return_values = FALSE,
#'   plot_theta = 200, plot_phi = 60)
#'  dcrit_equ_plot_td2pLL(des = des, theta = theta_0, return_values = FALSE,
#'   plot_theta = 150, plot_phi = 60)
#'   # check largest value:
#'   equ_values[which.max(equ_values$eq), ]
#' @export

dcrit_equ_plot_td2pLL <- function(des, theta, time_lim = c(1, 10), dose_lim = c(0, 1),
                                  n_grid = 75, return_values = T,
                                  title = NULL,
                                  plot_theta = 45, plot_phi = 15, plot_r = sqrt(3)){
  time_values <- seq(time_lim[1], time_lim[2], length = n_grid)
  dose_values <- seq(dose_lim[1], dose_lim[2], length = n_grid)
  df_values <- expand.grid(time = time_values, dose = dose_values)


  finfo_inv <- solve(finfo(w = des$w, M = des$supPoints, gradient = grad_td2pLL_pso,
                           theta = theta))


  df_values$eq <- apply(df_values, 1, function(x){dcrit_equ_td2pLL(M_x = x,
                                                                   des = des,
                                                                   theta = theta,
                                                                   finfo_inv = finfo_inv)
  })

  eq_values_matrix <- matrix(df_values$eq, nrow = n_grid, ncol = n_grid)

  if(is.null(title)) title <- paste0("Max: ", round(max(df_values$eq), 4)
                                    ," at t=", round(df_values$time[which.max(df_values$eq)], 4)
                                    ,", dose=", round(df_values$dose[which.max(df_values$eq)], 4)
                                    )

  graphics::persp(time_values, dose_values, eq_values_matrix,
        main = title,
                    theta = plot_theta, phi = plot_phi, r = plot_r,
                    ticktype = "detailed", xlab="time", ylab="dose",
                    zlab= "Equ. Thm. Value")

  if(return_values) return(list(values = df_values))

}
