
#' @title Information matrix of a td2pLL model with unknown E0
#' @description Used in [grid_opt]. E0, i.e. the effect at dose 0
#'  is assumed to be unknown as otherwise the optimal design will
#'  not propose any measurements at dose 0.
#' @param u_time (`numeric()`) \cr
#'   Contains all potential time points.
#' @param u_dose (`numeric()`) \cr
#'   Contains all potential dose points. Consider logarithmic scaling.
#' @param w (`array(length(u_time), length(u_dose))`) \cr
#'  Weight array. Each row are the weights of the doses for a given time.
#'  All weights must sum up to 1.
#' @inheritParams grad_td2pLL
#' @return (`list(2)`) \cr
#'   First element contains derivatives at all points. \cr
#'   Second element contains information matrix.

info <- function(u_time, u_dose, w, h, delta, gamma, c0) { # calculates information matrix,

  le_dose <- length(u_dose)
  le_time <- length(u_time)

  a <- array(0,c(le_time,le_dose, 5)) # "a" contains derivatives of dose response function for every design point


  # a contains derivates of td2pLL (inluding E0) for every design point
  # h, delta, gamma and c0 were passed to the funciton earlier

  # M will be the information matrix
  M <- array(0, c(5, 5))
  for(time_row in 1:le_time){
    for(dose_col in 1:le_dose){
      a[time_row, dose_col, ] <- as.vector(grad_td2pLL(time = u_time[time_row], dose = u_dose[dose_col],
                                                       h = h, delta = delta, gamma = gamma, c0=c0))

      M <- M + a[time_row, dose_col, ] %*% t(a[time_row, dose_col, ]) * w[time_row, dose_col]
    }
  }

  return(list(a, M))
}

################################################################################

#' @export
#' @title Optimized weights of td2pLL model (E0 unknown) using grid-search.
#' @description Optimized weights of a td2pLL model with assumed parameters
#'  `h`, `delta`, `gamma` and `c0`. Uses multiplicative optimization as
#'  originally described by Titterington (1976).\cr
#'  A less mathematical description is available in Holland-Letz &  Kopp-Schneider (2015).
#'  For the tuning parameter `beta_r`, we refer to Dette et al. (2008).
#' @inheritParams info
#' @param iter (`integer(1)`) \cr
#'   Number of iterations.
#' @param beta_r (`numeirc(1)`) \cr
#'  Tuning parameter in (0,1). Use 0 (default) for original Titterington.
#'  Some (Dette at al. (2008)) state that `beta_r=1` increases convergence speed.
#' @param gif_path (`character(1)`) \cr
#'  Path were a folder `plots` shall be generated and where all plots of
#'  the weights of each iterations shall be saved. No plots generated if NULL.
#' @param gif_name (`character(1)`) \cr
#'  Directory within `gif_path` where a gif from all the plots of the weights per
#'   iteration is saved. Is not generated if `iter`>1000 as this likely crashes the
#'   session. No plots or GIF generated if NULL.
#' @references
#' \itemize{
#'   \item Titterington, D. M. (1976). Algorithms for computing D-optimal design
#'    on finite design spaces. In Proc. of the 1976 Conf. on Information Science
#'     and Systems, John Hopkins University (Vol. 3, pp. 213-216).
#'   \item Holland-Letz, T., & Kopp-Schneider, A. (2015). Optimal experimental
#'   designs for dose-response studies with continuous endpoints.
#'   Archives of toxicology, 89(11), 2059-2068. \doi{10.1007/s00204-014-1335-2}
#'   \item Dette, H., Pepelyshev, A., & Zhigljavsky, A. (2008). Improving updating
#'    rules in multiplicative algorithms for computing D-optimal designs.
#'    Computational Statistics & Data Analysis, 53(2), 312-320. \doi{10.1016/j.csda.2008.10.002}
#' }


grid_opt <- function(u_dose, u_time, iter = 10,
                   h, delta, gamma, c0,
                   beta_r = 0,
                   gif_path = NULL,
                   gif_name = NULL){
  parnum <- 5

  le_dose <- length(u_dose)
  le_time <- length(u_time)
  le_tot <- le_dose * le_time

  w <- array(1/(le_time*le_dose),c(le_time, le_dose))

  # derivatives
  a <- info(u_time, u_dose, w, h, delta, gamma, c0)[[1]]

  # long_a is a, only re-arranged
  # first column has derivative of time 1, dose 1
  # second column has derivative of time 1, dose 2 ...
  long_a <- array(0, c(parnum, le_time*le_dose))
  for(time_row in 1:le_time){
    for(dose_col in 1:le_dose){
      long_a[, (time_row-1)*le_dose + dose_col] <- as.vector(a[time_row, dose_col, ])
    }
  }

  for(j in 1:iter){
    print(j)
    M <- info(u_time, u_dose, w, h, delta, gamma, c0)[[2]] # info matrix for current design
    m <- ginv(M)

    d <- (diag(t(long_a) %*% m %*% long_a) - beta_r) /
      (parnum - beta_r)

    d_mat <- matrix(d, nrow = le_time, ncol = le_dose, byrow = T)

    w <- w*d_mat # update weights
    if(!is.null(gif_path)){
      j_formatted <- formatC(j, width = 4, format = "d", flag = "0")
      dir.create(path = gif_path)
      dir.create(path = paste0(gif_path,"/plots"))
      jpeg(paste0(gif_path,"/plots/plot_",j_formatted,".jpg"))
      heatmap(w, Rowv = NA, Colv = NA, xlab = "Dose", ylab = "Time",
              labRow = u_time, labCol = round(u_dose, 3))
      dev.off()
    }
  }

  if(!is.null(gif_path)){
    if(is.null(gif_name)) gif_name <- "optDes"
    imgs <- list.files(paste0(gif_path,"/plots"), full.names = TRUE)
    img_list <- lapply(imgs, image_read)
    img_joined <- image_join(img_list)
    img_animated <- image_animate(img_joined, fps = 100)

    ## save to disk
    if(iter <=1000){
      image_write(image = img_animated,
                  path = paste0(gif_path,"/",gif_name,".gif"))
    } else {
      print("Iter > 1000. No GIF generated, as this might cause a crash.")
    }

  }

  critd <- det(M)^(1/parnum)
  return(list(w,critd,m))

}



# h<-1 # parameter settings
# delta<-0.5
# gamma<-2
# c0<-0.01
#
# iter<-1000 # number of iterations in the algorithm
# u_dose <-exp(seq(-8, 2, 0.1)) # design/dose level space (here from -5 to
# u_time <- seq(1,10, 0.1)


# grid_opt(u_dose = exp(seq(-8, 1, 0.1)),
#        u_time = seq(1, 10, 0.1),
#       iter = 500,
#       h = 1,
#       delta = 0.5,
#       gamma = 2,
#       c0 = 0.01,
#       gif_path = "./attics/optDesgif2",
#       gif_name = "optDes")

# opt2_2 <- grid_opt(u_dose = exp(seq(-8, 1, 0.1)),
#        u_time = seq(1, 10, 0.1),
#       iter = 500,
#       h = 1,
#       delta = 0.5,
#       gamma = 2,
#       c0 = 0.01,
#       beta_r = 1,
#       gif_path = "./attics/optDesgif2_beta_1",
#       gif_name = "optDes")


#
#
# opt3 <- grid_opt(u_dose = exp(seq(-8, 0, 0.1)),
#        u_time = seq(1, 10, 0.1),
#        iter = 500,
#        h = 1,
#        delta = 0.5,
#        gamma = 2,
#        c0 = 0.01,
#        gif_path = "./attics/optDesgif3",
#        gif_name = "optDes")

# opt4 <- grid_opt(u_dose = exp(seq(-8, 0, 0.1)),
#                        u_time = seq(1, 10, 0.1),
#                        iter = 500,
#                        h = 2,
#                        delta = 0.2,
#                        gamma = 1.3,
#                        c0 = 0.2,
#                        gif_path = "./attics/optDesgif4_compare",
#                        gif_name = "optDes")
# opt4[[1]]

# opt4_2 <- grid_opt(u_dose = exp(seq(-8, 0, 0.1)),
#                        u_time = seq(1, 10, 0.1),
#                        iter = 500,
#                        h = 2,
#                        delta = 0.2,
#                        gamma = 1.3,
#                        c0 = 0.2,
#                  beta_r = 1,
#                        gif_path = "./attics/optDesgif4_compare_beta_1",
#                        gif_name = "optDes")



#
#
# opt5 <- grid_opt(u_dose = exp(seq(-8, 0, 0.1)),
#                u_time = seq(1, 10, 0.1),
#                iter = 5000,
#                h = 2,
#                delta = 0.2,
#                gamma = 1.3,
#                c0 = 0.2,
#                gif_path = "./attics/optDesgif4_compare_n5000",
#                gif_name = "optDes")

