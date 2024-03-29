#' @title Fitting a 2pLL model
#'
#' @description `fit_joint_2pLL` fits a 2pLL model. The "joint" is in the
#' name because, despite being a regular 2pLL model, this function is used
#' in the anova-td2pLL pipeline [TDR()]: If the anova pre-test, that
#' checks if there is a difference in ED50 parameters between (exposure) times,
#' is not significant, the exposure times are ignored. In other words,
#' a single, 'joint', dose-response curve is fitted to the data,
#' ignoring the information of exposure time.
#'
#' @param data (`data.frame()`)\cr
#'  Contains the `dose` and `resp` variable..
#' @details The function is a wrapper function using the LL2.2() argument
#' from the [drc::drm()] function with the fixed asymptotes
#' upper=100 and lower=0.
#' As a first try, the Nelder-Mead method is used as optimization procedure.
#' If this yields an error, the BFGS method is used as opitmization
#' method.
#' @note When using the LL2.2() model from the [drc::drm()] function,
#' the ED50 parameter is parametrized as log(ED50). We do this for improved
#' stability.
#'
#' @return An object of class `drc`.

fit_joint_2pLL <- function(data) {

  msg_1 <- 'data must be a numeric data.frame with colnames "dose" and "resp"'
  if(is.null(data) | all(is.na(data)))
    stop(msg_1)
  if(!is.data.frame(data))
    stop(msg_1)
  if(ncol(data) < 2)
    stop(msg_1)
  if (!(all(c("dose", "resp") %in% colnames(data))))
    stop(msg_1)
 if(!(is.numeric(data$dose) & is.numeric(data$resp)))
    stop(msg_1)

  doses <- unique(data$dose)
  if(length(doses) < 2)
    stop("Data must contain at least two different doses.")

  tryCatch(
    {
      drc::drm(resp ~ dose,
               data = data, fct = drc::LL2.2(upper = 100),
               control = drc::drmc(method = "Nelder-Mead")
      )
    },
    error = function(cond) {
      return(drc::drm(resp ~ dose,
                      data = data, fct = drc::LL2.2(upper = 100)
      ))
    }
  )
}
