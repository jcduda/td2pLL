
#' @title Fitting seperate 2pLL models for (exposure) times
#'
#' @description `fit_sep_2pLL` is used in the anova pre-test of the
#'  time-dose-response analysis pipeline in [TDR()]. For each
#'  (exposure) time, a 2pLL dose-response model is fitted with
#'  upper limit 100 and lower limit 0. For each exposure time, an independent
#'  EC50 parameter is modeled. However, a common `h` parameter
#'  (or -b in the often used parametrization of the 4pLL model) is
#'  shared across the dose-response models of the (exposure) times.
#' @param data Data frame containing the dose, resp and time variable.
#' @details The model serves as the full model, denoted by Q_1 in duda et al.
#'  (2021). `fit_sep_2pLL` is a wrapper function that uses the
#' [drc::drm()] function. At first, the optimization method
#' BFGS is used to fit the model. If this fails, the Nelder-Mead method is used.
#' @return An object of class `drc`.
#' @note The LL.2 option from [drc::drm()] is used, as here,
#' opposed to [fit_joint_2pLL()], it seems to lead to more stable
#' fits.
#' @importFrom rlang .data


fit_sep_2pLL <- function(data) {

  data$time <- as.factor(data$time)
  tryCatch(
    {
      drc::drm(resp ~ dose,
        curveid = .data$time,
        data = data,
        fct = drc::LL.2(upper = 100),
        pmodels = list(
          ~1, # h
          ~time
        )
      ) #
    },
    error = function(cond) {
      return(
        drc::drm(resp ~ dose,
          curveid = .data$time,
          control = drc::drmc(method = "Nelder-Mead"),
          data = data,
          fct = drc::LL.2(upper = 100),
          pmodels = list(
            ~1, # h
            ~time
          )
        )
      )
    }
  )
}
