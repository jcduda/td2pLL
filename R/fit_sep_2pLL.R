
#' @title Fitting separate 2pLL models for (exposure) times
#'
#' @description `fit_sep_2pLL` is used in the anova pre-test of the
#'  time-dose-response analysis pipeline in [TDR()]. For each
#'  (exposure) time, a 2pLL dose-response model is fitted with
#'  upper limit 100 and lower limit 0. For each exposure time, an independent
#'  ED50 parameter is modeled. However, a common `h` parameter
#'  (or -b in the often used parametrization of the 4pLL model) is
#'  shared across the dose-response models of the (exposure) times.
#' @param data (`data.frame()`)\cr
#'  Data frame containing `dose` (numeric), `resp` (numeric) and `time` (factor) variables.
#' @details The model serves as the full model, denoted by Q_1 in duda et al.
#'  (2021). `fit_sep_2pLL` is a wrapper function that uses the
#'  [drc::drm()] function. At first, the optimization method
#'  BFGS is used to fit the model. If this fails, the Nelder-Mead method is used.
#' @return An object of class `drc`.
#' @note The LL.2 option from [drc::drm()] is used, as here,
#'  opposed to [fit_joint_2pLL()], it seems to lead to more stable
#'  fits.
#' @importFrom rlang .data


fit_sep_2pLL <- function(data) {

  time <- NULL

  msg_1 <- 'data must be a data.frame with colnames "time" (factor),
  "dose" (numeric) and "resp" (numeric)'
  if(is.null(data) | all(is.na(data)))
    stop(msg_1)
  if(!is.data.frame(data))
    stop(msg_1)
  if(ncol(data) < 3)
    stop(msg_1)
  if (!(all(c("time", "dose", "resp") %in% colnames(data))))
    stop(msg_1)

  if(!is.factor(data$time)){
    message("For the ANOVA pre-test, data$time was changed to a factor variable.")
    data$time <- as.factor(data$time)
  }

  if(!is.numeric(data$dose) | !is.numeric(data$dose) | !is.factor(data$time))
    stop(msg_1)


  doses <- unique(data$dose)
  if(length(doses) < 2)
    stop("Data must contain at least two different doses.")

  if(length(unique(data$time)) < 2)
    stop("Data must contain at least two different times.")



  tryCatch(
    {
      drc::drm(resp ~ dose,
               curveid = time,
               control = drc::drmc(method = "Nelder-Mead"),
               data = data,
               fct = drc::LL2.2(upper = 100),
               pmodels = list(
                 ~1, # h
                 ~time
               )
      ) #
    },
    error = function(cond) {
      return(
        drc::drm(resp ~ dose,
                 curveid = time,
                 data = data,
                 fct = drc::LL2.2(upper = 100),
                 pmodels = list(
                   ~1, # h
                   ~time
                 )
        )
      )
    }
  )
}
