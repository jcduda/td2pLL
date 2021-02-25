
# function: td2pLL
#' @export
#' @title Calculate response value for given td2pLL model
#' @description \code{td2pLL} returns the response value of a fully specified
#' time-dose two-parameter log-logistic model at a certain time and dose
#' value:
#' \deqn{f(d,t)=100-100\frac{d^h}{EC_{50}(t)^h + d^h}}
#'  with
#'  \deqn{EC_{50}(t) = \Delta \cdot t^{-\gamma} + C_0}
#' @param time numeric for the time value where the response shall be calculated
#' @param dose numeric for the dose (or concentration) value where the response
#' shall be calculated
#' @param h The \code{h} parameter of the model
#' @param gamma The \code{gamma} parameter of the model
#' @param c0 The \code{c0} parameter of the model
#' @param delta The \code{delta} parameter of the model
#' @return The response value of the model [in percent] at the given \code{time} and
#' \code{dose} value.
#' @examples
#' td2pLL(time=4, dose=0.1, h=2, gamma=2.5, c0=0.1, delta=0.3)
td2pLL <- function(time, dose, h, gamma, c0, delta){
  EC50 <- delta * time^(-gamma) + c0
  resp <- 100 - 100 * (dose^h) / (EC50^h + dose^h)
  return(resp = resp)
}



################
# FITTING
################


# function: interp_ED50
#' @title Starting value calculation for ED50 parameter
#' @description Calculates cheap starting parameter for the ED50 parameter of a
#'  2pLL model by interpolation. This function is used in get_starting_values.
#' @param  data A data.frame with columns dose (numeric) and mean_resp (numeric)
#'  of mean response data at given dose level.
#' @details The function assumes that the mean response values lay within 0 and
#'  100, i.e. have percent units.
#'  It also assumes a downward trend of the ED50 with dose.
#'  If the lowest mean_resp value is already above 50, the function
#'  returns the maximal dose.
#' @return Returns the numeric \code{ED50_interp}.
interp_ED50 <- function(data=NULL){
  if(min(data$mean_resp) > 50) {
    ED50_interp <- max(data$dose)
  } else {
    # dose-ids just above 50 and then below 50
    id_bef_after <- which.min(data$mean_resp > 50) - c(1, 0)
    dose_bef_after <- data$dose[id_bef_after]
    resp_bef_after <- data$mean_resp[id_bef_after]

    # intercept
    b <- (resp_bef_after[2] - resp_bef_after[1] *
            dose_bef_after[2]/dose_bef_after[1]) /
      (1 - dose_bef_after[2] / dose_bef_after[1])
    # slope
    m <- (resp_bef_after[2] - b) / dose_bef_after[2]

    ED50_interp <- (50-b)/m
  }
  return(ED50_interp)
}



# function: get_starting_values
#' @title Starting values for td2pLL model
#'
#' @description Calculates (cheap) default starting values for fitting a
#'  td2pLL model based on dose time response data.
#' @param data data frame with numeric columns named time, dose and resp
#' @param h_start optional starting value for the \code{h} parameter.
#' Default is 2.
#' @param c0_start optional starting value for the theshold parameter \code{c0}.
#' Default is 0.
#' @details As starting value for \code{delta} and \code{gamma}, a pair of
#' cheap starting values of (dose, ED50) at the lowest and the highest
#' (exposure) time are calculated via the \code{\link{interp_ED50}} function.
#' With these two pairs (dose_1, ED50_1) and (dose_2, ED50_2), as
#' well as the set starting values of \code{h} and \code{c0},
#' the model equation of the td2pLL model
#' \deqn{f(d,t)=100-100\frac{d^h}{EC_{50}(t)^h + d^h},}
#' \deqn{EC_{50}(t) = \Delta \cdot t^{-\gamma} + C_0}
#'  is solved to get the starting values \code{gamma_start} and
#'  \code{delta_start} for \code{gamma} and  \code{delta}.
#' @return List with starting values for h, delta, gamma and c0.

get_starting_values <- function(data, h_start = 2, c0_start = 0){
  data_low_time_mean <- data %>% filter(time == min(time)) %>%
    group_by(dose) %>%
    summarize(mean_resp = mean(resp), .groups = "drop")

  data_high_time_mean <- data %>% filter(time == max(time)) %>%
    group_by(dose) %>%
    summarize(mean_resp = mean(resp), .groups = "drop") %>%
    ungroup()

  ED50_start_low_time <- interp_ED50(data = data_low_time_mean)
  ED50_start_high_time <- interp_ED50(data = data_high_time_mean)


  time_low_high <- range(data$time)

  # get starting values for delta and gamma:
  gamma_start <- log(ED50_start_low_time / ED50_start_high_time,
                     base = time_low_high[1]*time_low_high[2])
  delta_start <- mean(c(ED50_start_low_time / time_low_high[1]^(gamma_start),
                        ED50_start_high_time / time_low_high[2]^(gamma_start)))
  c0_start <- c0_start

  return(list(h = h_start,
              delta = delta_start,
              gamma = gamma_start,
              c0 = c0_start)
  )
}



# function: fit_td2pLL (fit td2pLL model)
#' @export
#' @title Fit a td2pLL model
#'
#' @description \code{fit_td2pLL} is used to fit time-dose two-parameter
#'  log-logistic functions to time-dose-response data by the least-squares
#'  approach. This application of this model is tailored to dose-repsonse
#'  cytotoxicity data, where also the exposure time is varied in the experiments.
#'  The model formula is
#'  \deqn{f(d,t)=100-100\frac{d^h}{EC_{50}(t)^h + d^h}}
#'  with
#'  \deqn{EC_{50}(t) = \Delta \cdot t^{-\gamma} + C_0}
#'  where
#'  \itemize{
#'    \item \code{d} is the dose or concetration, \code{t} is the (exposure)
#'  time,
#'    \item  \code{h} is the hill- or slope parameter as known in the
#'  classical 4pLL model (there often parametrized as -b), a.k.a. sigmoid
#'  Emax model \code{\link[DoseFinding]{drmodels}},
#'    \item \code{gamma} represents the influence of (exposure) time on
#'     the dose-response relationship. Note that the model
#'     has identifiability issues if \code{gamma}=0. This is
#'     teh case when there is no infuence of (exposure) time \code{t}
#'     on the dose-response relationship. Hence, for such a
#'     situation, the model is not appropriate.
#'    \item \code{delta} is the maximum effect of (exposure) time on
#'     the EC50 parameter and
#'    \item \code{c0} is the threshold or minimal value of the
#'     EC50 value at all (exposure) times.
#'  }
#' @param data numeric data frame with columns named time, dose and resp.
#' @param start Optional listwith named numeric startig values for
#' \code{h}, \code{delta}, \code{gamma} and \code{c0}. When no starting values
#' are provided, the default is used which is 2 for \code{h}, 0 for \code{c0}
#' and a linear interpolation procedure that leads to starting values for
#' \code{delta} and \code{gamma}. For details, see
#' \code{\link{get_starting_values}}.
#' @param control Optional control argument for numerical optimization that will
#' be passed to the \code{nls} function that is used here for non-linear
#' fitting.
#' @param lower Optional named list or named numeric vector for lower
#' boundaries for the parameters \code{h}, \code{delta}, \code{gamma} and
#' \code{c0} (in this order). As default, 1, -3*max(dose), -10 and 0 are used.
#' @param upper Optional named list or named numeric vector for upper
#' boundaries for the parameters \code{h}, \code{delta}, \code{gamma} and
#' \code{c0} (in this order). As default, 10, 3*max(dose), 10 and 3*max(dose)
#'  are used.
#' @param trace Ooptinal argument passed to nls function to trace (print) the
#' optimization status at each iteration.
#' @details The non-linear fitting minimizes the sum of squared errors.
#' We use the \code{nls} function with the port algorithm.
#' Note that the fitting assumes the response data to be measured in percent,
#' i.e. ranges between 100 and 0 where it is assumed to be 100 at
#' dose=0 and decreases with increasing doses.
#' @return An object of class \code{c("td2pLL", "nls")}.

fit_td2pLL <- function(data, start = NULL, control = NULL, lower = NULL,
                       upper = NULL, trace = FALSE){
  if(!(all(colnames(data) %in% c("time", "dose", "resp")))){
    stop("Argument data must be a data frame with columns expo, dose, resp.")
  }


  doses <- unique(data$dose)


  if(is.null(start)){
    start <- get_starting_values(data)
  } else {
    if((length(start) != 4) | !is.numeric(start)) stop("Argument start must be numeric of length 4.")
  }

  if(is.null(control)){
    control <- list(maxiter = 100, warnOnly = TRUE, printEval = FALSE)
  }

  if(is.null(lower)){
    lower <- c(h = 1, delta = -3*max(doses), gamma = -10, c0 = 0)
  } else {
    if((length(lower) != 4) |  !is.numeric(lower)) stop("Argument lower must be numeric of length 4.")
  }

  if(is.null(upper)){
    upper <- c(h = 10, delta = max(doses)*3, gamma = 10, c0 = max(doses)*3)
  } else {
    if((length(upper) != 4) |  !is.numeric(upper)) stop("Argument upper must be numeric of length 4.")
  }

  if(is.null(trace)){
    trace <- FALSE
  }
  # use means and weights
  data_w <-
    data %>% group_by(time, dose) %>%
    summarize(n = n(),
              resp_m = mean(resp)
    ) %>%
    ungroup()

  fit <- nls(resp_m ~ 100 - 100*(dose^h) / ( ( delta * time^(-gamma) + c0)^h + dose^h),
             data = data_w,
             weights = data_w$n,
             algorithm = "port",
             lower = lower,
             upper = upper,
             start = start,
             trace = trace,
             control = control
  )

  attr(fit, "class") <- c("td2pLL", "nls") # to use plot.tdp2LL as method

  return(fit)

}



#' @title Fitting a 2pLL model
#'
#' @description \code{fit_joint_2pLL} fits a 2pLL model. The "joint" is in the
#' name because, despite being a regular 2pLL model, this function is used
#' in the anova-td2pLL pipeline \code{\link{TDR}}: If the anova pre-test, that
#' checks if there is a difference in EC50 parameters between (exposure) times,
#' is not significant, the exposure times are ignored. In other words,
#' a single, 'joint', dose-response curve is fitted to the data,
#' ignoring the information of exposure time.
#'
#' @param data data frame containing the dose and resp variable
#' @details The function is a wrapper function using the LL2.2() argument
#' from the \code{\link[drc]{drm}} function with the fixed asymptotes
#' upper=100 and lower=0.
#' As a first try, the BFGS method is used as optimization procedure.
#' If this yields an error, the Nelder-Mead method is used as opitmization
#' method.
#' @note When using the LL2.2() model from the \code{\link[drc]{drm}} function,
#' the EC50 parameter is parametrized as log(EC50). We do this for improved
#' stability.
#'
#' @return An object of class \code{drc}.

fit_joint_2pLL <- function(data){
  tryCatch(
    {
      drm(resp ~ dose,
          data = data, fct = LL2.2(upper = 100))
    },
    error = function(cond){
      return(drm(resp ~ dose,
                 data = data, fct = LL2.2(upper = 100),
                 control = drmc(method = "Nelder-Mead")))
    }
  )
}

#' @title Fitting seperate 2pLL models for (exposure) times
#'
#' @description \code{fit_sep_2pLL} is used in the anova pre-test of the
#'  time-dose-response analysis pipeline in \code{\link{TDR}}. For each
#'  (exposure) time, a 2pLL dose-response model is fitted with
#'  upper limit 100 and lower limit 0. For each exposure time, an independent
#'  EC50 parameter is modeled. However, a common \code{h} parameter
#'  (or -b in the often used parametrization of the 4pLL model) is
#'  shared across the dose-response models of the (exposure) times.
#' @param data Data frame containing the dose, resp and time variable.
#' @details The model serves as the full model, denoted by Q_1 in duda et al.
#'  (2021). \code{fit_sep_2pLL} is a wrapper function that uses the
#' \code{\link[drc]{drm}} function. At first, the optimization method
#' BFGS is used to fit the model. If this fails, the Nelder-Mead method is used.
#' @return An object of class \code{drc}.
#' @note The LL.2 option from \code{\link[drc]{drm}} is used, as here,
#' opposed to \code{\link{fit_joint_2pLL}}, it seems to lead to more stable
#' fits.


fit_sep_2pLL <- function(data){
  tryCatch(
    {
      drm(resp ~ dose,
          curveid = time,
          data = data %>% dplyr::mutate(time = factor(time)),
          fct = LL.2(upper = 100),
          pmodels = list(~1, # h
                         ~time)) # EC50
    },
    error = function(cond){
      return(drm(resp ~ dose,
                 curveid = expo,
                 control = drmc(method = "Nelder-Mead"),
                 data = data %>% dplyr::mutate(time = factor(time)),
                 fct = LL.2(upper = 100),
                 pmodels = list(~1, # h
                                ~time)) # EC50
      )
    }
  )
}




