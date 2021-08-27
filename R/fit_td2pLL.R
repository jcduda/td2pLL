# function: fit_td2pLL (fit td2pLL model)
#' @export
#' @title Fit a td2pLL model
#'
#' @description `fit_td2pLL` is used to fit time-dose two-parameter
#'  log-logistic functions to time-dose-response data by the least-squares
#'  approach. This application of this model is tailored to dose-response
#'  cytotoxicity data, where also the exposure time is varied in the experiments.
#'  The model formula is
#'  \deqn{f(t,d)=100-100\frac{d^h}{EC_{50}(t)^h + d^h}}
#'  with
#'  \deqn{EC_{50}(t) = \Delta \cdot t^{-\gamma} + C_0}
#'  where
#'  \itemize{
#'    \item `d` is the dose or concentration, `t` is the (exposure)
#'  time,
#'    \item  `h` is the hill- or slope parameter as known in the
#'  classical 4pLL model (there often parametrized as -b), a.k.a. sigmoid
#'  Emax model [DoseFinding::drmodels()],
#'    \item `gamma` represents the influence of (exposure) time on
#'     the dose-response relationship. Note that the model
#'     has identifiability issues if `gamma`=0. This is
#'     the case when there is no infuence of (exposure) time `t`
#'     on the dose-response relationship. Hence, for such a
#'     situation, the model is not appropriate.
#'     The default boundaries for gamma are -0.01 or 0.01, respectively.
#'     If this happens, there might be no time-dependency, the td2pLL
#'     model might not be appropriate and the parameters delta, h and c0
#'     are not interpretable. A warning will appear.
#'    \item `delta` is the maximum effect of (exposure) time on
#'     the ED50 parameter and
#'    \item `c0` is the threshold or minimal value of the
#'     ED50 value at all (exposure) times.
#'  }
#' @param data (numeric `data.frame()`)\cr
#'  Data frame with columns named `time`, `dose` and `resp` storing numeric
#'  data on time, dose and response measures, respectively.
#' @param start (`list(4)`)\cr
#'  Optional list with named numeric startig values for
#' `h`, `delta`, `gamma` and `c0`. When no starting values
#'  are provided, the default is used which is 2 for `h`, 0 for `c0`
#'  and a linear interpolation procedure that leads to starting values for
#'  `delta` and `gamma`. For details, see
#'  [get_starting_values()].
#' @param control (`list()`)\cr
#'  Optional control argument for numerical optimization that will
#'  be passed to the [nls] function for non-linear fitting.
#' @param lower (`list(4)` or `numeric(4)`)\cr
#'  Optional named list or named numeric vector for lower
#'  boundaries for the parameters `h`, `delta`, `gamma` and
#'  `c0` (in this order). As default, 1, -3*max(dose), -10 and 0 are used.
#' @param upper (`list(4)` or `numeric(4)`)\cr
#'  Optional named list or named numeric vector for upper
#'  boundaries for the parameters `h`, `delta`, `gamma` and
#' `c0` (in this order). As default, 10, 3\*max(dose), 10 and 3\*max(dose)
#'  are used.
#' @param trace (`logical(1)`)\cr
#'  Optinal argument passed to nls function to trace (print) the
#'  optimization status at each iteration.
#' @details The non-linear fitting minimizes the sum of squared errors.
#'  We use the `nlsLM` function from the [minpack.lm::nlsLM] package.
#'  Note that the fitting assumes the response data to be measured in percent,
#'  i.e. ranges between 100 and 0 where it is assumed to be 100 at
#'  dose=0 and decreases with increasing doses.
#' @return An object of class `c("td2pLL_mod", "nls")`.
#' @examples
#' data(cytotox)
#' data_subset <- cytotox[cytotox$compound == "ASP", c("expo", "dose", "resp")]
#' colnames(data_subset)[1] <- "time"
#' fit <- fit_td2pLL(data = data_subset)
#' plot(fit, add_data = data_subset)
#' # Note that you can also plot a td2pLL model where simply the parameters are
#' # specified using [plot_td2pLL].



fit_td2pLL <- function(data, start = NULL, control = NULL, lower = NULL,
                       upper = NULL, trace = FALSE) {

  msg_1 <- 'data must be a numeric data.frame with colnames "time", "dose" and "resp"'
  if(is.null(data) | all(is.na(data)))
    stop(msg_1)
  if(!is.data.frame(data))
    stop(msg_1)
  if(ncol(data) < 3 | !all(apply(data, 2, is.numeric)))
    stop(msg_1)
  if (!(all(c("time", "dose", "resp") %in% colnames(data))))
    stop(msg_1)
  if(any(is.na(data)))
    stop("There must not be missing values in data.")

  if(min(data$time) <= 0)
    stop("Time cannot be smaller or equal to 0.")

  if(min(data$dose) < 0)
    stop("Dose cannot be smaller than 0.")

  doses <- unique(data$dose)
  if(length(doses) < 2)
    stop("Data must contain at least two different doses.")

  if(length(unique(data$time)) < 3)
    stop("Data must contain at least 3 different times.")

  if (is.null(start)) {
    start <- get_starting_values(data)
  } else {
    if ((length(start) != 4) | !is.numeric(start))
      stop("Argument start must be numeric of length 4.")
  }

  if (is.null(control)) {
    control <- list(maxiter = 100, warnOnly = TRUE, printEval = FALSE)
  }

  if (is.null(lower)) {
    lower <- c(h = 1,
               delta = -3 * max(doses),
               # gamma = -10,
               # Do let let gamma get to close to zero to avoid non-identifiability problems
               # that would lead to convergence problems
               gamma = ifelse(start$gamma < 0, -10, 0.01),
               c0 = 0)
  } else {
    if ((length(lower) != 4) | !is.numeric(lower)) stop("Argument lower must be numeric of length 4.")
  }

  if (is.null(upper)) {
    upper <- c(h = 10,
               delta = max(doses) * 3,
               # gamma = 10,
               gamma = ifelse(start$gamma < 0, -0.01, 10),
               c0 = max(doses) * 3)
  } else {
    if ((length(upper) != 4) | !is.numeric(upper)) stop("Argument upper must be numeric of length 4.")
  }

  if (is.null(trace)) {
    trace <- FALSE
  }
  # use means and weights
  data_w <-
    data %>%
    dplyr::group_by(.data$time, .data$dose) %>%
    summarize(
      n = dplyr::n(),
      resp_m = mean(.data$resp)
    ) %>%
    dplyr::ungroup()

  # fit <- nls(resp_m ~ 100 - 100 * (dose^h) / ((delta * time^(-gamma) + c0)^h + dose^h),
  #            data = data_w,
  #            weights = data_w$n,
  #            algorithm = "port",
  #            lower = lower,
  #            upper = upper,
  #            start = start,
  #            trace = trace,
  #            control = control
  # )

  fit <- minpack.lm::nlsLM(resp_m ~ 100 - 100 * (dose^h) / ((delta * time^(-gamma) + c0)^h + dose^h),
                           data = data_w,
                           weights = data_w$n,
                           algorithm = "LM",
                           lower = lower,
                           upper = upper,
                           start = start,
                           trace = trace,
                           control = list(maxiter = control$maxiter,
                                          nprint = ifelse(control$printEval, control$maxiter, 0),
                                          ptol = 1e-04,
                                          ftol = 1e-05)
  )

  if(coef(fit)["gamma"] %in% c(0.01, -0.01)){
    warning(paste0("gamma fit = ", coef(fit)["gamma"]," is on the boundary and close to zero.
                   The model might not be appropriate, there might be no time-dependency.
                   The parameters delta and c0 might not be interpretable."))
  }



  fit$orig_data <- data

  attr(fit, "class") <- c("td2pLL_mod", "nls") # to use plot.tdp2LL as method

  return(fit)
}

