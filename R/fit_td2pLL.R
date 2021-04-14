# function: fit_td2pLL (fit td2pLL model)
#' @export
#' @title Fit a td2pLL model
#'
#' @description `fit_td2pLL` is used to fit time-dose two-parameter
#'  log-logistic functions to time-dose-response data by the least-squares
#'  approach. This application of this model is tailored to dose-repsonse
#'  cytotoxicity data, where also the exposure time is varied in the experiments.
#'  The model formula is
#'  \deqn{f(t,d)=100-100\frac{d^h}{EC_{50}(t)^h + d^h}}
#'  with
#'  \deqn{EC_{50}(t) = \Delta \cdot t^{-\gamma} + C_0}
#'  where
#'  \itemize{
#'    \item `d` is the dose or concetration, `t` is the (exposure)
#'  time,
#'    \item  `h` is the hill- or slope parameter as known in the
#'  classical 4pLL model (there often parametrized as -b), a.k.a. sigmoid
#'  Emax model [DoseFinding::drmodels()],
#'    \item `gamma` represents the influence of (exposure) time on
#'     the dose-response relationship. Note that the model
#'     has identifiability issues if `gamma`=0. This is
#'     teh case when there is no infuence of (exposure) time `t`
#'     on the dose-response relationship. Hence, for such a
#'     situation, the model is not appropriate.
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
#'  be passed to the [nls] function that for non-linear fitting.
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
#'  We use the `nls` function with the port algorithm.
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


fit_td2pLL <- function(data, start = NULL, control = NULL, lower = NULL,
                       upper = NULL, trace = FALSE) {
  if (!(all(colnames(data) %in% c("time", "dose", "resp")))) {
    stop("Argument data must be a data frame with columns expo, dose, resp.")
  }


  doses <- unique(data$dose)


  if (is.null(start)) {
    start <- get_starting_values(data)
  } else {
    if ((length(start) != 4) | !is.numeric(start)) stop("Argument start must be numeric of length 4.")
  }

  if (is.null(control)) {
    control <- list(maxiter = 100, warnOnly = TRUE, printEval = FALSE)
  }

  if (is.null(lower)) {
    lower <- c(h = 1, delta = -3 * max(doses), gamma = -10, c0 = 0)
  } else {
    if ((length(lower) != 4) | !is.numeric(lower)) stop("Argument lower must be numeric of length 4.")
  }

  if (is.null(upper)) {
    upper <- c(h = 10, delta = max(doses) * 3, gamma = 10, c0 = max(doses) * 3)
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

  fit <- nls(resp_m ~ 100 - 100 * (dose^h) / ((delta * time^(-gamma) + c0)^h + dose^h),
             data = data_w,
             weights = data_w$n,
             algorithm = "port",
             lower = lower,
             upper = upper,
             start = start,
             trace = trace,
             control = control
  )

  fit$orig_data <- data

  attr(fit, "class") <- c("td2pLL_mod", "nls") # to use plot.tdp2LL as method

  return(fit)
}

