
#' @export
#' @title Time-Dose-Response analysis pipeline
#'
#' @description `TDR` performs the time-dose-response analysis pipeline as
#'  presented in Duda et al. (2021). That is: For a dose-response or
#'  concentration-response data set where, additionally also the (exposure)
#'  time is varied, this procedure can be applied. The main aim of this
#'  procedure are cytotoxicity data.
#'  The approach is divided into two steps. In step one, it is tested in a
#'  nested ANOVA if the (exposure) time has an influence on the dose-response
#'  relationship. In case of a significant results, the time-dose two-parameter
#'  log-logistic model is fitted to the data:
#'  \deqn{f(d,t)=100-100\frac{d^h}{EC_{50}(t)^h + d^h}}
#'  with
#'  \deqn{EC_{50}(t) = \Delta \cdot t^{-\gamma} + C_0.}
#'  If no significant result is obtained, a dose-response 2pLL curve is fitted,
#'  ignoring the information on (exposure) time.
#' @details For further details on the td2pLL model, check [fit_td2pLL()].
#'  For details on the ANOVA used, see [td2pLL_anova()]. More over,
#'  the entire procedure is explained in duda et al. (2021).
#' @param data (Numeric `data.frame()`)\cr
#'  Data frame with columns named `time`, `dose` and `resp`.\cr
#'  Note that the data is expected to be on the percent scale with values
#'  (roughly) within 0 and 100.
#' @param alpha (`numeric(1)`in (0,1))\cr
#'  1- alpha is the confidence level for testing in step 1.
#' @param strict_stop (`logical(1)`)\cr
#'  Optional logical. When `FALSE`, the default, then
#'  in case of an error due to non-convergence in the pre-test, then in the
#'  second step a simple 2pLL model is fitted as if the pre-test was non-
#'  significant.
#'  If `strict_stop` is `TRUE` and there is an error due to
#'  non-convergence in the pre-test, the procedure stops and no model is fitted
#'  in step 2.
#' @param ... Further arguments that can be passed on to [fit_td2pLL()].
#' @return A list with entries `pretest` and `fit`.
#'  `pretest` is captures the anova based pre-test result as a list with
#'  entires `signif` (TRUE/FALSE or NA if no-convergence), `alpha`,
#'  `anova_res` (the anova result from function \cite{anova}) and
#'  `conv` (logical: If the pre-test converged).
#'  `fit` Is, depending on the pre-test, either an object of class
#'  `td2pLL` or a 2pLL fit, i.e. an object of class `drc`.
#' @examples
#' data(cytotox)
#' data_subset <- cytotox[cytotox$compound == "ASP", c("expo", "dose", "resp")]
#' colnames(data_subset)[1] <- "time"
#' TDR_res <- TDR(data = data_subset)
#' # Pre-test rejected time dependency, so a regular 2pLL model is the result
#' plot(TDR_res$fit)
#' data_subset <- cytotox[cytotox$compound == "CHL", c("expo", "dose", "resp")]
#' colnames(data_subset)[1] <- "time"
#' TDR_res <- TDR(data = data_subset)
#' # Pre-test did not reject time dependency, so a  td2pLL model is the result
#' # Note that the interactive Plot is in the Viewer panel, not in the Plots panel
#' plot(TDR_res$fit)

TDR <- function(data, alpha = 0.05, strict_stop = FALSE, ...) {
  stopifnot(is.data.frame(data))
  if(any(is.na(data)))
    stop("There must not be missing values in data.")
  stopifnot(all(colnames(data) %in% c("time", "dose", "resp")))
  stopifnot(is.numeric(data$time) &
    is.numeric(data$dose) &
    is.numeric(data$resp))

  if(length(unique(data$dose)) < 2)
    stop("There must be at least two different dose levels.")

  if(length(unique(data$time)) < 3)
    stop("There must be at least three different time levels.")

  if(min(data$time) <= 0)
    stop("Time cannot be smaller or equal to 0.")

  if(min(data$dose) < 0)
    stop("Dose cannot be smaller than 0.")

  stopifnot(is.numeric(alpha))
  stopifnot(length(alpha) == 1)
  stopifnot(alpha > 0 & alpha < 1)

  stopifnot(length(strict_stop) == 1)
  stopifnot(is.logical(strict_stop))
  # ANOVA pre-test
  res_pretest <- tryCatch(
    {
      td2pLL_anova(data = data, alpha = alpha)
    },
    error = function(cond) {
      list(
        signif = NA,
        alpha = alpha,
        anova_res = NA,
        conv = FALSE
      )
    }
  )

  if (res_pretest$conv == FALSE & strict_stop == TRUE) {
    warning("The ANOVA pre-test did not converge. Since strict_stop=TRUE
      was set, no model is fitted in the second step.")
    return(list(pretest = res_pretest, fit = NA))
  }

  if(res_pretest$conv == FALSE & strict_stop == FALSE)
    warning("The ANOVA pre-test did not converge. Since strict_stop=FALSE
            was set,this will be handled as if time-dependency was NOT
            rejected and the time componenet will be ignored in the fitting step.")


  if ((res_pretest$conv == FALSE & strict_stop == FALSE) |
    res_pretest$signif == FALSE) {
    res_fit <- tryCatch(
      {
        fit_joint_2pLL(data = data)
      },
      error = function(cond) NA
    )
  }


  if (res_pretest$signif == TRUE) {
    res_fit <- tryCatch(
      {
        fit_td2pLL(data = data, ...)
      },
      error = function(cond) NA
    )
  }

  return(list(pretest = res_pretest, fit = res_fit))
}
