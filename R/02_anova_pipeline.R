

#' @title ANOVA pre-test in the TDR analysis pipeline
#'
#' @description In the time-dose-response (TDR) analysis pipeline,
#' data where in addition to dose (or concentration) and response, also
#' time (e.g. exposure times) are measured, it is checked in a (nested) anova
#' based pre-test if the td2pLL model is appropriate.
#' That means, it is checked if the exposure time has an effect on the
#' dose-response relationship or not.
#' @param data data frame containing the response \code{resp}, \code{time} and
#' \code{dose}
#' @param alpha Optional. alpha is 1 - confidence level.
#' @details A nested anova is performed where a null model is tested against
#' a full model. The null model is a regular 2pLL model with upper and
#' lower limit 100 and 0, respectively.
#' The (exposure) time variable is ignored and a single dose-response
#' curve is fitted for the null model.
#' In the full model, an individual EC50 value is calculated for each (exposure)
#' time level. However, only one common \code{h} parameter is included in
#' the full model.
#' Hence, if the test rejects the null model, there seems to be significant
#' evidence that (exposure) time has an effect on the dose-response
#' relationship. Thus, the td2pLL model will be used during the fitting step
#' of the TDR analysis pipeline. Otherwise, if there is no significant evidence
#' of the influence of (exposure) time on the dose-response relationship,
#' a simple 2pLL model is used in the fitting step that ignores (exposure)
#' time.
#' For more details see Duda et al. (2021).


td2pLL_anova <- function(data, alpha = 0.05){

  stopifnot(is.numeric(alpha) & alpha > 0 & alpha < 1)
  stopifnot(is.na(setdiff(colnames(data), c("time", "dose", "resp"))))

  # dose-response curve with seperate EC50 parameters, but shared h
  drm_seperate <- fit_sep_2pLL(data = data)

  # single dose-response curves ignoring exposure time
  drm_joint <- fit_joint_2pLL(data = data)

  anova_res <- anova(drm_joint, drm_seperate, details = F)
  p_value <- anova_res$`p value`[2]
  signif <- ifelse(p_value < 0.05, TRUE, FALSE)

  res <- list(signif = signif,
              alpha = alpha,
              anova = anova_res,
              conv = TRUE)

  return(res)
}

#' @export
#' @title Time-Dose-Response analysis pipeline
#'
#' @description \code{TDR} performs the time-dose-response analysis pipeline as
#' presented in Duda et al. (2021). That is: For a dose-response or
#' concentration-response data set where, additionally also the (exposure)
#' time is varied, this procedure can be applied. The main aim of this
#' procedure are cytotoxicity data.
#' The approach is divided into two steps. In step one, it is tested in a
#' nested ANOVA if the (exposure) time has an influence on the dose-response
#' relationship. In case of a significant results, the time-dose two-parameter
#' log-logistic model is fitted to the data:
#' \deqn{f(d,t)=100-100\frac{d^h}{EC_{50}(t)^h + d^h}}
#'  with
#' \deqn{EC_{50}(t) = \Delta \cdot t^{-\gamma} + C_0.}
#' If no significant result is obtained, a dose-response 2pLL curve is fitted,
#' inoring the information on (exposure) time.
#' @details For further details on the td2pLL model, check \code{\link{fit_td2pLL}}.
#' For details on the ANOVA used, see \code{\link{td2pLL_anova}}. More over,
#' the entire procedure is explained in duda et al. (2021).
#' @param data numeric data frame with columns named time, dose and resp. Note
#' that the data is expected to be on the percent scale and there have values
#' (roughly) within 0 and 100.
#' @param alpha 1- alpha is the confidence level for testing in step 1.
#' @param strict_stop Optional logical. When \code{FALSE}, the default, then
#' in case of an error due to non-convergence in the pre-test, then in the
#' second step a simple 2pLL model is fitted as if the pre-test was non-
#' significant.
#' If \code{strict_stop} is \code{TRUE} and there is an error due to
#' non-convergence in the pre-test, the procedure stops and no model is fitted
#' in step 2.
#' @param ... Further aguments that can be passed on to \code{\link{fit_td2pLL}}.
#' @return A list with entries \code{pretest} and \code{fit}.
#' \code{pretest} is captures the anova based pre-test result as a list with
#' entires \code{signif} (TRUE/FALSE or NA if no-convergence), \code{alpha},
#' \code{anova_res} (the anova result from function \cite{anova}) and
#' \code{conv} (logical: If the pre-test converged).
#' \code{fit} Is, depending on the pre-test, either an object of class
#' \code{td2pLL} or a 2pLL fit, i.e. an object of class \code{drc}.

TDR <- function(data, alpha=0.05, strict_stop = FALSE, ...){

  stopifnot(is.data.frame(data))
  stopifnot(all(colnames(data) %in% c("time", "dose", "resp")))
  stopifnot(is.numeric(data$time) &
              is.numeric(data$dose) &
              is.numeric(data$resp))

  stopifnot(is.numeric(alpha))
  stopifnot(length(alpha) == 1)
  stopifnot(alpha > 0 & alpha < 1)

  stopifnot(length(strict_stop) == 1)
  stopifnot(is.logical(strict_stop))
  # ANOVA pre-test
  res_pretest <- tryCatch(
    {td2pLL_anova(data=data, alpha=alpha)},

    error = function(cond){
      list(signif=NA,
           alpha=alpha,
           anova_res=NA,
           conv=FALSE)
    }
  )

  if(res_pretest$conv == FALSE & strict_stop == TRUE){
    warning("The ANOVA pre-test did not converge. Since strict_stop=TRUE
      was set, no model is fitted in the second step.")
    return(list(pretest = res_pretest, fit = NA))
  }

  if((res_pretest$conv == FALSE & strict_stop == FALSE) |
      res_pretest$signif == FALSE){
    res_fit <- tryCatch(
      {fit_joint_2pLL(data=data)},

      error = function(cond) NA
    )
  }


  if(res_pretest$signif == TRUE){
    res_fit <- tryCatch(
      {fit_td2pLL(data=data, ...)},

      error = function(cond) NA
    )
  }

  return(list(pretest = res_pretest, fit = res_fit))

}


