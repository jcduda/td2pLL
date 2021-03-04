
#' @export
#' @title ANOVA pre-test in the TDR analysis pipeline
#'
#' @description In the time-dose-response (TDR) analysis pipeline,
#'  data where in addition to dose (or concentration) and response, also
#'  time (e.g. exposure times) are measured, it is checked in a (nested) anova
#'  based pre-test if the td2pLL model is appropriate.
#'  That means, it is checked if the exposure time has an effect on the
#'  dose-response relationship or not.
#' @param data (`data.frame()`)\cr
#'  Numeric data frame containing the response `resp`, `time` and `dose`.
#' @param alpha (`numeric(1)`)\cr
#'  The significance level.
#' @details A nested anova is performed where a null model is tested against
#'  a full model. The null model is a regular 2pLL model with upper and
#'  lower limit 100 and 0, respectively.
#'  The (exposure) time variable is ignored and a single dose-response
#'  curve is fitted for the null model.
#'  In the full model, an individual ED50 value is calculated for each (exposure)
#'  time level. However, only one common `h` parameter is included in
#'  the full model.
#'  Hence, if the test rejects the null model, there seems to be significant
#'  evidence that (exposure) time has an effect on the dose-response
#'  relationship. Thus, the td2pLL model will be used during the fitting step
#'  of the TDR analysis pipeline. Otherwise, if there is no significant evidence
#'  of the influence of (exposure) time on the dose-response relationship,
#'  a simple 2pLL model is used in the fitting step that ignores (exposure)
#'  time.
#'  For more details see Duda et al. (2021).
#' @examples
#' data(cytotox)
#' data_subset <- cytotox[cytotox$compound == "ASP", c("expo", "dose", "resp")]
#' colnames(data_subset)[1] <- "time"
#' td2pLL_anova(data = data_subset)


td2pLL_anova <- function(data, alpha = 0.05) {
  stopifnot(is.numeric(alpha) & alpha > 0 & alpha < 1)
  stopifnot(is.na(setdiff(colnames(data), c("time", "dose", "resp"))))

  if(!(is.factor(data$time))) data$time <- as.factor(data$time)

  # dose-response curve with separate ED50 parameters, but shared h
  drm_seperate <- fit_sep_2pLL(data = data)

  # single dose-response curves ignoring exposure time
  drm_joint <- fit_joint_2pLL(data = data)

  anova_res <- anova(drm_joint, drm_seperate, details = F)
  p_value <- anova_res$`p value`[2]
  signif <- ifelse(p_value < 0.05, TRUE, FALSE)

  res <- list(
    signif = signif,
    alpha = alpha,
    anova = anova_res,
    conv = TRUE
  )

  return(res)
}
