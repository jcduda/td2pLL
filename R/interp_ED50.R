
# function: interp_ED50
#' @title Starting value calculation for ED50 parameter
#' @description Calculates cheap starting parameter for the ED50 parameter of a
#'  2pLL model by interpolation. This function is used in [get_starting_values].
#' @param  data (numeric `data.frame()`)\cr
#'  A data frame with columns `dose` for the dose values and `mean_resp`
#'  as (mean) response for a given dose value.
#'  of mean response data at given dose level.
#' @details The function assumes that the mean response values lay within 0 and
#'  100, i.e. have percent units.
#'  It also assumes a downward trend of the ED50 with dose.
#'  If the lowest mean_resp value is already above 50, the function
#'  returns the maximal dose.
#' @return Returns the numeric `ED50_interp`.
interp_ED50 <- function(data = NULL) {
  if (min(data$mean_resp) > 50) {
    ED50_interp <- max(data$dose)
  } else {
    # dose-ids just above 50 and then below 50
    id_bef_after <- which.min(data$mean_resp > 50) - c(1, 0)
    dose_bef_after <- data$dose[id_bef_after]
    resp_bef_after <- data$mean_resp[id_bef_after]

    # # intercept
    # b <- (resp_bef_after[2] - resp_bef_after[1] *
    #         dose_bef_after[2] / dose_bef_after[1]) /
    #   (1 - dose_bef_after[2] / dose_bef_after[1])
    # # slope
    # m <- (resp_bef_after[2] - b) / dose_bef_after[2]

    # slope
    m <- (resp_bef_after[2] - resp_bef_after[1]) /
      (dose_bef_after[2] - dose_bef_after[1])

    # intercept
    b <- resp_bef_after[1] - m * dose_bef_after[1]


    ED50_interp <- (50 - b) / m
  }
  return(ED50_interp)
}
