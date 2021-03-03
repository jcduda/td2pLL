# function: get_starting_values
#' @title Starting values for td2pLL model
#'
#' @description Calculates (cheap) default starting values for fitting a
#'  td2pLL model based on dose time response data.
#' @param data (numeric `data.frame()`)\cr
#'  Data frame with numeric columns named `time`, `dose` and `resp`.
#' @param h_start (`numeric(1)`)\cr
#'  Optional starting value for the `h` parameter.
#'  Default is 2.
#' @param c0_start (`numeric(1)`)\cr
#'  Optional starting value for the trheshold parameter `c0`.
#'  Default is 0.
#' @details As starting value for `delta` and `gamma`, a pair of
#'  cheap starting values of (dose, ED50) at the lowest and the highest
#'  (exposure) time are calculated via the [interp_ED50()] function.
#'  With these two pairs (dose_1, ED50_1) and (dose_2, ED50_2), as
#'  well as the set starting values of `h` and `c0`,
#'  the model equation of the td2pLL model
#'  \deqn{f(d,t)=100-100\frac{d^h}{EC_{50}(t)^h + d^h},}
#'  \deqn{EC_{50}(t) = \Delta \cdot t^{-\gamma} + C_0}
#'  is solved to get the starting values `gamma_start` and
#'  `delta_start` for `gamma` and  `delta`.
#' @return (`list(4)`)\cr
#'  List with starting values for h, delta, gamma and c0.

get_starting_values <- function(data, h_start = 2, c0_start = 0) {
  data_low_time_mean <- data %>%
    dplyr::filter(.data$time == min(.data$time)) %>%
    dplyr::group_by(.data$dose) %>%
    dplyr::summarize(mean_resp = mean(.data$resp)) %>%
    dplyr::ungroup()

  data_high_time_mean <- data %>%
    dplyr::filter(.data$time == max(.data$time)) %>%
    dplyr::group_by(.data$dose) %>%
    dplyr::summarize(mean_resp = mean(.data$resp), .groups = "drop") %>%
    dplyr::ungroup()

  ED50_start_low_time <- interp_ED50(data = data_low_time_mean)
  ED50_start_high_time <- interp_ED50(data = data_high_time_mean)


  time_low_high <- range(data$time)

  # get starting values for delta and gamma:
  gamma_start <- log(ED50_start_low_time / ED50_start_high_time,
                     base = time_low_high[1] * time_low_high[2]
  )
  delta_start <- mean(c(
    ED50_start_low_time / time_low_high[1]^(gamma_start),
    ED50_start_high_time / time_low_high[2]^(gamma_start)
  ))
  c0_start <- c0_start

  return(list(
    h = h_start,
    delta = delta_start,
    gamma = gamma_start,
    c0 = c0_start
  ))
}
