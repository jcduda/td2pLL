#' @title Plot interactive td2pLL models
#'
#' @description `plot.td2pLL_mod` is the plot method for the S3 class
#' `td2pLL_mod`.
#'  Model fits generated with the [fit_td2pLL()] function are of class
#'  `c("td2pLL_mod", "nls")` and can therefore be used for this plot method.
#'  If no fitted model but a through parameters pre-specified td2pLL model
#'  shall be plotted, this can be done via the `td2pLL_coefs` argument.
#'  For details on the `td2pLL` model, see [fit_td2pLL()].
#'  If the [TDR()] function is used which performs the two-step
#'  modeling pipeline, one can apply `plot.td2pLL_mod()` to the `fit` list entry of the
#'  object returned by [TDR()],
#'  if fitting a `td2pLL` model was chosen in accordance
#'  to the anova pre-test (see [td2pLL_anova()]) calculated in [TDR()].
#' @details For further details on the td2pLL model, check [fit_td2pLL()].
#'  For details on the ANOVA used, see [td2pLL_anova()]. More over,
#'  the entire procedure is explained in duda et al. (2021).
#'  For plotting, the `plot_ly` function of package `plotly` is used.
#' @param x (`td2pLL_mod` object)\cr
#'  A `td2pLL_mod` object generatet via [fit_td2pLL]. If not
#'  provided, alternatively, `td2pLL_coefs` can be provided.
#' @param td2pLL_coefs (named `numeric(4)`)\cr
#'  If `td2pLL_model` is not provided, then `td2pLL_coefs` contains
#'  parameters for `h`, `delta`, `gamma` and `c0`
#'  of the [td2pLL] model.
#' @param dose_lim (`numeric(2)`)\cr
#'  Boundaries of the doses (xaxis) for plotting.
#'  Note: If `xaxis_scale = "log"` (default), then `dose_lim` cannot include 0.
#'  If `dose_lim` shall include the 0, set `xaxis_scale = "linear"`.
#' @param time_lim (`numeric(2)`)\cr
#'  Boundaries for the time (yaxis) for plotting.
#' @param add_data (`data.frame()`)\cr
#'  Optional `data.frame` to add data points to the
#'  surface plot. Must include columns `dose`, `time` and `resp`.
#' @param n_grid (`integer(1)`)\cr
#'  `n_grid*n_grid` is the `dose*time` grid for surface evaluations
#'  that will be interpolated. Increase `n_grid` for a smoother plot.
#' @param title (`character(1)`)\cr
#'  Optional plot title.
#' @param xaxis_scale (`character(1)` in `c("log", "linear", "-")`\cr
#'  Scale of x-axis (dose-axis).\cr
#'  If `"-"` is set, then [plot_ly()] tries to guess which scale to use.
#' @param yaxis_scale (`character(1)` in `c("log", "linear", "-")`\cr
#'  Scale of y-axis (time-axis).\cr
#'  If `"-"` is set, then [plot_ly()] tries to guess which scale to use.
#' @param xaxis_title,yaxis_title,zaxis_title (`character(1)`)\cr
#'  Title for dose-axis, time-axis and response-axis.
#' @param add_ED50_line (`logical(1)`)\cr
#'  Indicates if the line of ED50 values shall be annotated (=`TRUE`).
#' @param ED50_line_col (`character(1)`)\cr
#'  Color for optionally added ED50 line.
#' @param ED50_line_width (`numeric(1)`)\cr
#'  Width for optionally added ED50 line.
#' @param ... [any] \cr
#'   Not used.
#' @examples
#' data(cytotox)
#' data_subset <- cytotox[cytotox$compound == "ASP", c("expo", "dose", "resp")]
#' colnames(data_subset)[1] <- "time"
#' fit <- fit_td2pLL(data = data_subset)
#' plot(fit, add_data = data_subset)
#' plot.td2pLL_mod(x=NULL, td2pLL_coefs = c(h = 2, delta = 5, gamma = 2, c0=1),
#' dose_lim = c(0, 10), time_lim = c(1, 2), xaxis_scale = "linear", n_grid = 200)
#' @export
plot.td2pLL_mod <- function(x = NULL, td2pLL_coefs = NULL,
                        dose_lim = NULL,
                        time_lim = NULL,
                        add_data = NULL,
                        n_grid = 100,
                        title = NULL,
                        xaxis_scale = "log",
                        yaxis_scale = "-",

                        xaxis_title = "Dose",
                        yaxis_title = "Time",
                        zaxis_title = "Response",

                        add_ED50_line = TRUE,
                        ED50_line_col = "red",
                        ED50_line_width = 6,
                        ...
                        ) {

  td2pLL_model = x

  if(is.null(dose_lim)) {
    if(is.null(add_data)) {
      stop('dose_lim has to be specified. If xaxis_scale="log", then
      dose_lim cannot include 0. Set xaxis_scale="linear" if you want to include 0.')
    }
    if(xaxis_scale == "log") {
      dose_lim <- c(min(add_data$dose[add_data$dose != 0]), max(add_data$dose))
    }
    if(xaxis_scale == "linear") {
      dose_lim <- range(add_data$dose)
    }
  }

  if(is.null(time_lim)) {
    if(is.null(add_data)) {
      stop('time_lim has to be specified if no add_data is provided.')
    } else {
      time_lim <- range(as.numeric(as.character(add_data$time)))
    }

  }

  time_seq <- seq(time_lim[1], time_lim[2], length.out = n_grid)
  if(xaxis_scale == "log") {
    dose_seq <- c(0, exp(seq(log(dose_lim[1]), log(dose_lim[2]), length.out = n_grid)))
  } else {
    dose_seq <- seq(dose_lim[1], dose_lim[2], length.out = n_grid)
  }

  # seq(dose_lim[1], dose_lim[2], length.out = n_grid)

  input_grid <- expand.grid(time = time_seq, dose = dose_seq) %>%
    as.data.frame()

  if (is.null(td2pLL_model) & is.null(td2pLL_coefs)) {
    stop("Either the td2pLL_model argument or the td2pLL_coefs argument must be specified.")
  }

  if (!(is.null(td2pLL_model)) & !(is.null(td2pLL_coefs))) {
    stop("Only td2pLL_model or td2pLL_coefs can be specified, not both.")
  }

  # calculate the responses at each grid-point
  if (!is.null(td2pLL_coefs)) {
    input_grid$resp <- apply(input_grid, 1, function(x) {
      td2pLL(
        dose = x["dose"],
        time = x["time"],
        h = td2pLL_coefs["h"],
        delta = td2pLL_coefs["delta"],
        gamma = td2pLL_coefs["gamma"],
        c0 = td2pLL_coefs["c0"]
      )
    })
  } else {
    input_grid$resp <- predict(td2pLL_model, newdata = input_grid)
    coefs <- coef(td2pLL_model)
  }

  # make matrix-style input for plot_ly function
  input_grid <- input_grid %>%
    pivot_wider(names_from = .data$dose, values_from = .data$resp) %>%
    dplyr::select(-.data$time) %>%
    as.matrix()

  if (is.null(title)) {
    title <- paste0(
      "h = ", round(td2pLL_coefs["h"], 3),
      "; delta = ", round(td2pLL_coefs["delta"], 3),
      ";\n gamma = ", round(td2pLL_coefs["gamma"], 3),
      "; c0 = ", round(td2pLL_coefs["c0"], 3)
    )
  }

  # actual plotting:

  res_plot <- plotly::plot_ly(
    x = dose_seq, y = time_seq, z = input_grid,
    type = "surface"
  ) %>%
    plotly::layout(
      title = title,
      scene = list(
        xaxis = list(
          title = xaxis_title,
          tick0 = 0,
          type = xaxis_scale
        ),
        yaxis = list(
          title = yaxis_title,
          type = yaxis_scale
        ),
        zaxis = list(title = zaxis_title)
      )
    )

  # add single data points, if available
  if (!is.null(add_data)) {
    res_plot <- res_plot %>% plotly::add_markers(
      x = add_data$dose,
      y = add_data$time,
      z = add_data$resp,
      marker = list(size = 2),
      showlegend = F
    )
  }

  # add ED50 line, if wanted
  if (add_ED50_line) {
    add_ED50 <- data.frame(
      time = time_seq,
      ED50 = td2pLL_coefs["delta"] * time_seq^(-td2pLL_coefs["gamma"]) +
        td2pLL_coefs["c0"]
    ) %>% dplyr::filter(
      ED50 > dose_lim[1] & ED50 < dose_lim[2]
      )

    # add_ED50 <- add_ED50 %>% filter(ED50 <= 1.02)
    add_ED50$resp <- 50

    res_plot <- res_plot %>% plotly::add_trace(
      x = add_ED50$ED50, y = add_ED50$time, z = add_ED50$resp,
      type = "scatter3d", mode = "lines",
      showlegend = F,
      line = list(
        color = ED50_line_col,
        width = ED50_line_width
      )
    )
  }

  return(res_plot)
}
