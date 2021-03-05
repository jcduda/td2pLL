
#' @export
#' @title Plot interactive td2pLL models
#'
#' @description `plot.td2pLL_mod` is the plot method for the S3 class
#' `td2pLL_mod`.
#'  Model fits generated with the [fit_td2pLL()] function are of class
#'  `c("td2pLL_mod", "nls")` and can therefore be used for this plot method.
#'  If no fitted model but a through parameters pre-specified td2pLL model
#'  shall be plotted, this can be done via the [plot_td2pLL].
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
#' @param dose_lim (`numeric(2)`)\cr
#'  Boundaries of the doses (xaxis) for plotting.
#'  Note: If `xaxis_scale = "log"` (default), then `dose_lim` cannot include 0.
#'  If `dose_lim` shall include the 0, set `xaxis_scale = "linear"`.
#' @param time_lim (`numeric(2)`)\cr
#'  Boundaries for the time (yaxis) for plotting.
#' @param add_model_data (`logical(1)`)\cr
#'  By default the original data
#'  used for the fit are added to the plot.
#' @param add_ext_data (`data.frame()`)\cr
#'  Optional numeric `data.frame` to add data points to the
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
#' plot(fit)
#' plot(fit, xaxis_scale = "linear")
#' plot(fit, title = "td2pLL model of Compound ASP", dose_lim = c(0.01, 100))
#' plot(fit, xaxis_scale = "linear", dose_lim = c(0, 15))
#' # If you want to see how the model looks like for certain parameters,
#' # use [plot_td2pLL]
#' plot_td2pLL(td2pLL_coefs = c(h = 2, delta = 3, gamma = 1.5, c0 = 1),
#'  dose_lim = c(0.01, 10), time_lim = c(1, 5))


plot.td2pLL_mod <- function(x = NULL,
                        dose_lim = NULL,
                        time_lim = NULL,
                        add_model_data = TRUE,
                        add_ext_data = NULL,
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
    if(is.null(add_ext_data) & is.null(x)) {
      stop('dose_lim has to be specified if no model is provided through x and no
      external data is provided through add_ext_data.
      If xaxis_scale="log", then dose_lim cannot include 0.
      Set xaxis_scale="linear" if you want to include 0.')
    }
    if(xaxis_scale == "log") {
      if(!is.null(x)){
        dose_lim <- c(min(x$orig_data$dose[x$orig_data$dose != 0]), max(x$orig_data$dose))
      }
      if(is.null(x) & !(is.null(add_ext_data))) {
        dose_lim <- c(min(add_ext_data$dose[add_ext_data$dose != 0]), max(add_ext_data$dose))
      }

    }
    if(xaxis_scale == "linear") {
      if(!is.null(x)){
        dose_lim <- range(x$orig_data$dose)
      }
      if(is.null(x) & !(is.null(add_ext_data))) {
        dose_lim <- range(add_ext_data$dose)
      }
    }
  }

  if(is.null(time_lim)) {
    if(is.null(x) & is.null(add_ext_data)) {
      stop('time_lim has to be specified if no model is provided through x and no
      external data is provided through add_ext_data.')
    } else {
      if(!is.null(x)){
        time_lim <- range(as.numeric(as.character(x$orig_data$time)))
      }
      if(is.null(x) & !(is.null(add_ext_data))) {
        time_lim <- range(as.numeric(as.character(add_ext_data$time)))
      }
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

  if (is.null(td2pLL_model)) {
    stop("td2pLL_model argument must be specified.")
  }

  if(!(class(td2pLL_model)[1] == "td2pLL_mod")){
    stop("td2pLL_model has to be an td2pLL_mod object generated with fit_td2pLL.")
  }

  # calculate the responses at each grid-point
  input_grid$resp <- predict(td2pLL_model, newdata = input_grid)
  td2pLL_coefs <- coef(td2pLL_model)


  # make matrix-style input for plot_ly function
  input_grid <- input_grid %>%
    pivot_wider(names_from = .data$dose, values_from = .data$resp) %>%
    dplyr::select(-.data$time) %>%
    as.matrix()


  if (is.null(title)) title <- td2pLL_plot_title(td2pLL_coefs)

  # actual plotting:
  res_plot <- td2pLL_plot_basis(dose_seq, time_seq, input_grid, xaxis_title,
                                xaxis_scale, yaxis_title, yaxis_scale, zaxis_title,
                                title)

  # add original data used for model fit, if a model fit is provided, by default
  if (add_model_data == TRUE) {
    res_plot <- td2pLL_plot_add_data(plotly_plot = res_plot,
                                     data_to_add = td2pLL_model$orig_data)
    }

  if (!is.null(add_ext_data)) {
    res_plot <- td2pLL_plot_add_data(plotly_plot = res_plot,
                                     data_to_add = add_ext_data)
  }

  # add ED50 line, if wanted
  if (add_ED50_line) {
    res_plot <- td2pLL_plot_add_ED50(res_plot, td2pLL_coefs = td2pLL_coefs,
                                     time_seq, dose_lim, ED50_line_col, ED50_line_width)

  }

  return(res_plot)
}
