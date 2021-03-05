
#' @export
#' @title Plot interactive td2pLL model without underlying model fit
#'
#' @description `plot_td2pLL` plots an interactive time-dose 2 parameter
#'  log-logistic (td2pLL) model when parameters are specified but there is no
#'  underlying fit. See [fit_td2pLL] for the model function.
#'  If a model fit is available through [fit_td2pLL], please use [plot.td2pLL_mod].
#' @param td2pLL_coefs (named `numeric(4)`)\cr
#'  Contains parameters for `h`, `delta`, `gamma` and `c0`
#'  of the [td2pLL] model.
#' @inheritParams plot.td2pLL_mod
#' @examples
#' plot_td2pLL(td2pLL_coefs = c(h = 2, delta = 3, gamma = 1.5, c0 = 1),
#'  dose_lim = c(0.01, 10), time_lim = c(1, 5))

plot_td2pLL <- function(td2pLL_coefs = NULL,
                        dose_lim = NULL,
                        time_lim = NULL,
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

  if(is.null(dose_lim)) {
    if(is.null(add_ext_data)) {
      stop('dose_lim has to be specified if  no
      external data is provided through add_ext_data.
      If xaxis_scale="log", then dose_lim cannot include 0.
      Set xaxis_scale="linear" if you want to include 0.')
    }

    if(xaxis_scale == "log") {
      dose_lim <- c(min(add_ext_data$dose[add_ext_data$dose != 0]), max(add_ext_data$dose))
    }

    if(xaxis_scale == "linear") {
      dose_lim <- range(add_ext_data$dose)
    }
  }


  if(is.null(time_lim)) {
    if(is.null(add_ext_data)) {
      stop('time_lim has to be specified if no
      external data is provided through add_ext_data.')
    }

    if(!(is.null(add_ext_data))) {
      time_lim <- range(as.numeric(as.character(add_ext_data$time)))
    }
  }


  time_seq <- seq(time_lim[1], time_lim[2], length.out = n_grid)
  if(xaxis_scale == "log") {
    dose_seq <- c(0, exp(seq(log(dose_lim[1]), log(dose_lim[2]), length.out = n_grid)))
  } else {
    dose_seq <- seq(dose_lim[1], dose_lim[2], length.out = n_grid)
  }

  stopifnot(!is.null(td2pLL_coefs))
  stopifnot(all(is.numeric(td2pLL_coefs)))
  stopifnot(length(td2pLL_coefs) == 4)
  stopifnot(length(unique(names(td2pLL_coefs))) == 4)
  stopifnot(length(setdiff(names(td2pLL_coefs), c("h", "delta","gamma","c0"))) == 0)

  input_grid <- expand.grid(time = time_seq, dose = dose_seq) %>%
    as.data.frame()

  input_grid$resp <- apply(input_grid, 1, function(x) {
    td2pLL(
      dose = x["dose"],
      time = x["time"],
      h = td2pLL_coefs["h"],
      delta = td2pLL_coefs["delta"],
      gamma = td2pLL_coefs["gamma"],
      c0 = td2pLL_coefs["c0"]
    )

  }
  )

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

  # if provided, add external data
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

