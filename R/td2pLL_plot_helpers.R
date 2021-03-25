
#' @title helper function for plotting td2pLL model
#' @description Generates title of plot based on coefficients
#' @inheritParams plot_td2pLL
#' @return `character(1)`

td2pLL_plot_title <- function(td2pLL_coefs) {
  title <- paste0("h = ", round(td2pLL_coefs["h"], 3) ,"; delta = ",
                  round(td2pLL_coefs["delta"], 3),";\n gamma = ",
                  round(td2pLL_coefs["gamma"], 3),"; c0 = ",
                  round(td2pLL_coefs["c0"], 3)
  )

  return(title)
}

#' @title helper function for plotting td2pLL model
#' @description Generates basis of td2pLL plot in plotly
#' @inheritParams plot.td2pLL_mod
#' @param dose_seq (`numeric()`)
#' @param time_seq (`numeric()`)
#' @param input_grid (`data.frame()`)\cr
#'  Contains the function values. Number of rows is length of `dose_seq`,
#'  number of cols is length of `time_seq`. The (i,j)th element is the
#'  response value at dose `dose_seq`(i) and `time_seq`(j).

td2pLL_plot_basis <- function(dose_seq, time_seq, input_grid, xaxis_title, xaxis_scale,
                              yaxis_title, yaxis_scale, zaxis_title, title) {
  res <- plotly::plot_ly(
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
  return(res)
}


#' @title helper function for plotting td2pLL model
#' @description Adds data points to plot
#' @param plotly_plot \cr
#'  A plotly plot generated with [td2pLL_plot_basis]
#' @param data_to_add (`data.frame`)\cr
#'  Numeric data frame with columns `dose`, `time` and `resp`.

td2pLL_plot_add_data <- function(plotly_plot, data_to_add) {
  plotly_plot <- plotly_plot %>% plotly::add_markers(
    x = data_to_add$dose,
    y = data_to_add$time,
    z = data_to_add$resp,
    marker = list(size = 2),
    showlegend = F
  )

  return(plotly_plot)
}

#' @title helper function for plotting td2pLL models
#' @description Adds ED50 line to td2pLL model plot
#' @inheritParams td2pLL_plot_add_data
#' @inheritParams plot.td2pLL_mod
#' @inheritParams plot_td2pLL
#' @inheritParams td2pLL_plot_basis

td2pLL_plot_add_ED50 <- function(plotly_plot, td2pLL_coefs, time_seq, dose_lim,
                                 ED50_line_col, ED50_line_width) {

  plotly_plot <- tryCatch(
    {
  add_ED50 <- data.frame(
    time = time_seq,
    ED50 = td2pLL_coefs["delta"] * time_seq^(-td2pLL_coefs["gamma"]) +
      td2pLL_coefs["c0"]
  ) %>% dplyr::filter(
    .data$ED50 > dose_lim[1] & .data$ED50 < dose_lim[2]
  )
  # add_ED50 <- add_ED50 %>% filter(ED50 <= 1.02)
  add_ED50$resp <- 50

  plotly_plot <- plotly_plot %>% plotly::add_trace(
    x = add_ED50$ED50, y = add_ED50$time, z = add_ED50$resp,
    type = "scatter3d", mode = "lines",
    showlegend = F,
    line = list(
      color = ED50_line_col,
      width = ED50_line_width
    )
  )
  return(plotly_plot)
  },
  error = function(cond){
    warning("ED50 line can not be drawn as these values are not reached yet.
            Try setting a larger range in dose_lim")
    return(plotly_plot)
  }
  )

}





