% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_td2pLL.R
\name{plot_td2pLL}
\alias{plot_td2pLL}
\title{Plot interactive td2pLL model without underlying model fit}
\usage{
plot_td2pLL(
  td2pLL_coefs = NULL,
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
)
}
\arguments{
\item{td2pLL_coefs}{(named \code{numeric(4)})\cr
Contains parameters for \code{h}, \code{delta}, \code{gamma} and \code{c0}
of the \link{td2pLL} model.}

\item{dose_lim}{(\code{numeric(2)})\cr
Boundaries of the doses (xaxis) for plotting.
Note: If \code{xaxis_scale = "log"} (default), then \code{dose_lim} cannot include 0.
If \code{dose_lim} shall include the 0, set \code{xaxis_scale = "linear"}.}

\item{time_lim}{(\code{numeric(2)})\cr
Boundaries for the time (yaxis) for plotting.}

\item{add_ext_data}{(\code{data.frame()})\cr
Optional numeric \code{data.frame} to add data points to the
surface plot. Must include columns \code{dose}, \code{time} and \code{resp}.}

\item{n_grid}{(\code{integer(1)})\cr
\code{n_grid*n_grid} is the \code{dose*time} grid for surface evaluations
that will be interpolated. Increase \code{n_grid} for a smoother plot.}

\item{title}{(\code{character(1)})\cr
Optional plot title.}

\item{xaxis_scale}{(\code{character(1)} in \code{c("log", "linear", "-")}\cr
Scale of x-axis (dose-axis).\cr
If \code{"-"} is set, then \code{\link[=plot_ly]{plot_ly()}} tries to guess which scale to use.}

\item{yaxis_scale}{(\code{character(1)} in \code{c("log", "linear", "-")}\cr
Scale of y-axis (time-axis).\cr
If \code{"-"} is set, then \code{\link[=plot_ly]{plot_ly()}} tries to guess which scale to use.}

\item{xaxis_title}{(\code{character(1)})\cr
Title for dose-axis, time-axis and response-axis.}

\item{yaxis_title}{(\code{character(1)})\cr
Title for dose-axis, time-axis and response-axis.}

\item{zaxis_title}{(\code{character(1)})\cr
Title for dose-axis, time-axis and response-axis.}

\item{add_ED50_line}{(\code{logical(1)})\cr
Indicates if the line of ED50 values shall be annotated (=\code{TRUE}).}

\item{ED50_line_col}{(\code{character(1)})\cr
Color for optionally added ED50 line.}

\item{ED50_line_width}{(\code{numeric(1)})\cr
Width for optionally added ED50 line.}

\item{...}{\link{any} \cr
Not used.}
}
\description{
\code{plot_td2pLL} plots an interactive time-dose 2 parameter
log-logistic (td2pLL) model when parameters are specified but there is no
underlying fit. See \link{fit_td2pLL} for the model function.
If a model fit is available through \link{fit_td2pLL}, please use \link{plot.td2pLL_mod}.
}
\examples{
plot_td2pLL(td2pLL_coefs = c(h = 2, delta = 3, gamma = 1.5, c0 = 1),
 dose_lim = c(0.01, 10), time_lim = c(1, 5))
plot_td2pLL(td2pLL_coefs = c(h = 1.5, delta = 5, gamma = 1, c0 = 0.1),
 dose_lim = c(0.1, 10), time_lim = c(1, 10), xaxis_scale = "log",
 yaxis_scale = "linear")
 # If time and response is on log-scale and gamme = 1(=original Haber's law
 # with equal contribution of time and exposure), then one can clearly see
 # haber's law since the red ED50 is a line, if viewed from above.
plot_td2pLL(td2pLL_coefs = c(h = 1.5, delta = 5, gamma = 1, c0 = 0.1),
 dose_lim = c(0.1, 10), time_lim = c(1, 10), xaxis_scale = "log",
 yaxis_scale = "log")
}
