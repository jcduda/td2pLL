% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_td_opt_des.R
\name{get_td_opt_des}
\alias{get_td_opt_des}
\title{Filtering the design from a minimax object}
\usage{
get_td_opt_des(td_opt_des)
}
\arguments{
\item{td_opt_des}{An object of class \code{minimax} from the \link{ICAOD} package}
}
\value{
A data.frame with columns \code{time} and \code{dose} for the support points
and \code{weight} for the corresponding weight
}
\description{
From a minimax object generated with \code{\link[=locally]{locally()}}
or \code{\link[=td_opt]{td_opt()}}, extract only the design, i.e. the support
points and weights
}
