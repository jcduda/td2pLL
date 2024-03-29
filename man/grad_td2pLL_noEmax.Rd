% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grad_td2pLL.R
\name{grad_td2pLL_noEmax}
\alias{grad_td2pLL_noEmax}
\title{Gradient function of the td2pLL model when Emax=-100 is NOT known}
\usage{
grad_td2pLL_noEmax(time, dose, h, delta, gamma, c0, Emax)
}
\arguments{
\item{time}{(\code{numeric(1)})}

\item{dose}{(\code{numeric(1)})}

\item{h}{(\code{numeric(1)})}

\item{delta}{(\code{numeric(1)})}

\item{gamma}{(\code{numeric(1)})}

\item{c0}{(\code{numeric(1)})}

\item{Emax}{(\code{numeric(1)})}
}
\value{
\code{numeric(6)} Gradient vector of the td2pLL model at given
parameters. \cr
Note that the first entry is a 1 for the intercept. In the the \link{td2pLL}
model, it is assumed to be 100. However, for initial normalization,
it is important to consider this parameter as unknown for optimal
design calculations.
The last entry is not the derivate by \code{Emax} which is assumed fixed at -100
in the \link{td2pLL} model.
It is know assumed unknown so that in the optimal design, measurements
for larger doses are also recommended.
}
\description{
For details on the td2pLL model, see \link{fit_td2pLL}
}
