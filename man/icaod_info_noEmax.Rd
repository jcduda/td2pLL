% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/icaod_info.R
\name{icaod_info_noEmax}
\alias{icaod_info_noEmax}
\title{Information matrix for \code{\link[=locally]{locally()}} assuming unknown e0 and unknown eMax}
\usage{
icaod_info_noEmax(x, w, param)
}
\arguments{
\item{x}{(\code{numeric(2*8)}) \cr
Contains the current support points under consideration where first all
time coordinates are stored in x and then all dose coordinates.}

\item{w}{(\code{numeric(8)}) \cr
Weights of the support points. \code{sum(w)} must be 1.}

\item{param}{(\code{numeric(6)}) \cr
Contains assumed true parameters for \code{e0}, \code{h}, \code{delta}, \code{gamma}, \code{c0}, \code{Emax} where
the values of \code{e0} and \code{Emax} have no effect on the optimal design calculation
because they are linear parameter.}
}
\value{
Information matrix of dimension 6 times 6.
}
\description{
If e0 is assumed unknown and eMax is assumed unknown in the td2pLL model
(where usually we assume e0=100 and eMax=100 for dytotoxicity assays)
the optimal design will propose a support point in dose=0 and in the maximal dose.
To assume that eMax is known, use \link{icaod_info}.
}
