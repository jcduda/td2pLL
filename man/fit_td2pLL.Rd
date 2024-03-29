% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_td2pLL.R
\name{fit_td2pLL}
\alias{fit_td2pLL}
\title{Fit a td2pLL model}
\usage{
fit_td2pLL(
  data,
  start = NULL,
  control = NULL,
  lower = NULL,
  upper = NULL,
  trace = FALSE
)
}
\arguments{
\item{data}{(numeric \code{data.frame()})\cr
Data frame with columns named \code{time}, \code{dose} and \code{resp} storing numeric
data on time, dose and response measures, respectively.}

\item{start}{(\code{list(4)})\cr
Optional list with named numeric startig values for
\code{h}, \code{delta}, \code{gamma} and \code{c0}. When no starting values
are provided, the default is used which is 2 for \code{h}, 0 for \code{c0}
and a linear interpolation procedure that leads to starting values for
\code{delta} and \code{gamma}. For details, see
\code{\link[=get_starting_values]{get_starting_values()}}.}

\item{control}{(\code{list()})\cr
Optional control argument for numerical optimization that will
be passed to the \link{nls} function for non-linear fitting.}

\item{lower}{(\code{list(4)} or \code{numeric(4)})\cr
Optional named list or named numeric vector for lower
boundaries for the parameters \code{h}, \code{delta}, \code{gamma} and
\code{c0} (in this order). As default, 1, -3*max(dose), -10 and 0 are used.}

\item{upper}{(\code{list(4)} or \code{numeric(4)})\cr
Optional named list or named numeric vector for upper
boundaries for the parameters \code{h}, \code{delta}, \code{gamma} and
\code{c0} (in this order). As default, 10, 3\*max(dose), 10 and 3\*max(dose)
are used.}

\item{trace}{(\code{logical(1)})\cr
Optinal argument passed to nls function to trace (print) the
optimization status at each iteration.}
}
\value{
An object of class \code{c("td2pLL_mod", "nls")}.
}
\description{
\code{fit_td2pLL} is used to fit time-dose two-parameter
log-logistic functions to time-dose-response data by the least-squares
approach. This application of this model is tailored to dose-response
cytotoxicity data, where also the exposure time is varied in the experiments.
The model formula is
\deqn{f(t,d)=100-100\frac{d^h}{EC_{50}(t)^h + d^h}}
with
\deqn{EC_{50}(t) = \Delta \cdot t^{-\gamma} + C_0}
where
\itemize{
\item \code{d} is the dose or concentration, \code{t} is the (exposure)
time,
\item  \code{h} is the hill- or slope parameter as known in the
classical 4pLL model (there often parametrized as -b), a.k.a. sigmoid
Emax model \code{\link[DoseFinding:drmodels]{DoseFinding::drmodels()}},
\item \code{gamma} represents the influence of (exposure) time on
the dose-response relationship. Note that the model
has identifiability issues if \code{gamma}=0. This is
the case when there is no infuence of (exposure) time \code{t}
on the dose-response relationship. Hence, for such a
situation, the model is not appropriate.
The default boundaries for gamma are -0.01 or 0.01, respectively.
If this happens, there might be no time-dependency, the td2pLL
model might not be appropriate and the parameters delta, h and c0
are not interpretable. A warning will appear.
\item \code{delta} is the maximum effect of (exposure) time on
the ED50 parameter and
\item \code{c0} is the threshold or minimal value of the
ED50 value at all (exposure) times.
}
}
\details{
The non-linear fitting minimizes the sum of squared errors.
We use the \code{nlsLM} function from the \link[minpack.lm:nlsLM]{minpack.lm::nlsLM} package.
Note that the fitting assumes the response data to be measured in percent,
i.e. ranges between 100 and 0 where it is assumed to be 100 at
dose=0 and decreases with increasing doses.
}
\examples{
data(cytotox)
data_subset <- cytotox[cytotox$compound == "ASP", c("expo", "dose", "resp")]
colnames(data_subset)[1] <- "time"
fit <- fit_td2pLL(data = data_subset)
plot(fit, add_data = data_subset)
# Note that you can also plot a td2pLL model where simply the parameters are
# specified using [plot_td2pLL].
}
