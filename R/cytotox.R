#' Exposure time-concentration-response cytotoxicity data
#'
#' A dataset containing cytotoxicity viability assay data of 30 compounds where
#' the exposure time `expo` and concentration `dose` were varied in the experiment.
#' For each compound, there are between 214 and 276 measurements.\cr
#' For each compound, there are measurements of 3 donors.\cr
#' For each compound and donor, there are 3 different exposure times (1, 2, and 7).\cr
#' For each compound, donor and exposure time, there are
#' mostly 6 (240 cases), rarely 7 (27 cases) or 8 (3 cases) different concentrations
#' in milli molar (mM).\cr
#' For each compound, donor, exposure time and concentration, there are mostly
#' 3 measurements. For the control (concentration=0) there
#' are often 8 measurements.
#'
#' @format A data frame with 7064 rows and 12 variables:
#' \describe{
#'   \item{compound}{compound name}
#'   \item{donor}{donor within compound, factor}
#'   \item{expo}{exposure time in days}
#'   \item{dose}{concentration in mM}
#'   \item{resp}{normalized response: resp = (raw_resp / mean_resp) * 100}
#'   \item{sample_id}{sample within compound, donor and exposure time}
#'   \item{donor_name}{donor name}
#'   \item{control_mean}{Arithmetic mean of raw control (dose=0) measurements
#'     `raw_resp` for a compound, donor and exposure time. Sometimes,
#'     there are two sets of conrol measurements for such a setting.
#'     See `by_control`.}
#'  \item{raw_resp}{Raw measurements in mM}
#'  \item{by_control}{If "0", there was just one set of control measurements
#'     for a compound, donor and exposure time combination, which
#'     leas to `control_mean` and was always used for the normalization in
#'     `resp`.\cr
#'     If "1", there are 2 sets of control measurements and the first
#'     one (control_1 in original excel files) is used for this observation for
#'     normalization.\cr
#'     If "2", there are 2 sets of control measurements and the second
#'     one (control_2 in original excel files) is used for this observation for
#'     normalization.}
#'  \item{left_asymp}{Left asymptote calculated in when as additional
#'     pre-processing step (refitting) is applied.}
#'  \item{resp_refit}{New normalizes responses when as aditional
#'     pre-processing, refitting is applied. That means that for a given compound,
#'     donor and exposure time, a dose-response curve (4pLL) is
#'     fitted. The corresponding measurements are divided by the
#'     resulting left (upper) asymptote `left_asymp` and multiplied by 100.
#'     This way, the data is expected to better follow the assumption of
#'     a left asymptote of 100 percent.}
#'   ...
#' }
#' @source \url{https://doi.org/10.1007/s00204-018-2302-0}
"cytotox"
