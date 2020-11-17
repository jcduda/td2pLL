########################################
# functions for the ext3pLL simulation #
########################################

library(dplyr)
library(DoseFinding)
library(tidyr)
library(plotly)
library(scales)
library(drc)


# function: ext4pLL
#           calculates value of extended 4pLL model for parameters dose, time, h, gamma, c0 and delta
ext3pLL <- function(dose, expo, h, gamma, c0, delta){
  EC50 <- delta * expo^(-gamma) + c0
  resp <- 100 - 100 * (dose^h) / (EC50^h + dose^h)
  return(resp = resp)
}



################
# FITTING
################


# function: interp_ED50
#' @title Starting value calculation for ED50 parameter
#' @description Calculates cheap starting parameter for the ED50 parameter of a
#'  2pLL model by interpolation. This function is used in get_starting_values.
#' @param  data A data.frame with columns dose (numeric) and mean_resp (numeric)
#'  of mean response data at given dose level.
#' @details The function assumes that the mean response values lay within 0 and
#'  100, i.e. have percent units.
#'  It also assumes a downward trend of the ED50 with dose.
#'  If the lowest mean_resp value is already above 50, the function
#'  returns the maximal dose.
#' @return Returns the numeric \code{ED50_interp}.
interp_ED50 <- function(data=NULL){
  if(min(data$mean_resp) > 50) {
    ED50_interp <- max(data$dose)
  } else {
    # dose-ids just above 50 and then below 50
    id_bef_after <- which.min(data$mean_resp > 50) - c(1, 0)
    dose_bef_after <- data$dose[id_bef_after]
    resp_bef_after <- data$mean_resp[id_bef_after]

    # intercept
    b <- (resp_bef_after[2] - resp_bef_after[1] *
            dose_bef_after[2]/dose_bef_after[1]) /
      (1 - dose_bef_after[2] / dose_bef_after[1])
    # slope
    m <- (resp_bef_after[2] - b) / dose_bef_after[2]

    ED50_interp <- (50-b)/m
  }
  return(ED50_interp)
}



# function: get_starting_values
#   calcultaes default starting values for exposure dose response data when ext3pLL model is fitted.
#
# Input:
#     - data:     data.frame(numerical) with columns expo, dose, resp
#     - h_start:  Starting value used for h in ext3pLL model. Default is 2.
#     - c0_start: Starting value used for threshold c0 in ext3pLL mpdel. Default is 0.
#
# Output:
#     - list with starting values for h, delta, gamma and c0
get_starting_values <- function(data, h_start = 2, c0_start = 0){
  # for starting value of "delta" consider ED50 estimate at lowest exposure time
  # idea: with linear interpolation: get two pairs (expo, ED50) and then starting values for gamma and delta
  data_low_exp_mean <- data %>% filter(expo == min(expo)) %>%
    group_by(dose) %>%
    summarize(mean_resp = mean(resp), .groups = "drop")

  data_high_exp_mean <- data %>% filter(expo == max(expo)) %>%
    group_by(dose) %>%
    summarize(mean_resp = mean(resp), .groups = "drop") %>%
    ungroup()

  ED50_start_low_exp <- interp_ED50(data = data_low_exp_mean)
  ED50_start_high_exp <- interp_ED50(data = data_high_exp_mean)


  expo_low_high <- range(data$expo)

  # get starting values for delta and gamma:
  gamma_start <- log(ED50_start_low_exp / ED50_start_high_exp, base = expo_low_high[1]*expo_low_high[2])
  delta_start <- mean(c(ED50_start_low_exp / expo_low_high[1]^(gamma_start),
                        ED50_start_high_exp / expo_low_high[2]^(gamma_start)))
  c0_start <- c0_start

  return(list(h = h_start,
              delta = delta_start,
              gamma = gamma_start,
              c0 = c0_start)
  )
}



# function: fitDERmod (fit Dose-Exposure-Response-Model)
#   Fits the extended 3pLL model to the data and returns a object using
#   the port algorithm in the nls function.
#   The bject has the specific class DERmod and general class nls.
# Input:
#   - data:     data.frame (numeric) with columns expo, dose and resp
#   - start:    list with named (numeric) starting values for h, delta, gamma and c0.
#               When the default (NULL) is used, h=2 and c0=0 are used and for
#               delta and gamma a linear interpolation procedure is used.
#   - control:  List of control arguments for nls function. For default (NULL),
#               maxiter = 100 and warnOnly = TRUE is used instead of default
#               maxiter = 50 and warnOnly = FALSE.
#   - lower:    named list or numeric of lower boundaries for parameters for
#               h, delta, gamma and c0 (in this order). For default (NULL),
#               1, -3*max(dose), -10 and 0 is used.
#   - upper:    Named list or numeric of upper boundaries for parameters for
#               h, delta, gamma and c0 (in this order). For default (NULL),
#               10, max(dose)*3, 10 and max(dose)*3 is used.
#   - trace:    trace argument for nls function. Default is FALSE.
#
# Output:
#   - fit:  an nls object of the fitted extended 3pLL model


fitDERmod <- function(data, start = NULL, control = NULL, lower = NULL, upper = NULL, trace = FALSE){

  if(!(all(colnames(data) %in% c("expo", "dose", "resp")))){
    stop("Argument data must be data.frame with columns expo, dose, resp.")
  }


  doses <- unique(data$dose)


  if(is.null(start)){
    start <- get_starting_values(data)
  } else {
    if((length(start) != 4) | !is.numeric(start)) stop("Argument start must be numeric of length 4.")
  }

  if(is.null(control)){
    control <- list(maxiter = 100, warnOnly = TRUE, printEval = FALSE)
  }

  if(is.null(lower)){
    lower <- c(h = 1, delta = -3*max(doses), gamma = -10, c0 = 0)
  } else {
    if((length(lower) != 4) |  !is.numeric(lower)) stop("Argument lower must be numeric of length 4.")
  }

  if(is.null(upper)){
    upper <- c(h = 10, delta = max(doses)*3, gamma = 10, c0 = max(doses)*3)
  } else {
    if((length(upper) != 4) |  !is.numeric(upper)) stop("Argument upper must be numeric of length 4.")
  }

  if(is.null(trace)){
    trace <- FALSE
  }
  # approach using means and weights
  data_w <-
    data %>% group_by(expo, dose) %>%
    summarize(n = n(),
           resp_m = mean(resp)
    ) %>%
      ungroup()

  # original. unweighted
  # fit <- nls(resp ~ 100 - 100*(dose^h) / ( ( delta * expo^(-gamma) + c0)^h + dose^h),
  #            data = data,
  #            algorithm = "port",
  #            lower = lower,
  #            upper = upper,
  #            start = start, trace = trace,
  #            control = control
  # )

  # new: weighted (hopefully more robust)
  fit <- nls(resp_m ~ 100 - 100*(dose^h) / ( ( delta * expo^(-gamma) + c0)^h + dose^h),
              data = data_w,
              weights = data_w$n,
              algorithm = "port",
              lower = lower,
              upper = upper,
              start = start,
             trace = trace,
              control = control
  )

  attr(fit, "class") <- c("DERmod", "nls") # to use plot.DERmod as method

  return(fit)

}

# function: fitJointmod
#           Fits a 2pLL model (upper limit = 100, lower limit = 0) to dose response data
# Input:
# - data: data.frame containing numeric columns resp and dose
#
# Output:
# drc-object containing the model
#
# Details:
# First, the method BFGS (default) is used. It his throws an error, the Nelder-Mead
# method is used.
# Note:
# We use the LL2.2 function. That means that the logEC50 parameterization is
# used.

fitJointmod <- function(data){
  tryCatch(
    {
      drm(resp ~ dose,
          data = data, fct = LL2.2(upper = 100))
    },
    error = function(cond){
      return(drm(resp ~ dose,
                 data = data, fct = LL2.2(upper = 100),
                 control = drmc(method = "Nelder-Mead")))
    }
  )
}

# function: fitSepmod
#           Fits seperate 2pLL model (upper limit = 100, lower limit = 0) to dose response data
# Input:
# - data: data.frame containing numeric columns resp and dose and expo
#
# Output:
# drc-object containing the model
#
# Details:
# First, the method BFGS (default) is used. It his throws an error, the Nelder-Mead
# method is used. Both the hilld parameter h and the ec50 are seperate for each exposure time.
# Note: We DO NO MORE use the LL2.2 function which uses log(EC50) as parameter for enhanced
# stability. When we want to get the true EC50s value of each, seperate curve,
# in case of LL2.2 we have to back-transform via exp("e:(Intercept)" + "e:expo2"), to get the
# EC50 of the second curve, for example.

fitSepmod <- function(data){
  tryCatch(
    {
      drm(resp ~ dose,
          curveid = expo,
          data = data %>% dplyr::mutate(expo = factor(expo)),
          fct = LL.2(upper = 100),
          pmodels = list(~expo, # h
                         ~expo) # EC50
      )
    },
    error = function(cond){
      drm(resp ~ dose,
          curveid = expo,
          control = drmc(method = "Nelder-Mead"),
          data = data %>% dplyr::mutate(expo = factor(expo)),
          fct = LL.2(upper = 100),
          pmodels = list(~expo, # h
                         ~expo) # EC50
      )
    }
  )
}

##############
# PLOTTING
##############


# function: plot.DERmod
#           Plot method for object of class DERmod as generated by fitDERmod OR
#           for vector of coefficients for ext3pLL model, when no DERmod object is generated
#           Uses the plotly function to create a three-dimensional, interactive plot
#           of the fitted model with data included.
# Input:
#   - DERmod:         DERmod object as genrated by fitDERmod function
#   - coefs:          Named vector with parameters for h, delta, gmma and c0
#   - dose_lim:       ranges of dose used in the plot: use first value above 0, because
#                     dose is on a logarithmic scale. Default is 1e-04
#   - expo_lim:       see dose_lim. Default is (1, 7)
#   - add_data:       data.frame (numeric, columns expo, dose and resp) with data.point that can
#                     be added to the plot
#   -n_grid:          Model is evaluated at equidistant grid of size n_grid^2 in dose_lim x expo_lim. Default is 100.
#   - title:          Plot title as character. If default (NULL), the parameters will be printed in the title.
#   - add_ED50_line:  Logical. Should a red ED50 line be added to the plot? Default is TRUE.
# Output: A plotly plot


plot.DERmod <- function(DERmod = NULL, coefs = NULL, dose_lim = c(1e-04,1), expo_lim = c(1,7),
                        add_data = NULL, n_grid = 100, title = NULL, add_ED50_line = T){

  expo_seq = seq(expo_lim[1], expo_lim[2], length.out = n_grid)
  dose_seq = c(0, exp(seq(log(dose_lim[1]), log(dose_lim[2]), length.out = n_grid)))
  # seq(dose_lim[1], dose_lim[2], length.out = n_grid)

  input_grid <- expand.grid(expo = expo_seq, dose = dose_seq) %>%
    as.data.frame()

  if(is.null(DERmod) & is.null(coefs)) stop("Either the DERmod argument or the coefs argument must be specified.")
  if(!(is.null(DERmod)) & !(is.null(coefs))) stop("Only DERmod or coefs can be specified.")

  # calculate the responses at each grid-point
  if(!is.null(coefs)){
    input_grid$resp <- apply(input_grid, 1, function(x){
      ext3pLL(dose = x["dose"],
              expo = x["expo"],
              h = coefs["h"],
              delta = coefs["delta"],
              gamma = coefs["gamma"],
              c0 = coefs["c0"])
    }
    )
  } else {
    input_grid$resp <- predict(DERmod, newdata = input_grid)
    coefs <- coef(DERmod)
  }

  # make matrix-style input for plot_ly function
  input_grid <- input_grid %>% pivot_wider(names_from = dose, values_from = resp) %>%
    dplyr::select(-expo) %>%
    as.matrix()

  if(is.null(title)){
    title <- paste0("h = ", round(coefs["h"], 3) ,"; delta = ", round(coefs["delta"], 3) ,";\n gamma = ",
                    round(coefs["gamma"],3) ,"; c0 = ", round(coefs["c0"], 3))
  }

  res_plot <- plot_ly(x = dose_seq, y = expo_seq, z = input_grid,
                      type="surface") %>%
    layout(title = title,
           scene = list(
             xaxis = list(title = "Dose",
                          tick0 = 0,
                          type = "log"),
             yaxis = list(title = "Exposure Time"),
             zaxis = list(title = "Response")
           ))

  # add single data points, if available
  if(!is.null(add_data)){
    res_plot <- res_plot %>% add_markers(x = add_data$dose, y = add_data$expo, z = add_data$resp,
                                         marker = list(size = 2), showlegend = F)
  }

  # add ED50 line, if wanted
  if(add_ED50_line){
    add_ED50 <- data.frame(
      expo = expo_seq,
      ED50 = coefs["delta"] * expo_seq^(-coefs["gamma"]) + coefs["c0"]
    )
    add_ED50 <- add_ED50 %>% filter(ED50 <= 1.02)
    add_ED50$resp <- apply(add_ED50, 1, function(x){
      ext3pLL(dose = x["ED50"],
              expo = x["expo"],
              h = coefs["h"],
              delta = coefs["delta"],
              gamma = coefs["gamma"],
              c0 = coefs["c0"])
    }
    )

    res_plot <- res_plot %>% add_trace(x = add_ED50$ED50, y = add_ED50$expo, z = add_ED50$resp,
                                       type = "scatter3d", mode = "lines",
                                       showlegend = F,
                                       line = list(color = "red", width = 6))
  }

  return(res_plot)
}




# function: plot_designs
#           plots the considered experimental designs as a ggplot.
#
# Input:
# - designs:  data.frame with columns expo, dose, n and design
#
# Details:
#   Right now, inceased doses in a log scale with base sqrt(10) is assumed
#   for plotting the dose-axes accordingly.

plot_designs <- function(designs){
  suppressWarnings(
  ggplot(data = designs, aes(x = dose, y = expo, size = as.factor(n))) +
    geom_point(alpha = 0.7, color = "steelblue") +
    geom_text(aes(label=n), size = 3, hjust=+0.5, vjust=0.5) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.0001, base = sqrt(10)),
                       breaks = unique(designs$dose),
                       labels = round(unique(designs$dose), 4)) +
    labs(y = "exposure time") +
    guides(size = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(facets = "design", labeller = "label_both")
  )
}



# function: my_pseudo_log
#           Creates a pseudo_log of the doses, so that dose=0
#           does not hrow an error
# Input:
# - x:  positive numeric
my_pseudo_log <- function(x) asinh(x/(2 * 0.0001))/log(sqrt(10))



#######################
# Data generation
#######################



# function: generate_data
#           generates data of an extended 3pLL model for specified noise and experimental design (expDes)
# Input:
# - model:    named, numeric character with model parameters h, delta, gamma and c0
# - noise:    character "N1", "N2, or "N3" for low, baseline or high noise level respectively
#             Details: For the baseline N2, we use a linear model that uses dose and expo as covariates,
#             which has to be tranformed via my_pseudo_log first, to define the variance of the mean-zero, normal distributed
#             noise at each dose-expo point of expDes. For N1, the variance is halfed. For N3, it is multiplied by 1.5.
# - expDes:   data.frame with columns expo, dose and n where n specifies how many replicates are measured at given expo-dose point
#
# Output:
# - res:      data.frame containing generated responses. Columns are expo, dose, resp

generate_data <- function(model, noise_id, expDes){
  # grid for expo, dose, n, h, delta, gamma, c0
  inputs <- cbind.data.frame(expDes %>% dplyr::select(-design), t(as.data.frame(model)) %>% `rownames<-`(NULL))

  # calculate mean_resp (true mean) at each expo dose point
  inputs$mean_resp <- apply(inputs, 1, function(x){
    ext3pLL(dose = x["dose"],
            expo = x["expo"],
            h = x["h"],
            gamma = x["gamma"],
            c0 = x["c0"],
            delta = x["delta"])
  }
  )

  # Calculate noise sd first
  inputs$noise_sd <- predict(lm_sd_noise, newdata = inputs[, c("expo", "dose")] %>% dplyr::mutate(dose = my_pseudo_log(dose)))

  # check via noise parameter, if sd will be halfed (N1)

  stopifnot(noise_id %in% c("N1", "N2", "N3"))

  if(noise_id == "N1"){
    inputs$noise_sd <- inputs$noise_sd * 0.5
  }

  if(noise_id == "N3"){
    inputs$noise_sd <- inputs$noise_sd * 1.5
  }

  # if noise_id == "N2": Do nothing, this is the baseline noise sd.


  # generate random noise values
  help_add_noise <- apply(inputs[, c("expo", "dose", "n", "noise_sd")], 1, function(x){
    suppressWarnings(
      cbind.data.frame(expo = x["expo"], dose = x["dose"], noise_value = rnorm(n = x["n"], sd = x["noise_sd"]))
    )
  }) %>% do.call(rbind.data.frame, .)

  # put together
  res_pre <- left_join(inputs, help_add_noise, by = c("expo", "dose")) %>%
    dplyr::mutate(resp = mean_resp + noise_value) %>%
    dplyr::select(expo, dose, resp)

  # apply regular pre-processing

  # First, divide by mean response at dose 0

  mean_resp_0 <- res_pre %>% dplyr::filter(dose == 0) %>% dplyr::select(resp) %>% unlist %>% mean
  res_pre <- dplyr::mutate(res_pre, resp = (resp / mean_resp_0) * 100 )

  # refit procedure: seperately fit 4pLL dose-response curves at each exposure.
  # Divide by resulting left (upper) asymptote
  all_expos <- unique(res_pre$expo)

  # get left asymptote (e0) or ach exposure time
  help_left_asymp <- sapply(all_expos, function(curr_expo){
    c(expo = curr_expo,
      fitMod(dose, resp,
             data = res_pre %>% dplyr::filter(expo == curr_expo) %>% dplyr::select(dose, resp),
             model = "sigEmax") %>%
        coef() %>%
        .["e0"]
    )
  }) %>% t() %>% as.data.frame()

  # divide by left asymptote stratified by exposure time
  res <- left_join(res_pre, help_left_asymp, by = "expo") %>%
    dplyr::mutate(resp = (resp / e0) * 100) %>%
    dplyr::select(expo, dose, resp)


  return(res)

}

##############################
# Data analysis
##############################



###############
# Pre-tests
###############

# We consider two pre-tests to decide if there is a dependency: anova and profileLL

# function: ext3pLL_anova_check
#           Checks with an anova-based approach, if the EC50 estimates depend on exposure time
# Input:
# - data:   data.frame with columns expo, dose, resp containing the data
# - alpha:  numeric in (0,1) representing the significance level. Default is 0.05.
#
# Output:
# - vector with two elements: p_value (numeric in (0,1)) and signif (logical). If signif = TRUE,
#           we reject the likelihood ratio null-hypothesis of EC50 being a common parameter.
#           Hence, if we reject, we will model the ext3pLL model, otherwise, we pool the data into
#           one dose-response curve (we ignrore the expo values).
#
# Details: We use an approximate ANOVA approach by comparing a model where the EC50 parameter
#          is shared between dose-response curves (i.e.: only one dose-response curve is fitted;
#          the exposure is neglected)
#          and a model where individual ec50 paraemters are fitted for each exposure time
#          (but the hill paraeter h is still shared over exposure times).
#
#          We use the LL2.2 instead of the LL.2 function for joint modeling for stability:
#          In LL2.2 the parameter log(EC50) instead of EC50 is used.
#          For the seperate model however, the LL.2 is used, as here, the LL2.2
#          seems to generate highly varying EC50 estimates

ext3pLL_anova_check <- function(data, alpha = 0.05){

  stopifnot(is.numeric(alpha) & alpha > 0 & alpha < 1)
  stopifnot(is.na(setdiff(colnames(data), c("expo", "dose", "resp"))))

  # dose-response curve with seperate EC50 parameters: Try default method BFGS first, then Nelder-Mead
  drm_seperate <- tryCatch(
    {
      drm(resp ~ dose,
         curveid = expo,
         data = data %>% dplyr::mutate(expo = factor(expo)),
         fct = LL.2(upper = 100),
         pmodels = list(~1, # h
                        ~expo)) # EC50
      },
    error = function(cond){
      return(drm(resp ~ dose,
          curveid = expo,
          control = drmc(method = "Nelder-Mead"),
          data = data %>% dplyr::mutate(expo = factor(expo)),
          fct = LL.2(upper = 100),
          pmodels = list(~1, # h
                         ~expo)) # EC50
      )
    }
  )

  # dose-response curves with shared EC50 parameters: Try default method BFGS first, then Nelder-Mead
  drm_joint <- fitJointmod(data = data)


  anova_res <- anova(drm_joint, drm_seperate)
  p_value <- anova_res$`p value`[2]
  signif <- ifelse(p_value < 0.05, TRUE, FALSE)

  return(c(p_value = p_value, signif = signif))
}


# function: get_EC50s
#           For parameters of a DERmod object (accesible via coef()) and given expo values,
#           the EC50 at said expo values is calculated.
# Input:
# - coefs:  named, numeric coefficients vector with paraemters h, delta, gamma and c0
# - expos:  numeric containing the exposure time values for which the EC50 is to be caculated
#
# Output:
# res:     data.frame with columns expo and EC50 containing the result

get_EC50s <- function(coefs, expos){

  stopifnot(length(setdiff(names(coefs), c("h", "delta", "gamma", "c0"))) == 0)
  stopifnot(is.numeric(expos) & length(expos) >= 1)
  stopifnot(all(expos >= 1 & expos <= 7))

  EC50s <- sapply(expos, function(expo) coefs["delta"] * expo^(-coefs["gamma"]) + coefs["c0"])

  res <- data.frame(expo = expos, EC50 = EC50s)

  return(res)
}

