#' @title Filtering the design from a minimax object
#' @description From a minimax object generated with [ICAOD::locally()]
#' or [td2pLL::td_opt()], extract only the design, i.e. the support
#' points and weights
#' @param td_top_des An objetc of class `minimax`
#' @return A data.frame with columns `time` and `dose` for the support points
#'  and `weight` for the corresponding weight

get_td_opt_des <- function(td_opt_des){
  des_names <- names(td_opt_des$design)
  times <- unlist(td_opt_des$design[grep("x1\\d", des_names, value = T)])
  doses <- unlist(td_opt_des$design[grep("x2\\d", des_names, value = T)])
  weights <- unlist(td_opt_des$design[grep("w\\d", des_names, value = T)])

  return(data.frame(time = times,
                    dose = doses,
                    weight = weights, row.names = NULL))
}
