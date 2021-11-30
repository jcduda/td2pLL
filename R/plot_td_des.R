
#' @export
#' @title Visualize design calculated with [td_opt]
#' @param td_opt_des Optimal design for td2pLL model calculated with [td_opt]
#' @param label_weights (`logical(1)`) \cr
#'  Add labels with weights to the support points.
#' @param log10_dose (`logical(1)`) \cr
#'  If `TRUE` (default), the dose is displayed on a log10 scale.


plot_td_des <- function(td_opt_des, label_weights = TRUE, log10_dose = TRUE){
  des <- get_td_opt_des(td_opt_des)

  p <- ggplot(des, aes(x = dose, y = time, color = weight, size = weight,
                       label = round(weight, 3))) +
    geom_point() +
    scale_color_binned() +
    scale_size_binned() +
    guides(color = "none") +
    theme_bw()

  if(label_weights)  p <- p + geom_text_repel(color = "black", size = rel(3))

  if(log10_dose){
    p <- p + scale_x_log10()
  }

  return(p)

}
