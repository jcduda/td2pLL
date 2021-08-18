


# mytheta = c(h=2, delta=0.2, gamma=1.3, c0=0.2)
# Lb0 = c(1, 0);  Ub0 = c(10, 1)



opt_des_td2pLL = psoOptDesign(crit=Dcrit, nPoints=5, dimension=2, control = list(numIt = 100, numPart = 300, setProgressBar= TRUE),
                       wFixed=NULL, Lb= Lb0, Ub =Ub0,
                       gradient= grad_td2pLL_pso, theta= mytheta)


opt_des_td2pLL_n4 = psoOptDesign(crit=Dcrit, nPoints=4, dimension=2, control = list(numIt = 300, numPart = 300, setProgressBar= TRUE),
                              wFixed=NULL, Lb= Lb0, Ub =Ub0,
                              gradient= grad_td2pLL_pso, theta= mytheta)
# highest point in equivalence plot: 1.00 0.59 0.1863050603

opt_des_td2pLL_n6 = psoOptDesign(crit=Dcrit, nPoints=6, dimension=2, control = list(numIt = 300, numPart = 300, setProgressBar= TRUE),
                              wFixed=NULL, Lb= Lb0, Ub =Ub0,
                              gradient= grad_td2pLL_pso, theta= mytheta)
# highest point in equivalence plot: 1.00 0.27 0.181398373

opt_des_td2pLL_n6_2 = psoOptDesign(crit=Dcrit, nPoints=6, dimension=2, control = list(numIt = 200, numPart = 300, setProgressBar= TRUE),
                                 wFixed=NULL, Lb= Lb0, Ub =Ub0,
                                 gradient= grad_td2pLL_pso, theta= mytheta)

opt_des_td2pLL_n8 = psoOptDesign(crit=Dcrit, nPoints=8, dimension=2, control = list(numIt = 300, numPart = 300, setProgressBar= TRUE),
                              wFixed=NULL, Lb= Lb0, Ub =Ub0,
                              gradient= grad_td2pLL_pso, theta= mytheta)
# highest point in equivalence plot: 1.00 0.27 0.186544725

opt_des_td2pLL_n9 = psoOptDesign(crit=Dcrit, nPoints=9, dimension=2, control = list(numIt = 300, numPart = 300, setProgressBar= TRUE),
                                 wFixed=NULL, Lb= Lb0, Ub =Ub0,
                                 gradient= grad_td2pLL_pso, theta= mytheta)
# highest point in equivalence plot: 10.00 0.15 0.0264534054
# 2. run: highest point:              2.26 0.39 0.0204721413

opt_des_td2pLL_n9_2 = psoOptDesign(crit=Dcrit, nPoints=9, dimension=2, control = list(numIt = 500, numPart = 300, setProgressBar= TRUE),
                                 wFixed=NULL, Lb= Lb0, Ub =Ub0,
                                 gradient= grad_td2pLL_pso, theta= mytheta)
# highest point in equivalence plot: 1.00 0.59 0.1863050603

opt_des_td2pLL_n10 = psoOptDesign(crit=Dcrit, nPoints=10, dimension=2, control = list(numIt = 300, numPart = 300, setProgressBar= TRUE),
                                 wFixed=NULL, Lb= Lb0, Ub =Ub0,
                                 gradient= grad_td2pLL_pso, theta= mytheta)
# highest point in equivalence plot: 10.00 0.31 0.3493505109

opt_des_td2pLL_n11 = psoOptDesign(crit=Dcrit, nPoints=11, dimension=2, control = list(numIt = 300, numPart = 300, setProgressBar= TRUE),
                              wFixed=NULL, Lb= Lb0, Ub =Ub0,
                              gradient= grad_td2pLL_pso, theta= mytheta)
# highest point in equivalence plot: 10.00 0.30 0.069471944

opt_des_td2pLL_n12 = psoOptDesign(crit=Dcrit, nPoints=12, dimension=2, control = list(numIt = 300, numPart = 300, setProgressBar= TRUE),
                              wFixed=NULL, Lb= Lb0, Ub =Ub0,
                              gradient= grad_td2pLL_pso, theta= mytheta)
# highest point in equivalence plot: 1.00 0.27 1.962578e-01

opt_des_td2pLL = psoOptDesign(crit=Dcrit, nPoints=12, dimension=2, control = list(numIt = 100, numPart = 300, setProgressBar= TRUE),
                              wFixed=NULL, Lb= Lb0, Ub =Ub0,
                              gradient= grad_td2pLL_pso, theta= mytheta)




dcrit_equ_plot_td2pLL <- function(des, theta, time_lim = c(1, 10), dose_lim = c(0, 1),
                                  n_grid = 101, return_values = T){
  time_values <- seq(time_lim[1], time_lim[2], length = n_grid)
  dose_values <- seq(dose_lim[1], dose_lim[2], length = n_grid)
  df_values <- expand.grid(time = time_values, dose = dose_values)

  finfo_inv <- solve(finfo(w = des$w, M = des$supPoints, gradient = grad_td2pLL_pso,
              theta = theta))


  df_values$eq <- apply(df_values, 1, function(x){dcrit_equ_td2pLL(M_x = x,
                                                                des = des,
                                                                theta = theta,
                                                                finfo_inv = finfo_inv)
    })

  eq_values_matrix <- matrix(df_values$eq, nrow = n_grid, ncol = n_grid)

  equ_plot <- graphics::persp(time_values, dose_values, eq_values_matrix,
                theta = 45, ticktype = "detailed", xlab="time", ylab="dose",
                zlab= "equ_thm_val")

  if(return_values) return(list(values = df_values))

}


res <- dcrit_equ_plot_td2pLL(des = opt_des_td2pLL_n9,
                             theta = mytheta)

res$equ_plot

dcrit_aeq()
### ?berpr?fung, ob Design D-optimal ist
d1_values = seq(1.00, 10, length= 101) # Werte im Intervall der ersten Komponente
d2_values = seq(0.00, 1, length= 101) # Werte im Intervall der zweiten Komponente
d_values = expand.grid(d1_values, d2_values )
## Berechnung der Werte der ?quivalenzsatz-Ungleichung
z_values = apply(d_values, MARGIN = 1, FUN = function(x){dcrit_aeq(M_x=x, des = opt_des_td2pLL_n9, gradient = grad_td2pLL_pso, theta = mytheta)})
z_values = matrix(z_values, nrow=length(d1_values), ncol=length(d2_values))
### Plot
persp(d1_values, d2_values, z_values, theta = 45, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = 30, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = 10, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = -10, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = -30, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = -50, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = -75, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = -100, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = -125, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = -150, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")


persp(d1_values, d2_values, z_values, theta = 45, phi = 20, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = 45, phi = 40, ticktype = "detailed", xlab="time", ylab="dose", zlab= "")
persp(d1_values, d2_values, z_values, theta = 45, phi = 60, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = 45, phi = 80, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
persp(d1_values, d2_values, z_values, theta = 45, phi = 20, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")

df_sens_plot <- data.frame(d_values, z = as.vector(z_values))
df_sens_plot %>% filter(z >= 0) %>% arrange(z)

as.vector(z_values) %>% head
z_values %>% head

persp(d1_values, d2_values, z_values, theta = 30,
      phi = 45,
      r = sqrt(3),
      ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")
library(car)
scatter3d(x = d_values[,1], z = d_values[,2], y = as.vector((z_values)),
          grid = FALSE, residuals = "FALSE", fit = NULL, surface = TRUE, surface.alpha = 0,
          point.col = "black", fogtype = "linear")

