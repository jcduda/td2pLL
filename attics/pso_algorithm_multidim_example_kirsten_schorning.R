# Example Kirsten Schorning for multidimensional PSO



### Betrachte nun konkret: f(x, y)= (1, x, y, x*y) als Regressionsfunktion
### Informationsmatrix des Beispiels
info  = function(w, M)
{
  k = length(w)
  result = matrix(0, ncol= 4, nrow= 4)
  for(l in 1:k)
  {
    result = result + w[l]*c(1, M[l,1], M[l, 2], M[l,1]*M[l,2])%*%t(c(1, M[l,1], M[l, 2], M[l,1]*M[l,2]))
  }
  return(result)
}

### D-Informationskriterium
dcrit = function(w, M)
{
  det(info(w, M))
}


### Aequivalenzsatz-Ungleichung f?r die vierdimensionale Regressionsfunktion f(x, y) = (1, x, y, x*y)
dcrit_aeq = function(d, des)
{
  w= des$weights
  M = des$supPoints
  info0 = info(w, d)
  minfo0= solve(info0)
  result = t(c(1, d[1], d[2], d[1]*d[2])) %*% minfo0 %*% c(1, d[1], d[2], d[1]*d[2]) -nrow(info0)
  return(result)
}


### Berechnung des optimalen Versuchsplans mit 4 Punkten
###(Wenn Anzahl Punkte = Anzahl Parameter, sind die Gewichte immer gleich: hier 1/4, siehe wFixed)
### Design Space [0, 1]x[0.5, 10] -> Lb =c(0, 0.5) (beide lower bounds), Ub = c(1, 10) (beide upper bounds)
### dimension des Design Space : 2
opt_des = psoOptDesign(dcrit,nPoints=4,dimension=2, Lb=c(0,0.5),Ub=c(1,10), wFixed = rep(1/4, 4))

### ?berpr?fung, ob Design D-optimal ist
d1_values = seq(0, 1, length= 10) # Werte im Intervall der ersten Komponente
d2_values = seq(0.5, 10, length= 10) # Werte im Intervall der zweiten Komponente
d_values = expand.grid(d1_values, d2_values )
## Berechnung der Werte der ?quivalenzsatz-Ungleichung
z_values = apply(d_values, MARGIN = 1, FUN = function(x){dcrit_aeq(d=x, des = opt_des)})
z_values = matrix(z_values, nrow=length(d1_values), ncol=length(d2_values))
### Plot
persp(d1_values, d2_values, z_values, theta = 45, ticktype = "detailed", xlab=" ", ylab=" ", zlab= "")

