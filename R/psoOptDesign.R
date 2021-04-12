

#' @title PSO algorithm for multidimensional optimal designs
#' @description Calculate an optimal design for a multi-dimensional
#'  model, e.g. the td2pLL [fit_td2pLL] model, using particle swarm optimization
#'  (PSO) algorithm.
#'  Code kindly provided by Kirsten Schorning.
#'
#' @param crit (`function`)cr
#'   Function that calcultion the optimality criterion of interest.
#'   I.e. for D-optimality the determinant.
#'   The `crit` funciton must have parameters:
#'   `M` (matrix) with support points as rows,
#'   `w` (vector) with weights,
#'   `gradient` a function for the model's gradient depending on
#'     `x` (numeric(`dimension`)) which is a single point of the design and
#'     `theta` (numeric()) which contains the model parameters at which the
#'       gradient shall be evaluated.
#' @param control (named `list()`) \cr
#'   Named list with further arguments for the PSO algorithm: \cr
#'   `numit` (`integer(1)`) Number of iterations, default is 100 \cr
#'   `numPart` (`integer(1)`) Number of particles, default is 500 \cr
#'   `beta` (`numeric(1)`) ???, default is 0.5 \cr
#'   `gamma` (`numeric(1)`) ???, default is 0.7 \cr
#'   `setRsol` (`logical(1)`) ??, default is FALSE \cr
#'   `setProgressBar` (`logical(1)`) If a progress bar shall be printed,
#'    default is FALSE \cr
#'    `OutCritValue` (`logical(1)`) If the final value of the criterion shall be
#'      added to the output, default is FALSE. \cr
#'    `intmRes` (logical(1)) If intermediate results shall be ???, default is FALSE \cr
#' @param nPoints (`integer()`) \cr
#'  Number of design support points (number of rows of argument `M` for other functions).
#' @param dimension (`integer()`) \cr
#'  Number of dimensions for one support point (number of columns of argument `M`
#'  for other functions). Default is 1.
#' @param Lb (`numeric()`) \cr
#'  Lower bounds for the support points (e.g. for time and dose in td2pLL).
#' @param Ub (`numeric()`) \cr
#'  Upper bound for the support points (used e.g. for time and dose in td2pLL).
#' @param xFixed (Optional numeric `matrix`) \cr
#'  (Fixed) support Points, default is NULL. One point is one row. \cr
#' @param wFixed (Optional numeric(`nPoints`)) \cr
#'   (Fixed) weights. Must have length `nPoints` and weights must sum up to 1.
#'    Giving equal weights might improve the performance of the algorithm.
#' @param xold ???
#' @param nold ???
#' @param nold ???
#' @param n2 ???
#' @param ... ???



psoOptDesign<-function(crit,control=list(),nPoints=3,dimension=1, Lb=0,Ub=150, xFixed=NULL,wFixed=NULL,xold=NULL,nold=rep(0,length(xold)),n2=rep(0,length(old)),...){
  ###############################################################################

  con <- list(numIt = 100, numPart = 500, beta=0.5, gamma=0.7,setRsol=FALSE,setProgressBar=FALSE,OutCritValue=TRUE,intmRes=FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))


  if(!is.null(xFixed)&!is.null(wFixed))
    warning("weights and support points cannot both be predetermined")
  #if weights or support points are predetermined, so is the value of 'nPoints'
  if(!is.null(xFixed)){nPoints<-length(xFixed)}
  if(!is.null(wFixed)){nPoints<-length(wFixed)}
  Blow<-Lb
  Bup<-Ub
  timesteps<-con$numIt
  n<-con$numPart
  betas<-con$beta
  gammas<-con$gamma
  best<-array(0,dim=c(timesteps,nPoints, 1+dimension))

  Lb<-rep(c(0,Lb),each=nPoints)
  Ub<-rep(Ub,2*nPoints)



  #Initializing particles


  arg<-array(0,dim=c(nPoints*2,n, dimension));

  if(is.null(wFixed)) {
    arg[1:nPoints, ,dimension]<-matrix(runif(n*nPoints,0,1),nrow=nPoints, ncol=n)
    h<-colSums(arg[1:nPoints, , dimension])


    for(i in 1:nPoints){
      arg[i, ,dimension]<-arg[i, ,dimension]/h;		# Normierung der Gewichte durch h
    }
  }
  else{
    arg[1:nPoints, ,dimension]<-matrix(wFixed, nrow= nPoints, ncol=n, byrow=FALSE)
  }

  if(is.null(xFixed)) {
    supportvector=c()
    for (i in 1: (nPoints*n*dimension)){supportvector[i]=runif(1, Blow[(i%%dimension)+1], Bup[(i%%dimension)+1])}
    arg[(nPoints+1):(2*nPoints), ,]<-array(supportvector, dim=c(nPoints,n, dimension) )
  }
  else			# vielleicht dimension und npoints vertauschen, damit startvektor passt
  {
    arg[(nPoints+1):(2*nPoints), ,]<-array(xFixed,dim=c(nPoints,n, dimension))
  }


  zn<-1:n

  #==================================================================
  #This is where the actual alogorithm begins
  if(con$setProgressBar) pb<-txtProgressBar(min=0,max=timesteps,style=3)
  for(i in 1:timesteps){

    alphas<-gammas^(i/10)

    for(k in 1:n){
      if(sum(nold)==0) zn[k]<--crit(w=arg[1:nPoints,k,dimension],M=arg[(nPoints+1):(2*nPoints),k, ],...)
      else{
        wts<-c(1/(sum(nold)+n2)*nold,(1-sum(nold)/(sum(nold)+n2))*arg[1:(length(arg[,k, dimension])/2),k, ])
        dos<-c(xold,arg[(length(arg[,k, dimension])/2+1):(length(arg[,k, ])),k])
        zn[k]<--crit(wts,dos,...)
      }
    }
    #zn_min<-min(zn)

    wo<-1:nPoints
    xo<-matrix(0, nrow=nPoints, ncol=dimension)
    for(k in 1:nPoints){
      wo[k]<-arg[k,which.min(zn),dimension] ## hier stand noch min, wieso?
    }
    for(k in 1:nPoints){
      xo[k,]<-arg[(k+nPoints),which.min(zn), ] ##hier stand noch min, wieso?
    }

    zo<-min(zn)
    best[i,,]<-cbind(wo,xo) # Matrix: enth?lt beste Crit-Werte der i-ten Iteration zeilenweise erst Gewicht, dann zugeh?rigen Tr?gerpunkt (dim=dimension)

    #Move the particles
    #the weights...
    if(is.null(wFixed)){
      for(l in 1:n){
        arg[1:nPoints,l,dimension]<-arg[1:nPoints,l,dimension]*(1-betas)+wo*betas+alphas*rnorm(nPoints,0,0.25)
      }
    }
    # ... and the support points
    if(is.null(xFixed)){
      for(l in 1:n){
        arg[(nPoints+1):(2*nPoints),l,]<-arg[(nPoints+1):(2*nPoints),l,]*(1-betas)+xo*betas+alphas*matrix(rnorm(nPoints*dimension,0,(Bup-Blow)/4),
                                                                                                          nrow=nPoints, ncol=dimension, byrow=TRUE)
      }
    }

    #Make sure that the particles are still in the search room after having moved
    if(is.null(wFixed)){for (l in 1:n){arg[1:nPoints,l,dimension]<-pmax(arg[1:nPoints,l,dimension], .01)}}
    #if(is.null(wFixed)){for (k in 1:dimension) {arg[1:nPoints,,k]<-pmax(arg[1:nPoints,k,k],.005)}} ## Schranke fuer Gewicht eventuell zu hoch?
    if(is.null(xFixed)){for (k in 1:dimension) {arg[(nPoints+1):(2*nPoints),,k]<-pmax(arg[(nPoints+1):(2*nPoints),,k],Blow[k])}}
    if(is.null(wFixed)){for (l in 1:n){arg[1:nPoints,l,dimension]<-pmin(arg[1:nPoints, l, dimension], 1)}}
    #if(is.null(wFixed)){for (k in 1:dimension) {arg[1:nPoints,,k]<-pmin(arg[1:nPoints,,k],1)}}
    if(is.null(xFixed)){for (k in 1:dimension) {arg[(nPoints+1):(2*nPoints),,k]<-pmin(arg[(nPoints+1):(2*nPoints),,k],Bup[k])}}

    if(is.null(wFixed)){
      h<-colSums(arg[1:nPoints,,dimension])
      for(t in 1:nPoints){
        arg[t,,dimension]<-arg[t,,dimension]/h;
      }
    }
    if(con$setProgressBar) setTxtProgressBar(pb,i)
    if(con$intmRes)print(best[i,,])
  }
  if(con$setProgressBar) close(pb)
  if(con$setRsol){
    print("PSO-Algorithmus has finished. Intermediate results:")
    print(list(weights=round(wo,digits=5),supPoints=round(xo,digits=5)))
    print("start using solnp from package Rsol, please wait...") #Lagrange-Optimierung
    require(Rsolnp, quietly = TRUE)
    eqfun<-function(w) sum(w[1:nPoints])
    fac<-1/crit(wo,xo,...)
    fun<-function(w) -fac*crit(w[1:nPoints],w[(nPoints+1):(2*nPoints)],...)
    #{ sink("NUL"); z<-round(solnp(pmin(pmax(c(wo,xo),Blow+1e-15),Bup-1e-15),fun, eqfun, eqB = 1,LB=Lb,UB=Ub,control=list(tol=1e-15))$pars,digits = 5); sink(); }
    z<-round(solnp(pmin(pmax(c(wo,xo),Blow+1e-15),Bup-1e-15),fun, eqfun, eqB = 1,LB=Lb,UB=Ub)$pars,digits = 5)
    wo<-z[1:nPoints]
    xo<-z[(nPoints+1):(2*nPoints)]
  }

  #hin<-function(w) c(w[1:nPoints]+.1,w[(nPoints+1):(2*nPoints)]-Lb,Ub+.1-w[(nPoints+1):(2*nPoints)])
  #constrOptim.nl(c(wo,xo), fun, heq=eqfun, hin=hin)
  #Output the results
  #cat("\n")
  if(con$OutCritValue)
    return(list(value=crit(round(wo, digits=5),round(xo, digits=5),...),weights=round(wo,digits=5),supPoints=round(xo,digits=5)))
  else return(list(weights=round(wo,digits=5),supPoints=round(xo,digits=5)))
}
