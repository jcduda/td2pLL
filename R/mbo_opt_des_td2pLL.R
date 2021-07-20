


theta =  c(h=2, delta=0.2, gamma=1.3, c0=0.2)

n_points = 7
weights = NULL

# to optimize:
# time
# doses
# weights

# for defining constraints:
# time_bnds
# dose_bnds
# weight_bnds

# preparing objective function first:

# get n_points

get_n_points <- function(n_points, weights){
  if(!is.null(weights) & !is.null(n_points)) {
    stopifnot(round(n_points) == n_points)
    stopifnot(length(weights) != length(n_points))
  }
  if(!is.null(n_points)){
    stopifnot(round(n_points) == n_points)
    if(n_points < 7) {
      warning("Choosing n_points < 7 is not advised. n_points.
              n_points will be set to 7.")
      return(7)
    } else {
      return(n_points)
    }
  }
  if(is.null(n_points) & !(is.null(weights))){
    warning(paste0("weights provided but not n_points. Set n_points to
            length of weights:", length(weights)))
    return(7)
  }
  if(is.null(n_points)) return(7)

}

n_points <- get_n_points(n_points, weights)

# get gradient

get_td4pLL_grad <- function(e0_fix, Emax_fix){
  if(is.null(e0_fix) & !(is.null(Emax_fix))){
    grad <- function(time, dose, theta){
      h <- theta["h"]; delta <- theta["delta"]; gamma <- theta["gamma"]
      c0 <- theta["c0"]
      # To avoid log(0) throwing an error
      lg2 <- function(x) ifelse(x == 0, 0, log(x))

      B <- delta * time^(-gamma) + c0
      const1 <- Emax_fix*dose^h / (B^h + dose^h)^2

      # h
      g1 <- const1 * B^h * lg2(dose / B)
      # delta
      g2 <- const1 * sign(Emax_fix) * h * B^(h-1) * time^(-gamma)
      # gamma
      g3 <- const1 * B^(h-1) * h * delta * lg2(time) * time^(-gamma)
      # c0
      g4 <- const1 * sign(Emax_fix) * B^(h-1) * h

      cbind(e0 = 1, h = g1, delta = g2, gamma = g3, c0 = g4)
    }
    return(grad)
  }
}

e0_fix = NULL; Emax_fix = -100
grad <- get_td4pLL_grad(e0_fix, Emax_fix)


# compare to previous implementation:
# i <- 1
# grad(time =  times[i], dose = doses[i], theta = theta)
# grad_td2pLL_pso(x = cbind(times[i], doses [i]), theta = theta)
# fits!



times <- c(3, 1, 1, 3, 3, 7, 7)
doses <- c(0, 0.4, 0.6, 0.3, 0.5, 0.1, 0.7)
weights <- c(1/5, rep(2/15, 6))

get_finfo_td4pLL <- function(e0_fix, Emax_fix, theta, grad){
  if(is.null(e0_fix) & !is.null(Emax_fix)){
    finfo <- function(times, doses, weights){

      # check dimensions
      tmp <- weights[1] * t(grad(times[1], doses[1], theta)) %*% grad(times[1], doses[1], theta)
      for(i in 2:length(times)){
        tmp <- tmp + weights[i] * t(grad(times[i], doses[i], theta)) %*% grad(times[i], doses[i], theta)
      }
      return(tmp)
    }
    return(finfo)
  }
}

# finfo is a function
finfo <- get_finfo_td4pLL(e0_fix, Emax_fix, theta, grad)

# ex:
finfo1 <- finfo(times, doses, weights)
finfo1

get_Dcrit <- function(times, doses, weights, finfo){
  function(times, doses, weights){
    return(-log(det(finfo(times, doses, weights))))
  }
}

Dcrit <- get_Dcrit(times, doses, weights, finfo)

# Dcrit works:
times;doses;weights
Dcrit(times, doses, weights)



# still missing
# init_des <- function(n_points, theta, e0_fix, Emax_fix, time_bnds, dose_bnds, weight_bnds){}
#
# des <- list(times = times, doses = doses, weights = weights,
#             y = Dcrit(times, doses, weights, finfo))

# make initial design des


time_bnds <- c(1, 10)
dose_bnds <- c(0, 1)

obj.fun <- makeSingleObjectiveFunction(
  name = "td2pLL_odes_nU",
  fn = Dcrit,
  has.simple.signature = F,
  constraint.fn = function(weights) {round(sum(weights), 4) == 1},
  par.set = makeParamSet(
    makeNumericVectorParam("times", len = n_points,
                           lower = time_bnds[1], upper = time_bnds[2]),
    makeNumericVectorParam( "doses", len = n_points,
                      lower = dose_bnds[1], upper = dose_bnds[2]),
    makeNumericVectorParam("weights", len = n_points,
                     lower = 0, upper = 1)
  ),
  minimize = TRUE
)
des <- generateDesign(n = 1, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
des$y = apply(des, 1, obj.fun)



mbo(obj.fun, )

#
#
# run = mbo(obj.fun, design = des, learner = surr.km, control = control, show.info = TRUE)
#
#
# makeSingleObjectiveFunction(
#   name = "my_sphere",
#   fn = function(x) {
#     sum(x*x) + 7
#   },
#   par.set = makeParamSet(
#     makeNumericVectorParam("x", len = 2L, lower = -5, upper = 5)
#   ),
#   minimize = TRUE
# )







# compare to previous implementation:
#
#
# finfo <- function(w, M, gradient, theta){
#   stopifnot(length(w) == nrow(M))
#   # Add linear parameter e0
#   dimTheta <- length(theta) + 1
#   tmp <- matrix(0, dimTheta, dimTheta)
#   for(i in 1:nrow(M)){
#     tmp = tmp + w[i] * (gradient(M[i, ], theta) %*% t(gradient(M[i, ], theta)))
#   }
#
#   return(tmp)
# }
# w = weights
# M = cbind(times, doses)
# gradient = grad_td2pLL_pso
# finfo2 <- finfo(w = weights, M = M, gradient = gradient, theta = theta)
# finfo2
# det(finfo2)






