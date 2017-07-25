### Bayesian on the error rates of methylation in reads
##
#


## Simulations, with specific number of cytosines and error rates
#

error.rates <- 0.02  ### Conversion of 98% 

num.cytosines <- 4


likelihoodFunction <- function(p, mC, tC){
  a <- choose(n = tC,k = mC)
  b <- p ^ mC
  c <- (1-p) ^ (tC - mC)
  return(a * b * c)
}

likelihoodFunction(5, 2, 0.99999)

getPoints <- function(mC, tC){
  x <- seq(0, 1, length.out = 200)
  y <- sapply(x, FUN = function(x){likelihoodFunction(p = x, mC = mC, tC = tC)})
  return(cbind(x,y))
}


plot(getPoints(1,5))   ### 2/5
points(getPoints(4,10), col="blue")   ## 4/10
points(getPoints(40,100), col="red")   ## 4/10
points(getPoints(1,3), col = "orange")



optimise(interval = c(0,1),f = likelihoodFunction, mC = 1, tC = 3, maximum = T, tol = 0.001)












