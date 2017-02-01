## bayesian stat
## 

library(ggplot2)
theta_true <- runif(1,0,1)

N = sample(seq(100000,400000),1)
A = round(theta_true*N)
B = N - A
Zpop <- sample(c(rep(1,A),rep(0,B)))



N_samp <- 500
Zsamp <- sample(Zpop,N_samp)


mean(Zsamp)


betaplot <- function(a,b){
  theta = seq(0,1,0.005)
  p_theta = dbeta(theta, a, b)
  p <- qplot(theta, p_theta, geom='line')
  p <- p + theme_bw()
  p <- p + ylab(expression(paste('p(',theta,')', sep = '')))
  p <- p + xlab(expression(theta))
  return(p)
}

dbeta(theta, 1, 1)

betaplot(10,10)

## Essentially, higher values of the ratio of α to β weights higher values of θ, lower values of that ratio place greater weight on lower θ values, and higher value of α+β indicates higher certainty.


prior <- function(m,n){
  a = n * m
  b = n * (1 - m)
  dom <- seq(0,1,0.005)
  val <- dbeta(dom,a,b)
  return(data.frame('x'=dom, 'y'=val))
}
likelihood <- function(N,Y){
  a <- Y + 1
  b <- N - Y + 1
  dom <- seq(0,1,0.005)
  val <- dbeta(dom,a,b)
  return(data.frame('x'=dom, 'y'=val))
}

getprior <- prior(0.5, 5)


hist(dbeta(c(rep(0.1, 10), rep(0.6, 5)), 10, 10), breaks = 100)
hist(dbeta(c(rep(0.1, 10), rep(0.6, 5)), 1, 10), breaks = 100, add = T)

#qplot(getprior$x, getprior$y)



