## "heart" dataset is used from "kmed" package
data = kmed::heart

## data pre-processing
data_X = data[, 1:13]
data_X[, 2] = as.factor(as.numeric(data_X[, 2]))
data_X[, 6] = as.factor(as.numeric(data_X[, 6]))
data_X[, 9] = as.factor(as.numeric(data_X[, 9]))
data_X = data.matrix(data_X)
#data_X[,c(4,5,8)] = data_X[,c(4,5,8)]/1e3
X = data_X
Y = data.matrix(ifelse(data[, 14]==0, 0, 1))

## defining target (posterior) distribution
f <- function(beta, s=20) sum(dbinom(Y, size=1, prob=1/(1+exp(-X%*%beta))), log=TRUE)+
  sum(dnorm(beta, 0, s, log = TRUE))

## posterior mean calculated from 1e6 iterations is taken as the true mean
post.beta <-c( -0.027515432,  0.271487958,  0.477965324,  0.012726262,  0.001501614,  0.008701202,
               0.208373899, -0.034807646,  0.424081323,  0.462574606, -0.025540309,   0.560421448,
               0.436959322)

## NON ADAPTIVE MCMC
# plain RWM with Gaussian(0,sigma^2*I_d)

run.rw <- function(N = 1e3, d = 13){
  beta = matrix(0, nrow = N, ncol = d)
  acc.prob <- numeric(N)
  for(i in 2:N){
    # if(i%%1e3==0) print(i/1e3)
    y = rnorm(d, 0, 0.008) + beta[i-1, ]
    alpha = f(y) - f(beta[i-1, ])
    if(log(runif(1)) < alpha){
      beta[i, ] = y
      acc.prob[i] = 1
    }else{
      beta[i, ] = beta[i-1, ]
    }
  }
  return(list(mc = beta, acc.prob = acc.prob))
}
#mean(run.rw()$acc.prob)

checkcoverage<-function(alpha = 0.05, iterations = 5e2, N = 1e3, d = 13){
  in.or.out <- numeric(iterations)
  for (j in 1:iterations){
    if(j %% 1e2 == 0) cat("Iteration-->", j/1e2, "\n")
    simpleMC <- run.rw(N = N, d = d)
    MC <- simpleMC$mc
    Gamma.n.sq.h <- mcmcse::mcse.multi(MC, method = "bartlett", r = 1, size = "sqroot")$cov
    T <- (N-d)/d/(N-1) * t(colMeans(MC)-post.beta) %*% solve(Gamma.n.sq.h) %*% (colMeans(MC)-post.beta)
    in.or.out[j] <- as.numeric(T >= qf(alpha/2, d, N-d) & T <= qf(1-alpha/2, d, N-d))
  }
  return (mean(in.or.out))
}

# ADAPTIVE MCMC
run.adaprw <- function(N = 1e3, d = 13){
  beta = matrix(1, nrow = N, ncol = d)
  beta[1,] = rnorm( d , 0 , 0.1 / sqrt(d) )
  Sigma = 2.38^2/d*diag(d)
  acc.prob <- numeric(N)
  for(i in 2:N){
    # if(i%%1e3==0) print(i/1e3)
    if( i < 2*d ) {
      y = rnorm( d , beta[i-1, ] , 0.1 / sqrt(d) )
    }else {
      b = rbinom(1, size=1, prob = 0.15)
      Sigma = cov(beta[1:(i-1), ])
      y = b * rnorm(d, beta[i-1, ], 0.1/sqrt(d)) + (1 - b) *
        as.vector(mvtnorm::rmvnorm(1, mean = beta[i-1, ],sigma = 2.38^2*Sigma/d))
    }
    alpha = f(y) - f(beta[i-1, ])
    if(log(runif(1)) < alpha){
      beta[i, ] = y
      acc.prob[i] = 1
    }else{
      beta[i,] = beta[i-1, ]
    }
  }
  return(list(mc = beta, acc.prob = acc.prob))
}
#mean(run.adaprw()$acc.prob)

checkcoverage.adap<-function(alpha = 0.05, iterations = 5e2, N = 1e3, d = 13){
  in.or.out <- numeric(iterations)
  for (j in 1:iterations){
    if(j %% 1e2 == 0) cat("Iteration-->", j/1e2, "\n")
    adapMC <- run.adaprw(N = N, d = d)
    MC <- adapMC$mc
    Gamma.n.sq.h <- mcmcse::mcse.multi(MC, method = "bartlett", r = 1, size = "sqroot")$cov
    T <- (N-d)/d/(N-1) * t(colMeans(MC)-post.beta) %*% solve(Gamma.n.sq.h) %*% (colMeans(MC)-post.beta)
    in.or.out[j] <- as.numeric(T >= qf(alpha/2, d, N-d) & T <= qf(1-alpha/2, d, N-d))
  }
  return (mean(in.or.out))
}
# checkcoverage.adap()

set.seed(8)
coverage <- matrix(0, nrow = 1, ncol = 3)
coverage[1] <- checkcoverage(N = 1e3)
coverage
coverage[2] <- checkcoverage(N = 1e4)
coverage
coverage[3] <- checkcoverage(N = 1e5)
colnames(coverage) <- c("1e3 iterations","1e4 iterations","1e5 iterations")
rownames(coverage) <- "coverage probability"
coverage
# coverage[3,] <- checkcoverage(N = 1e5)
# coverage
# 1     2     3     4     5     6     7     8     9    10    11    12    13
# 1e3 iterations 0.058 0.032 0.006 0.040 0.054 0.130 0.070 0.036 0.004 0.016 0.124 0.002 0.006
# 1e4 iterations 0.060 0.102 0.102 0.006 0.024 0.094 0.116 0.016 0.084 0.082 0.092 0.102 0.100

coverage.adap <- matrix(0, nrow = 1 , ncol = 3)
coverage.adap[1] <- checkcoverage(N = 1e3)
coverage.adap
coverage.adap[2] <- checkcoverage(N = 1e4)
coverage.adap
coverage.adap[3] <- checkcoverage(N = 1e5)
colnames(coverage.adap) <- c("1e3 iterations", "1e4 iterations", "1e5 iterations")
rownames(coverage.adap) <- "coverage probability"
coverage.adap
