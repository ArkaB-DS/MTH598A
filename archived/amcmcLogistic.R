## PROGRAM TO COMPARE COVERAGE PROBABILITIES IN ADAPTIVE AND NON-ADAPTIVE MCMC USING "HEART" DATASET

# IMPORTING NECESSARY LIBRARIES

library(kmed)
library(mvtnorm)
library(mcmcse)

# DEFINING DATA, X AND Y

data = heart

Y = ifelse(data[, 14] == 0, 0, 1)

model = glm(Y ~ . - 1 - class, family = "binomial", data = data)
X = model.matrix(model)[,-2]

model2 = glm(Y ~ X - 1, family = "binomial")

# DEFINING TARGET DISTRIBUTION

f <- function(beta, s = 20) {
  beta = as.matrix(beta, ncol = 1)
  sum(dbinom(Y, size=1, prob=1/(1+exp(-X%*%beta)  ), log = TRUE  )  ) +
    sum(dnorm(beta, 0, s, log = TRUE))
}

# LOADING POSTERIOR MEAN OF BETA

load("C://Users//ARKAJYOTI//Desktop//MTH598A//betaPosteriorMean.Rdata")

# post.beta = colMeans(sample$mc)

# NON ADAPTIVE MCMC

run.rw <- function(N = 1e3, d = 18){

  beta = matrix(0, nrow = N, ncol = d)
  beta[1,] = coef(model2)

  acc.prob <- numeric(N)

  for(i in 2:N){
    # if(i%%1e3==0) print(i/1e3)

    y = rnorm(d, 0, 0.0026) + beta[i-1, ]

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
# set.seed(8)
# sample = run.rw(N=1e6)
# mean(run.rw(N=1e5)$acc.prob)
# plot(cumsum(run.rw(N=1e5)$acc.prob)/1:1e5, col="blue", type = "l")

# ADAPTIVE MCMC

run.adaprw <- function(N = 1e3, d = 18){

  beta <- matrix(0, nrow = N, ncol = d)
  beta[1,] <- coef(model2)

  Sigma <- 2.38^2/d*diag(d)

  acc.prob <- numeric(N)

  for(i in 2:N){

    # if(i %% 1e3 == 0) print(i/1e3)

    if( i < 2*d ) {
      y <- rnorm( d , beta[i-1, ] , 0.1 / sqrt(d) )
    }else {
      b <- rbinom(1, size=1, prob = 0.35)
      Sigma <- cov(beta[1:(i-1), ])
      y <- b * rnorm(d, beta[i-1, ], 0.1/sqrt(d)) + (1 - b) *
        as.vector(rmvnorm(1, mean = beta[i-1, ],sigma = 2.38^2*Sigma/d))
    }

    alpha <- f(y) - f(beta[i-1, ])

    if(log(runif(1)) < alpha){
      beta[i, ] <- y
      acc.prob[i] <- 1
    }else{
      beta[i,] <- beta[i-1, ]
    }
  }

  return(list(mc = beta, acc.prob = acc.prob))
}
# set.seed(8)
# MC = run.adaprw(N=1e5)
# mean(MC$acc.prob)
# plot(cumsum(MC$acc.prob)/1:length(MC$acc.prob), col="blue", type = "l")


# FUNCTION FOR FINDING COVERAGE POBABILITY IN NON-ADAPTIVE MCMC

checkcoverage<-function(alpha = 0.05, iterations = 5e2, N = 1e3, d = 18){

  in.or.out <- numeric(iterations)

  for (j in 1:iterations){

    if(j %% 1e2 == 0) cat("Iteration-->", j/1e2, "\n")

    simpleMC <- run.rw(N = N, d = d)
    MC <- simpleMC$mc

    Gamma.n.sq.h <- mcmcse::mcse.multi(MC, method = "bartlett", r = 1, size = "sqroot")$cov

    T <- N*t(colMeans(MC)-coef(model2)) %*% solve(Gamma.n.sq.h) %*% (colMeans(MC)-coef(model2))

    in.or.out[j] <- as.numeric(T >= qchisq(alpha/2, d) & T <= qchisq(1-alpha/2, d))
  }
  return (mean(in.or.out))
}


checkcoverage.adap<-function(alpha = 0.05, iterations = 5e2, N = 1e3, d = 18){

  in.or.out <- numeric(iterations)

  for (j in 1:iterations){

    if(j %% 1e2 == 0) cat("Iteration-->", j/1e2, "\n")

    adapMC <- run.adaprw(N = N, d = d)
    MC <- adapMC$mc

    Gamma.n.sq.h <- mcse.multi(MC, method = "bartlett", r = 1, size = "sqroot")$cov

    T <- N*t(colMeans(MC)-post.beta) %*% solve(Gamma.n.sq.h,tol = 1e-20) %*% (colMeans(MC)-post.beta)

    in.or.out[j] <- as.numeric(T >= qchisq(alpha/2, d) & T <= qchisq(1-alpha/2, d))
  }
  return (mean(in.or.out))
}


set.seed(8)
coverage <- matrix(0, nrow = 1, ncol = 3)
coverage[,1] <- checkcoverage(N = 1e3)
coverage
coverage[,2] <- checkcoverage(N = 1e4)
coverage
coverage[,3] <- checkcoverage(N = 1e5)
colnames(coverage) <- c("1e3 iterations","1e4 iterations","1e5 iterations")
rownames(coverage) <- "coverage probability"
coverage

# coverage
#                       1e3 iterations 1e4 iterations 1e5 iterations
# coverage probability          0.006          0.054              0

coverage.adap <- matrix(0, nrow = 1 , ncol = 3)
coverage.adap[,1] <- checkcoverage(N = 1e3)
coverage.adap
coverage.adap[,2] <- checkcoverage(N = 1e4)
coverage.adap
coverage.adap[,3] <- checkcoverage(N = 1e5)
colnames(coverage.adap) <- c("1e3 iterations", "1e4 iterations", "1e5 iterations")
rownames(coverage.adap) <- "coverage probability"
coverage.adap
# > coverage.adap
#                       1e3 iterations 1e4 iterations 1e5 iterations
# coverage probability          0.004          0.038              0
