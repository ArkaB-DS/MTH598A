## CODE TO CALCULATE THE POSTERIOR MEAN

## IMPORTING NECESSARY LIBRARIES
library(kmed) # REQUIRED FOR "heart" DATASET

data = heart

Y = ifelse(data[, 14] == 0, 0, 1)
model = glm(Y ~ . - 1 - class, family = "binomial", data = data)
X = model.matrix(model)[,-2]

model2 = glm(Y ~ X - 1, family = "binomial")

## defining target (posterior) distribution
f <- function(beta, s = 20) {
  beta = as.matrix(beta, ncol = 1)
  sum(dbinom(Y, size=1, prob=1/(1+exp(-X%*%beta)  ), log = TRUE  )  ) +
    sum(dnorm(beta, 0, s, log = TRUE))
}

# ADAPTIVE MCMC
run.adaprw <- function(N = 1e3, d = 18){

  beta <- matrix(0, nrow = N, ncol = d)
  beta[1,] <- coef(model2)

  Sigma <- 2.38^2/d*diag(d)

  acc.prob <- numeric(N)

  for(i in 2:N){

    if(i %% 1e3 == 0) print(i/1e3)

    if( i < 2*d ) {
      y <- rnorm( d , beta[i-1, ] , 0.1 / sqrt(d) )
    }else {
      b <- rbinom(1, size=1, prob = 0.35)
      Sigma <- cov(beta[1:(i-1), ])
      y <- b * rnorm(d, beta[i-1, ], 0.1/sqrt(d)) + (1 - b) *
        as.vector(mvtnorm::rmvnorm(1, mean = beta[i-1, ],sigma = 2.38^2*Sigma/d))
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
set.seed(8)
#mean(run.adaprw(N=1e4)$acc.prob)
#set.seed(8)
sample = run.adaprw(N = 1e5)
post.beta = colMeans(sample$mc)
#mean(sample$acc.prob)
#max(abs(colMeans(sample$mc)-coef(model2)))
#plot(cumsum(sample$acc.prob)/1:1e5,type="l",col="blue")

save(sample, post.beta, file = "betaPosteriorMean.Rdata")
