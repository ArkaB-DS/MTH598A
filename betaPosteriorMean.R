library(kmed)
data <- heart

glmModel <- glm(Y ~ . - class - 1, data = data)
X <- model.matrix(glmModel)[, -2]
Y <- ifelse(data[, 14] == 0, 0, 1)

## defining target (posterior) distribution
f <- function(beta, s=20) sum(dbinom(Y, size = 1, prob = 1/(1 + exp(-X %*% beta))), log = TRUE) +
  sum(dnorm(beta, 0, s, log = TRUE))

run.rw <- function(N = 1e3, d = 18){

  beta = matrix(0, nrow = N, ncol = d)
  acc.prob <- numeric(N)

  for(i in 2:N){
    if(i%%1e3==0) print(i/1e3)
    y = rnorm(d, 0, 0.5) + beta[i-1, ]
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

set.seed(8)
sigmaHat <- cov(run.rw(N = 1e5)$mc)

N <- 1e5
d <- 18

z <- matrix(rnorm(N*d), nrow = N, ncol = d)

beta = matrix(0, nrow = N, ncol = d)
beta[1, ] = coef(glmModel)[-2]

acc.prob <- numeric(N)

c <- 1

for(i in 2:N){
  if(i%%1e3==0) print(i/1e3)

  y <- as.vector(beta[i-1, ] + c * z[i, ] %*% chol(sigmaHat))

  alpha = f(y) - f(beta[i-1, ])

  if(log(runif(1)) < alpha){
    beta[i, ] = y
    acc.prob[i] = 1
  }else{
    beta[i,] = beta[i-1, ]
  }

  c = c + 1/i * ( exp(alpha) - 0.23 )
}
