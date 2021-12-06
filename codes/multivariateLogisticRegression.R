##############################################
## Muktivariate Logistic Regression Example ##
##############################################

data(logit, package = "mcmc") # LOADING "LOGIT" DATASET FROM "MCMC" PACKAGE
data <- logit

X <- as.matrix(data[, -1]) # DEFINING DESIGN MATRIX
Y <- data[, 1] # DEFINING Y VARIABLE

# DEFINING TARGET DISTRIBUTION
f <- function(beta, s = 100) {
  beta = as.matrix(beta, ncol = 1)
  sum(dbinom(Y, size=1, prob=1/(1+exp(-X%*%beta)  ), log = TRUE  )  ) +
                                 sum(dnorm(beta, 0, s, log = TRUE))
}

## NON-ADAPTIVE MCMC
# Random-Walk Metropolis with Gaussian(0,sigma^2*I_d)

run.rw <- function(N = 1e3, d = 4){

  beta <- matrix(0, nrow = N, ncol = d)

  acc.prob <- numeric(N)

  for(i in 2:N){

    # if(i%%1e3==0) print(i/1e3)

    y <- rnorm(d, 0, 0.5) + beta[i-1, ]

    alpha <- f(y) - f(beta[i-1, ])

    if(log(runif(1)) < alpha){

      beta[i, ] <- y
      acc.prob[i] <- 1

      }else{

        beta[i, ] <- beta[i-1, ]
    }
  }

  return(list(mc = beta, acc.prob = acc.prob))
}

# set.seed(8)
# mean(run.rw(N = 1e5)$acc.prob)


# ADAPTIVE MCMC
run.adaprw <- function(N = 1e3, d = 4){

  beta <- matrix(0, nrow = N, ncol = d)
  beta[1,] <- coef(glm(Y ~ X - 1, family = "binomial"))

  Sigma <- 2.38^2/d*diag(d)

  acc.prob <- numeric(N)

  for(i in 2:N){

    # if(i %% 1e3 == 0) print(i/1e3)

    if( i < 2*d ) {

      y <- rnorm( d , beta[i-1, ] , 0.1 / sqrt(d) )

      }else {

        b <- rbinom(1, size=1, prob = 0)
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
# set.seed(8)
# mean(run.adaprw(N = 1e5)$acc.prob)

# sample <- run.adaprw(N = 1e5)$mc
# par(mfrow = c(2,2))
# plot.ts(sample)

## CALCULATION OF POSTERIOR MEAN
load("C:/Users/ARKAJYOTI/Desktop/MTH598A/betaPosteriorMean.Rdata")
#post.beta <- c(0.7781330, 0.9457948, 0.4558997, 0.6298492)

## CALCULATION OF COVERAGE PROBABILITY FOR NON-ADAPTIVE MC
checkcoverage <- function(alpha = 0.05, iterations = 5e2, N = 1e3, d = 4){

  in.or.out <- numeric(iterations)

  for (j in 1:iterations){
    if(j %% 1e2 == 0) cat("Iteration-->", j/1e2, "\n")

    simpleMC <- run.rw(N = N, d = d)
    MC <- simpleMC$mc

    Gamma.n.sq.h <- mcmcse::mcse.multi(MC, method = "bartlett", r = 1, size = "sqroot")$cov

    T <- N*t(colMeans(MC)-post.beta) %*% solve(Gamma.n.sq.h) %*% (colMeans(MC)-post.beta)

    in.or.out[j] <- as.numeric(T >= qchisq(alpha/2, d) & T <= qchisq(1-alpha/2, d))
  }
  return (mean(in.or.out))
}

## CALCULATION OF COVERAGE PROBABILITY FOR ADAPTIVE MC
checkcoverage.adap <- function(alpha = 0.05, iterations = 5e2, N = 1e3, d = 4){

  in.or.out <- numeric(iterations)

  for (j in 1:iterations){

    if(j %% 1e2 == 0) cat("Iteration-->", j/1e1, "\n")

    adapMC <- run.adaprw(N = N, d = d)
    MC <- adapMC$mc

    Gamma.n.sq.h <- mcmcse::mcse.multi(MC, method = "bartlett", r = 1, size = "sqroot")$cov

    T <- N*t(colMeans(MC)-post.beta) %*% solve(Gamma.n.sq.h) %*% (colMeans(MC)-post.beta)

    in.or.out[j] <- as.numeric(T >= qchisq(alpha/2, d) & T <= qchisq(1-alpha/2, d))
  }
  return (mean(in.or.out))
}
# checkcoverage.adap()

set.seed(8) #  SETTING SEED FOR REPRODUCIBILITY OF RESULTS

coverage <- matrix(0, nrow = 1, ncol = 3)
coverage[,1] <- checkcoverage(N = 1e3)
coverage
coverage[,2] <- checkcoverage(N = 1e4)
coverage
coverage[,3] <- checkcoverage(N = 1e5)
coverage

colnames(coverage) <- c("1e3 iterations","1e4 iterations","1e5 iterations")
rownames(coverage) <- "coverage probability"
coverage

# coverage
#                       1e3 iterations 1e4 iterations 1e5 iterations
# coverage probability          0.844           0.928          0.666

coverage.adap <- matrix(0, nrow = 1 , ncol = 3)
coverage.adap[,1] <- checkcoverage.adap(N = 1e3)
coverage.adap
coverage.adap[,2] <- checkcoverage.adap(N = 1e4)
coverage.adap
coverage.adap[,3] <- checkcoverage.adap(N = 1e5)
colnames(coverage.adap) <- c("1e3 iterations", "1e4 iterations", "1e5 iterations")
rownames(coverage.adap) <- "coverage probability"
coverage.adap
# coverage.adap
#                       1e3 iterations 1e4 iterations 1e5 iterations
# coverage probability          0.746           0.916         0.676
