#############################################################
## Random-Walk Metropolis-Hastings using Gaussian proposal ##
## to sample from N(0,1)                                   ##
#############################################################

## simple MCMC
run.rw <- function(sample.size=1e3,s=3.1){

  mc <- numeric(sample.size)
  acc.prob <- numeric(sample.size)
  z <- rnorm(sample.size)
  mc[1] <- z[1]

  for (i in 2:sample.size){

    y <- mc[i-1] + s* z[i]

    alpha <- min(1,dnorm(y)/dnorm(mc[i-1]))

    if(runif(1)<alpha) {

      mc[i] <- y
      acc.prob[i] <- 1

      }else{

        mc[i] <- mc[i-1]
    }
  }

  return(list(mc=mc,acc.prob=acc.prob))

}

# adaptive MCMC
run.adaprw <- function(sample.size=1e3,s=3.1){

  mc <- numeric(sample.size)
  acc.prob <- numeric(sample.size)
  z <- rnorm(sample.size)
  mc[1] <- z[1]

  for (i in 2:sample.size){

    y <- mc[i-1] + s* z[i]

    alpha <- min(1,dnorm(y)/dnorm(mc[i-1]))

    if(runif(1)<alpha) {

      mc[i] <- y
      acc.prob[i] <- 1

      }else{

        mc[i] <- mc[i-1]
    }

    s = exp(log(s) + (1/i)*(alpha-0.44))

  }

  return(list(mc=mc,acc.prob=acc.prob))

}
# MC <- run.adaprw()
# plot(cumsum(MC$acc.prob)/(1:length(MC$acc.prob)),type="l",col="red",
#      ylab="Acceptance Probability",xlab="Iterations")
# abline(h=0.44,lty=2,col="blue")

## checking coverage probability for AMCMC
checkcoverage.adap<-function(alpha=0.05,iterations=5e2,sample.size=1e3,prop.sd=3){

  in.or.out <- numeric(iterations)

  for (j in 1:iterations){

    if(j%%1e2==0) print(j/1e2)

    adapMC <- run.adaprw(sample.size=sample.size,s=prop.sd)
    MC <- adapMC$mc

    Gamma.n.sq.h <- as.numeric(mcmcse::mcse.multi(MC, method = "bartlett", r = 1, size = "sqroot")$cov)

    in.or.out[j] <- as.numeric( (mean(MC)-qnorm(alpha/2,lower.tail=FALSE)*sqrt(Gamma.n.sq.h/length(MC))) < 0 &
                           (mean(MC)+qnorm(alpha/2,lower.tail=FALSE)*sqrt(Gamma.n.sq.h/length(MC))) > 0 )

  }

  return (mean(in.or.out))

}

## checking coverage probability for simple MCMC
checkcoverage<-function(alpha=0.05,iterations=5e2,sample.size=1e3,prop.sd=3.1){

  in.or.out <- numeric(iterations)

  for (j in 1:iterations){

    if(j%%1e2==0) print(j/1e2)

    simpleMC <- run.rw(sample.size=sample.size,s=prop.sd)
    MC <- simpleMC$mc

    Gamma.n.sq.h <- as.numeric(mcmcse::mcse.multi(MC, method = "bartlett", r = 1, size = "sqroot")$cov)

    in.or.out[j] <- as.numeric( (mean(MC)-qnorm(alpha/2,lower.tail=FALSE)*sqrt(Gamma.n.sq.h/length(MC))) < 0 &
                                  (mean(MC)+qnorm(alpha/2,lower.tail=FALSE)*sqrt(Gamma.n.sq.h/length(MC))) > 0 )
  }

  return (mean(in.or.out))

}

set.seed(8)

coverage.adap <- matrix(0,nrow=1,ncol=3)
coverage.adap[1] <- checkcoverage.adap(sample.size = 1e3)
coverage.adap[2] <- checkcoverage.adap(sample.size = 1e4)
coverage.adap[3] <- checkcoverage.adap(sample.size = 1e5)
colnames(coverage.adap) <- c("1e3 iterations","1e4 iterations","1e5 iterations")
rownames(coverage.adap) <- "coverage probability"
coverage.adap
#                      1e3 iterations 1e4 iterations 1e5 iterations
# coverage probability          0.922           0.94         0.94

coverage <- matrix(0,nrow=1,ncol=3)
coverage[1] <- checkcoverage(sample.size = 1e3)
coverage[2] <- checkcoverage(sample.size = 1e4)
coverage[3] <- checkcoverage(sample.size = 1e5)
colnames(coverage) <- c("1e3 iterations","1e4 iterations","1e5 iterations")
rownames(coverage) <- "coverage probability"
coverage
#                      1e3 iterations 1e4 iterations 1e5 iterations
# coverage probability           0.918        0.95          0.95
