run.adaprw <- function(sample.size=1e3,s=3){
  mc <- numeric(sample.size)
  #lag.est <- numeric(sample.size)
  acc.prob <- numeric(sample.size)
  z <- rnorm(sample.size)
  mc[1] <- z[1]
  for (i in 2:sample.size){
    y <- mc[i-1] + s* z[i]
    alpha <- min(1,dnorm(y)/dnorm(mc[i-1],sd=s))
    if(runif(1)<alpha) {
        mc[i] <- y
        acc.prob[i] <- 1
    }else{
        mc[i] <- mc[i-1]
    }
    s = exp(log(s) + (1/i)*(alpha-0.44))
    #lag.est[i] = as.numeric(mcmcse::mcse.multi(mc[1:i], method = "bartlett", r = 1, size = "sqroot")$cov)
  }
  #return(list(mc=mc,acc.prob=acc.prob,lag.est=lag.est))
  return(list(mc=mc,acc.prob=acc.prob))
}
# MC <- run.adaprw()
# plot(cumsum(MC$acc.prob)/(1:length(MC$acc.prob)),type="l",col="red",
#      ylab="Acceptance Probability",xlab="Iterations")
# abline(h=0.44,lty=2,col="blue")

checkcoverage<-function(alpha=0.05,iterations=5e2,sample.size=1e3,prop.sd=3){
  in.or.out <- numeric(iterations)
  for (j in 1:iterations){
    if(j%%1e2==0) print(j/1e2)
    adapMC <- run.adaprw(sample.size=sample.size,s=prop.sd)
    MC <- adapMC$mc
    Gamma.n.sq.h <- as.numeric(mcmcse::mcse.multi(MC, method = "bartlett", r = 1, size = "sqroot")$cov)
    # Gamma.n.sq.h <- tail(adapMC$lag.est,1)
    in.or.out[j] <- as.numeric( (mean(MC)-qnorm(alpha/2,lower.tail=FALSE)*sqrt(Gamma.n.sq.h/length(MC))) < 0 &
                           (mean(MC)+qnorm(alpha/2,lower.tail=FALSE)*sqrt(Gamma.n.sq.h/length(MC))) > 0 )
  }
  return (mean(in.or.out))
}
display.coverage <- matrix(0,nrow=1,ncol=3)
display.coverage[1] <- checkcoverage(sample.size = 1e3)
display.coverage[2] <- checkcoverage(sample.size = 1e4)
display.coverage[3] <- checkcoverage(sample.size = 1e5)
colnames(display.coverage) <- c("1e4 iterations","1e5 iterations","1e6 iterations")
rownames(display.coverage) <- "coverage probability"
display.coverage
#                      1e4 iterations 1e5 iterations 1e6 iterations
# coverage probability          0.902           0.94          0.948
