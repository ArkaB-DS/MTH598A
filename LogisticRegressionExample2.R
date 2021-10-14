rm(list=ls())
data(logit,package="mcmc")

Y <- logit$y
X <- data.matrix(logit[,2:5])
N=1e5
d=4
s=20
f <- function(beta) prod(dbinom(Y,size=1,prob=1/(1+exp(-X%*%beta))))*prod(dnorm(beta,0,s))

# markov chain
beta = matrix(0,nrow=N,ncol=d)
acc.prob = numeric(N)

## NON-ADAPTIVE MCMC
for(i in 2:N){
  if(i%%1e3==0) print(i/1e3)
  y =rnorm(d,0,0.48) + beta[i-1,]
  alpha = min(1, f(y)/f(beta[i-1,]))
  if(log(runif(1))<log(alpha)){
    beta[i,] = y
    acc.prob[i] = 1
  }else{
    beta[i,] = beta[i-1,]
  }
}
mean(acc.prob)
plot(cumsum(acc.prob)/(1:N),col="red",type="l",xlab="",ylab="")
abline(h=0.23,col="blue",lty=2)

#posterior mean of beta(3)
true.nada = mean(beta[,3])

## ADAPTIVE MCMC
rm(list=ls())

data(logit,package="mcmc")
Y <- logit$y
X <- data.matrix(logit[,2:5])
N=1e6
d=4
s=20

beta = matrix(0,nrow=N,ncol=d)
beta[1,]<-rep(1,d)
acc.prob = numeric(N)

c = 0
Sigma=exp(c)*diag(d)
mu = rep(0,d)

f <- function(beta) prod(dbinom(Y,size=1,prob=1/(1+exp(-X%*%beta))))*prod(dnorm(beta,0,s))

for(i in 2:N){
  if(i%%1e3==0) print(i/1e3)
  y =   as.vector(mvtnorm::rmvnorm(1,mean=beta[i-1,],sigma=Sigma))
  alpha = min(1, f(y)/f(beta[i-1,]))
  if(runif(1)<alpha){
    beta[i,] = y
    acc.prob[i] = 1
  }else{
    beta[i,] = beta[i-1,]
  }
  if(isSymmetric(Sigma + 1/i*(tcrossprod(beta[i,]-mu)-Sigma)) &
   min(eigen(Sigma + 1/i*(tcrossprod(beta[i,]-mu)-Sigma))$values)>0){
    Sigma = Sigma + 1/i*(tcrossprod(beta[i,]-mu)-Sigma)
   }
  mu = mu + 1/i*(beta[i,]-mu)
  if ((c + 1/i*(alpha-0.23))>=-10 & (c + 1/i*(alpha-0.23))<=10){
   c = c + 1/i*(alpha-0.23)
  }
  Sigma = exp(c)*Sigma
}
mean(acc.prob)
plot(cumsum(acc.prob)/(1:N),col="red",type="l",xlab="",ylab="")
abline(h=0.23,col="blue",lty=2)

#posterior mean of beta(3)
true.ada= .4571507


run.adaprw <- function(sample.size=1e3,d=4){
  beta <- matrix(0, nrow=sample.size, ncol=d)
  beta[1,] <- rep(1,d)
  acc.prob <- numeric(sample.size)

  c = 0
  Sigma=exp(c)*diag(d)
  mu = rep(0,d)

  for(i in 2:sample.size){
    # if(i%%1e2 == 0) print(i/1e2)
    y = as.vector(mvtnorm::rmvnorm(1, mean = beta[i-1,], sigma = Sigma))
    alpha = min(1, f(y)/f(beta[i-1,]))
    if( runif(1) < alpha ){
      beta[i,] = y
      acc.prob[i] = 1
    }else{
      beta[i,] = beta[i-1,]
    }
    if(isSymmetric(Sigma + 1/i*(tcrossprod(beta[i,]-mu)-Sigma)) &
       min(eigen(Sigma + 1/i*(tcrossprod(beta[i,]-mu)-Sigma))$values)>0){
      Sigma = Sigma + 1/i*(tcrossprod(beta[i,]-mu)-Sigma)
    }
    mu = mu + 1/i*(beta[i,]-mu)
    if ((c + 1/i*(alpha-0.23))>=-10 & (c + 1/i*(alpha-0.23))<=10){
      c = c + 1/i*(alpha-0.23)
    }
    Sigma = exp(c)*Sigma
  }
  return(list(mc=beta,acc.prob=acc.prob))
}

checkcoverage.adap<-function(alpha=0.05,iterations=5e2,sample.size=1e3,d=4){
  in.or.out <- numeric(iterations)
  for (j in 1:iterations){
    if(j%%1e1==0) print(j/1e1)
    adapMC <- run.adaprw(sample.size=sample.size,d=d)
    MC <- adapMC$mc
    Gamma.n.sq.h <- as.numeric(mcmcse::mcse.multi(MC[,3], method = "bartlett", r = 1, size = "sqroot")$cov)
    # Gamma.n.sq.h <- tail(adapMC$lag.est,1)
    in.or.out[j] <- as.numeric( (mean(MC[,3])-qnorm(alpha/2,lower.tail=FALSE)*sqrt(Gamma.n.sq.h/length(MC[,3]))) < .4571507 &
                                  (mean(MC[,3])+qnorm(alpha/2,lower.tail=FALSE)*sqrt(Gamma.n.sq.h/length(MC[,3]))) > .4571507 )
  }
  return (mean(in.or.out))
}
coverage.adap <- matrix(0,nrow=1,ncol=3)
coverage.adap[1] <- checkcoverage.adap(sample.size = 1e3)
coverage.adap[2] <- checkcoverage.adap(sample.size = 1e4)
coverage.adap[3] <- checkcoverage.adap(sample.size = 1e5)
colnames(coverage.adap) <- c("1e4 iterations","1e5 iterations","1e6 iterations")
rownames(coverage.adap) <- "coverage probability"
coverage.adap














