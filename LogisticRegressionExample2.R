rm(list=ls())
data(logit,package="mcmc")

Y <- logit$y
X <- data.matrix(logit[,2:5])
N=1e5
d=4
s=20
f <- function(beta) prod(dbinom(Y,size=1,prob=1/(1+exp(-X%*%beta))))*prod(dnorm(beta,0,s))
#w1 <- function(x) ifelse(x>=0 & x<=1,1-abs(x),0)
#w2 <- function(x) ifelse(x>=0 & x<=1,1-x^2,0)

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
true.nada=0.4531783

# K = 1:200
# n=3e4
# cov.prob1 = numeric(K)
# cov.prob1 = numeric(K)
# half.width1 = numeric(K)
# half.width2 = numeric(K)
#
# gamma_nk <- function(beta,n,k) {
#   sum((beta[1:(n-k),3] - mean(beta[,3]))*(beta[(k+1):n,3] - mean(beta[,3])))
# }
#
# Gamma_nh1 = function(c,beta) sum(w1(k/c)*gamma_nk(beta,n,k))
#
#
# for (i in K ){
#  beta = matrix(0,nrow=n,ncol=d)
#  for(i in 2:n) {
#   if(i%%1e3==0) print(i/1e3)
#   y =rnorm(d,0,0.48) + beta[i-1,]
#   alpha = min(1, f(y)/f(beta[i-1,]))
#   if(log(runif(1))<log(alpha)){
#     beta[i,] = y
#    }else{
#     beta[i,] = beta[i-1,]
#    }
#   }
#  # k = seq(-n+1,n-1)
#  # delta = seq(0,1,length=10)
#  # cn = n^delta
#  # half.width1 =
#  #
# }

## ADAPTIVE MCMC
rm(list=ls())

data(logit,package="mcmc")
Y <- logit$y
X <- data.matrix(logit[,2:5])
N=1e5
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
true.ada=0.4433487
