rm(list=ls())
data = kmed::heart
data_X = data[,1:13]
data_X[,2] = as.factor(as.numeric(data_X[,2]))
data_X[,6] = as.factor(as.numeric(data_X[,6]))
data_X[,9] = as.factor(as.numeric(data_X[,9]))
data_X = data.matrix(data_X)
#data_X[,c(4,5,8)] = data_X[,c(4,5,8)]/1e3
X = data_X
Y = data.matrix(ifelse(data[,14]==0,0,1))

## NON ADAPTIVE MCMC
# plain RWM with Gaussian(0,sigma^2*I_d)

N = 1e6
d=13
s=20

f <- function(beta) prod(dbinom(Y,size=1,prob=1/(1+exp(-X%*%beta))))*prod(dnorm(beta,0,s))

# markov chain
beta = matrix(0,nrow=N,ncol=d)
acc.prob = numeric(N)
for(i in 2:N){
  if(i%%1e3==0) print(i/1e3)
  y = rnorm(d,0,0.0022) + beta[i-1,]
  alpha = min(1, f(y)/f(beta[i-1,]))
  if(runif(1)<alpha){
    beta[i,] = y
    acc.prob[i] = 1
  }else{
    beta[i,] = beta[i-1,]
  }
}
mean(acc.prob)
plot(cumsum(acc.prob)/(1:N),col="red",type="l",xlab="",ylab="")
abline(h=0.23,col="blue",lty=2)
# > colMeans(beta)
# [1] -0.027515432  0.271487958  0.477965324  0.012726262  0.001501614  0.008701202  0.208373899 -0.034807646  0.424081323  0.462574606 -0.025540309
# [12]  0.560421448  0.436959322


# ADAPTIVE MCMC
rm(list=ls())
data = kmed::heart
data_X = data[,1:13]
data_X[,2] = as.factor(as.numeric(data_X[,2]))
data_X[,6] = as.factor(as.numeric(data_X[,6]))
data_X[,9] = as.factor(as.numeric(data_X[,9]))
data_X = data.matrix(data_X)
#data_X[,c(4,5,8)] = data_X[,c(4,5,8)]/1e3
X = data_X
Y = data.matrix(ifelse(data[,14]==0,0,1))

N = 1e6
d = 13
s = 20

beta = matrix(0, nrow=N,ncol=d)
beta[1,] = rnorm( d , 0 , 0.1 / sqrt(d) )
Sigma = 2.38^2/d*diag(d)

acc.prob = numeric(N)

# c = -2.3
# Sigma = exp(c)*diag(d)
#f <- function(beta) sum(dbinom(Y,size=1,prob=1/(1+exp(-X%*%beta)),log=TRUE))+sum(dnorm(beta,0,s,log=TRUE))

f <- function(beta) prod(dbinom(Y,size=1,prob=1/(1+exp(-X%*%beta))))*prod(dnorm(beta,0,s))

for(i in 2:N){
  if( i %% 1e3 == 0) {
    print( i / 1e3)
    cat("AP= ",mean(acc.prob[1:i]),"\n")
  }
  if( i < 2*d ) {
    y = rnorm( d , beta[i-1,] , 0.1 / sqrt(d) )
  }else {
    b = rbinom(1,size=1,prob = 0.05)
    Sigma = cov(beta[1:(i-1),])
    y = b * rnorm(d,beta[i-1,],0.1/sqrt(d)) + (1 - b) * as.vector(mvtnorm::rmvnorm(1,mean = beta[i-1,],sigma = 2.38^2*Sigma/d))
  }
  alpha = f(y) / f(beta[i-1,])

  if(runif(1) < alpha){
    beta[i,] = y
    acc.prob[i] = 1
  }else{
    beta[i,] = beta[i-1,]
  }
 # Sigma = cov(beta[1:i,]) + diag(1e-3,d)
  # if ( (c + 1 / i * (alpha-0.23)) >= -5 & (c + 1 / i * (alpha-0.23)) <= 5 ){
  #   c = c + 1 / i * (alpha-0.23)
  # }
  # Sigma=exp(c)*Sigma

}
mean(acc.prob)
plot(cumsum(acc.prob)/(1:N),col="red",type="l",xlab="",ylab="")
abline(h=0.23,col="blue",lty=2)
# colMeans(beta)
# [1] -0.04753335  0.55069270  0.46754834  0.01575564  0.00168761 -0.40223771  0.27813489 -0.03983079  0.57816903
# [10]  0.38083418  0.11345734  1.23574338  0.75533312
