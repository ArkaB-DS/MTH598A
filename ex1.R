## Example 1 from Yves F. Atchade's "ADAPTIVE MARKOV CHAIN MONTE CARLO CONFIDENCE INTERVALS"

# define beta
beta <- 0.75

# define target distribution
pi_x <- function(x){
  return (ifelse( (x >= -beta-1 & x <= -beta) | (x >= beta & x <= beta+1), 0.5, 0))
}

# define proposal distribution
q_xy <- function(x,y,t){
   return(ifelse( (y >= x-t) & (y <= x + t), 1/(2*t), 0))
}

# define w
# w<-function(x){
#   return(ifelse(abs(x)<1,1-abs(x),0))
# }

# define g_nl
# g<-function(X,n,l){
#   j <- 1:(n-l)
#   return(1/n*sum((X[j]-mean(X[1:n]))*(X[j+l]-mean(X[1:n]))))
# }

sample.size <- 5e5
# create the markov chain
x <- numeric(sample.size)
# create theta_n's
theta <- numeric(sample.size)
# create lag-window estimator
#lag.est <- numeric(250)
# the acceptance probability
alpha <- numeric(sample.size)
# setting up a and A required, from D. Hastie's PhD thesis (Eq 4.1 and Eq 4.2)
a <- 2*beta+1
A <- a + 1

#starting values
x[1] <- .75
theta[1] <- 1.56

for (i in 2:sample.size){
  if(i%%1e5==0) print(i/1e5)
  y <- runif(n=1, min=x[i-1]-theta[i-1], max=x[i-1]+theta[i-1])
  acc.prob <- min(1, pi_x(y)*q_xy(y,x[i-1],theta[i-1])/(pi_x(x[i-1])*q_xy(x[i-1],y,theta[i-1])))
  if(runif(1) < acc.prob){
    x[i] <- y
    alpha[i] <- 1
    theta[i] <- max(a, min(theta[i-1] + (theta[1]/(i-1))*(1-0.30), A))

  }else{
    x[i] <- x[i-1]
    alpha[i] <- 0
    theta[i] <- max(a, min(theta[i-1] + (theta[1]/(i-1))*(0-0.30), A))
  }
  # if(i<=250)
  # k <- seq(-i+1, i-1)
  # lag.est[i] <- sum(sapply(k/i,w)*sapply(k, function(y) g(x,i,abs(y))))
}
#tictoc::toc()
#par(mfcol=c(1,2))
plot(theta,type="l",ylab=expression(theta[n]),xlab ="iterations",col="blue")
# plot(lag.est[1:250],type="l",xlab="iterations",ylab="Lag-windows Est.")

## self tests
## see if alpha is converging to 0.23
#plot(cumsum(alpha)/(1:sample.size),type="l",col="red")
## see if alpha is close to 0.23
#mean(alpha)
## see if the distribution of pi is as it should look
#hist(x)
