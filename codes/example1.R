##########################################################################################
## Example 1 from YF Atchade's "ADAPTIVE MARKOV CHAIN MONTE CARLO CONFIDENCE INTERVALS" ##
##########################################################################################

# DEFINE beta
beta <- 0.75

# DEFINE TARGET DISTRIBUTION
pi_x <- function(x){
  return (ifelse( (x >= -beta-1 & x <= -beta) | (x >= beta & x <= beta+1), 0.5, 0))
}

# DEFINE PROPOSAL DISTRIBTUION
q_xy <- function(x,y,t){
  return(ifelse( (y >= x-t) & (y <= x + t), 1/(2*t), 0))
}

sample.size <- 1e6

chains<-function(shift=0){

  # create the markov chain
  x <- numeric(sample.size)

  # create theta_n's
  theta <- numeric(sample.size)

  # the acceptance probability
  alpha <- numeric(sample.size)

  # setting up a and A required, from D. Hastie's PhD thesis (Eq 4.1 and Eq 4.2)
  a <- 2*beta+shift
  A <- a+1

  #starting values
  x[1] <- .75
  theta[1] <- (a+A)/2

  for (i in 2:sample.size){

    if(i%%1e5==0) print(i/1e5)

    y <- runif(n=1, min=x[i-1]-theta[i-1], max=x[i-1]+theta[i-1])

    acc.prob <- min(1, pi_x(y)*q_xy(y,x[i-1],theta[i-1])/(pi_x(x[i-1])*q_xy(x[i-1],y,theta[i-1])))

    if(runif(1) < acc.prob){

      x[i] <- y
      alpha[i] <- 1
      theta[i] <- max(a, min(theta[i-1] + (10/(i-1))*(1-0.3), A))

    }else{

      x[i] <- x[i-1]
      alpha[i] <- 0
      theta[i] <- max(a, min(theta[i-1] + (10/(i-1))*(0-0.3), A))

      }
  }

  return(list(alpha=alpha,theta=theta,x=x))

}

set.seed(8)

par(mfcol=c(2,2))

for(i in c(0,1.5)){

  chain <-chains(i)
  # plot of theta_n
  plot(chain$theta,type="l",ylab=expression(theta[n]),xlab ="iterations",col="blue")
  # plot of running acceptance probability
  plot(cumsum(chain$alpha)/(1:sample.size),type="l",col="blue",
       ylab=expression(alpha[n]),xlab="iterations")
  abline(h=0.30,lty=2,col="red")
}

set.seed(8)

par(mfcol=c(1,2))

for(i in c(0,1.5)){

  chain <-chains(i)

  x = NULL

  for (j in seq(5e3,1e6,5e3)){

    if (j%%5e2==0) cat("j=",j/5e2,"\n")

    x = c(x,as.numeric(mcmcse::mcse.multi(chain$x[1:j], method = "bartlett", r = 1, size = "sqroot")$cov))

  }
  # plot of sample path of lag-window estimator
  plot(x,type="l",xlab="iterations",ylab="lag-window estimator",col="blue",xaxt="n")

}
