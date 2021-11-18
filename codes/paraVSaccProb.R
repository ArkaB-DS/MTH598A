## Example 1 from Yves F. Atchade's "ADAPTIVE MARKOV CHAIN MONTE CARLO CONFIDENCE INTERVALS"

# define beta
beta <- 0.75

# define target distribution
pi_x <- function(x){
  return ( ifelse( (x >= -beta - 1 & x <= -beta) | (x >= beta & x <= beta + 1), 0.5, 0) )
}

# define proposal distribution
q_xy <- function(x, y, t){
  return( ifelse( (y >= x - t) & (y <= x + t), 1 / (2 * t), 0) )
}

sample.size <- 1e5
x <- numeric(sample.size)
theta <- seq(1.5, 6, length = 1e3)
alpha <- numeric(sample.size)
a<-numeric(1e3)
x[1] <- .75
for ( j in 1 : 1e3 ){
for (i in 2 : sample.size){
  if(i %% 1e4 == 0) cat(j, "/", i/1e4, "\n")
  y <- runif(n = 1, min = x[i-1] - theta[j], max = x[i-1] + theta[j] )
  acc.prob <- min(1, pi_x(y) * q_xy(y,x[i-1], theta[j]) / (pi_x(x[i-1]) * q_xy(x[i-1], y, theta[j])))
  if(runif(1) < acc.prob){
    x[i] <- y
    alpha[i] <- 1
  }else{
    x[i] <- x[i-1]
    alpha[i] <- 0
  }
}
a[j] <- mean(alpha)
}
plot(x=theta,y=a,type="l",col="blue",xlab=expression(theta),ylab="Acceptance Probability")
abline(h=0.30,col = "red", lty = 2)
