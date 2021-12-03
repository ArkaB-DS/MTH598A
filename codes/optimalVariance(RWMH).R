### optimal scaling and motivation for adaptive mcmc
set.seed(8)

# target - N(0,1)
# proposal - N(x,h) where h is to be chosen

N <- 1e5  # length of the chain

rwmh <- function(N=1e3,h){
X <- numeric(N) # creating MC
X <- -1 # starting value
p <- numeric(N) # for calculating acc. prob.
# implementing RWMH
for(i in 2:N)
{
 Y <- rnorm(1,mean=X[i-1],sd=sqrt(h))
 if(runif(1)<dnorm(Y)/dnorm(X[i-1])) 
 {
  X[i] <- Y
  p[i] <- 1
  } else {
  X[i] <- X[i-1]
  }
} 
return( list(acc.prob=p, chain=X) )
}

# plots

#layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=TRUE), heights=c(4,4,1))
MC1 = rwmh(h=1000)
MC2 = rwmh(h=0.0001)
MC3 = rwmh(h=6)

par(mfrow=c(3,3))

plot.ts(MC1$chain, col = "blue",
xlab="iterations", ylab = "x")
ts.plot(MC2$chain, col = "blue",
xlab="iterations", ylab = "x")
ts.plot(MC3$chain, col = "blue",
xlab="iterations", ylab = "x")

plot(cumsum(MC1$acc.prob)/1:length(MC1$acc.prob),
type = "l", col = "blue", xlab ="iterations", ylab="Acceptance Probability")
plot(cumsum(MC2$acc.prob)/1:length(MC2$acc.prob),
type = "l", col = "blue", xlab ="iterations", ylab="Acceptance Probability")
plot(cumsum(MC3$acc.prob)/1:length(MC3$acc.prob),
type = "l", col = "blue", xlab ="iterations", ylab="Acceptance Probability")

acf(MC1$chain, col = "blue", main = "", sub = "(a)")
acf(MC2$chain, col = "blue", main = "",  sub = "(b)")
acf(MC3$chain, col = "blue", main = "",  sub = "(c)")

#title("Fig. 2. Simple Metropolis algorithm with (a) too-large variance (left plots), (b) too-small variance (middle) and (c) appropriate variance
#(right). Trace plots (top) and autocorrelation plots (below) are shown for each case.",
#line = -45, outer = TRUE,font.main=3)
