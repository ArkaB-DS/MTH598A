beta <- 0.5

pi_x<-function(x){
  ifelse((x >= -beta-1 & x <=-beta) | (x >= beta & x <= beta+1),0.5,0)
}
