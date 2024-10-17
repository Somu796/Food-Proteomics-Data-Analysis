paretoScaling <- function(data){
  x<- data
  x.centered <- apply(x, 2, function(x) x - mean(x))
  # Then we perform scaling on the mean-centered matrix
  x.sc <- apply(x.centered, 2, function(x) x/sqrt(sd(x)))
  
  return(as.data.frame(x.sc))
}