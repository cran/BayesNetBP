
SymmetricKLD <- function(post1, post2, discrete) {
  
  if (!discrete) {
    return(SymKLD.continuous(post1, post2))
  } else {
    return(SymKLD.discrete(post1, post2))
  }
  
}

SymKLD.continuous <- function(post1, post2) {
  step <- 0.01
  x <- seq(-20,20,by=step)
  
  y1 <- y2 <- rep(0, length(x))
  for (i in 1:nrow(post1)){
    y1 <- y1 + post1[i,1]*dnorm(x, 
                                mean=post1[i,2], 
                                sd=sqrt(post1[i,3]))
  }
  
  for (i in 1:nrow(post2)){
    y2 <- y2 + post2[i,1]*dnorm(x, 
                                mean=post2[i,2], 
                                sd=sqrt(post2[i,3]))
  }
  
  m1 <- sum(step*x*y1)
  m2 <- sum(step*x*y2)
  
  kld1 <- sum(step*y1*log(y1/y2))
  kld2 <- sum(step*y2*log(y2/y1))
  return(0.5*(kld1+kld2)*sign(m2-m1))
}


SymKLD.discrete <- function(p1, p2) {
  kld1 <- sum(p1*log(p1/p2))
  kld2 <- sum(p2*log(p2/p1))
  return(0.5*(kld1+kld2))
}