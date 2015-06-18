#========================================================================================================================
# Author: Francisco Javier Arceo
# Contact: francisco.arceo@cba.com.au
# Date Last Update: 06/18/2015
# Purpose: Below is a pathological example of why SVD should not be used for supervised learning
#========================================================================================================================
n <- 1e3 ; k <- 10      # These are the parameters of the simulation
xs <- matrix(rnorm(n*k),nrow=n,ncol=k)
# This function creates the new matrix that is represented in a lower dimension, except it maintains the shape for 
# ease of creating the epsilon matrix
createMatrix <- function(xs,k=5){ 
  #  X = U D V'
  s <- svd(xs) ; l <- dim(xs)[2]
  s$d[(k+1):l] <- 0
  xsn <- s$u %*% diag(s$d) %*% t(s$v)
  return(xsn)
}
xsn <- createMatrix(xs,k=4)
epsilon <- xsn-xs
betas <- rnorm(10)
prd <- epsilon%*%betas
summary(lm(prd~xsn))
summary(lm(prd~xs))
#========================================================================================================================
# End 
#========================================================================================================================
