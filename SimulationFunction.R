SimulateData<- function(n,k1,k2,seedval=1023){
  # n is the sample size of the simulation
  # k1 is the number of continous variables
  # k2 is the number of categorical variables
  # note k2 can only be at most 26 categories
  # -- I can update this later
  library(Matrix)
  set.seed(seedval)
  vars <- matrix(rnorm(n*k1),n,k1)
  cats <- sample(letters[1:k2],n,replace=T)
  res <- rnorm(n)
  tmp <- data.frame(vars,cats)
  xs <- sparse.model.matrix(~.-1,data=tmp)
  betas <- rnorm(k1+k2)
  yhat <- as.numeric(xs%*%betas)
  actual <- yhat + res
  return(list(Actual=actual,
              Predicted=yhat,
              Xs=xs))  
}
sim <- SimulateData(1e4,20,10,1023)
cor(sim$Actual,sim$Predicted)**2
