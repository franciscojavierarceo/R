# Quick way to simulate data
require(gains)
SimulateData<- function(n,k1,k2,seedval=1023,resmult=1){
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
  actual <- yhat + res*resmult
  prob_a <- 1/(1+exp(-actual))
  prob_p <- 1/(1+exp(-yhat))
  classp <- ifelse(prob_p>=0.5,1,0)
  classa <- ifelse(prob_a>=0.5,1,0)
  return(list(Actual=classa,
              ActualProb=prob_a,
              ActualLog = actual,
              PredClass=classp,
              PredProb=prob_p,
              PredClass=classp,
              PredLog=yhat,
              Xs=xs))  
}
sim <- SimulateData(1e4,20,10,1023,2)
plot(sim$PredProb,sim$ActualProb,col='blue',main='Scatter Plot',ylab='Actual',xlab='Predicted')
plot(density(sim$PredProb-sim$ActualProb),col='red',main='Error')
print(paste('The accuracy is',sum(diag(table(sim$PredClass,sim$Actual) / length(sim$Actual)))))
print(paste('The error is',1-sum(diag(table(sim$PredClass,sim$Actual) / length(sim$Actual)))))
print(table(-sim$PredClass,-sim$Actual) / length(sim$Actual))
