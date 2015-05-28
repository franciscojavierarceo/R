#========================================================================================================================
require(Hmisc)
require(devtools)
require(caret)
# install_github('tqchen/xgboost',subdir='R-package')
require(pROC)
require(verification)
require(xgboost)
require(tm)
require(gbm)
require(foreach)
# require(rgl)
require(mgcv)
require(data.table)
require(glmnet)
require(Matrix)
require(coefplot)
require(stringi)
require(ROCR)
require(gdata)
require(reshape2)
require(ggplot2)
require(gridExtra)
require(scales)
require(doParallel)
require(randomForest)
options(gsubfn.engine = "R")
require(sqldf)
#========================================================================================================================
# Author: Francisco Javier Arceo
# Contact: francisco.arceo@cba.com.au
# Date Last Update: 01/16/2015
# Purpose: The code below contains a variety of different functions for a variety of different things; I've used these
# for general modeling, and for things that I simply got tired of manually coding.
#========================================================================================================================
# Easy things that are annoying to constantly repeat
#========================================================================================================================
# This cleans up data
ReplaceMissingFunc<- function(mydf,start_col=2){
  # You have to specify the starting column
  print(paste('The input data has',ncol(mydf),'columns.'))
  j <- ncol(mydf)
  for(i in start_col:j){
    if(is.factor(mydf[,i])==TRUE){
      mydf[,i] <- replacefactor(mydf[,i])
    }
    else if(is.numeric(mydf[,i])==TRUE){
      if(sum(is.na(mydf[,i]))>0){
        mydf[,paste('missing_',names(mydf)[i],sep='')] <- ifelse(is.na(mydf[,i])==T,1,0)
        mydf[,i] <- replacecont(mydf[,i])        
      }
    }
  }
  print(paste('The output data has',ncol(mydf),'columns.'))
  return(mydf)
}
# Summary of continuous variables in a data frame
# I prefer this because it allows me to refer to these summaries in a dataframe
cvarsum <-function(df){
  n <- dim(df)[1]
  k <- dim(df)[2]
  out <- data.frame(Variable_Name=rep('.',k),
                    Prcnt_0=rep(0,k),
                    Prcnt_25=rep(0,k),
                    Prcnt_75=rep(0,k),
                    Prcnt_100=rep(0,k),
                    Mean=rep(0,k),
                    Median=rep(0,k),
                    Length=rep(0,k),
                    Missing=rep(0,k),
                    N_Distinct=rep(0,k))
  for(i in 1:k){
    out$Prcnt_0[i]   	<- quantile(df[,i],probs=0.00)[[1]]
    out$Prcnt_25[i]  	<- quantile(df[,i],probs=0.25)[[1]]
    out$Prcnt_75[i]		<- quantile(df[,i],probs=0.75)[[1]]
    out$Prcnt_100[i]	<- quantile(df[,i],probs=1.00)[[1]]
    out$Mean[i]			<- mean(df[,i])
    out$Median[i]		<- median(df[,i])
    out$Length[i]		<- length(df[,i])
    if(table(is.na(df[,i]))[[1]]==n){
      out$Missing[i] <- 0  
    } else {
      out$Missing[i]  <- table(is.na(df[,i]))[[2]]  
    }
    out$N_Distinct[i]	<- length(unique(df[,i]))
  }
  out$Variable_Name <- names(df)
  return(out)
}
# Use this to replace missing categorical variables with a level that says "MISSING"
replacefactor <- function(x,replacement_value="MISSING"){
  x <- as.character(x)
  x[x==''] <- replacement_value
  x[is.na(x)==TRUE] <- replacement_value
  return(factor(x))
}
# This does mean replacement for continuous variables
replacecont <- function(x){
  x[is.na(x)] <- mean(x,na.rm=T)
  return(x)
}
num <- function(x){
  x <- as.numeric(x)
  return(x)
}
# Decile function, and can be extended to percentiles using the n parameter
deciles <- function(x,n=10,print=F){
  out <- data.frame(quantile(x,probs=1:n/n))
  out <- data.frame(Decile=row.names(out),round(out,2))
  row.names(out) <- NULL ; names(out)[2] <- 'Value'
  if(print==T){
    print(out,row.names=FALSE)
  }
  return(out)
}
# I like to print the names of a dataset by sort order
snames    <-function(x,mat=F){
  if(mat==FALSE){
    print(sort(names(x)))
  }
  if(mat==TRUE){
    print(matrix(sort(names(x))))
  }
}
# This tells how many distinct values of a variable I have
ndistinct <- function(x){
  print(paste("There are",length(unique(x)),'unique values.'))    
}
# The mode--No idea why this isn't standard
Mode     <-function(x){
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# Another summary function
fsummarize <- function(x){
  print('Deciles:')
  print(quantile(x,probs=1:10/10))
  print(paste('Variance:',var(x)))
  print(paste('Standard Deviation:', sd(x)))
  print(paste('Mean:', mean(x)))
  print(paste('Median:', median(x)))
}
# For consistincy between Python
len <- function(x){
  return(length(x))
}
# Capping function
my_cap<-function(x, cap_l=NULL, cap_h=NULL) {
  if(is.null(cap_l)==TRUE){cap_l <- min(x)}
  if(is.null(cap_h)==TRUE){cap_h <- max(x)}
  x1<-x
  x1[x<cap_l]<-cap_l
  x1[x>cap_h]<-cap_h
  return(x1)
}
# This produces the natural log of a variable and sets 0 as the lower bound
my_ln_cap<-function(x, cap_l=1, cap_h=NULL) {
  if(is.null(cap_h)==TRUE){cap_h <- max(x)}
  x1<-x
  x1[x<cap_l]<-cap_l
  x1[x>cap_h]<-cap_h  
  return(log(x1,exp(1)))
}
# Takes x variable, exports a K level ranking
myrank <- function(x,k=10,factor=F){
  df <- data.frame(x)
  df$indx <- 1:length(x)
  df <- df[order(df$x),]
  df$rank <- 1
  df$cmsm <- cumsum(df$rank)
  df$cmsmrank <- round(df$cmsm*k/length(df$rank)+0.5,0)
  df$cmsmrank[df$cmsmrank<1] <- 1
  df$cmsmrank[df$cmsmrank>k] <- k
  df <- df[order(df$indx),]
  if(factor==T){
    df$cmsmrank <- factor(df$cmsmrank)
  }
  return(df$cmsmrank)
}
# A lag function I found that's very useful
mylag<-function(x,shift_by){
  stopifnot(is.numeric(shift_by))
  stopifnot(is.numeric(x))
  
  if (length(shift_by)>1)
    return(sapply(shift_by,shift, x=x))
  
  out<-NULL
  abs_shift_by=abs(shift_by)
  if (shift_by > 0 )
    out<-c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
  else if (shift_by < 0 )
    out<-c(rep(NA,abs_shift_by), head(x,-abs_shift_by))
  else
    out<-x
  return(out)
}
# Trim Whitespace
trimws <- function (x){
  gsub("^\\s+|\\s+$", "", x)
} 
# Quickly generates an indicator for simple string search
word_search <- function(word,string_vector){
  out <- c(regexpr(toupper(word),toupper(string_vector)))
  out <- ifelse(out>0,1,0)
  return(out)
}
# I like code that's more explicit than "sub", so I made this 
# simple code
word_replace <- function(word,string_vector){
  out <- gsub(string_vector,word,replacement='')
  return(out)
}
# Makes the words properly capitalized
propercase <- function(x) {
  stri_trans_totitle(x)
}
month_stamp <- function(x){
  return(format(Sys.time(), "%Y%m"))
}
day_stamp <- function(x){
  return(format(Sys.time(), "%Y%m%d"))
}
date_stamp <- function(x){
  return(format(Sys.time(), "%Y%m%d_%H%M"))
}
export <- function(mydf,filename){
  outfilename <- paste(date_stamp(),'_',filename,'.csv',sep='')
  write.csv(x=mydf,
            file=outfilename,
            row.names=F)
  tmpout <- paste(outfilename,'has been exported.')
  print(tmpout)
}
ptable <- function(var,sort=T){
  x <- data.frame(table(var))
  x <- x[order(x$Freq,decreasing=sort),]
  return(x)
}
na_table<- function(var){
  return(table(is.na(var)))
}
ReduceCat <- function(df,xvar,dep,minval=10,depminval=NULL){
  tmpdf <- df[,c(xvar,dep)]
  names(tmpdf) <- c('xvar','yvar')
  out <- sqldf('select xvar,sum(yvar) as yvar_tot, sum(1) as cnt
               from tmpdf group by 1')
  out$xvar[out$cnt < minval] <- 'ALL OTHERS'
  if(len(depminval)>0){
   out$xvar[out$yvar_tot < depminval] <- 'ALL OTHERS' 
  } else {
    out$xvar[out$cnt < minval] <- 'ALL OTHERS' 
  }
  out2 <- sqldf("select xvar, sum(yvar_tot) as yvar_tot, sum(cnt) as cnt
                from out group by 1")
  final <- sqldf('select b.xvar as newxvar, b.yvar_tot as ytot,
                 b.cnt from tmpdf a 
                 left join out2 b
                 on a.xvar = b.xvar')
  final$newxvar[is.na(final$newxvar)] <- 'ALL OTHERS'
  final$ytot[final$newxvar=='ALL OTHERS'] <- out2$yvar_tot[out2$xvar=='ALL OTHERS']
  final$cnt[final$newxvar=='ALL OTHERS'] <- out2$cnt[out2$xvar=='ALL OTHERS']
  return(final)
}
# A function to summarize discrete variables
cat_summary <- function(df,x,percent=T){
  df <- data.frame('var'=df[,x])
  out<- sqldf('SELECT var,
              COUNT(*) as Count
              FROM df
              GROUP BY 1
              ORDER BY 2 DESC')
  temp <- out$Count/sum(out$Count)*100
  if(percent==TRUE){
    temp <- sprintf("%.2f", round(temp,2))
    out$Percent <- paste(temp,'%',sep='')  
  } 
  if(percent==FALSE){
    out$Proportion <- out$Count/sum(out$Count)
  }
  names(out)[1] <- x
  return(out)
}
# Tabulates the mean, median, or ratio of dependent
cat_dep_summary <- function(df,x,y,fun='mean',indx=T,plt=FALSE){
  df <- data.frame(df[,c(x,y)])
  names(df) <-c('var','dep')
  out<- sqldf('SELECT var,
              COUNT(*) as Count,
              avg(dep) as Average_Y
              FROM df
              GROUP BY 1
              ORDER BY 2 DESC')
  temp <- out$Count/sum(out$Count)*100
  temp <- sprintf("%.2f", round(temp,2))
  out$Percent <- paste(temp,'%',sep='')  
  if(indx==T){
    if(fun=='mean'){
      out$Average_Y <- round(out$Average_Y/ mean(df$dep),2)
      names(out)[3] <- paste('Average',y)    
    }
    if(fun=='median'){
      out$Average_Y <- round(out$Average_Y/ median(df$dep),2)
      names(out)[3] <- paste('Median',y)    
    }
    names(out)[1] <- x
    names(out)[3] <- paste("Index",y)    
  }
  return(out)
  if(plt==TRUE & indx==FALSE){
    outplt <- with(out,plot(var,Average_Y,col='blue',pch=16,type='b'))
    return(outplt)
  }
}
# This creates a dot plot of the dependent variable by a 
# particular categorical variable
prop_plot <- function(x,y,xlab=NULL,ylab=NULL,main='',flip=FALSE){
  if(is.null(xlab)==TRUE){
    xlab <- 'Category'
  }
  if(is.null(ylab)==TRUE){
    ylab <- 'Dependent Variable'
  }
  df <- data.frame(x,y)
  plt <- ggplot(df, aes(x, y)) + 
    stat_summary(fun.y = mean,colour='#0099FF',geom = 'point',size=5) + 
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", mult = 1)+
    theme_bw()+coord_cartesian(ylim=c(-0.1,1.1))+ 
    scale_y_continuous(labels=percent_format())+
    xlab(xlab)+ylab(ylab)+ggtitle(main)
  if(flip==TRUE){
    plt <- plt+coord_flip()
  }
  return(plt)
}
# This takes variables with  p-values below some specified threshold
# and maps them to an "other" category and outputs the modified 
# categorical variable
cat_signf_cr8tr <- function(df,x,y,signif=0.10){
  tmp_df <- df[,c(x,y)]; names(tmp_df) <- c('x','y')
  if(length(unique(tmp_df$y))==2){
    dist <- 'binomial'
  }
  else {
    dist <- 'guassian'
  }
  xlab <- x
  x <- summary(glm(y~x-1,data=tmp_df,family=dist))
  sdf <- data.frame(x[12])
  names(sdf) <- c('Betas','STD_ERR','Z_VALUE','P_VALUE')
  sdf$Names <- row.names(sdf)
  sdf$P_VALUE <- round(sdf$P_VALUE,4)
  row.names(sdf) <- NULL
  sdf <- sdf[c(5,1,2,3,4)]
  sdf <- sdf[sdf$P_VALUE <= signif,]
  sdf$level <- gsub('x',sdf$Names,rep='')
  lvls <- unique(sdf$level)
  tmp_df$new <- as.character(tmp_df$x)
  tmp_df$new[!tmp_df$x %in% lvls] <- 'OTHER'
  return(tmp_df$new)
}
my_exp1<-function(d1, Xvar1, Xvar2, y, AverageY, filter, cred_k=10, r_k=0){
  # d1 is the data frame
  # Xvar1 is the  first factor variable
  # Xvar2 is the second factor variable
  # y is the actual
  # AverageY is the mean of the dependent variable (training sample)
  # filter is the trianing filter (can be true/false)
  # r_k is the random jitter size
  # cred_k is the credibility adjustment size 
  # To do the one way just create a dummy variable with all ones
  d2<-d1[, c(Xvar1, Xvar2, y, AverageY)]
  names(d2)<-c("f1", "f2", "a", "p")
  d2$filter<-filter
  sum1<-sqldf("select f1, f2, sum(1) as cnt, sum(p) as sump, sum(a) as suma from d2 where filter=1 group by 1,2")
  tmp1<-sqldf("select a.p, b.cnt, b.sump, b.suma from d2 a left join sum1 b on a.f1=b.f1 and a.f2=b.f2")
  tmp1$cnt[is.na(tmp1$cnt)]<-0
  tmp1$avgp<-with(tmp1, sump/cnt)
  tmp1$avgp[is.na(tmp1$avgp)]<-0
  tmp1$suma[is.na(tmp1$suma)]<-0
  tmp1$cnt[filter]<-tmp1$cnt[filter]-1
  tmp1$suma[filter]<-tmp1$suma[filter]-d2$a[filter]
  tmp1$exp_a<-with(tmp1, suma/cnt)
  tmp1$adj_a<-with(tmp1, (suma+p*cred_k)/(cnt+cred_k))
  tmp1$exp_a[is.na(tmp1$exp_a)]<-tmp1$p[is.na(tmp1$exp_a)]
  tmp1$adj_a[is.na(tmp1$adj_a)]<-tmp1$p[is.na(tmp1$adj_a)]
  set.seed(10)
  tmp1$adj_a[filter]<-tmp1$adj_a[filter]*(1+(runif(sum(filter))-0.5)*r_k)
  return(tmp1)
}
# Getting ggplot's color scale
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
myformula <- function(features,y=NULL){
  if(len(y)==0){
    out <- as.formula(paste('',paste(paste(features,collapse='+'),'-1',sep=''),sep='~'))
  }
  if(len(y)>0){
    out <- as.formula(paste(y,paste(paste(features,collapse='+'),'-1',sep=''),sep='~'))
  }
  return(out)
}
# This is to get coefficients
glmnet_betas <- function(glmnet_model){
  bs <- as.matrix(coef(glmnet_model,s='lambda.min'))
  bs <- data.frame(bs) ; names(bs) <- 'Coefficients'
  prcnt0 <- round(len(bs$Coefficients[bs$Coefficients==0])/len(bs$Coefficients),4)
  out <- paste(prcnt0*100,'% of the variables were zero.',sep='')
  Nonzero <- len(bs$Coefficients[bs$Coefficients!=0])
  out2<- paste(Nonzero,'of the variables were non-zero.')
  print(out)
  print(out2)
  return(bs)
}
# This will just nicely automate stuff that's annoying
glmnetout <- function(model,
                      xmatrix,
                      ydep,
                      trainfilter,
                      validfilter,
                      filename='GLMNET',
                      owd=NA,
                      ModelName='My GLMNET Model',
                      frmla=NA,
                      df=NA){
  time <- format(Sys.time(), "%Y%m%d_%H%M") 
  beta_filename <- paste(time,'_',filename,'_Betas_','.csv',sep='')
  pred_filename <- paste(time,'_',filename,'_Preds_','.csv',sep='')
  setwd(owd)
  if(is.na(owd)==TRUE){
    owd <- getwd()
  }
  if(is.na(c(frmla))==FALSE){ # If this is populated, execute
    frmla_filename <- paste(time,'_',filename,'_ModelDetails_','.txt',sep='')
    frmla_out <- sink(frmla_filename)
    cat(paste('#',paste(rep('=',60),collapse=''),sep=''),sep='\n')
    cat('Start of file.',sep='\n')
    cat(paste('#',paste(rep('=',60),collapse=''),sep=''),sep='\n')
    cat('',sep='\n')# Blank Space
    cat('Model Author: Francisco Javier Arceo',frmla_out,sep='\n')
    cat(paste('Model Execution Date:',time),frmla_out,sep='\n')
    cat(paste('Location of Model Output:',owd),sep='\n')
    cat(paste('File Name of Predictions:',pred_filename),frmla_out,sep='\n')
    cat(paste('File Name of Coefficents:',beta_filename),frmla_out,sep='\n')
    cat('',sep='\n')# Blank Space
    cat(paste(ModelName,'Equation:'),frmla_out,sep='\n')
    cat(paste(frmla[1]),frmla_out,sep='\n')
    cat(paste(frmla[2]),frmla_out,sep='\n')
    cat('',sep='\n')# Blank Space
    cat(paste('#',paste(rep('=',60),collapse=''),sep=''),sep='\n')
    cat('End of file.',sep='\n')
    cat(paste('#',paste(rep('=',60),collapse=''),sep=''),sep='\n')
    sink(frmla_out)
  }
  pred <- c(predict(model,newx=xs,type='response',s='lambda.min'))
  outdf <- data.frame(Actual=ydep,Predicted=pred)
  outdf$Training <- trainfilter
  outdf$Validation <- validfilter
  if(is.null(dim(df))==FALSE){  
    outdf <- cbind(outdf,df)
  }
  coefs <- coef(model,s='lambda.min',exact=T)
  betas <- data.frame(Names=rownames(coefs),as.matrix(coefs)); names(betas)[2] <- 'Coefficients'
  betas$ModelName <- ModelName; betas$TimeStamp <- time
  write.csv(betas,beta_filename,row.names=T)
  write.csv(outdf,pred_filename,row.names=F)
}
xgboost_out <- function(model,
                        xmatrix,
                        actual,
                        trainfilter,
                        validfilter,
                        filename='XGBoost',
                        ModelName='My XGBoost Model',
                        frmla=NA,
                        df=NA){
  time <- format(Sys.time(), "%Y%m%d_%H%M") 
  pred_filename <- paste(time,'_',filename,'_Preds','.csv',sep='')
  setwd(owd)
  if(is.na(owd)==TRUE){
    owd <- getwd()
  }
  pred <- predict(model,xs)
  outdf <- data.frame(Actual=actual,Predicted=pred)
  outdf$Training <- trainfilter
  outdf$Validation <- validfilter
  if(is.null(dim(df))==FALSE){  
    outdf <- cbind(outdf,df)
  }
  write.csv(outdf,pred_filename,row.names=F)  
  if(is.na(c(frmla))==FALSE){ # If this is populated, execute
    frmla_filename <- paste(time,'_',filename,'_ModelDetails','.txt',sep='')
    frmla_out <- sink(frmla_filename)
    cat(paste('#',paste(rep('=',60),collapse=''),sep=''),sep='\n')
    cat('Start of file.',sep='\n')
    cat(paste('#',paste(rep('=',60),collapse=''),sep=''),sep='\n')
    cat('',sep='\n')# Blank Space
    cat('Model Author: Francisco Javier Arceo',frmla_out,sep='\n')
    cat(paste('Model Execution Date:',time),frmla_out,sep='\n')
    cat(paste('Location of Model Output:',owd),sep='\n')
    cat(paste('File Name of Predictions:',pred_filename),frmla_out,sep='\n')
    cat('',sep='\n')# Blank Space
    cat(paste(ModelName,'Equation:'),frmla_out,sep='\n')
    cat(paste(frmla[1]),frmla_out,sep='\n')
    cat(paste(frmla[2]),frmla_out,sep='\n')
    cat('',sep='\n')# Blank Space
    cat('The Training AUC is:',auc(actual[trainfilter==1],pred[trainfilter==1]),sep='\n')
    cat('',sep='\n')# Blank Space
    cat('The Validation AUC is:',auc(actual[validfilter==1],pred[validfilter==1]),sep='\n')
    cat('',sep='\n')# Blank Space
    cat(paste('#',paste(rep('=',60),collapse=''),sep=''),sep='\n')
    cat('End of file.',sep='\n')
    cat(paste('#',paste(rep('=',60),collapse=''),sep=''),sep='\n')
    sink(frmla_out)
  }
}
#========================================================================================================================
# Code for validating models
#========================================================================================================================
# This is a psuedo-R-squared. Note, it's not a real R-squared because it's not adjusted for the 
# actual degrees of freedom, but I like this to understand the percent of predicted variance
rsq <- function(pred,actual,printrsq=T){
  rsqt <- summary(lm(actual~pred))$r.squared
  if(printrsq==T){
    print(paste('The Model explains ',round(rsqt,2)*100,'% of the variation.',sep=''))
  }
  if(printrsq==F){
    return(rsqt)
  }
}
# Root Mean Squared Error
rmse <- function(pred,actual){
  return(sqrt(mean((pred-actual)^2)))
}
# Mean Absolute Percent Error
mape <- function(pred,actual){
  mape <- sum(abs(pred-actual))/sum(actual)
  return(mape)
}
# This calculates the gini
mygini<-function(pred, actual) {
  w<-rep(1, length(pred))
  v<-data.frame(o=pred, p=pred, a=actual, w=w)
  v<-v[order(v$o),]
  v$cumm_w<-cumsum(v$w)
  v$cumm_y<-cumsum(v$w*v$a)
  total_w<-sum(v$w)
  gini<-with(v,1-2* sum(cumm_y*w)/(sum(a*w)*total_w))
  return(gini)
}
# True Positive Rate
tpr <- function(pred,actual){
  smry_table <- data.frame(table(pred,actual))
  num <-smry_table$Freq[smry_table$pred==1 & smry_table$actual==1]
  den <- sum(smry_table$Freq[smry_table$actual==1])
  out <- num/den
  if(length(out)==0){
    out <- NA
  }
  return(out)
}
# True Negative Rate
tnr <- function(pred,actual){
  smry_table <- data.frame(table(pred,actual))
  num <-smry_table$Freq[smry_table$pred==0 & smry_table$actual==0]
  den <- sum(smry_table$Freq[smry_table$actual==0])
  out <- num/den
  if(length(out)==0){
    out <- NA
  }
  return(out)
}
# False Positive Rate
fpr <- function(pred,actual){
  smry_table <- data.frame(table(pred,actual))
  num <-smry_table$Freq[smry_table$pred==1 & smry_table$actual==0]
  den <- sum(smry_table$Freq[smry_table$actual==0])
  out <- num/den
  if(length(out)==0){
    out <- NA
  }
  return(out)
}
# False Negative Rate
fnr <- function(pred,actual){
  smry_table <- data.frame(table(pred,actual))
  num <-smry_table$Freq[smry_table$pred==0 & smry_table$actual==1]
  den <- sum(smry_table$Freq[smry_table$actual==1])
  out <- num/den
  if(length(out)==0){
    out <- NA
  }
  return(out)
}
BOW2SparseMatrix <- function(var,sparseval=NULL,tf=FALSE,idf=FALSE){
  # Bag of Words to Sparse Matrix Creator
  var <- as.character(var)
  MyCorpus <- VCorpus(VectorSource(var))
  CorpusMatrix <- DocumentTermMatrix(MyCorpus)
  if(tf==TRUE){
    CorpusMatrix <- DocumentTermMatrix(MyCorpus,control = 
                                         list(weighting = function(x) weightTfIdf(x, normalize = idf)))    
  }
  print(CorpusMatrix)
  if(length(sparseval)>0){
    CorpusMatrix <- removeSparseTerms(CorpusMatrix,sparseval)
    print(paste(1-sparseval,'Sparse Terms Removed from the Document Matrix'))
    print("Processed Document Matrix is:")
    print(CorpusMatrix)    
  }
  words <- colnames(CorpusMatrix)
  out <- Matrix(CorpusMatrix,sparse=T,nrow=nrow(CorpusMatrix))
  colnames(out) <-words
  return(out)
}
# ROC function for one plot
my_roc <- function(pred,actual,thresh=0,add_p=0,colr='red',ttl=''){
  mydf <- data.frame(cont_pred=pred,actual_pred=actual)
  mydf$true_bin <- with(mydf,ifelse(actual_pred>thresh,1,0))
  mydf$sort <-with(mydf,cont_pred + min(cont_pred)+1)
  mydf <- sqldf('select * from mydf order by sort')
  n <- length(mydf$sort)
  mydf$pred_rank <- seq(from=1,to=n,length=n) / n
  index <- seq(from=-.1,to=1.2,length=n) 
  roc_pred <- with(mydf,prediction(pred_rank,true_bin))
  perf <-  performance(roc_pred,'tpr','fpr')
  main_title <- paste(ttl,'\n','ROC Curve \n','AUC =',round(auc(actual,pred),4))
  plot(perf,col=colr,main=main_title,add=add_p);grid(5,5,'gray44')
  lines(index,index,col="black",lty=1)
}
# The function below looks at the fit statistics of various predictions
predcompare <- function(pred_df,actual){
  preds <- data.frame(pred_df)
  out <- data.frame(Predictions=names(preds),
                    Rsq=0,
                    RMSE=0,
                    MAPE=0,
                    GINI=0,
                    NormGINI=0)
  for(j in 1:ncol(preds)){
    out$Rsq[j] <- rsq(preds[,j],actual,F)
    out$RMSE[j] <- rmse(preds[,j],actual)
    out$MAPE[j] <- mape(preds[,j],actual)
    out$GINI[j] <- mygini(preds[,j],actual)
  } 
  out$NormGINI <- out$GINI / mygini(actual,actual)
  return(out)
}
# Classification accuracy
cls_acc <- function(pred,actual){
  out <- sum(abs(pred - actual)) / length(actual)
  out <- round((1- out)*100,2)
  #print(paste("The classification accuracy is ",out,"%.",sep=''))
  return(out)
}
myrocplot <- function(yvar, predprob){
  out <- paste("The AUC is",round(auc(yvar,predprob),4))
  print(out)
  roc.plot(yvar,predprob,main=out)
}
# This function creates a spline summary
# X is the continuous input
# Y is the dependent variable
# K is the number of categories
# I admit that this code isn't terribly elegant, but it works.
mysplinefunc <- function(x,y,k=10,plt=F,out=F,print_ranks=F,regline=T,main='Spline Plot'){
  if(k>length(unique(x))){k <- length(unique(x))}
  df <- data.frame(x,y)
  df$indx <- 1:dim(df)[1]
  df <- df[order(df$x),]
  df$rank <- 1
  df$cmsm <- cumsum(df$rank)
  df$cmsmrank <- round(df$cmsm*k/length(df$rank)+0.5,0)
  df$cmsmrank[df$cmsmrank<1] <- 1
  df$cmsmrank[df$cmsmrank>k] <- k
  df <- df[order(df$indx),]
  if(plt==TRUE | print_ranks==TRUE | plt=='rank'){
    out_df <- sqldf('select  cmsmrank,
                    min(x) as x_min,
                    max(x) as x_max,
                    sum(y) as y_sum,
                    count(y) as count_y
                    from df
                    group by 1')
    out_df$prop_y_by_x <- with(out_df,y_sum/count_y)
    if(print_ranks==TRUE){ print(out_df) }
    if(plt==TRUE){
      plt <- ggplot(out_df,aes(x=x_max,y=prop_y_by_x))+ggtitle(main)+
        geom_line(colour='blue')+theme_bw()+geom_point(colour='blue')
      if(regline==TRUE){
        plt <- plt + geom_smooth(colour='red',se=F,linetype=2)
      }
      return(plt)
    }
    if(plt=='rank'){
      plt <- ggplot(out_df,aes(x=cmsmrank,y=prop_y_by_x))+ggtitle(main)+
        geom_line(colour='blue')+theme_bw()+geom_point(colour='blue')					 
      if(regline==TRUE){
        plt <- plt + geom_smooth(colour='red',se=F,linetype=2)
      }
      return(plt)
    }
  }
  if(out==TRUE){
    return(df$cmsmrank)    
  }
}
# Owen's Lift
	ownslift<-function(pred, orderby,actual, w=NULL, n=10) {
  if (length(w)==0) {
    w<-rep(1.0, length(pred))
  }
  v<-data.frame(o=orderby, p=pred, a=actual, w=w)
  v<-v[order(v$o),]
  #print(head(v, 100))
  v$cumm_w<-cumsum(v$w)
  v$cumm_y<-cumsum(v$w*v$a)
  total_w<-sum(v$w)
  gini_val <- paste("gini=", with(v,1-2* sum(cumm_y*w)/(sum(a*w)*total_w)))
  print(gini_val)
  v$pidx<-round(v$cumm_w*n/total_w+0.5)
  v$pidx[v$pidx>n]<-n
  v$pidx[v$pidx<1]<-1
  
  sum1<-sqldf("select pidx, sum(w) as w, min(o) as min, max(o) as max, sum(o*w)/sum(w) as o, sum(p*w)/sum(w) as p, sum(a*w)/sum(w) as a, sum(a*w) as a_cnt from v group by pidx")
  sum1$err<-with(sum1, a/p)
  print(sum1)
  with(sum1,plot(p, type='b', ylim=c(min(c(p, a)), max(c(p, a))),col='blue',pch=16,
                 main=gini_val,xlab='Decile',ylab='Actual versus Predicted'))
  lines(sum1$a,col='red',type='b',pch=16);grid(5,5,'gray44')
  legend('topleft',c('Predicted','Actual'),col=c('blue','red'),pch=c(16,16),lty=c(1,1))
}

# A prettier version of Owen's original Lift Function
mylift<-function(pred, orderby=NULL, actual, w=NULL, n=10,yscl='dollar',metric=NULL){   	
  if (length(w)==0) {
    w<-rep(1.0, length(pred))
  }
  if (is.null(orderby)==TRUE){
  	orderby <- pred
  }
  pred <- as.numeric(pred)
  v<-data.frame(o=orderby, p=pred, a=actual, w=w)
  v<-v[order(v$o),]
  v$cumm_w<-cumsum(v$w)
  v$cumm_y<-cumsum(v$w*v$a)
  total_w<-sum(v$w)
  ginival <- round(with(v,1-2* sum(cumm_y*w)/(sum(a*w)*total_w)),8)
  gini<-paste("gini=",ginival)
  Normgini<-paste("Normalized gini=", ginival/mygini(v$a,v$a))
  v$pidx<-round(v$cumm_w*n/total_w+0.5)
  v$pidx[v$pidx>n]<-n
  v$pidx[v$pidx<1]<-1
  sum1<-sqldf("select pidx, 
              sum(w) as w, 
              min(o) as min, 
              max(o) as max, 
              sum(o*w)/sum(w) as o, 
              sum(p*w)/sum(w) as p, 
              sum(a*w)/sum(w) as a, 
              sum(a*w) as a_cnt 
              from v group by pidx")
  sum1$err<-with(sum1, a/p)
  sum1$gini <-with(v,1-2* sum(cumm_y*w)/(sum(a*w)*total_w))
  x2 <- melt(sum1[,c('pidx','p','a')],id='pidx')
  names(x2) <- c('Decile','Variable','Value')
  x2$Variable <- as.character(x2$Variable)
  x2$Variable[x2$Variable=='a'] <- 'Actual'
  x2$Variable[x2$Variable=='p'] <- 'Predicted'
  if(length(unique(actual))==2){
    main_title <- paste('Deciles Sorted by Prediction \n','AUC =',auc(actual,pred))
    plt <- ggplot(x2)+geom_line(aes(x=Decile,y=Value,colour=Variable,linetype=Variable),size=0.75)+
      geom_point(aes(x=Decile,y=Value,colour=Variable),size=4)+
      theme_bw()+scale_linetype_manual(values=c("longdash", "solid"))+
      scale_colour_manual(values=c("blue","red"))+
      theme(legend.position='bottom',plot.title=element_text(size=10))+
      scale_x_continuous(labels=comma,breaks=1:n)+scale_y_continuous(labels=percent,limits=c(0,1))+
      ggtitle(main_title)+ylab('Actual and Predicted')
    return(plt)
  }
  if(length(unique(actual))>2 & yscl=='dollar'){
    main_title <- paste('Deciles Sorted by Prediction \n',gini,'\n',Normgini)
    plt <- ggplot(x2)+geom_line(aes(x=Decile,y=Value,colour=Variable,linetype=Variable),size=0.75)+
      geom_point(aes(x=Decile,y=Value,colour=Variable),size=4)+
      theme_bw()+scale_linetype_manual(values=c("longdash", "solid"))+
      scale_colour_manual(values=c("blue","red"))+
      theme(legend.position='bottom',plot.title=element_text(size=10))+
      scale_x_continuous(labels=comma,breaks=1:n)+scale_y_continuous(labels=dollar)+
      ggtitle(main_title)+ylab('Actual and Predicted')
    return(plt)
  }
  if(length(unique(actual))>2 & yscl=='percent'){
    main_title <- paste('Deciles Sorted by Prediction \n',gini,'\n',Normgini)
    plt <- ggplot(x2)+geom_line(aes(x=Decile,y=Value,colour=Variable,linetype=Variable),size=0.75)+
      geom_point(aes(x=Decile,y=Value,colour=Variable),size=4)+
      theme_bw()+scale_linetype_manual(values=c("longdash", "solid"))+
      scale_colour_manual(values=c("blue","red"))+
      theme(legend.position='bottom',plot.title=element_text(size=10))+
      scale_x_continuous(labels=comma,breaks=1:n)+
      scale_y_continuous(labels=percent_format())+
      ggtitle(main_title)+ylab('Actual and Predicted')
    return(plt)
  }
  if(length(unique(actual))>2 & yscl=='other'){
    main_title <- paste('Deciles Sorted by Prediction \n',gini,'\n',Normgini)
    plt <- ggplot(x2)+geom_line(aes(x=Decile,y=Value,colour=Variable,linetype=Variable),size=0.75)+
      geom_point(aes(x=Decile,y=Value,colour=Variable),size=4)+
      theme_bw()+scale_linetype_manual(values=c("longdash", "solid"))+
      scale_colour_manual(values=c("blue","red"))+
      theme(legend.position='bottom',plot.title=element_text(size=10))+
      scale_x_continuous(labels=comma,breaks=1:n)+
      ggtitle(main_title)+ylab('Actual and Predicted')
    return(plt)
  }
}
# This is a really long piece of code, and it was really annoying to write
# but it creates pretty graphs with minimal effort
my_gplot <- function(var1,var2=NULL,var_name1='X Variable',var_name2='Y Variable',TYPE='HISTOGRAM',
                     byvar=NULL,main='',xlabl='X Variable',ylabl='Y Variable',regline=T,lfit='auto',
                     lsize=1,ltype=1,dsize=3,dshape=16,size_g=14,xscal=comma,xlims=NULL,ylims=NULL,
                     yscal=comma,lincol='plum3',dotcol='deepskyblue'){
  my_df <- data.frame(var1)
  if(length(var1)>999){dsize <- 1}
  if(is.null(xlims)==TRUE){xlims<- c(min(var1),max(var1))}
  if(is.null(var2)==FALSE){
    my_df <- data.frame(var1,var2)
    if(is.null(byvar)==FALSE){
      my_df$byvar <- factor(byvar)
    }
    TYPE <- 'SCATTER'
  }
  if(toupper(TYPE)=='SCATTER'){
    if(is.null(byvar)==TRUE){
      plt <- ggplot(my_df,aes(x=var1,y=var2))+
        geom_point(size=dsize,shape=dshape)+theme_bw()+
        ggtitle(main)+xlab(xlabl)+ylab(ylabl)+
        theme(legend.position='bottom',plot.title=element_text(size=size_g))+
        scale_x_continuous(labels=xscal)+scale_y_continuous(labels=yscal)
    } else {
      plt <- ggplot(my_df,aes(x=var1,y=var2,colour='blue'))+
        geom_point(size=dsize,shape=dshape)+theme_bw()+
        ggtitle(main)+xlab(xlabl)+ylab(ylabl)+
        theme(legend.position='bottom',plot.title=element_text(size=size_g))+
        scale_x_continuous(labels=xscal)+scale_y_continuous(labels=yscal)  
    }
    if(regline==TRUE){
      plt <- plt+geom_smooth(method=lfit,colour=lincol,linetype=3,size=1)
    }
  }
  if(toupper(TYPE)=='HISTOGRAM'){
    # if(is.null(ylabl)==TRUE){ylabl <- 'Percent'}
    if(is.null(ylabl)==TRUE){ylabl <- 'Count'}
    # plt <- ggplot(my_df,aes(x=var1,y=..density..))+geom_histogram(fill='blue')+
    # scale_y_continuous(labels=percent_format())+
    # scale_x_continuous(labels=xscal)+theme_bw()+
    # xlab(xlabl)+ylab(ylabl)+
    # ggtitle(main)+
    # theme(legend.position='bottom',plot.title=element_text(size=size_g))
    plt <- ggplot(my_df,aes(x=var1))+geom_histogram(fill='blue')+
      scale_y_continuous(labels=comma)+
      scale_x_continuous(labels=xscal)+theme_bw()+
      xlab(xlabl)+ylab(ylabl)+
      ggtitle(main)+theme(legend.position='bottom',plot.title=element_text(size=size_g))
    plt <- plt +coord_cartesian(xlim=xlims)
  }
  if(toupper(TYPE)=='ECDF'){
    if(is.null(my_df$byvar)==TRUE){
      plt <- ggplot(my_df,aes(x=var1))+stat_ecdf(colour='blue')+
        scale_y_continuous(labels=percent_format(),limits=c(0,1))+
        scale_x_continuous(labels=xscal,limits=xlims)+theme_bw()+
        xlab(xlabl)+ylab(ylabl)
      ggtitle(main)+theme(legend.position='bottom',plot.title=element_text(size=size_g))
    }
    if(ylabl=='Y Variable'){
      plt <- plt+ylab('Cumulative Percent')
    }
    if(is.null(my_df$byvar)==FALSE){
      plt <- ggplot(my_df,aes(x=var1,colour=byvar))+stat_ecdf()+
        scale_y_continuous(labels=percent_format())+
        scale_x_continuous(labels=xscal)+theme_bw()+
        xlab(xlabl)+ylab(ylabl)
      ggtitle(main)+theme(legend.position='bottom',plot.title=element_text(size=size_g))
    }
  }
  if(toupper(TYPE)=='DENSITY'){
    plt <- ggplot(my_df,aes(x=var1))+geom_density(color='blue')+
      scale_y_continuous(labels=percent_format())+
      scale_x_continuous(labels=xscal)+theme_bw()+
      xlab(xlabl)+ylab('Percent')+
      ggtitle(main)+
      theme(legend.position='bottom',plot.title=element_text(size=size_g))
  }
  if(toupper(TYPE)=='LINE'){
    plt <- ggplot(my_df,aes(x=var1,y=var2,byvar))+geom_line(colour=dotcol,size=lsize,lty=ltype)+theme_bw()+
      ggtitle(main)+xlab(xlabl)+ylab(ylabl)+
      theme(legend.position='bottom',plot.title=element_text(size=size_g))+
      scale_x_continuous(labels=xscal)+scale_y_continuous(labels=yscal)+
      geom_point(colour=dotcol,size=dsize,shape=dshape)
  }
  return(plt)
  # if(toupper(TYPE)=='LINE'){
  # plt <- ggplot(my_df,aes(reorder(xvar,by)),by+geom_bar(fill="0000FF")+coord_flip()+theme_bw()+
  # ggtitle(main)+xlab(xlabl)+ylab(ylabl)+
  # theme(legend.position='bottom',plot.title=element_text(size=size_g))+
  # +scale_y_continuous(labels=percent_format(),limits=cc(0,1))+
  # geom_text(label=paste(round(xvar*100,0),'%',sep=''),size=3,hjust=-0.25)
  # }
}
# Bootstrapped Standard Errors in GLMNET
glmnet_boot <- function(ydep,xs,trainfilter,model,b=1e3,myalpha){
  lamb <- model$lambda.min
  train_xs <- xs[trainfilter==1,]
  n_obs <- dim(train_xs)[1]
  bs_out<-data.frame(as.matrix(t(coef(model,s='lambda.min'))))
  for(i in 1:dim(bs_out)[2]){bs_out[,i] <- 0 }
  # one thousand bootstrap estimates
  for(i in 1:b){
    indx <- sample(1:n_obs,1e4,replace=TRUE)
    bs_xs <- train_xs[indx,]
    bs_ys <- ydep[indx]
    bs_mod <- glmnet(y=bs_ys,x=bs_xs,
                     family='poisson',
                     alpha=myalpha,
                     lambda=lamb)
    tmp <- data.frame(as.matrix(t(coef(bs_mod))))
    bs_out <- rbind(bs_out,tmp)
    if((i%%100)==0){
    	print(paste("Iteration",i,"complete."))
    }
  }
  names(bs_out)[1] <- '(Intercept)'
  n <- dim(bs_out)[1]
  bs_out <- bs_out[2:n,]
  final <- data.frame(mean=sapply(bs_out,FUN=mean))
  final$SEs <- sapply(bs_out,FUN=sd)
  tmp <- data.frame(var=as.matrix(coef(model,s='lambda.min')))
  final$mean <- tmp[,1]
  final$CI_025 <- final$mean-1.94*final$SEs
  final$CI_975 <- final$mean+1.94*final$SEs
  names(final)[1] <- 'Estimate'
  final$T_stat <- final$Estimate / final$SEs
  return(final)
}
glmnet_preds <- function(glmnet_model,xs){
	preds <- c(predict(glmnet_model,xs,s='lambda.min',type='class'))
	probs <- c(predict(glmnet_model,xs,s='lambda.min',type='response'))
	return(data.frame(Class=preds,Probs=probs))
}
# Bootstrap AUC
BootStrapAUC <- function(preds,actual,BS_Reps=100,CI_Lower=0.025,CI_Upper=0.975){
  auc_out <- c()
  for(i in 1:BS_Reps){
    vals <- 1:length(actual)
    indx <- sample(vals,size=500,replace=T)
    BSPred <- preds[indx]
    BSActual <- actual[indx]
    auc_out[i] <- auc(BSActual,BSPred)
  }  
  print("The Expected AUC is:")
  print(mean(auc_out))
  tmpval <- paste((CI_Upper - CI_Lower)*100,'%',sep='')
  print(paste("The",tmpval,"Confidence Interval is:"))
  print(quantile(auc_out,probs=c(CI_Lower,CI_Upper)))
  par(mfrow=c(1,2))
  auc_val <- round(mean(auc_out),4)
  auc_val2<- round(median(auc_out),4)
  ttle <- paste("Boot Strapped AUC\n",'Mean =',auc_val)
  ttle2<- paste("Boot Strapped AUC\n",'Median =',auc_val2)
  plot(density(auc_out),main=ttle,col='blue',lty=2);grid(5,5,'gray44')
  plot(ecdf(auc_out),main=ttle2,col='red',lty=1,
       xlab='Area Under the Curve (AUC)',
       ylab='Cumulative Percent');grid(5,5,'gray44')
  par(mfrow=c(1,1))
}
PartialDependence <- function(xs,varname,model_name,upper=90){
	# This is a partial depedence plot for a given feature
	# this takes in a column, zeros out all of the other columns
	# then populates the feature with the min value to the 
	# highest specified percentile value (for max, set upper=100)
	# then gives the partial prediction.
  n <- 100
  k <- dim(xs)[2]
  mfx <- xs[1:n,]
  varlist <- colnames(mfx)
  x_j <- which(varlist==varname)
  if(is.numeric(xs[,x_j])==T){
    # Using the quantile makes it cleaner than using the min to max
    xvar <- array(quantile(xs[,x_j],1:upper/100))
    xvar <- c(rep(min(xs[,x_j]),100-length(xvar)),xvar)
    mfx[,1:k] <- 0
    mfx[,x_j] <-xvar
    y_dx <- predict(model_name,mfx)
    mfx_df <- data.frame(y_dx,xvar)
    plt <- ggplot(mfx_df,aes(x=xvar,y=y_dx))+
      geom_line(colour='red ')+
      scale_y_continuous(labels=comma)+
      theme_bw()+xlab(varname)+ylab('Partial Prediction')+
      ggtitle(paste('Estimated Partial Prediciton of',varname))
    return(plt)    
  } else {
    return(NULL)
  }
}
EmpiricalDependence <- function(xs,varname,model_name,RankNum=10){
	# This is a partial dependence plot for a given feature
	# this takes in a column, zeros out all of the other columns,
	# then sorts by the features value and takes the average 
	# prediction for that bucket rank -- the default is 10
  n <- dim(xs)[1]
  k <- dim(xs)[2]
  mfx <- xs[1:n,]
  varlist <- substr(colnames(mfx),1,nchar(varname))
  x_j <- which(varlist==varname)
  indx <- 1:k
  newvars <- indx[(!indx %in% x_j)==T]
  mfx[,newvars] <- 0 
  xvar <- predict(model_name,xs)
  y_dx <- predict(model_name,mfx)
  if(len(xvar)==len(y_dx)){
    mfx_df <- data.frame(y_dx,xvar)
    mfx_df$xrank <- myrank(mfx_df$xvar,k=RankNum)
    tmp <- sqldf('select xrank, max(xvar) as MaxValue, avg(y_dx) as Prediction from mfx_df group by 1')
    plt <- ggplot(tmp,aes(x=xrank,y=Prediction))+
      geom_line(colour='blue')+
      scale_y_continuous(labels=comma)+
      theme_bw()+xlab(varname)+ylab('Partial Prediction')+
      ggtitle(paste('Empirical Partial Prediciton of',varname))
    return(plt)
    } else {
    return(NULL)
  }
}
EmpVarDep <- function(xs,varlist,GBM_Model,NRanks=100){
  # Empirical Dependence for All Variables
  for(i in varlist){
    EmpiricalDependence(xs,i,GBM_Model,NRanks)  
  }  
}
ParVarDep <- function(xs,varlist,GBM_Model,NRanks=100){
  # Empirical Dependence for All Variables
  for(i in varlist){
    PartialDependence(xs,i,GBM_Model,NRanks)  
  }  
}
AllVariableDependence<-function(xs,varlist,GBM_Model,NRanks){
  for(i in varlist){
    print(i)
    VariableDependence(xs,i,GBM_Model,NRanks)
  }
}
VariableDependence<- function(xs,varname,GBM_Model,NRanks=100){
  p1 <- EmpiricalDependence(xs,varname,GBM_Model,NRanks)
  p2 <- PartialDependence(xs,varname,GBM_Model,NRanks)
  if(is.null(p1)==T & is.null(p1)==T){
    out <- NULL
  } else
  if(is.null(p1)==F & is.null(p2)==F){
    out <- grid.arrange(p1,p2)        
  } else
  if(is.null(p1)==T & is.null(p2)==F){
    out <- p2
  } else
  if(is.null(p1)==F & is.null(p2)==T){
    out <- p1    
  }
  return(out)
}
MyExperVar <- function(tmpdf,tflt,y,xvar1,xvar2=NULL,alphaval=0.8,fldsval=10){
  ydep <- tmpdf[,y]
  features <- data.frame(tmpdf[,c(xvar1)]) ;names(features) <- xvar1
  if(len(xvar2)>0){
    features <- data.frame(tmpdf[,c(xvar1,xvar2)])
    features$xvar12 <- paste(features[,1],features[,2],sep='_')
    names(features) <- c(xvar1,xvar2,paste(substr(xvar1,1,3),substr(xvar2,1,3),sep=''))
  }
  for(i in 1:ncol(features)){
    features[,i] <- factor(features[,i])
  }
  features <-ReplaceMissingFunc(features,start_col=1)
  mymodelformula  <- myformula(colnames(features))  
  xs <- sparse.model.matrix(mymodelformula,data=features)
  print(dim(xs))
  registerDoParallel(cores=detectCores(all.test=TRUE))
  getDoParWorkers()
  set.seed(60)
  TmpModel <- cv.glmnet(
    x=xs[tflt,],
    y=ydep[tflt],
    nfolds=fldsval,
    type.measure='auc',
    parallel=TRUE,
    family='binomial',
    alpha=alphaval)
  TmpPred <- c(predict(TmpModel,newx=xs,s='lambda.min',type='response'))
  return(TmpPred)  
}
# My Count Experience Variable
MyCExperVar <- function(tmpdf,tflt,y,idvar,xvar1,xvar2=NULL){
  if(len(xvar2)==0){
    features <- data.frame(tmpdf[,c(idvar,xvar1)]); 
    names(features)[1] <- 'idvar'
    names(features)[2] <- 'xvar1'
    features$ydep <- tmpdf[,y]
    features$TmpFilter <- ifelse(tflt==1,1,0)
    out <- sqldf('select xvar1, sum(ydep) as XSum, count(*) as XCount
                 from features where TmpFilter=1 group by 1')
    out2<- sqldf('select a.idvar ,b.xvar1,b.XSum,b.XCount from features a left join out b
                 on a.xvar1 = b.xvar1 order by a.idvar')
    out2$XSum[is.na(out2$XSum)==T]<-0
    out2$XCount[is.na(out2$XCount)==T]<-0
    names(out2) <- c(idvar,xvar1,paste(xvar1,'Cnt',sep=''),paste(xvar1,'Sum',sep=''))    
  }
  if(len(xvar2)>0){
    features <- data.frame(tmpdf[,c(idvar,xvar1,xvar2)])
    names(features) <- c('idvar','xvar1','xvar2')
    names(features)[1] <- 'idvar'
    names(features)[2] <- 'xvar1'
    names(features)[3] <- 'xvar2'
    features$ydep <- tmpdf[,y]
    features$TmpFilter <- ifelse(tflt==1,1,0)
    out <- sqldf('select xvar1, xvar2,
                 sum(ydep) as XSum, 
                 count(*) as XCount
                 from features where TmpFilter=1 group by 1,2')
    out2<- sqldf('select a.idvar ,a.xvar1,a.xvar2, b.XSum,b.XCount from features a left join out b
                 on a.xvar1 = b.xvar1  and a.xvar2 = b.xvar2 order by a.idvar')
    out2$XSum[is.na(out2$XSum)==T]<-0
    out2$XCount[is.na(out2$XCount)==T]<-0
    names(out2) <- c(idvar,xvar1,xvar2,paste(xvar1,xvar2,'Cnt',sep=''),paste(xvar1,xvar2,'Sum',sep=''))
  }
  return(out2)
}
GBMExpVar <- function(tmpdf,tflt,y,xvar1,xvar2=NULL,prnt=T){
  ydep <- tmpdf[,y]
  if(len(xvar2)==0){
    features <- data.frame(tmpdf[,c(xvar1)])
    names(features) <- xvar1
  }
  if(len(xvar2)>0){
    features <- data.frame(tmpdf[,c(xvar1,xvar2)])
    features$xvar12 <- paste(features[,1],features[,2],sep='_')
    names(features) <- c(xvar1,xvar2,paste(substr(xvar1,1,3),substr(xvar2,1,3),sep=''))
  }
  for(i in 1:dim(features)[2]){
    # Turning them into integers
    features[,i] <- as.integer(factor(features[,i]))
  }
  xs <- model.matrix(~.,data=features)
  set.seed(504)
  TmpModel<- xgboost(data = xs[tflt,],
                        label = ydep[tflt],
                        max.depth = 30,
                        eta = 1,
                        nround = 200,
                        objective = "binary:logistic",
                        verbose=0)
  out<- predict(TmpModel,xs)
  if(prnt==T){
    print(auc(ydep[tflt],out[tflt]))    
  }
  return(out)
}
#========================================================================================================================
# ROC databuild and use.
#========================================================================================================================
cmplot <- function(pred,actual,colrs=c('2','2','4','4'),ymax=1,ttl){
  # This is a confusion matrix plot
  print(table(pred,actual))
  sdf <- data.frame(table(pred,actual))
  sdf$Names <- factor(paste('Pred=',sdf$pred,' | Actual=',sdf$actual,sep=''))
  sdf$Names <- factor(sdf$Names,levels=sdf$Names[order(sdf$Freq)])
  sdf$Percent <- sdf$Freq/sum(sdf$Freq)
  accuracy <- round(sum(diag(table(pred,actual)))/length(pred),4)*100
  ttl <- paste("The Classification Accuracy is ", accuracy,'% \n',ttl,sep='')
  sdf$Plabels <- paste(round(sdf$Percent,3)*100,'%',sep='')
  pout <- ggplot(sdf,aes(x=Names,y=Percent,fill=Names))+geom_bar(stat='identity')+
  theme_bw()+coord_flip()+geom_text(aes(label=Plabels),size=5,hjust = -0.05)+
  scale_y_continuous(labels=percent_format())+xlab('')+ggtitle(ttl)+
  scale_fill_manual(values=colrs)+ylim(c(0,ymax))+ylab('Percent of Data')+
  theme(legend.position='bottom',legend.title=element_blank(),
  axis.text=element_text(size=12,face='bold'))
  print("Ignore any warning signs")
  pout
}
rocdata <- function(grp, pred){
  # Produces x and y coordinates for ROC curve plot
  # Arguments: grp - labels classifying subject status
  #            pred - values of each observation
  # Output: List with 2 components:
  #         roc = data.frame with x and y coordinates of plot
  #         stats = data.frame containing: area under ROC curve, p value, upper and lower 95% confidence interval
  grp <- as.factor(grp)
  if (length(pred) != length(grp)) {
    stop("The number of classifiers must match the number of data points")
  } 
  
  if (length(levels(grp)) != 2) {
    stop("There must only be 2 values for the classifier")
  }
  
  cut <- unique(pred)
  tp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[2])))
  fn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[2])))
  fp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[1])))
  tn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[1])))
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  roc = data.frame(x = fpr, y = tpr)
  roc <- roc[order(roc$x, roc$y),]
  
  i <- 2:nrow(roc)
  auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
  
  pos <- pred[grp == levels(grp)[2]]
  neg <- pred[grp == levels(grp)[1]]
  q1 <- auc/(2-auc)
  q2 <- (2*auc^2)/(1+auc)
  se.auc <- sqrt(((auc * (1 - auc)) + ((length(pos) -1)*(q1 - auc^2)) + ((length(neg) -1)*(q2 - auc^2)))/(length(pos)*length(neg)))
  ci.upper <- auc + (se.auc * 0.96)
  ci.lower <- auc - (se.auc * 0.96)
  
  se.auc.null <- sqrt((1 + length(pos) + length(neg))/(12*length(pos)*length(neg)))
  z <- (auc - 0.5)/se.auc.null
  p <- 2*pnorm(-abs(z))
  
  stats <- data.frame (auc = auc,
                       p.value = p,
                       ci.upper = ci.upper,
                       ci.lower = ci.lower
  )
  
  return (list(roc = roc, stats = stats))
}
rocplot.single <- function(grp, pred, p.value = FALSE){
  require(ggplot2)
  plotdata <- rocdata(grp, pred)
  
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 3), " (P=", signif(p.value, 3), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 3), " (95%CI ", signif(ci.upper, 3), " - ", signif(ci.lower, 2), ")", sep=""))
  }
  title <- 'ROC Plot'
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
    geom_line(aes(colour = "red")) +
    geom_abline (intercept = 0, slope = 1) +
    theme_bw() +
    scale_x_continuous("False Positive Rate (1-Specificity)") +
    scale_y_continuous("True Positive Rate (Sensitivity)") +
    scale_colour_manual(labels = annotation, values = "red")  +
    ggtitle(title)+
    theme(legend.position='bottom',
          plot.title = element_text(face="bold", size=14), 
          axis.title.x = element_text(face="bold", size=12),
          axis.title.y = element_text(face="bold", size=12, angle=90),
          legend.justification=c(1,0), 
          legend.position=c(1,0),
          legend.title=element_blank(),
          legend.key = element_blank())
  return(p)
}
rocplot.multiple <- function(test.data.list, groupName = "grp", predName = "res", title = "ROC Plot", p.value = FALSE){
  require(plyr)
  require(ggplot2)
  plotdata <- llply(test.data.list, function(x) with(x, rocdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)))))
  plotdata <- list(roc = ldply(plotdata, function(x) x$roc),
                   stats = ldply(plotdata, function(x) x$stats)
  )
  
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 3), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 3), " (95%CI ", signif(ci.lower, 2), " - ", signif(ci.upper, 2), ")", sep=""))
  }
  
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
    geom_line(aes(colour = .id)) +
    geom_abline (intercept = 0, slope = 1) +
    theme_bw() +
    scale_x_continuous("False Positive Rate (1-Specificity)") +
    scale_y_continuous("True Positive Rate (Sensitivity)") +
    scale_colour_brewer(palette="Set1", breaks = names(test.data.list), labels = paste(names(test.data.list), ": ", annotation, sep = "")) +
    theme(legend.position='bottom',
          plot.title = element_text(face="bold", size=14), 
          axis.title.x = element_text(face="bold", size=12),
          axis.title.y = element_text(face="bold", size=12, angle=90),
          legend.justification=c(1,0), 
          legend.position=c(1,0),
          legend.title=element_blank(),
          legend.key = element_blank())
  return(p)
}
# Cross Fold Validation for XGBoost
xgboost_cv <- function(xs,yvar,tflt,loss='binary:logistic',nfolds=10,nrounds=200,Verbose=0,maxseq=NULL){
  # This is a way to do cross fold validation for the Maximum Depth
  Xs <- xs[tflt,]
  yvar <- yvar[tflt]
  set.seed(510)
  folds <- cut(seq(1,nrow(Xs)),breaks=nfolds,labels=FALSE)
  flds <- unique(folds)
  if(len(maxseq)==0){
    max_depth_var <- seq(1,21,length=11); 
    max_depth_var<- round(max_depth_var,0)    
  }
  if(len(maxseq)>0){
    max_depth_var <- maxseq
  }
  out <- expand.grid(FoldIteration=flds,MaxValueIterated=max_depth_var)
  out$TrainAUC <- NA
  out$TestAUC <- NA
  for(i in 1:length(flds)){
    xtst_tmp <- Xs[folds %in% flds[i],]
    xtrn_tmp <- Xs[!folds %in% flds[i],]
    ytst_tmp <- yvar[folds %in% flds[i]]
    ytrn_tmp <- yvar[!folds %in% flds[i]]
    for(j in 1:length(max_depth_var)){
      print(paste('Fold',i,'Iteration',j))
      fold_filter <- out$FoldIteration==i & out$MaxValueIterated==max_depth_var[j]
      mod <- xgboost(data = xtrn_tmp,
                     label = ytrn_tmp,
                     max.depth = max_depth_var[j],
                     eta = 0.05, 
                     nround = nrounds,
                     objective = loss,
                     verbose=Verbose)
      trn_pred <- predict(mod,xtrn_tmp); trn_clss <- ifelse(trn_pred>=0.5,1,0)
      tst_pred <- predict(mod,xtst_tmp); tst_clss <- ifelse(tst_pred>=0.5,1,0)
      out$TestAUC[fold_filter] <- auc(ytst_tmp,tst_pred)
      out$TrainAUC[fold_filter]<- auc(ytrn_tmp,trn_pred)
      out$TestFPR[fold_filter] <- fpr(ytst_tmp,tst_clss)
      out$TrainFPR[fold_filter]<- fpr(ytrn_tmp,trn_clss)
      out$TestFNR[fold_filter] <- fnr(ytst_tmp,tst_clss)
      out$TrainFNR[fold_filter]<- fnr(ytrn_tmp,trn_clss)
      out$Iteration[fold_filter] <- j
    }
  }
  out <- out[complete.cases(out),]
  smry<-sqldf('select Iteration, MaxValueIterated,
              avg(TestAUC) as AVG_Test_AUC,
              stdev(TestAUC) as SD_Test_AUC,
              avg(TestFNR) as AVG_Test_FNR,
              avg(TestFPR) as AVG_Test_FPR,
              avg(TrainAUC) as AVG_Train_AUC,
              stdev(TrainAUC) as SD_Train_AUC,
              avg(TrainFNR) as AVG_Train_FNR,
              avg(TrainFPR) as AVG_Train_FPR,
              count(*) as nrows
              from out group by 1,2')
  with(smry,plot(Iteration,AVG_Test_AUC,ylim=c(0,1),col='red',pch=16,type='b'))
  with(smry,lines(Iteration,AVG_Test_AUC+2*SD_Test_AUC,ylim=c(0,1),col='blue',pch=16,type='b'))
  with(smry,lines(Iteration,AVG_Test_AUC-2*SD_Test_AUC,ylim=c(0,1),col='blue',pch=16,type='b'))
  grid(5,5,'gray44')
  smry <- smry[order(smry$AVG_Test_AUC,decreasing=F),]
  print(smry)
  return(smry$MaxValueIterated[dim(smry)[1]])
}
# Useful for modeling
PredictorImportance<- function(xs,var_list,ydep,model,filename='',newvarnames=NA,ToCSV=F){
  # This finds the predictor importance for GLMNET and Other types of models
  Sys.setlocale('LC_ALL','C')  # Takes care of foreign characters
  mfx <- xs
  varnames <- colnames(xs)
  varbuckets <- var_list
  labels <- c()
  labeldf <- data.frame(varbuckets,TEST=NA)
  for(i in 1:length(varnames)){
    for(j in 1:length(labeldf$varbuckets)){
      labeldf$TEST[j] <- grepl(labeldf$varbuckets[j],varnames[i])    
    }
    # Output only the first appropriate class
    labels[i] <- as.character(labeldf$varbuckets[labeldf$TEST==T][1])
  }
  mfx_details <-data.frame(labels,varnames)
  if((class(model)=='cv.glmnet')==TRUE){
    pred_total <- predict(model,mfx,s='lambda.min',type='response')    
  } else {
    pred_total <- predict(model,mfx)    
  }
  actual <- ydep
  outdf <-data.frame(varbuckets,RMSE_Without=0)
  StartTime <- Sys.time()
  for(i in 1:length(varbuckets)){
    print(i)
    indx <- 1:length(labels)
    colfilter <- indx[word_search(varbuckets[i],labels)==1]
    Xs_without <- mfx
    if(length(colfilter)>0){
      Xs_without[,colfilter] <- 0    
      if((class(model)=='cv.glmnet')==TRUE){
        pred_without <- predict(model,Xs_without,s='lambda.min',type='response')
      } else {
        pred_without <- predict(model,Xs_without)
      }
      fitstats <- predcompare(pred_without,actual)  
      outdf$RMSE_Without[i] <- fitstats$RMSE
    }
  }
  totfitstats <- predcompare(pred_total,actual)  
  outdf$RMSE_Total <- totfitstats$RMSE
  outdf$RMSE_Reduction <- with(outdf,RMSE_Without-RMSE_Total)
  outdf$RMSE_Reduction[outdf$RMSE_Reduction<0] <- 0
  outdf$RelativeImportance <- with(outdf,RMSE_Reduction/max(RMSE_Reduction))
  FinishTime <- Sys.time()
  RunTime <- FinishTime - StartTime
  print(RunTime)
  outfilename <- paste(date_stamp(),filename,'.csv',sep='')
  if(ToCSV==T){
    write.csv(outdf,outfilename,row.names=FALSE)    
  }
  if(length(newvarnames)>1){
    outdf$varbuckets <- newvarnames
  }
  outdf <- outdf[order(outdf$RelativeImportance,decreasing=T),]
  outdf$varbuckets <- factor(outdf$varbuckets,
                             levels=outdf$varbuckets[order(outdf$RelativeImportance)])
  myplt <- ggplot(outdf,aes(x=varbuckets,y=RelativeImportance))+
    geom_bar(stat='identity',fill='blue')+coord_flip()+
    theme_bw()+xlab('')+
    ylab('Model Degradation from Exclusion of Feature')+
    ggtitle("Relative Predictor Importance")
  outlist <- list()
  outlist[[1]] <- myplt
  outdf <- outdf[,c('varbuckets','RMSE_Reduction','RelativeImportance')]
  outdf$RMSE_Reduction <- round(outdf$RMSE_Reduction,4)
  names(outdf) <-c('Variable','Model Improvement','Relative Importance')
  outlist[[2]] <- outdf
  return(outlist)
}
#========================================================================================================================
print("Your Functions are now loaded, dawg.")
#========================================================================================================================
######################################################## END ############################################################
#========================================================================================================================
