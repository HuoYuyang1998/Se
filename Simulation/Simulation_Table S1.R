library(MASS)
library(kernlab)
library(foreach)
library(doParallel)
require(tidyverse)
require(ggplot2)
library(glmnet)
library(caret)
library(reshape2)
library(latex2exp)
library(ggpubr)

library(randomForest)
library(quantreg)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

###Source functions and algorithm classes
source("Functions-scop.R")
source("AlgorithmClass-scop.R")


### Define functions for quantile regression score
QRscore<-function(algo,model_low,model_high,X_cal,Y_cal,islength=FALSE)
{if(islength)
{  Y_low<-Pred(algo,model_low,t(X_cal))
Y_high<-Pred(algo,model_high,t(X_cal))
ylen<-length(Y_cal)
return(pmax(rep(Y_low,ylen)-Y_cal,Y_cal-rep(Y_high,ylen)))
}else
{  Y_low<-Pred(algo,model_low,X_cal)
Y_high<-Pred(algo,model_high,X_cal)
return(pmax(Y_low-Y_cal,Y_cal-Y_high))
}
  
}


QRlength<-function(algo,model_low,model_high,X_test)
{
  Y_low<-Pred(algo,model_low,X_cal)
  Y_high<-Pred(algo,model_high,X_cal)
  return(Y_high-Y_low)
}

setClass("QR",slots = list(name="character"), prototype=list(name="QR") )
setMethod("fitting","QR",function(obj,X,Y,lambda){
  datawork=data.frame(X,y=Y)
  return(rq(y~.,lambda,datawork))
  
})
setMethod("Pred","QR",function(obj,model,X_test,lens=0){
  predict(model,as.data.frame(x=X_test))
  
})
setMethod("DataGen","QR",function(obj,X,beta,sigma){
  uX=X%*%beta
  Y=uX+sapply(1:dim(X)[1],function(i){rnorm(1,0,sigma+abs(uX[i]))})#生成Y
  return(data.frame(x=X,y=Y))
})


### Parameter setting
d<-10 ### dimension of covariates
n<-200 ### size of calibration set
n_train<-200 ### size of training set
n_rest<-n+n_train ### size of labeled set
m<-200 ### size of unlabeled set
N<-n_train+n+m ### size of the overall data
alpha<-0.1 ### target FCR
sigma<-1 ### noise strength
lambda_low=0.05 ### lower-percentile of quantile regression
lambda_high=0.95 ### higher-percentile of quantile regression
ns<-1000 ### replication times



Thresholding_type_array=c("cons","POS","top")
FDR_level=0.3 ### POS selection rule
K=60  ### Tok-K selection rule

algo_array=c(new('QR')) ### algorithms for predictions and corresponding scenarios



###并行运算
cl = makeCluster(10)
registerDoParallel(cl)
Result<-foreach(iter=1:ns,.combine="rbind",.packages = c('MASS',"randomForest","kernlab","quantreg"))%dopar% {
  info<-data.frame()
  
  for (algo in algo_array) {
    
    b_0=-1
    
    
    X<-matrix(runif(N*d,-1,1),nrow=N,ncol=d)
    beta<-runif(d,-2,2)
    data<-DataGen(algo,X,beta,sigma)
    
    datawork<-DataSplit(data,N,m,n,n_rest)
    data_train<-datawork$data_train
    data_cal<-datawork$data_cal
    data_test<-datawork$data_test
    
    X_train<-as.matrix(data_train[colnames(data_train)[-d-1]])
    Y_train<-as.matrix(data_train$y)
    X_cal<-as.matrix(data_cal[colnames(data_cal)[-d-1]])
    Y_cal<-as.matrix(data_cal$y)
    X_test<-as.matrix(data_test[colnames(data_test)[-d-1]])
    Y_test<-as.matrix(data_test$y)
    
    
    model_train_low<-fitting(algo,X_train,Y_train,lambda = lambda_low)
    model_train_high<-fitting(algo,X_train,Y_train,lambda = lambda_high)
    
    model_train<-fitting(algo,X_train,Y_train,lambda = 0.5)
    
    
    Y_cal_hat<-predict(model_train,as.data.frame(x=X_cal))
    Y_test_hat<-predict(model_train,as.data.frame(x=X_test))
    
    
    R_cal<-QRscore(algo,model_train_low,model_train_high,X_cal,Y_cal)
    R_test<-QRscore(algo,model_train_low,model_train_high,X_test,Y_test)
    
    ### Selection procedure
    T_test<-Y_test_hat
    T_cal<-Y_cal_hat
    
    
    for (Thresholding_type in Thresholding_type_array) {
      
      Selection_rule<-list(Thresholding_type=Thresholding_type,K=K,FDR_level=FDR_level,b_0=b_0)
      
      L<-Selection_thresholding(Selection_rule,T_cal,T_test,Y_cal)
      Select_test<-which(T_test<=L)
      
      
      
      Result_SCOP<-SCOP(R_cal,R_test,T_cal,T_test,Select_test,L,alpha,Selection_rule,exchangeable_adapted=TRUE)
      Result_OCP<-OCP(R_cal,R_test,Select_test,alpha)
      Result_ACP<-ACP(R_cal,R_test,T_test,Select_test,alpha,Selection_rule,simplified=TRUE)
      
      Len=mean(QRlength(algo,model_train_low,model_train_high,X_test[Reject_test]))
      
      info<-rbind(info,list(FCR=Result_OCP$FCP,method="OCP",Length=Result_OCP$Length+Len,Thresholding_type=Thresholding_type,Setting=algo@name,n=n,m=m))
      info<-rbind(info,list(FCR=Result_ACP$FCP,method="ACP",Length=Result_ACP$Length+Len,Thresholding_type=Thresholding_type,Setting=algo@name,n=n,m=m))
      info<-rbind(info,list(FCR=Result_SCOP$FCP,method="SCOP",Length=Result_SCOP$Length+Len,Thresholding_type=Thresholding_type,Setting=algo@name,n=n,m=m))
      
    }
  }
  info
  return(info)
}
stopCluster(cl)

## processing the results


Result$method<-factor(Result$method,levels=c("SCOP","OCP","ACP"))


Result$FCR<-100*Result$FCR








### obtain the contents of the table

Table<-Result%>%
  group_by(method,Thresholding_type)%>%
  dplyr::summarize(FCR=mean(FCR),Length=mean(Length,na.rm=TRUE))
Table
