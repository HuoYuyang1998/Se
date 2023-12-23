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


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

###Source functions and algorithm classes
source("Functions-scop.R")
source("AlgorithmClass-scop.R")


### Parameter setting
d<-10 ### dimension of covariates
n<-200 ### size of calibration set
n_train<-200 ### size of training set
n_rest<-n+n_train ### size of labeled set
m<-200 ### size of unlabeled set
N<-n_train+n+m ### size of the overall data
alpha<-0.1 ### target FCR
sigma<-1 ### noise strength
ns<-1000 ### replication times



Thresholding_type_array=c("cons","POS","top")
FDR_level=0.2 ### POS selection rule
K=60  ### Tok-K selection rule

algo_array=c(new('OLS'),new('SVM')) ### algorithms for predictions and corresponding scenarios



###²¢ÐÐÔËËã
cl = makeCluster(10)
registerDoParallel(cl)
Result<-foreach(iter=1:ns,.combine="rbind",.packages = c('MASS',"randomForest","kernlab"))%dopar% {
  info<-data.frame()
  
  for (algo in algo_array) {
    
    b_0=if(algo@name=="OLS"){-1}else if(algo@name=="SVM"){-8}
    
    
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
    
    
    
    model_train<-fitting(algo,X_train,Y_train)
    
    
    Y_cal_hat<-predict(model_train,as.data.frame(x=X_cal))
    Y_test_hat<-predict(model_train,as.data.frame(x=X_test))
    
    
    R_cal<-abs(Y_cal-Y_cal_hat)
    R_test<-abs(Y_test-Y_test_hat)
    
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
        
        
        info<-rbind(info,list(FCR=Result_OCP$FCP,method="OCP",Length=Result_OCP$Length,Thresholding_type=Thresholding_type,Setting=algo@name,n=n,m=m))
        info<-rbind(info,list(FCR=Result_ACP$FCP,method="ACP",Length=Result_ACP$Length,Thresholding_type=Thresholding_type,Setting=algo@name,n=n,m=m))
        info<-rbind(info,list(FCR=Result_SCOP$FCP,method="SCOP",Length=Result_SCOP$Length,Thresholding_type=Thresholding_type,Setting=algo@name,n=n,m=m))
        
      }
  }
  info
  return(info)
}
stopCluster(cl)

## processing the results


Result$method<-factor(Result$method,levels=c("SCOP","OCP","ACP"))

Result$Setting[Result$Setting=="OLS"]<-"Scenario A"
Result$Setting[Result$Setting=="SVM"]<-"Scenario B"
Result$Setting<-factor(Result$Setting,levels=c("Scenario A","Scenario B"))

Result$FCR<-100*Result$FCR








### obtain the contents of the table

Table<-Result%>%
  group_by(method,Thresholding_type,Setting)%>%
  dplyr::summarize(FCR=mean(FCR),Length=mean(Length,na.rm=TRUE))
Table
