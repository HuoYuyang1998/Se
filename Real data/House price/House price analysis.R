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
source("Functions-scop-house.R")
source("AlgorithmClass-scop.R")


### Loading and processing data

data<-read.table("data.csv",header=T,sep = ",")
data=na.omit(data)
colnames(data)[2]<-"y"
data$yr_built=2022-data$yr_built
data=data[data$y>0,]
data$y=data$y/1000000
continvar=subset(data,select=c(bedrooms,bathrooms,sqft_living,sqft_lot,floors,waterfront,view,condition,sqft_above,sqft_basement,yr_built))

data=cbind(continvar,y=data$y)


### parameter setting
d<-dim(data)[2]-1 ### dimension of covariates
n<-500 ### size of calibration set
n_train<-500 ### size of training set
n_rest<-n+n_train ### size of labeled set
m<-500 ### size of unlabeled set
N<-n_train+n+m ### size of the overall data
alpha<-0.1 ### target FCR
sigma<-1 ### noise strength
ns<-500 ### replication times


algo=new('RF')
Thresholding_type_array=c("test","POS","clu")
Thresholding_quantile=0.7
b_0=0.6
FDR_level=0.2

### parallel computing
cl = makeCluster(10) ### number of parallel threads
registerDoParallel(cl)
Result<-foreach(iter=1:ns,.combine="rbind",.packages = c('MASS',"randomForest","kernlab"))%dopar% {
  info<-data.frame()
  
    
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
    T_test<- Y_test_hat
    T_cal<- Y_cal_hat
    
    
    for (Thresholding_type in Thresholding_type_array) {

        Selection_rule<-list(Thresholding_type=Thresholding_type,Thresholding_quantile=Thresholding_quantile,b_0=b_0,FDR_level=FDR_level)
        
        L<-Selection_thresholding(Selection_rule,T_cal,T_test,Y_cal)
        Select_test<-which(T_test>=L)
        
        
        
        Result_SCOP<-SCOP(R_cal,R_test,T_cal,T_test,Select_test,L,alpha,Selection_rule,exchangeable_adapted=TRUE)
        Result_OCP<-OCP(R_cal,R_test,Select_test,alpha)
        Result_ACP<-ACP(R_cal,R_test,T_test,Select_test,alpha,Selection_rule,simplified=TRUE)
        
        
        info<-rbind(info,list(FCR=Result_OCP$FCP,method="OCP",Length=Result_OCP$Length,Thresholding_type=Thresholding_type,Setting=algo@name,n=n,m=m))
        info<-rbind(info,list(FCR=Result_ACP$FCP,method="ACP",Length=Result_ACP$Length,Thresholding_type=Thresholding_type,Setting=algo@name,n=n,m=m))
        info<-rbind(info,list(FCR=Result_SCOP$FCP,method="SCOP",Length=Result_SCOP$Length,Thresholding_type=Thresholding_type,Setting=algo@name,n=n,m=m))
        
  }
  info
  return(info)
}
stopCluster(cl)

## processing the results

Result2=Result[!is.na(Result$Length),]
pp=Result2%>%
  group_by(method,Thresholding_type)%>%
  dplyr::summarize(FCR=mean(FCR),Length=mean(Length))

pp$FCR=pp$FCR*100

pp$Thresholding_type=factor(pp$Thresholding_type,levels = c("test","POS","clu"))
pp$method=factor(pp$method,levels=c("SCOP","OCP","ACP"))

pp2=arrange(pp,method,Thresholding_type)
pp2

