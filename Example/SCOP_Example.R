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
m<-600 ### size of unlabeled set
N<-n_train+n+m ### size of the overall data
alpha<-0.1 ### target FCR
sigma<-1 ### noise strength



### Choose a selection rule

Selection_rule<-list(Thresholding_type="test",Thresholding_quantile=0.2)
#Selection_rule<-list(Thresholding_type="cal",Thresholding_quantile=0.2)
#Selection_rule<-list(Thresholding_type="exch",Thresholding_quantile=0.2)
#Selection_rule<-list(Thresholding_type="POS",FDR_level=0.2,b_0=-8)
#Selection_rule<-list(Thresholding_type="cons",b_0=-1)
#Selection_rule<-list(Thresholding_type="top",K=60)
#Selection_rule<-list(Thresholding_type="clu")



### Choose a learning algorithm

#algo<-new('OLS')
algo<-new('SVM')

### Generating data
X<-matrix(runif(N*d,-1,1),nrow=N,ncol=d)
beta<-runif(d,-2,2)
data<-DataGen(algo,X,beta,sigma) ## The data generation function is related to Scenario A (OLS) and Scenario B (SVM)


### Splitting data
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


### Train the model
model_train<-fitting(algo,X_train,Y_train)

Y_cal_hat<-predict(model_train,as.data.frame(x=X_cal))
Y_test_hat<-predict(model_train,as.data.frame(x=X_test))

R_cal<-abs(Y_cal-Y_cal_hat)
R_test<-abs(Y_test-Y_test_hat)

### Selection procedure
T_test<-Y_test_hat
T_cal<-Y_cal_hat
L<-Selection_thresholding(Selection_rule,T_cal,T_test,Y_cal)
Select_test<-which(T_test<=L)





### Results for our proposed method, including FCP and average length for constructed prediction intervals
SCOP(R_cal,R_test,T_cal,T_test,Select_test,L,alpha,Selection_rule,exchangeable_adapted=TRUE)
### If the length of SCOP is NA, 
### it means there is no test data or calibration data being selected.

### Results for ordinary conformal prediction
OCP(R_cal,R_test,Select_test,alpha)
  

### Results for adjusted conformal prediction
ACP(R_cal,R_test,T_test,Select_test,alpha,Selection_rule,simplified=TRUE)
  




