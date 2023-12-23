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

z_true_compute<-function(Y_test,X_test){
  z_true=rep(0,length(Y_test))
  for (i in 1:length(Y_test)) {
    sd=1+2*abs(X_test[i,1])
    z_true[i]=qnorm(1-alpha/2,0,sd)
  }
  return(z_true)
}



d<-1 ### dimension of covariates
n<-500 ### size of calibration set
n_train<-200 ### size of training set
n_rest<-n+n_train ### size of labeled set
m<-500 ### size of unlabeled set
N<-n_train+n+m ### size of the overall data
alpha<-0.1 ### target FCR
sigma<-1 ### noise strength
ns<-1000 ### replication times



Selection_rule<-list(Thresholding_type="test",Thresholding_quantile=0.3)
algo<-new('OLS')



beta=2
X<-matrix(runif(N*d,-1,1),nrow=N,ncol=d)
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

Result_SCOP<-SCOP(R_cal,R_test,T_cal,T_test,Select_test,L,alpha,Selection_rule,exchangeable_adapted=TRUE)
Result_OCP<-OCP(R_cal,R_test,Select_test,alpha)
Result_ACP<-ACP(R_cal,R_test,T_test,Select_test,alpha,Selection_rule,simplified=FALSE)




### determine the length of intervals
z_convention=Result_OCP$Length/2
z_adjusted=Result_ACP$Length/2
z_conditional=Result_SCOP$Length/2
z_true=z_true_compute(Y_test,X_test)

data_plot=data.frame(X=X_test[,1],Y=Y_test,hY=Y_test_hat,z_true=z_true)

### determine the plotting parameters
plotalpha=0
linewidth=1.3


### plotting

P2<-ggplot(data = data_plot[T_test<=L,], aes(x = X)) +
  #geom_smooth(aes(y=2*X+z_true),colour="#B15928",se=FALSE,linewidth=1,linetype="longdash")+geom_smooth(aes(y=2*X-z_true),colour="#B15928",se=FALSE,linewidth=1,linetype="longdash")+
  geom_ribbon(aes(ymin=2*X-predict(loess(z_true ~ X)),ymax=2*X+predict(loess(z_true ~ X) )),colour=NA,fill="#B15928",alpha=0.15)+
  
  geom_ribbon(aes(ymin=hY-z_convention,ymax=hY+z_convention),alpha=plotalpha,colour="#33A02C",fill="green",linewidth=linewidth)+
  geom_ribbon(aes(ymin=hY-z_adjusted,ymax=hY+z_adjusted),alpha=plotalpha,colour="#1F78B4",fill="blue",linewidth=linewidth)+
  geom_ribbon(aes(ymin=hY-z_conditional,ymax=hY+z_conditional),alpha=plotalpha,colour="#E31A1C",fill="red",linewidth=linewidth)+
  geom_point(aes(y=Y))+theme_bw()+ylab("Y")

P2



