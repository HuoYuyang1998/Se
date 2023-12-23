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

Selection_rule<-list(Thresholding_type="cons",b_0=-1)
algo<-new('OLS')

beta=rep(1,d)
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
Select_cal<-which(T_cal<=L)


R_cal_sel=R_cal[Select_cal]
R_test_sel=R_test[Select_test]


### Summarized the information of residuals
REwork3=data.frame(residual=c(R_cal,R_cal_sel,R_test_sel),type=c(rep("Marginal",length(R_cal)),rep("Conditional",length(R_cal_sel)),rep("Test",length(R_test_sel))))
REwork3$type=factor(REwork3$type,levels=c("Marginal","Conditional","Test"))

### Plotting

Pd<-REwork3 %>% 
  ggplot(aes(x = residual,y=after_stat(density),fill=type,color=type)) +geom_density(alpha=0.35,linewidth=0.7)+
  scale_x_continuous(name = TeX("$R_i$")) + scale_fill_manual(values=c("#1F78B4","#33A02C","#E31A1C"))+
  scale_color_manual(values=c("#1F78B4","#33A02C","#E31A1C"))+ theme_bw()+
  theme(plot.title = element_text(size = 14, face =  "bold",hjust=0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11),legend.position="",
        legend.text = element_text(size=9),legend.title = element_blank())#+ggtitle("(b)")
Pd
