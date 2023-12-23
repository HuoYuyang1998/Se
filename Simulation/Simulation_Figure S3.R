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

### Functions for LORD-CI
LORD_CI<-function(T_test,R_test,R_cal,L,Lord_gamma,W0)
{n_test=length(T_test)
Selection=rep(0,n_test)
z=rep(0,n_test)
for (i in 1:n_test) {
  Selection[i]=T_test[i]<=L
  if(Selection[i])
  {
    selected_seq=which(Selection==1)
    if(length(selected_seq)==1){alpha_lord=alpha*W0}else
    {alpha_lord=Lord_gamma[i]*W0+(alpha-W0)*Lord_gamma[i-which.min(Selection==1)]+alpha*(sum(Lord_gamma[i-selected_seq[2:length(selected_seq)]]))}
    z[i]=quantile(R_cal,min((1-alpha_lord)*(1+1/length(R_cal)),1))
  }
}
FCP=1-sum(R_test[Selection==1]<=z[Selection==1])/max(1,sum(Selection))
length=mean(2*z[Selection==1])

return(list(FCP=FCP,Length=length,method="LORD-CI"))
}

d<-10 ### dimension of covariates
n<-200 ### size of calibration set
n_train<-200 ### size of training set
n_rest<-n+n_train ### size of labeled set
m<-200 ### size of unlabeled set
N<-n_train+n+m ### size of the overall data
alpha<-0.1 ### target FCR
sigma<-1 ### noise strength
ns<-1000 ### replication times
Lord_gamma=0.0722*sapply(1:m,function(j){log(max(j,2))/(j*exp(sqrt(log(j))))}) ### used in LORD-CI
W0=alpha/2  ### used in LORD-CI


Thresholding_type_array=c("cal","test","exch") ### the selection rules
Thresholding_value_array=seq(0.2,1,0.1) ### the quantile level
algo_array=c(new('OLS'),new('SVM')) ### algorithms for predictions and corresponding scenarios



### parallel computing
cl = makeCluster(10) ### number of parallel threads
registerDoParallel(cl)
Result<-foreach(iter=1:ns,.combine="rbind",.packages = c('MASS',"randomForest","kernlab"))%dopar% {
  info<-data.frame()
  
  for (algo in algo_array) {
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
      for (Thresholding_value in Thresholding_value_array) {
        
        Selection_rule<-list(Thresholding_type=Thresholding_type,Thresholding_quantile=Thresholding_value)
        
        L<-Selection_thresholding(Selection_rule,T_cal,T_test,Y_cal)
        Select_test<-which(T_test<=L)
        
        
        
        Result_SCOP<-SCOP(R_cal,R_test,T_cal,T_test,Select_test,L,alpha,Selection_rule,exchangeable_adapted=TRUE)
        Result_OCP<-OCP(R_cal,R_test,Select_test,alpha)
        Result_ACP<-ACP(R_cal,R_test,T_test,Select_test,alpha,Selection_rule,simplified=TRUE)
        Result_LORDCI<-LORD_CI(T_test,R_test,R_cal,L,Lord_gamma,W0)

        
        info<-rbind(info,list(FCR=Result_OCP$FCP,method="OCP",Length=Result_OCP$Length,Thresholding_type=Thresholding_type,Thresholding_value=Thresholding_value,Setting=algo@name,n=n,m=m))
        info<-rbind(info,list(FCR=Result_ACP$FCP,method="ACP",Length=Result_ACP$Length,Thresholding_type=Thresholding_type,Thresholding_value=Thresholding_value,Setting=algo@name,n=n,m=m))
        info<-rbind(info,list(FCR=Result_SCOP$FCP,method="SCOP",Length=Result_SCOP$Length,Thresholding_type=Thresholding_type,Thresholding_value=Thresholding_value,Setting=algo@name,n=n,m=m))
        info<-rbind(info,list(FCR=Result_LORDCI$FCP,method="LORDCI",Length=Result_LORDCI$Length,Thresholding_type=Thresholding_type,Thresholding_value=Thresholding_value,Setting=algo@name,n=n,m=m))
        
      }
    }
  }
  info
  return(info)
}
stopCluster(cl)

## processing the results


Result$method<-factor(Result$method,levels=c("SCOP","OCP","ACP","LORDCI"))
Result$Thresholding_type[Result$Thresholding_type=="cal"]<-"T-cal"
Result$Thresholding_type[Result$Thresholding_type=="test"]<-"T-test"
Result$Thresholding_type[Result$Thresholding_type=="exch"]<-"T-exch"
Result$Thresholding_type<-factor(Result$Thresholding_type,levels=c("T-cal","T-test","T-exch"))

Result$Setting[Result$Setting=="OLS"]<-"Scenario A"
Result$Setting[Result$Setting=="SVM"]<-"Scenario B"
Result$Setting<-factor(Result$Setting,levels=c("Scenario A","Scenario B"))

Result$FCR<-100*Result$FCR


pp<-Result%>%
  group_by(method,Thresholding_type,Thresholding_value,Setting)%>%
  dplyr::summarize(FCR=mean(FCR),Length=mean(Length,na.rm=TRUE))


pp$Thresholding_value<-pp$Thresholding_value*100

### plotting

P<-pp %>% 
  ggplot(aes(x = Thresholding_value, y = FCR, group = method)) +theme_bw()+
  geom_line(aes(linetype=method,color=method),linewidth=0.4) + geom_point(aes(shape=method,color=method),size=1.2)+
  scale_linetype_manual(values=c(1,3,4,5))+
  scale_y_continuous(name="FCR(%)",position="left")+scale_x_continuous(name = "q%",breaks = c(30,60,90))+
  geom_hline(yintercept = alpha*100,color="black",linetype="dashed",size=0.4) +
  facet_grid(Setting~Thresholding_type ,scales="free")+theme(legend.position="")


P2<-pp %>% 
  ggplot(aes(x = Thresholding_value, y = Length, group = method))+theme_bw() +
  geom_line(aes(linetype=method,color=method),linewidth=0.4) + geom_point(aes(shape=method,color=method),size=1.2)+
  scale_linetype_manual(values=c(1,3,4,5))+
  scale_y_continuous(name="Length",position="left")+scale_x_continuous(name = "q%",breaks = c(30,60,90))+
  facet_grid(Setting~Thresholding_type ,scales="free")+theme(legend.position="")


PP<-ggarrange(P,P2,ncol = 2)
PP