#### Main functions 

SCOP<-function(R_cal,R_test,T_cal,T_test,Select_test,L,alpha,Selection_rule,exchangeable_adapted=TRUE){
###  Using selective conditional conformal prediction to give prediction intervals after selection
  Select_cal<-which(T_cal<=L)
  if(length(Select_test)==0|length(Select_cal)==0){
    return(data.frame(FCP=0,Length=NA,method="SCOP"))
  }
  Thresholding_type=Selection_rule$Thresholding_type
  
  if(Thresholding_type=="cons"|Thresholding_type=="exch"|!exchangeable_adapted){### The constant threshold and exch threshold does not require exchangeable-adaptation
    R_cal_rej<-R_cal[Select_cal]
    z_conditional<-quantile(R_cal_rej,min((1-alpha)*(1+1/length(R_cal_rej)),1))
    Len<-2*z_conditional
    FCP<-1-(sum(R_test[Select_test]<=z_conditional))/length(Select_test)
  }else{
    T_ordered<-T_test[order(T_test)]
    kappa<-max(which(T_ordered<=L))
    
    m<-length(R_test)
    kappa_index<-order(T_test)[kappa]
    cali_j<-R_cal[which(T_cal<=T_ordered[min(kappa+1,m)])]
    z_conditional=quantile(cali_j,min((1-alpha)*(1+1/length(cali_j)),1))
    
    Len<-2*z_conditional
    FCP<-1-(sum(R_test[Select_test]<=z_conditional))/length(Select_test)
    
    
    
    
    
  }
  
  
  return(data.frame(FCP=FCP,Length=Len[[1]],Method="SCOP"))
}







#### Other functions

DataSplit<-function(data,N,m,n,n_rest)
{### splitting the generated data into training, calibrating and test datasets
  if(m>0)
  {  index_test<-sample(1:N,m,replace=FALSE)
  data_test<-data[index_test,]
  data_rest2<-data[-index_test,]
  data_rest<-data_rest2[sample(1:dim(data_rest2)[1],n_rest),]
  index_cal<-sample(1:dim(data_rest)[1],n)
  data_train<-data_rest[-index_cal,]
  data_cal<-data_rest[index_cal,]
  return(list(data_train=data_train,data_cal=data_cal,data_test=data_test,data_rest=data_rest))}else
  {
    data_rest<-data[sample(1:dim(data)[1],n_rest,replace=FALSE),]
    index_cal<-sample(1:dim(data_rest)[1],n)
    data_train<-data_rest[-index_cal,]
    data_cal<-data_rest[index_cal,]
    return(list(data_train=data_train,data_cal=data_cal,data_test=0,data_rest=data_rest))
  }
  
}



POS_threshold<-function(T_test,T_cal,Y_cal,b_0,beta)
{# computing the threshold by selection rule T-pos(b_0,FDRlevel)
  Null_cal<-which(Y_cal>b_0)
  n0<-length(Null_cal)
  pvalues<-sapply(T_test, function(t){
    (sum(T_cal[Null_cal]<=t)+1)/(n0+1)
  })
  m0<-n0/length(Y_cal)*length(T_test)
  L_seq<-seq(0,0.7,0.7/5000)
  FDP_select<-sapply(L_seq, function(t){
    m0*t/max(1,sum(pvalues<=t))
  })
  L<-L_seq[max(which(FDP_select<=beta))]
  return(  max(T_test[pvalues<=L]))
}


Clustering_threshold<-function(T_test,T_cal)
{# computing the threshold by selection rule T-clu
  T_sort<-sort(c(T_cal,T_test))
  Tlen<-length(T_sort)
  vaT<-sapply(2:(Tlen-2), function(t){
    T1<-var(T_sort[1:t])
    T2<-var(T_sort[(t+1):Tlen])
    p<-t/Tlen
    return(p*T1+(1-p)*T2)
  })
  return(T_sort[which.min(vaT)])
}





Selection_thresholding<-function(Selection_rule,T_cal=0,T_test=0,Y_cal=0){
# computing the threshold by assigning a selection rule
  Thresholding_type<-Selection_rule$Thresholding_type
  
  
  if(Thresholding_type=="cons")
  {L<-Selection_rule$b_0}else if(Thresholding_type=="test")
  {L<-quantile(T_test,Selection_rule$Thresholding_quantile)}else if(Thresholding_type=="cal")
  {L<-quantile(Y_cal,Selection_rule$Thresholding_quantile)}else if(Thresholding_type=="exch")
  {L<-quantile(c(T_cal,T_test),Selection_rule$Thresholding_quantile)}else if (Thresholding_type=="POS")
  {L<-POS_threshold(T_test,T_cal,Y_cal,Selection_rule$b_0,Selection_rule$FDR_level)}else if (Thresholding_type=="clu")
  {L<-Clustering_threshold(T_test,T_cal)}else if (Thresholding_type=="top")
  {L<-T_test[order(T_test)][Selection_rule$K]}
    
  return(L)
}



OCP<-function(R_cal,R_test,Select_test,alpha){
# using ordinary conformal prediction to give prediction intervals after selection
  z_convention<-quantile(R_cal,min((1-alpha)*(1+1/length(R_cal)),1))
  Len<-2*z_convention
  if(length(Select_test)==0){
    FCP<-0
  }else{
    FCP<-1-(sum(R_test[Select_test]<=z_convention))/length(Select_test)
    
  }
  return(data.frame(FCP=FCP,Length=Len[[1]],Method="OCP"))
  
}

ACP<-function(R_cal,R_test,T_test,Select_test,alpha,Selection_rule,simplified=TRUE){
# using adjusted conformal prediction to give prediction intervals after selection
#  If the parameter "simplified=TRUE", we will use |\hat{\mathcal{S}}_u| as approximation of \mathcal{M}^j_{\min}.
  Thresholding_type<-Selection_rule$Thresholding_type
  
  if(Thresholding_type=="test"&!simplified){
    minimal_value<-min(T_test)-100
    Threshold<-sapply(1:length(Select_test), function(i){
      set_size<-sum(T_test<=quantile(c(minimal_value,T_test[-Select_test[i]]),Selection_rule$Thresholding_quantile))
      z_adjusted<-quantile(R_cal,min((1-alpha*set_size/length(T_test))*(1+1/length(R_cal)),1))
      return(z_adjusted)
    })
    Len<-2*mean(Threshold)
    if(length(Select_test)==0){FCP<-0}else{
      FCP<-1-(sum(R_test[Select_test]<=Threshold))/length(Select_test)
    }
  }else{
    
    z_adjusted<-quantile(R_cal,min((1-alpha*length(Select_test)/length(T_test))*(1+1/length(R_cal)),1))
    Len<-2*z_adjusted
    if(length(Select_test)==0){FCP<-0}else{
      FCP<-1-(sum(R_test[Select_test]<=z_adjusted))/length(Select_test)
    }
  }
  
  return(data.frame(FCP=FCP,Length=Len[[1]],Method="ACP"))
}











