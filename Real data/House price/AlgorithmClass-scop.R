library(methods)


setGeneric("fitting",function(obj,...) standardGeneric("fitting"))
setGeneric("Pred",function(obj,...) standardGeneric("Pred"))
setGeneric("DataGen",function(obj,...) standardGeneric("DataGen"))


###ridge
setClass("ridge",slots = list(name="character",alpha="numeric",family="character"), prototype=list(name="RR",alpha=0,family="gaussian") )
setMethod("fitting","ridge",function(obj,X,Y){
  glmnet(x=X,y=Y,family =obj@family,alpha=obj@alpha)
})
setMethod("Pred","ridge",function(obj,model,X_test,lens=0){
  predict(model,X_test)
})
setMethod("DataGen","ridge",function(obj,X,beta,sigma){
  Y=X%*%beta+rnorm(dim(X)[1],0,sigma)
  return(data.frame(x=X,y=Y))
})

###lasso
setClass("lasso",slots = list(name="character",alpha="numeric",family="character"), prototype=list(name="Lasso",alpha=1,family="gaussian") )
setMethod("fitting","lasso",function(obj,X,Y){
  glmnet(x=X,y=Y,family =obj@family,alpha=obj@alpha)
})
setMethod("Pred","lasso",function(obj,model,X_test,lens=0){
  predict(model,X_test)
})
setMethod("DataGen","lasso",function(obj,X,beta,sigma){
  b=c(1,-1,2,-2,rep(0,dim(X)[2]-4))
  Y=X%*%b+rnorm(dim(X)[1],0,sigma)
  return(data.frame(x=X,y=Y))
})


###random forest
setClass("RF",slots = list(name="character"),prototype = list(name="RF"))
setMethod("fitting","RF",function(obj,X,Y){
  datawork=data.frame(X,y=Y)
  randomForest(y~.,data=datawork,mtry=round(dim(X)[2]/3),ntree=500)
})
setMethod("Pred","RF",function(obj,model,X_test,lens=0){
  predict(model,as.data.frame(x=X_test))})
setMethod("DataGen","RF",function(obj,X,beta,sigma){
  uX=4*(X[,1]+1)*abs(X[,3])*(X[,2]>-0.4)+4*(X[,1]-1)*(X[,2]<=-0.4)
  Y=uX+rnorm(dim(X)[1],0,sigma)
  return(data.frame(x=X,y=Y))
})

###MLP
setClass("NN",slots = list(name="character",stepmax="numeric",rep="numeric",hidden="numeric"),
         prototype = list(name="NN",stepmax=10000,rep=1,hidden=20))
setMethod("fitting","NN",function(obj,X,Y){
  datawork=data.frame(X,y=Y)
  neuralnet(y~.,data=datawork,hidden=obj@hidden,stepmax =obj@stepmax,rep=obj@rep,learningrate=0.005,act.fct = "tanh")
})
setMethod("Pred","NN",function(obj,model,X_test,lens=0){
  predict(model,X_test)})
setMethod("DataGen","NN",function(obj,X,beta,sigma){
  Y=X%*%beta+rnorm(dim(X)[1],sigma)
  return(data.frame(x=X,y=Y))
})

###SVM
setClass("SVM",slots = list(name="character"),
         prototype = list(name="SVM"))
setMethod("fitting","SVM",function(obj,X,Y){
  datawork=data.frame(X,y=Y)
  ksvm(y~.,data=datawork,C=1)
})
setMethod("Pred","SVM",function(obj,model,X_test,lens,type="decision"){
  predraw=predict(model,X_test)
  return(predraw)
  })
setMethod("DataGen","SVM",function(obj,X,beta,sigma){
  uX=X[,1]*X[,2]+X[,3]-2*exp(X[,4]+1)
  Y=uX+rnorm(dim(X)[1],0,1)
  return(data.frame(x=X,y=Y))
})

###ordinary least square regression
setClass("OLS",slots = list(name="character"), prototype=list(name="OLS") )
setMethod("fitting","OLS",function(obj,X,Y){
  if(dim(X)[1]>=dim(X)[2])
  {datawork=data.frame(X,y=Y)
  return(lm(y~.,datawork))}else
  {X_a1=cbind(rep(1,dim(X)[1]),X)
  return(ginv(t(X_a1)%*%X_a1)%*%t(X_a1)%*%Y)}
})
setMethod("Pred","OLS",function(obj,model,X_test,lens=0){

  if(class(model)=="lm")
  {predict(model,as.data.frame(x=X_test))}else
  {X_a1=cbind(rep(1,dim(X_test)[1]),X_test)
  X_a1%*%model}
})
setMethod("DataGen","OLS",function(obj,X,beta,sigma){
  uX=X%*%beta
  Y=uX+sapply(1:dim(X)[1],function(i){rnorm(1,0,sigma+abs(uX[i]))})#Éú³ÉY
  return(data.frame(x=X,y=Y))
})


###Mean

setClass("Cons",slots = list(name="character"), prototype=list(name="Constant-max") )
setMethod("fitting","Cons",function(obj,X,Y){
return(max(Y))
})
setMethod("Pred","Cons",function(obj,model,X_test,lens=0){
return(rep(model,dim(X_test)[1]))
})
setMethod("DataGen","Cons",function(obj,X,beta,sigma){
  Y=X%*%beta+rnorm(dim(X)[1],sigma)
  return(data.frame(x=X,y=Y))
})