#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import random
import math
import DeepPurpose.DTI as models
from DeepPurpose.dataset import *
from DeepPurpose import utils, dataset
from DeepPurpose import DTI as models


### loading the data and log-scale the response
X_drugs, X_targets, y = dataset.read_file_training_dataset_drug_target_pairs('Davis3.txt')

for i in range(len(y)):
    y[i]=math.log(y[i])





### defining some functions

def list_select(mylist,selection):
    newlist=list()
    for a in range(len(mylist)):
        if selection[a]:
            newlist.append(mylist[a])
    return(newlist)

def PIresults(res_cal,res_test,alpha):
    res_thresh=np.percentile(res_cal,(1-alpha)*100)
    FCR=1-np.mean(res_test<=res_thresh)
    Length=2*res_thresh
    return([FCR,Length])

def POS(pred_cal,pred_test,beta,b0):
    pred_cal_null=list_select(pred_cal,pred_cal>=b0)
    cp=[0]*len(pred_test)
    for i in range(len(pred_test)):
        for j in range(len(pred_cal_null)):
            cp[i]=cp[i]+(pred_cal_null[j]<=pred_test[i])
        cp[i]=(cp[i]+1)/(len(pred_cal_null)+1)
    cp_sort=sorted(cp)
    sort_index=np.argsort(cp)
    index=np.min(np.where(np.array(cp_sort)/(np.array(range(len(cp)))+1)*len(cp)<=beta))
    return(pred_test[sort_index[index]])

def NewSCOP(thresh,pred_test,pred_cal,res_cal,res_test_ordered):
    selection_cal=pred_cal>=thresh
    sort_index=np.argsort(pred_test)
    pred_test_ordered=np.array(pred_test)[sort_index]
    m=len(pred_test)

    kappa=len(pred_test)-1
    while(pred_test_ordered[kappa]>=thresh):
        kappa=kappa-1
    cali_j=list_select(res_cal,np.array(pred_cal)>=pred_test_ordered[max(kappa-1,0)])
    z_j=np.percentile(cali_j,min((1-alpha)*(1+1/len(cali_j))*100,100))

    res_cal_select=list_select(res_cal,selection_cal)
    res_thresh=np.percentile(res_cal_select,(1-alpha)*100)



    FCR=1-(sum(res_test_ordered[range(kappa+1,m)]<=z_j)+np.array(res_test_ordered[kappa]<=z_j))/(m-kappa)
    Length=2*(res_thresh*((m-kappa)-1)+z_j)/(m-kappa)
    return([FCR,Length])



### determine the parameters

alpha=0.1
beta=0.2
b0=np.float64(9.210340371976184)
df=pd.DataFrame()


drug_encoding, target_encoding = 'Morgan', 'Conjoint_triad'

config = utils.generate_config(drug_encoding = drug_encoding, 
                         target_encoding = target_encoding, 
                         cls_hidden_dims = [1024,1024,512], 
                         train_epoch = 5, 
                         LR = 0.001, 
                         batch_size = 128,
                         hidden_dim_drug = 128,
                         mpnn_hidden_size = 128,
                         mpnn_depth = 3, 
                         cnn_target_filters = [32,64,96],
                         cnn_target_kernels = [4,8,12]
                        )


### 100 repetitions for three methods and three selection rules


for ol in range(100):

    train, val, test = utils.data_process(X_drugs, X_targets, y, 
                                drug_encoding, target_encoding, 
                                split_method='random',frac=[0.7,0.1448,0.1552])

    model = models.model_initialize(**config)
    model.train(train, val, test)

    y_pred = model.predict(test)
    caldex=random.sample(range(4000),2000)
    testdex=list(set(range(4000))-set(caldex))
    calset=test.iloc[caldex]
    testset=test.iloc[testdex]
    calset.index=range(2000)
    testset.index=range(2000)

    pred_cal=model.predict(calset)
    pred_test=model.predict(testset)
    sort_index=np.argsort(pred_test)
    
    
    
    thresh=np.percentile(pred_test,70)
    selection_test=pred_test>=thresh
    selection_cal=pred_cal>=thresh
    res_cal=abs(pred_cal-calset["Label"])
    res_test=abs(pred_test-testset["Label"])

    res_cal_select=list_select(res_cal,selection_cal)
    res_test_select=list_select(res_test,selection_test)

    res_test_ordered=np.array(res_test)[sort_index]
    
    
    R=NewSCOP(thresh,pred_test,pred_cal,res_cal,res_test_ordered)
    df=df.append(pd.DataFrame({"FCR":[R[0]],"Length":[R[1]],"method":["SCOP"],"Thresholding_type":["T-test"]}),ignore_index=True)
    R=PIresults(res_cal,res_test_select,alpha)
    df=df.append(pd.DataFrame({"FCR":[R[0]],"Length":[R[1]],"method":["OCP"],"Thresholding_type":["T-test"]}),ignore_index=True)
    R=PIresults(res_cal,res_test_select,alpha*sum(selection_test)/len(pred_test))
    df=df.append(pd.DataFrame({"FCR":[R[0]],"Length":[R[1]],"method":["ACP"],"Thresholding_type":["T-test"]}),ignore_index=True)

    thresh=np.percentile(calset["Label"],70)
    selection_test=pred_test>=thresh
    selection_cal=pred_cal>=thresh
    res_cal=abs(pred_cal-calset["Label"])
    res_test=abs(pred_test-testset["Label"])


    res_cal_select=list_select(res_cal,selection_cal)
    res_test_select=list_select(res_test,selection_test)

    R=NewSCOP(thresh,pred_test,pred_cal,res_cal,res_test_ordered)
    df=df.append(pd.DataFrame({"FCR":[R[0]],"Length":[R[1]],"method":["SCOP"],"Thresholding_type":["T-cal"]}),ignore_index=True)
    R=PIresults(res_cal,res_test_select,alpha)
    df=df.append(pd.DataFrame({"FCR":[R[0]],"Length":[R[1]],"method":["OCP"],"Thresholding_type":["T-cal"]}),ignore_index=True)
    R=PIresults(res_cal,res_test_select,alpha*sum(selection_test)/len(pred_test))
    df=df.append(pd.DataFrame({"FCR":[R[0]],"Length":[R[1]],"method":["ACP"],"Thresholding_type":["T-cal"]}),ignore_index=True)

    thresh=np.float64(POS(pred_cal,pred_test,beta,b0))
    selection_test=pred_test>=thresh
    selection_cal=pred_cal>=thresh
    res_cal=abs(pred_cal-calset["Label"])
    res_test=abs(pred_test-testset["Label"])


    res_cal_select=list_select(res_cal,selection_cal)
    res_test_select=list_select(res_test,selection_test)

    R=NewSCOP(thresh,pred_test,pred_cal,res_cal,res_test_ordered)
    df=df.append(pd.DataFrame({"FCR":[R[0]],"Length":[R[1]],"method":["SCOP"],"Thresholding_type":["T-pos"]}),ignore_index=True)
    R=PIresults(res_cal,res_test_select,alpha)
    df=df.append(pd.DataFrame({"FCR":[R[0]],"Length":[R[1]],"method":["OCP"],"Thresholding_type":["T-pos"]}),ignore_index=True)
    R=PIresults(res_cal,res_test_select,alpha*sum(selection_test)/len(pred_test))
    df=df.append(pd.DataFrame({"FCR":[R[0]],"Length":[R[1]],"method":["ACP"],"Thresholding_type":["T-pos"]}),ignore_index=True)



### the results are summarized in df

df


#df.to_csv("plotdata.csv")






