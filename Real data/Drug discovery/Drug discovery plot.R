require(tidyverse)
require(ggplot2)
library(reshape2)
library(latex2exp)
library(ggpubr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

### processing the data
data=read.csv("plotdata.csv")
data$FCR=data$FCR*100
alpha=0.1

data$Thresholding_type[data$Thresholding_type=="T-cal"]="T-cal(70%)" 
data$Thresholding_type[data$Thresholding_type=="T-test"]="T-test(70%)" 
data$Thresholding_type[data$Thresholding_type=="T-pos"]="T-pos(9.21,20%)" 


data$Thresholding_type=factor(data$Thresholding_type,levels = c("T-cal(70%)","T-test(70%)","T-pos(9.21,20%)"))
data$method=factor(data$method,levels = c("SCOP","OCP","ACP"))


### plotting

P1<-ggplot(data = data, aes(x = Thresholding_type, y = FCR,fill=method)) +
  geom_boxplot(alpha=0.7,width=0.5) +
  ylab("FCP (%)")+ylim(0,35)+
  scale_x_discrete(name = "") +
  theme_bw() +                
  geom_hline(yintercept = alpha*100,color="black",linetype="dashed",linewidth=0.8)+
  stat_summary(mapping=aes(group=method),                    
               fun="mean",                                   
               geom="point",shape=23,size=1.1,fill="red",    
               position=position_dodge(0.5))+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        legend.position="",axis.text.x = element_text(angle = 0,hjust=0.5))




P2<-ggplot(data = data, aes(x = Thresholding_type, y = Length,fill=method)) +
  geom_boxplot(alpha=0.7,width=0.5) +
  scale_y_continuous(name = "Length")+
  scale_x_discrete(name = "") +
  theme_bw() +                
  stat_summary(mapping=aes(group=method),                    
               fun="mean",                                   
               geom="point",shape=23,size=1.1,fill="red",    
               position=position_dodge(0.5))+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        legend.position="",axis.text.x = element_text(angle = 0,hjust=0.5))


PP=ggarrange(P1,P2)


PP