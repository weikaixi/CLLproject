setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("allmodels.RData")
library(survival)
library(caret)
library("glmnet")
library(survminer)

#below is using fPM model to do prediction
rt=read.table("demo_data.csv",header=T,sep=",",row.names = 1)
fix(rt)
rt$Vandetanib=ifelse(rt$Vandetanib>7,1,0) #the median value of Vandetanib DSS values is 7. 
#detailed calculation process of DSS values in drug response curve can be found in https://github.com/IanevskiAleksandr/Breeze/blob/master/dss_scores.R

#in Bcl.2.Btk_pY223_Itk_pY180 column, if the level of Bcl2 is greater than Btk_pY223_Itk_pY180
#then the values are 1,  otherwise the values are 0.

#in p90RSK_pS380.PLC.gamma2_pY759 column, if the level of p90RSK_pS380 is greater than PLC.gamma2_pY759
#then the values are 1,  otherwise the values are 0.
rt$Ruxolitinib_Venetoclax=ifelse(rt$Ruxolitinib_Venetoclax>17.35,1,0) #the median value of Ruxolitinib_Venetoclax DSS values in training set is 17.35. 
rt$Acalabrutinib_ZSTK474=ifelse(rt$Acalabrutinib_ZSTK474>10.6,1,0) #the median value of Ruxolitinib_Venetoclax DSS values in training set is 10.6. 
rt$CD8plus4CD3_T_cells=ifelse(rt$CD8plus4CD3_T_cells>49.9385,1,0) #the median percentage value of CD8plus4CD3_T_cells in training set is 49.9385. 
#this column means the proportion of CD8+ T cells in  CD3+ T cells
fPM=predict(ridge_model, s = optimal_lambda,newx = as.matrix(rt), type = "response")
riskgroup=ifelse(fPM>0.36173,"highrisk","lowrisk")
#predict(ridge_model, s = optimal_lambda,newx = as.matrix(rt), type = "link")
#you can use either  type = "response" or  type = "link", if you choose "response", the median value is 0.36173
#if you choose type = "link", the median value is -1.016857.The different choice will not influence the ranking of risk score.
outputtable=cbind(rt,fPM,riskgroup)
colnames(outputtable)[c(7,8)]=c("fPM","riskgroup")
outputtable
write.table(outputtable, file = "prediction_result.csv", quote = FALSE, sep = ",", row.names = T)


