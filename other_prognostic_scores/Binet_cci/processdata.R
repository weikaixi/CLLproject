setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
####below is the binet performance in combined test set####################
###########################################################################
library(survival)
library(survminer)
library(timeROC)
library(rms)
library(regplot)
riskFile="riskTestRIDGE_97.txt"   
cliFile="binetCCI_score_151patients.txt"  
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
fix(risk)
risk=risk[,c("futime", "fustat", "riskScore")]

#
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#cli$age=as.numeric(cli$age)

#
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
fix(rt)

table(rt$binet_group)
#do below process manually to get p values and adj.p value
#rt=rt[which(rt[,"binet_group"]%in%c("A","C")),]
#rt[which(rt[,"binet_group"]%in%c("A","B")),"binet_group"]="A"
#rt$binetsurvival=ifelse(rt$binet>2.5,"high","low")
###
diff=survdiff(Surv(futime, fustat) ~binet_group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(futime, fustat) ~ binet_group, data = rt)
summary(fit)    #²é¿´ÎåÄêÉú´æÂÊ
pdf(file="survivalcombinedTest_binet.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","green","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Binet stage in combined test set (adj.p=0.909)",sep=""),
     mark.time=T)
legend("bottomright", 
       c("Stage C", "Stage B","Stage A"),
       lwd=2,
       col=c("red","green","blue"))
dev.off()





ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=-rt$binet, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="binet_ROCtest_combined_test_66patients.pdf", width=5, height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))
       ),
       col=c("green","blue",'red'),lwd=2,bty = 'n')
dev.off()

library(Hmisc)
test=rt
rcorr.cens(test$binet,Surv(test$futime,test$fustat))
##############below is the performance of binet in the test set 1.#############
##############################################################################
riskFile="77riskTestRIDGE.txt"   
cliFile="binetCCI_score_151patients.txt"  
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
fix(risk)
risk=risk[,c("futime", "fustat", "riskScore")]

#
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#cli$age=as.numeric(cli$age)

#?ϲ?????
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
fix(rt)




ROC_rt1=timeROC(T=rt$futime, delta=rt$fustat,
               marker=-rt$binet, cause=1,
               weighting='aalen',
               times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=-rt$binet, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=-rt$binet, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="binet_ROCtest__test1_49patients.pdf", width=5, height=5)
plot(ROC_rt1,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt2,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt3,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt1$AUC[2])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt2$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt3$AUC[2]))
       ),
       col=c("green","blue",'red'),lwd=2,bty = 'n')
dev.off()

library(Hmisc)
test=rt
rcorr.cens(test$binet,Surv(test$futime,test$fustat))


##############below is the performance of binet in the test set 2.#############
##############################################################################
riskFile="20patient_riskTestRIDGE.txt"   
cliFile="binetCCI_score_151patients.txt"   
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
fix(risk)
risk=risk[,c("futime", "fustat", "riskScore")]

#
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#cli$age=as.numeric(cli$age)

#
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
fix(rt)




ROC_rt1=timeROC(T=rt$futime, delta=rt$fustat,
                marker=-rt$binet, cause=1,
                weighting='aalen',
                times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=-rt$binet, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=-rt$binet, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="binet_ROCtest__test2_17patients.pdf", width=5, height=5)
plot(ROC_rt1,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt2,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt3,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt1$AUC[2])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt2$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt3$AUC[2]))
       ),
       col=c("green","blue",'red'),lwd=2,bty = 'n')
dev.off()

library(Hmisc)
test=rt
rcorr.cens(test$binet,Surv(test$futime,test$fustat))

####below is the CCI performance in combined test set####################
###########################################################################
riskFile="riskTestRIDGE_97.txt"   
cliFile="binetCCI_score_151patients.txt"  
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
fix(risk)
risk=risk[,c("futime", "fustat", "risk","riskScore")]

#
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#cli$age=as.numeric(cli$age)

#
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
fix(rt)
table(rt$risk)
table(rt$cci_group)
#do below process manually to get p values and adj.p value
# rt=rt[which(rt[,"cci_group"]%in%c("mild","moderate")),]
# rt=rt[which(rt[,"cci_group"]%in%c("mild","severe")),]
# rt[which(rt[,"cci_group"]%in%c("moderate","severe")),"cci_group"]="moderate_severe"
#rt$ccisurvival=ifelse(rt$cci>1.5,"high","low")
###36 binet 3 stage as high risk group, 38 binet 1+2 stage as low risk group
diff=survdiff(Surv(futime, fustat) ~cci_group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(futime, fustat) ~ cci_group, data = rt)
summary(fit)    #²é¿´ÎåÄêÉú´æÂÊ
pdf(file="survivalcombinedTest_cci.pdf",width=5.5,height=5)
# Plot the survival curve
plot(fit, 
     lwd = 2,  # Line width
     col = c("blue", "#5fa460","red"),  # Colors for the groups
     xlab = "Time (year)",  # X-axis label
     ylab = "Survival rate",  # Y-axis label
     main = paste("Survival curve (p=", pValue, ")", sep = ""),  # Main title with dynamic p-value
     mark.time = TRUE)  # Show marks at censored times

# Add legend to the plot
legend("bottomright", 
       c("mild group", "moderate group","severe group"),
       lwd = 2,  # Line width
       col = c("blue", "#5fa460","red"),
       cex = 0.7)  # Adjusting cex for legend font size
dev.off()

table(rt$risk)



ROC_rt1=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$cci, cause=1,
                weighting='aalen',
                times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$cci, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$cci, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="CCI_ROCtest__combined_66patients.pdf", width=5, height=5)
plot(ROC_rt1,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt2,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt3,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt1$AUC[2])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt2$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt3$AUC[2]))
       ),
       col=c("green","blue",'red'),lwd=2,bty = 'n')
dev.off()

library(Hmisc)
test=rt
1-rcorr.cens(test$cci,Surv(test$futime,test$fustat))
##############below is the performance of CCI in the test set 1.#############
##############################################################################
riskFile="77riskTestRIDGE.txt"   
cliFile="binetCCI_score_151patients.txt"  
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
fix(risk)
risk=risk[,c("futime", "fustat", "riskScore")]


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#cli$age=as.numeric(cli$age)

#?ϲ?????
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
fix(rt)


ROC_rt1=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$cci, cause=1,
                weighting='aalen',
                times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$cci, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$cci, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="CCI_ROCtest__test1_49patients.pdf", width=5, height=5)
plot(ROC_rt1,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt2,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt3,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt1$AUC[2])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt2$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt3$AUC[2]))
       ),
       col=c("green","blue",'red'),lwd=2,bty = 'n')
dev.off()

library(Hmisc)
test=rt
1-rcorr.cens(test$cci,Surv(test$futime,test$fustat))


##############below is the performance of binet in the test set 2.#############
##############################################################################
riskFile="20patient_riskTestRIDGE.txt"   
cliFile="binetCCI_score_151patients.txt"   
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
fix(risk)
risk=risk[,c("futime", "fustat", "riskScore")]


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#cli$age=as.numeric(cli$age)


samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
fix(rt)


ROC_rt1=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$cci, cause=1,
                weighting='aalen',
                times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$cci, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$cci, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="CCI_ROCtest__test2_17patients.pdf", width=5, height=5)
plot(ROC_rt1,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt2,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt3,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt1$AUC[2])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt2$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt3$AUC[2]))
       ),
       col=c("green","blue",'red'),lwd=2,bty = 'n')
dev.off()

library(Hmisc)
test=rt
1-rcorr.cens(test$cci,Surv(test$futime,test$fustat))


