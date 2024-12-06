setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(survival)  
####below is the 4factor performance in combined test set####################
###########################################################################
riskFile="riskTestRIDGE_97.txt"   
cliFile="4_factor_score_151patients.txt"  
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

table(rt$fourriskscore_group)

#do below process manually to get p values and adj.p value
# rt=rt[which(rt[,"fourriskscore_group"]%in%c("low","high")),]
# 
# rt[which(rt[,"fourriskscore_group"]%in%c("intermediate","high")),"fourriskscore_group"]="intermediate_high"


#rt$binetsurvival=ifelse(rt$fourriskscore>1.5,"high","low")
###36 binet 3 stage as high risk group, 38 binet 1+2 stage as low risk group
diff=survdiff(Surv(futime, fustat) ~fourriskscore_group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(futime, fustat) ~ fourriskscore_group, data = rt)
summary(fit)    #²é¿´ÎåÄêÉú´æÂÊ
pdf(file="survivalcombinedTest_fourriskfactor.pdf",width=5.5,height=5)
# Plot the survival curve
plot(fit, 
     lwd = 2,  # Line width
     col = c("red", "#5fa460","blue"),  # Colors for the groups
     xlab = "Time (year)",  # X-axis label
     ylab = "Survival rate",  # Y-axis label
     main = paste("Survival curve (p=", pValue, ")", sep = ""),  # Main title with dynamic p-value
     mark.time = TRUE)  # Show marks at censored times

# Add legend to the plot
legend("bottomright", 
       c("high risk", "intermediate risk","low risk"),
       lwd = 2,  # Line width
       col = c("red", "#5fa460","blue"),
       cex = 0.7)  # Adjusting cex for legend font size
dev.off()







ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$fourriskscore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="4factor_ROCtest_combined_test_66patients.pdf", width=5, height=5)
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
1-rcorr.cens(test$fourriskscore,Surv(test$futime,test$fustat))
##############below is the performance of 4factor in the test set 1.#############
##############################################################################
riskFile="77riskTestRIDGE.txt"   
cliFile="4_factor_score_151patients.txt"  
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
               marker=rt$fourriskscore, cause=1,
               weighting='aalen',
               times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$fourriskscore, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$fourriskscore, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="4factor_ROCtest__test1_49patients.pdf", width=5, height=5)
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
1-rcorr.cens(test$fourriskscore,Surv(test$futime,test$fustat))

##############below is the performance of 4 factor score in the test set 2.#############
##############################################################################
riskFile="20patient_riskTestRIDGE.txt"   
cliFile="4_factor_score_151patients.txt"  
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
                marker=rt$fourriskscore, cause=1,
                weighting='aalen',
                times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$fourriskscore, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$fourriskscore, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="4factor_ROCtest__test2_17patients.pdf", width=5, height=5)
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
1-rcorr.cens(test$fourriskscore,Surv(test$futime,test$fustat))