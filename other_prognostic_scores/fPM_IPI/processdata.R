setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
####below is the IPI performance in combined test set####################
###########################################################################
library(survival)
library(survminer)
library(timeROC)
library(rms)
library(regplot)
riskFile="riskTestRIDGE_97.txt"   
cliFile="ipi_score_151patients.txt"  
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
fix(risk)
risk=risk[,c("futime", "fustat", "riskScore")]

#??ȡ?ٴ???????Score??
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#cli$age=as.numeric(cli$age)

#?ϲ?????
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
fix(rt)

table(rt$IPI_group)
median(rt$IPI)

#using below codes to calculate p values between different groups manually, and then calculate the adj.p
# rt=rt[which(rt[,"IPI_group"]%in%c("low","very_high")),]
# rt[which(rt[,"IPI_group"]%in%c("high","very_high")),"IPI_group"]="very_high_high"
# rt[which(rt[,"IPI_group"]%in%c("intermediate","low")),"IPI_group"]="low_intermediate"
# #rt$IPIsurvival=ifelse(rt$IPI>3,"high","low")
diff=survdiff(Surv(futime, fustat) ~IPI_group,data = rt)
fit <- survfit(Surv(futime, fustat) ~ IPI_group, data = rt)
###below is new size
# Convert dimensions from mm to inches (1 inch = 25.4 mm)
width_in_inches <- 67 / 25.4
height_in_inches <- 67 / 25.4

# Set up the PDF device with specified width and height in inches
pdf(file = "newsize_survivalcombinedTest_IPI.pdf", width = width_in_inches, height = height_in_inches)

# Set font family to Helvetica and font size to 7
par(family = "Helvetica", cex = 0.7)  # Adjusting cex and ps to set font size appropriately

# Plot the survival curve
plot(fit, 
     lwd = 2,  # Line width
     col = c("red","#5fa460","blue","#a45a79"),  # Colors for the groups
     xlab = "Time (year)",  # X-axis label
     ylab = "Survival rate",  # Y-axis label
     main = paste("Survival curve (adj.p=0.679)", sep = ""),  # Main title with dynamic p-value
     mark.time = TRUE)  # Show marks at censored times

# Add legend to the plot
legend("bottomright", 
       c("high risk","intermediate risk", "low risk","very high risk"),
       lwd = 2,  # Line width
       col = c("red","#5fa460","blue","#a45a79"),
       cex = 0.7)  # Adjusting cex for legend font size

# Close the PDF device
dev.off()



ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$IPI, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="IPI_ROCtest_combined_test_66patients.pdf", width=5, height=5)
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



library(survival)
# Assuming you have a dataframe 'data' with 'futime', 'fustat', and 'riskScore'

library(Hmisc)
test=rt
1-rcorr.cens(test$IPI,Surv(test$futime,test$fustat))


# Print the uno C-index in test set
# library(survAUC)
# lpnew <- test$riskScore
# Surv.rsp <- survival::Surv(train$futime, train$fustat)
# Surv.rsp.new <- survival::Surv(test$futime, test$fustat)             
# Cstat <- UnoC(Surv.rsp, Surv.rsp.new, lpnew)
# Cstat

##############below is the performance of IPI in the test set 1.#############
##############################################################################
riskFile="77riskTestRIDGE.txt"   
cliFile="ipi_score_151patients.txt"  
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
               marker=rt$IPI, cause=1,
               weighting='aalen',
               times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$IPI, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$IPI, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="IPI_ROCtest__test1_49patients.pdf", width=5, height=5)
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

library(survival)
# Assuming you have a dataframe 'data' with 'futime', 'fustat', and 'riskScore'

library(Hmisc)
test=rt
1-rcorr.cens(test$IPI,Surv(test$futime,test$fustat))

##############below is the performance of IPI in the test set 2.#############
##############################################################################
riskFile="20patient_riskTestRIDGE.txt"   
cliFile="ipi_score_151patients.txt"  
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
                marker=rt$IPI, cause=1,
                weighting='aalen',
                times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$IPI, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$IPI, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="IPI_ROCtest__test2_17patients.pdf", width=5, height=5)
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
1-rcorr.cens(test$IPI,Surv(test$futime,test$fustat))



####below is the fPM performance in combined test set####################
###########################################################################
riskFile="riskTestRIDGE_97.txt"   
cliFile="ipi_score_151patients.txt"  
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




ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="fPM_ROCtest_combined_test_66patients.pdf", width=5, height=5)
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
1-rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))
##############below is the performance of MFP in the test set 1.#############
##############################################################################
riskFile="77riskTestRIDGE.txt"   
cliFile="ipi_score_151patients.txt"  
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
                marker=rt$riskScore, cause=1,
                weighting='aalen',
                times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$riskScore, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$riskScore, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="fPM_ROCtest__test1_49patients.pdf", width=5, height=5)
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
1-rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))

##############below is the performance of fPM in the test set 2.#############
##############################################################################
riskFile="20patient_riskTestRIDGE.txt"   
cliFile="ipi_score_151patients.txt"  
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
                marker=rt$riskScore, cause=1,
                weighting='aalen',
                times=c(1), ROC=TRUE)
ROC_rt2=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$riskScore, cause=1,
                weighting='aalen',
                times=c(2), ROC=TRUE)
ROC_rt3=timeROC(T=rt$futime, delta=rt$fustat,
                marker=rt$riskScore, cause=1,
                weighting='aalen',
                times=c(3), ROC=TRUE)
pdf(file="fPM_ROCtest__test2_17patients.pdf", width=5, height=5)
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
1-rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))

