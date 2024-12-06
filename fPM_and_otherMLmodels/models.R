setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
load("allmodels.RData")
library(survival)
library(caret)
library("glmnet")
library(survminer)
library("survivalROC")
library(timeROC)
library(randomForestSRC)

#below is fPM model doing prediction
rt=read.table("20patient_riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
fix(rt)
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    #²é¿´ÎåÄêÉú´æÂÊ
pdf(file="test2_survivalTestRIDGE_link.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("bottomright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()




###below is the new figure size of fPM in validation set#######################################
# Convert dimensions from mm to inches (1 inch = 25.4 mm)
width_in_inches <- 84 / 25.4
height_in_inches <- 84 / 25.4

# Set up the PDF device with specified width and height in inches
pdf(file = "newsize_test2_survivalTestRIDGE_link.pdf", width = width_in_inches, height = height_in_inches)

# Set font family to Helvetica and font size to 7
par(family = "Helvetica", cex = 0.7)  # Adjusting cex to scale font size appropriately

# Plot the survival curve
plot(fit, 
     lwd = 2,  # Line width
     col = c("red", "blue"),  # Colors for the groups
     xlab = "Time (year)",  # X-axis label
     ylab = "Survival rate",  # Y-axis label
     main = paste("Survival curve (p=", pValue, ")", sep = ""),  # Main title with dynamic p-value
     mark.time = TRUE)  # Show marks at censored times

# Add legend to the plot
legend("bottomright", 
       c("high risk", "low risk"),
       lwd = 2,  # Line width
       col = c("red", "blue"))  # Colors corresponding to the plot lines

# Close the PDF device
dev.off()

##below is to draw ROC curves of fPM score in validaiton set.
rt=read.table("20patient_riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
fix(rt)

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="test2ROC.pdf", width=5, height=5)

# First plot call without the title parameter
plot(ROC_rt, time=1, col='green', lwd=2, title=FALSE)

# Adding a title using the title function
title("MFP score in test set 2")

# The rest of the plotting code remains the same
plot(ROC_rt, time=2, col='blue', add=TRUE, title=FALSE, lwd=2)
plot(ROC_rt, time=3, col='red', add=TRUE, title=FALSE, lwd=2)
legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]),"*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]),"*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]),"*")
       ),
       col=c("green", "blue", 'red'), lwd=2, bty = 'n')

dev.off()

####below is new size figure of fPM in validation set######ROC curves######################
# Set the width and height in millimeters
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file = "newsize_test2ROC_fPM.pdf", width = 84/25.4, height = 84/25.4) ##test 2 set is validation set in the mauscript

# Set font to Helvetica with size 7
par(family = "Helvetica", cex = 0.7)

# First plot call without the title parameter
plot(ROC_rt, time = 1, col = 'green', lwd = 2, title = FALSE)

# Adding a title using the title function
title("fPM score in validation set")

# The rest of the plotting code remains the same
plot(ROC_rt, time = 2, col = 'blue', add = TRUE, title = FALSE, lwd = 2)
plot(ROC_rt, time = 3, col = 'red', add = TRUE, title = FALSE, lwd = 2)
legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]), "*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]), "*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]), "*")
       ),
       col = c("green", "blue", 'red'), lwd = 2, bty = 'n')

# Close the PDF device
dev.off()


###below to calculate Cindex
library(survival)
# Assuming  have a dataframe 'data' with 'futime', 'fustat', and 'riskScore'
library(Hmisc)
test=rt
1-rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))


# Print the uno C-index in test set
library(survAUC)
lpnew <- test$riskScore
Surv.rsp <- survival::Surv(train$futime, train$fustat)
Surv.rsp.new <- survival::Surv(test$futime, test$fustat)             
Cstat <- UnoC(Surv.rsp, Surv.rsp.new, lpnew)
Cstat



##################below is the result of test set.
########################################################
rt=read.table("77riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
fix(rt)
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    
pdf(file="test1_survivalTestRIDGE.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("bottomright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()



# Convert dimensions from mm to inches (1 inch = 25.4 mm)
width_in_inches <- 100 / 25.4
height_in_inches <- 100 / 25.4

# Set up the PDF device with specified width and height in inches
pdf(file = "test1_survivalTestRIDGEfont7.pdf", width = width_in_inches, height = height_in_inches)

# Set font family to Helvetica and font size to 7
par(family = "Helvetica", cex = 0.7)  # Adjusting cex to scale font size appropriately

# Plot the survival curve
plot(fit, 
     lwd = 2,  # Line width
     col = c("red", "blue"),  # Colors for the groups
     xlab = "Time (year)",  # X-axis label
     ylab = "Survival rate",  # Y-axis label
     main = paste("Survival curve (p=", pValue, ")", sep = ""),  # Main title with dynamic p-value
     mark.time = TRUE)  # Show marks at censored times

# Add legend to the plot
legend("bottomright", 
       c("high risk", "low risk"),
       lwd = 2,  # Line width
       col = c("red", "blue"))  # Colors corresponding to the plot lines

# Close the PDF device
dev.off()

# Convert dimensions from mm to inches (1 inch = 25.4 mm)
width_in_inches <- 95 / 25.4
height_in_inches <- 60 / 25.4

# Set up the PDF device with specified width and height in inches
pdf(file = "newsize_test1_survivalTestRIDGE_newyrange.pdf", width = width_in_inches, height = height_in_inches)

# Set font family to Helvetica and font size to 7
par(family = "Helvetica", cex = 0.7)  # Adjusting cex to scale font size appropriately

# Plot the survival curve
plot(fit, 
     lwd = 2,  # Line width
     col = c("red", "blue"),  # Colors for the groups
     xlab = "Time (year)",  # X-axis label
     ylab = "Survival rate",  # Y-axis label
     main = paste("Survival curve (p=", pValue, ")", sep = ""),  # Main title with dynamic p-value
     mark.time = TRUE,
     ylim = c(0.4, 1))  # Show marks at censored times

# Add legend to the plot
legend("bottomright", 
       c("high risk", "low risk"),
       lwd = 2,  # Line width
       col = c("red", "blue"))  # Colors corresponding to the plot lines

# Close the PDF device
dev.off()



rt=read.table("77riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
fix(rt)

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="test1ROC.pdf", width=5, height=5)

# First plot call without the title parameter
plot(ROC_rt, time=1, col='green', lwd=2, title=FALSE)

# Adding a title using the title function
title("MFP score in test set 1")

# The rest of the plotting code remains the same
plot(ROC_rt, time=2, col='blue', add=TRUE, title=FALSE, lwd=2)
plot(ROC_rt, time=3, col='red', add=TRUE, title=FALSE, lwd=2)
legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]),"*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]),"*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]),"*")
       ),
       col=c("green", "blue", 'red'), lwd=2, bty = 'n')

dev.off()



#####new size below

rt=read.table("77riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
fix(rt)

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
# Set the width and height in millimeters
pdf(file = "newsize_fPMtest1ROC.pdf", width = 70/25.4, height = 70/25.4)

# Set font to Helvetica with size 7
par(family = "Helvetica", cex = 0.7)

# First plot call without the title parameter
plot(ROC_rt, time = 1, col = 'green', lwd = 2, title = FALSE)

# Adding a title using the title function
title("MFP score in test set 1")

# The rest of the plotting code remains the same
plot(ROC_rt, time = 2, col = 'blue', add = TRUE, title = FALSE, lwd = 2)
plot(ROC_rt, time = 3, col = 'red', add = TRUE, title = FALSE, lwd = 2)
legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]), "*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]), "*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]), "*")
       ),
       col = c("green", "blue", 'red'), lwd = 2, bty = 'n')

# Close the PDF device
dev.off()


library(survival)
# Assuming you have a dataframe 'data' with 'futime', 'fustat', and 'riskScore'
library(Hmisc)
test=rt
1-rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))


# Print the uno C-index in test set
library(survAUC)
lpnew <- test$riskScore
Surv.rsp <- survival::Surv(train$futime, train$fustat)
Surv.rsp.new <- survival::Surv(test$futime, test$fustat)             
Cstat <- UnoC(Surv.rsp, Surv.rsp.new, lpnew)
Cstat

##################below is the result of fPM score in training set.
########################################################
rt=read.table("riskTrainRIDGE.txt",header=T,sep="\t",row.names = 1)
fix(rt)
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    #²é¿´ÎåÄêÉú´æÂÊ
pdf(file="training_survivalTestRIDGE.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("bottomright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()


width_in_inches <- 95 / 25.4
height_in_inches <- 60 / 25.4
# Set up the PDF device with specified width and height in inches
pdf(file = "training_survivalTestRIDGEfont7_newyrange.pdf", width = width_in_inches, height = height_in_inches)

# Set font family to Helvetica and font size to 7
par(family = "Helvetica", cex = 0.7)  # Adjusting cex to scale font size appropriately

# Plot the survival curve
plot(fit, 
     lwd = 2,  # Line width
     col = c("red", "blue"),  # Colors for the groups
     xlab = "Time (year)",  # X-axis label
     ylab = "Survival rate",  # Y-axis label
     main = paste("Survival curve (p=", pValue, ")", sep = ""),  # Main title with dynamic p-value
     mark.time = TRUE,
     ylim = c(0.4, 1))  # Show marks at censored times

# Add legend to the plot
legend("bottomright", 
       c("high risk", "low risk"),
       lwd = 2,  # Line width
       col = c("red", "blue"))  # Colors corresponding to the plot lines

# Close the PDF device
dev.off()



rt=read.table("riskTrainRIDGE.txt",header=T,sep="\t",row.names = 1)
fix(rt)

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="trainingROC.pdf", width=5, height=5)

# First plot call without the title parameter
plot(ROC_rt, time=1, col='green', lwd=2, title=FALSE)

# Adding a title using the title function
title("MFP score in training set")

# The rest of the plotting code remains the same
plot(ROC_rt, time=2, col='blue', add=TRUE, title=FALSE, lwd=2)
plot(ROC_rt, time=3, col='red', add=TRUE, title=FALSE, lwd=2)
legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]),"*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]),"*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]),"*")
       ),
       col=c("green", "blue", 'red'), lwd=2, bty = 'n')

dev.off()

library(survival)
# Assuming you have a dataframe 'data' with 'futime', 'fustat', and 'riskScore'
library(Hmisc)
test=rt
1-rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))


# Print the uno C-index in test set
library(survAUC)
lpnew <- test$riskScore
Surv.rsp <- survival::Surv(train$futime, train$fustat)
Surv.rsp.new <- survival::Surv(test$futime, test$fustat)             
Cstat <- UnoC(Surv.rsp, Surv.rsp.new, lpnew)
Cstat



##################below is the result of COX model.
########################################################
load("allmodels.RData")
rt2=read.table("20patient_riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
ncol(rt2)
fix(rt2)
rt2=rt2[,-c(9,10)]

rt1=read.table("77riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
ncol(rt1)
rt1=rt1[,-c(9,10)]
trainingdata=read.table("riskTrainRIDGE.txt",header=T,sep="\t",row.names = 1)
trainingdata=trainingdata[,-c(9,10)]
riskScore=predict(multiCox,type="risk",newdata=train) 
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
riskScoreTest1=predict(multiCox,type="risk",newdata=rt1)      
riskTest1=as.vector(ifelse(riskScoreTest1>medianTrainRisk,"high","low"))
testRiskOut1=cbind(id=rownames(cbind(rt1,riskScore=riskScoreTest1,risk=riskTest1)),cbind(rt1,riskScore=riskScoreTest1,risk=riskTest1))
write.table(testRiskOut1,file="riskTest1COX.txt",sep="\t",quote=F,row.names=F)
rt=read.table("riskTest1COX.txt",header=T,sep="\t",row.names = 1)
#rt=testRiskOut1[,-1]
fix(rt)
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    #²é¿´ÎåÄêÉú´æÂÊ
pdf(file="test1_survivalTestCOX.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Test set 1 (p=", pValue ,")",sep=""),
     mark.time=T)
legend("bottomright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

# Set the width and height in millimeters
pdf(file = "newsize_test1_survivalTestCOX2.pdf", width = 105/25.4, height = 105/25.4)

# Set font to Helvetica with size 7
par(family = "Helvetica", cex = 0.7)

# Plot survival curve
plot(fit, 
     lwd = 2,
     col = c("red", "blue"),
     xlab = "Time (year)",
     ylab = "Survival rate",
     main = paste("Test set 1 (p=", pValue ,")", sep = ""),
     mark.time = TRUE)

# Add legend
legend("bottomright", 
       c("high risk", "low risk"),
       lwd = 2,
       col = c("red", "blue"))

# Close the PDF device
dev.off()






rt=testRiskOut1[,-1]

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="test1ROCCOX.pdf", width=5, height=5)

# First plot call without the title parameter
plot(ROC_rt, time=1, col='green', lwd=2, title=FALSE)

# Adding a title using the title function
title("Test set 1")

# The rest of the plotting code remains the same
plot(ROC_rt, time=2, col='blue', add=TRUE, title=FALSE, lwd=2)
plot(ROC_rt, time=3, col='red', add=TRUE, title=FALSE, lwd=2)
legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]),"*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]),"*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]),"*")
       ),
       col=c("green", "blue", 'red'), lwd=2, bty = 'n')

dev.off()




####new size ROC figures of test set, of COX model 
rt=read.table("riskTest1COX.txt",header=T,sep="\t",row.names = 1)
fix(rt)
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
# Set the width and height in millimeters
pdf(file = "newsize_test1ROCCOX.pdf", width = 70/25.4, height = 70/25.4)

# Set font to Helvetica with size 7
par(family = "Helvetica", cex = 0.7)

# First plot call without the title parameter
plot(ROC_rt, time = 1, col = 'green', lwd = 2, title = FALSE)

# Adding a title using the title function
title("Test set 1")

# The rest of the plotting code remains the same
plot(ROC_rt, time = 2, col = 'blue', add = TRUE, title = FALSE, lwd = 2)
plot(ROC_rt, time = 3, col = 'red', add = TRUE, title = FALSE, lwd = 2)

legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]), "*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]), "*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]), "*")
       ),
       col = c("green", "blue", 'red'), lwd = 2, bty = 'n')

# Close the PDF device
dev.off()





library(survAUC)
library(survival)
# Assuming have a dataframe 'data' with 'futime', 'fustat', and 'riskScore'
library(Hmisc)

test=testRiskOut1[,-1]
1-rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))
# Print the uno C-index in test set
lpnew <- test$riskScore
Surv.rsp <- survival::Surv(train$futime, train$fustat)
Surv.rsp.new <- survival::Surv(test$futime, test$fustat)             
Cstat <- UnoC(Surv.rsp, Surv.rsp.new, lpnew)
Cstat





##################below is the result of RSF model
########################################################
load("allmodels.RData")
rt2=read.table("20patient_riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
ncol(rt2)
rt2=rt2[,-c(9,10)]

rt1=read.table("77riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
ncol(rt1)
rt1=rt1[,-c(9,10)]

trainingdata=read.table("riskTrainRIDGE.txt",header=T,sep="\t",row.names = 1)
trainingdata=trainingdata[,-c(9,10)]
riskScore=predict(object = rsfmodel,newdata=train)$predicted
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
riskScoreTest1=predict(object = rsfmodel,newdata=rt1)$predicted       #????train?õ?ģ??Ԥ??test??Ʒ????
riskTest1=as.vector(ifelse(riskScoreTest1>medianTrainRisk,"high","low"))
testRiskOut1=cbind(id=rownames(cbind(rt1,riskScore=riskScoreTest1,risk=riskTest1)),cbind(rt1,riskScore=riskScoreTest1,risk=riskTest1))
write.table(testRiskOut1,file="riskTest1RSF.txt",sep="\t",quote=F,row.names=F)
rt=testRiskOut1[,-1]
rt=read.table("riskTest1RSF.txt",header=T,sep="\t",row.names = 1)
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)   
pdf(file="test1_survivalTestRSF.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Test set 1 (p=", pValue ,")",sep=""),
     mark.time=T)
legend("bottomright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

####new size 
# Set the width and height in millimeters
pdf(file = "newsize_test1_survivalTestRSF.pdf", width = 105/25.4, height = 105/25.4)

# Set font to Helvetica with size 7
par(family = "Helvetica", cex = 0.7)

# Plot survival curve
plot(fit, 
     lwd = 2,
     col = c("red", "blue"),
     xlab = "Time (year)",
     ylab = "Survival rate",
     main = paste("Test set 1 (p=", pValue ,")", sep = ""),
     mark.time = TRUE)

# Add legend
legend("bottomright", 
       c("high risk", "low risk"),
       lwd = 2,
       col = c("red", "blue"))

# Close the PDF device
dev.off()




rt=testRiskOut1[,-1]

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="test1ROCRSF.pdf", width=5, height=5)

# First plot call without the title parameter
plot(ROC_rt, time=1, col='green', lwd=2, title=FALSE)

# Adding a title using the title function
title("Test set 1")

# The rest of the plotting code remains the same
plot(ROC_rt, time=2, col='blue', add=TRUE, title=FALSE, lwd=2)
plot(ROC_rt, time=3, col='red', add=TRUE, title=FALSE, lwd=2)
legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]),"*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]),"*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]),"*")
       ),
       col=c("green", "blue", 'red'), lwd=2, bty = 'n')

dev.off()



####new size ROC figures of test set, of RSF model 
rt=read.table("riskTest1RSF.txt",header=T,sep="\t",row.names = 1)
fix(rt)
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)



ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)

# Set the width and height in millimeters
pdf(file = "newsize_test1ROCRSF.pdf", width = 70/25.4, height = 70/25.4)

# Set font to Helvetica with size 7
par(family = "Helvetica", cex = 0.7)

# First plot call without the title parameter
plot(ROC_rt, time = 1, col = 'green', lwd = 2, title = FALSE)

# Adding a title using the title function
title("Test set 1")

# The rest of the plotting code remains the same
plot(ROC_rt, time = 2, col = 'blue', add = TRUE, title = FALSE, lwd = 2)
plot(ROC_rt, time = 3, col = 'red', add = TRUE, title = FALSE, lwd = 2)

legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]), "*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]), "*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]), "*")
       ),
       col = c("green", "blue", 'red'), lwd = 2, bty = 'n')

# Close the PDF device
dev.off()









library(survAUC)
library(survival)
# Assuming you have a dataframe 'data' with 'futime', 'fustat', and 'riskScore'
library(Hmisc)


test=testRiskOut1[,-1]
1-rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))
# Print the uno C-index in test set
lpnew <- test$riskScore
Surv.rsp <- survival::Surv(train$futime, train$fustat)
Surv.rsp.new <- survival::Surv(test$futime, test$fustat)             
Cstat <- UnoC(Surv.rsp, Surv.rsp.new, lpnew)
Cstat




##################below is the result of GBM model.
########################################################
library(gbm)
load("allmodels.RData")
rt2=read.table("20patient_riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
ncol(rt2)
rt2=rt2[,-c(9,10)]

rt1=read.table("77riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
ncol(rt1)
rt1=rt1[,-c(9,10)]

trainingdata=read.table("riskTrainRIDGE.txt",header=T,sep="\t",row.names = 1)
trainingdata=trainingdata[,-c(9,10)]
riskScore=predict(gbmmodel,train)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
riskScoreTest1=predict(gbmmodel,rt1)       #????train?õ?ģ??Ԥ??test??Ʒ????
riskTest1=as.vector(ifelse(riskScoreTest1>medianTrainRisk,"high","low"))
testRiskOut1=cbind(id=rownames(cbind(rt1,riskScore=riskScoreTest1,risk=riskTest1)),cbind(rt1,riskScore=riskScoreTest1,risk=riskTest1))
write.table(testRiskOut1,file="riskTest1GBM.txt",sep="\t",quote=F,row.names=F)
rt=testRiskOut1[,-1]
rt=read.table("riskTest1GBM.txt",header=T,sep="\t",row.names = 1)
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    #²é¿´ÎåÄêÉú´æÂÊ
pdf(file="test1_survivalTestGBM.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Test set 1 (p=", pValue ,")",sep=""),
     mark.time=T)
legend("bottomright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()
###new size############################################
# Set the width and height in millimeters
pdf(file = "newsize_test1_survivalTestGBM.pdf", width = 105/25.4, height = 105/25.4)

# Set font to Helvetica with size 7
par(family = "Helvetica", cex = 0.7)

# Plot survival curve
plot(fit, 
     lwd = 2,
     col = c("red", "blue"),
     xlab = "Time (year)",
     ylab = "Survival rate",
     main = paste("Test set 1 (p=", pValue ,")", sep = ""),
     mark.time = TRUE)

# Add legend
legend("bottomright", 
       c("high risk", "low risk"),
       lwd = 2,
       col = c("red", "blue"))

# Close the PDF device
dev.off()

rt=testRiskOut1[,-1]

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="test1ROCGBM.pdf", width=5, height=5)

# First plot call without the title parameter
plot(ROC_rt, time=1, col='green', lwd=2, title=FALSE)

# Adding a title using the title function
title("Test set 1")

# The rest of the plotting code remains the same
plot(ROC_rt, time=2, col='blue', add=TRUE, title=FALSE, lwd=2)
plot(ROC_rt, time=3, col='red', add=TRUE, title=FALSE, lwd=2)
legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]),"*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]),"*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]),"*")
       ),
       col=c("green", "blue", 'red'), lwd=2, bty = 'n')

dev.off()


####new size ROC figures of test set, of gbm model 
rt=read.table("riskTest1GBM.txt",header=T,sep="\t",row.names = 1)
fix(rt)
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)

# Set the width and height in millimeters
pdf(file = "newsize_test1ROCGBM.pdf", width = 70/25.4, height = 70/25.4)

# Set font to Helvetica with size 7
par(family = "Helvetica", cex = 0.7)

# First plot call without the title parameter
plot(ROC_rt, time = 1, col = 'green', lwd = 2, title = FALSE)

# Adding a title using the title function
title("Test set 1")

# The rest of the plotting code remains the same
plot(ROC_rt, time = 2, col = 'blue', add = TRUE, title = FALSE, lwd = 2)
plot(ROC_rt, time = 3, col = 'red', add = TRUE, title = FALSE, lwd = 2)

legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]), "*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]), "*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]), "*")
       ),
       col = c("green", "blue", 'red'), lwd = 2, bty = 'n')

# Close the PDF device
dev.off()







library(survAUC)
library(survival)
# Assuming you have a dataframe 'data' with 'futime', 'fustat', and 'riskScore'
library(Hmisc)

test=testRiskOut1[,-1]
1-rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))
# Print the uno C-index in test set
lpnew <- test$riskScore
Surv.rsp <- survival::Surv(train$futime, train$fustat)
Surv.rsp.new <- survival::Surv(test$futime, test$fustat)             
Cstat <- UnoC(Surv.rsp, Surv.rsp.new, lpnew)
Cstat





##################below is the result of SVM model
########################################################
library(survivalsvm)
load("allmodels.RData")
rt2=read.table("20patient_riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
ncol(rt2)
rt2=rt2[,-c(9,10)]

rt1=read.table("77riskTestRIDGE.txt",header=T,sep="\t",row.names = 1)
ncol(rt1)
rt1=rt1[,-c(9,10)]

trainingdata=read.table("riskTrainRIDGE.txt",header=T,sep="\t",row.names = 1)
trainingdata=trainingdata[,-c(9,10)]
riskScore=predict(svmmodel,train)$predicted[1,]
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
riskScoreTest1=predict(svmmodel,rt1)$predicted[1,]      
riskTest1=as.vector(ifelse(riskScoreTest1>medianTrainRisk,"high","low"))
testRiskOut1=cbind(id=rownames(cbind(rt1,riskScore=riskScoreTest1,risk=riskTest1)),cbind(rt1,riskScore=riskScoreTest1,risk=riskTest1))
write.table(testRiskOut1,file="riskTest1SVM.txt",sep="\t",quote=F,row.names=F)
rt=testRiskOut1[,-1]
rt=read.table("riskTest1SVM.txt",header=T,sep="\t",row.names = 1)
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = F)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    #²é¿´ÎåÄêÉú´æÂÊ
pdf(file="test1_survivalTestSVM.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Test set 1 (p=", pValue ,")",sep=""),
     mark.time=T)
legend("bottomright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()
####newsize##############################################
# Set the width and height in millimeters
pdf(file = "newsize_test1_survivalTestSVM2.pdf", width = 105/25.4, height = 105/25.4)

# Set font to Helvetica with size 7
par(family = "Helvetica", cex = 0.7)

# Plot survival curve
plot(fit, 
     lwd = 2,
     col = c("red", "blue"),
     xlab = "Time (year)",
     ylab = "Survival rate",
     main = paste("Test set 1 (p=", pValue ,")", sep = ""),
     mark.time = TRUE)

# Add legend
legend("bottomright", 
       c("high risk", "low risk"),
       lwd = 2,
       col = c("red", "blue"))

# Close the PDF device
dev.off()

rt=testRiskOut1[,-1]
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="test1ROCSVM.pdf", width=5, height=5)
# First plot call without the title parameter
plot(ROC_rt, time=1, col='green', lwd=2, title=FALSE)
# Adding a title using the title function
title("Test set 1")
# The rest of the plotting code remains the same
plot(ROC_rt, time=2, col='blue', add=TRUE, title=FALSE, lwd=2)
plot(ROC_rt, time=3, col='red', add=TRUE, title=FALSE, lwd=2)
legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]),"*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]),"*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]),"*")
       ),
       col=c("green", "blue", 'red'), lwd=2, bty = 'n')

dev.off()


####new size ROC figures of test set, of SVM model 
rt=read.table("riskTest1SVM.txt",header=T,sep="\t",row.names = 1)
fix(rt)
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
# Set the width and height in millimeters
pdf(file = "newsize_test1ROCSVM.pdf", width = 70/25.4, height = 70/25.4)

# Set font to Helvetica with size 7
par(family = "Helvetica", cex = 0.7)

# First plot call without the title parameter
plot(ROC_rt, time = 1, col = 'green', lwd = 2, title = FALSE)

# Adding a title using the title function
title("Test set 1")

# The rest of the plotting code remains the same
plot(ROC_rt, time = 2, col = 'blue', add = TRUE, title = FALSE, lwd = 2)
plot(ROC_rt, time = 3, col = 'red', add = TRUE, title = FALSE, lwd = 2)

legend('bottomright',
       c(paste0('AUC at 1st year: ', sprintf("%.03f", ROC_rt$AUC[1]), "*"),
         paste0('AUC at 2nd year: ', sprintf("%.03f", ROC_rt$AUC[2]), "*"),
         paste0('AUC at 3rd year: ', sprintf("%.03f", ROC_rt$AUC[3]), "*")
       ),
       col = c("green", "blue", 'red'), lwd = 2, bty = 'n')

# Close the PDF device
dev.off()


library(survAUC)
library(survival)
# Assuming you have a dataframe 'data' with 'futime', 'fustat', and 'riskScore'
library(Hmisc)


test=testRiskOut1[,-1]
1-rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))
rcorr.cens(test$riskScore,Surv(test$futime,test$fustat))
# Print the uno C-index in test set
lpnew <- test$riskScore
Surv.rsp <- survival::Surv(train$futime, train$fustat)
Surv.rsp.new <- survival::Surv(test$futime, test$fustat)             
Cstat <- UnoC(Surv.rsp, Surv.rsp.new, lpnew)
1-Cstat

