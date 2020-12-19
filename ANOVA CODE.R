
#######################################################################
########################### packages needed ###########################
#######################################################################
install.packages("PerformanceAnalytics")
install.packages("stargazer")
install.packages("lmtest")
install.packages("car")
install.packages("nortest")
install.packages("MASS")
install.packages("Rmisc")
install.packages("texreg")
library(texreg)
library(MASS)
library(car)
library(dae)
library(lmerTest)
library(stargazer)
library(xtable)
library(car)
library(nortest)
library(Rmisc)
library("PerformanceAnalytics")

#######################################################################
########################### data input ################################
#######################################################################
maindata<-read.csv(file.choose(),header=T)
head(maindata)
nrow(maindata)
attach(maindata)

adhesive<-factor(adhesive)
primer<-factor(primer)
thickness<-factor(thickness)
rm(Primer)
rm(Thickness)
rm(Adhesive)

######################################################################################
########################### preliminary data analysis ################################
######################################################################################
tapply(peelstrength,adhesive,mean)
tapply(peelstrength,primer,mean)
tapply(peelstrength,thickness,mean)
tapply(peelstrength,adhesive,var)
tapply(peelstrength,primer,var)
tapply(peelstrength,thickness,var)
tapply(peelstrength,adhesive,length)
tapply(peelstrength,primer,length)
tapply(peelstrength,thickness,length)
summarySE(maindata,measurevar="peelstrength",groupvars=c("adhesive","primer","thickness"))
mean(peelstrength)
sum(peelstrength)
mytable<-xtabs(peelstrength~adhesive+primer+thickness)

# Plots

boxplot(peelstrength~adhesive,ylab="Peel Strength",xlab="Adhesive")
boxplot(peelstrength~primer,ylab="Peel Strength",xlab="Primer")
boxplot(peelstrength~thickness,ylab="Peel Strength",xlab="Thickness")
boxplot(peelstrength~adhesive+primer+thickness,ylab="Peel Strength",xlab="Trt",data=maindata)
abline(52.77,0,lty=2)

# Interaction plots
# 2-way plots

par(mfrow = c(2,2))
interaction.plot(x.factor = adhesive, trace.factor = primer, response = peelstrength,
                 fun = mean, legend = F,  type = "b", lwd = 3,axes=F,
                 trace.label="Primer", xlab = "Adhesive 1 2 3 4", ylab = 
                         "Peel Strength", main = "Two-Way Plot (Adhesive vs Primer)",las=1,
                 col=1:2, pch=1:2 )
axis(1, at=1:2:3:4)
axis(2, at=seq(5,30,5))
text(1.5, 22, "cold room")
text(1.5, 9, "room temp")

interaction.plot(x.factor = potato, trace.factor = temp, response = leak,
                 fun = mean, legend = F,  type = "b", lwd = 3,axes=F,
                 trace.label="temp", xlab = "Species 1 vs 2", ylab = 
                         "Damaged score for ion leakage", main = "Two-Way Plot (potato vs cold temp)",las=1,
                 col=1:2, pch=1:2 )
axis(1, at=1:2)
axis(2, at=seq(5,30,5))
text(1.5, 22, "-4 degree c")
text(1.5, 6, "-8 degree c")

interaction.plot(x.factor = temp, trace.factor = regime, response = leak,
                 fun = mean, legend = F,  type = "b", lwd = 3,axes=F,
                 trace.label="regime", xlab = "-4 deg C vs -8 deg C", ylab = 
                         "Damaged score for ion leakage", main = "Two-Way Plot (cold temp vs regime)",las=1,
                 col=1:2, pch=1:2 )
axis(1,at=1:2)
axis(2, at=seq(5,30,5))
text(1.5, 19, "cold room")
text(1.5, 10, "room temp")

# split the data for 3-way plot
maindata1 <- subset(maindata, maindata$thickness == 1)
maindata2 <- subset(maindata, maindata$thickness == 2)
nrow(maindata1)
par(mfrow = c(1,1))

interaction.plot(x.factor = maindata1$adhesive, trace.factor = maindata1$primer, response = maindata1$peelstrength,
                 fun = mean, legend = F,  type = "b", lwd = 3,axes=F,
                 trace.label="Primer", xlab = "Adhesive", ylab = 
                         "Peel Strength",las=1,
                 col=1:2, pch=1:2 )
title (main = "Thickness A (Adhesive vs Primer)")
axis(1, at=1:4)
axis(2, at=10:80)
text(1.5, 55, "With Primer")
text(1.5, 47, "Without Primer")
interaction.plot(x.factor = maindata2$adhesive, trace.factor = maindata2$primer, response = maindata2$peelstrength,
                 fun = mean, legend = F,  type = "b", lwd = 3,axes=F,
                 trace.label="Primer", xlab = "Adhesive", ylab = 
                         "Peel Strength",las=1,
                 col=1:2, pch=1:2 )
title (main = "Thickness B (Adhesive vs Primer)")
axis(1, at=1:4)
axis(2, at=10:80)
text(2, 75, "With Primer")
text(1.5, 47, "Without Primer")

######################################################################################
########################### ANOVA model & diagnostics ################################
######################################################################################
# Full model
fit1<-lm(peelstrength~adhesive*primer*thickness,contrasts=c(adhesive=contr.sum,
                                             primer=contr.sum,thickness=contr.sum))
Anova(fit1,type=3) ##Type III SS
anova(fit1) ##Type I SS
summary(fit1)
stargazer(Anova(fit1),summary=F)
stargazer(anova(fit1),summary=F)

# Residual analysis

par(mfrow = c(2,2))
fits<-fitted.values(fit1)
resid<-residuals(fit1)
plot(fits,resid,las=1,main="Residuals vs Fits", xlab="Fitted values"
     , ylab="Residuals")
plot(fit1)
summary(fit1)

# input data again
trtdata<-read.csv(file.choose(),header=T)
trtdata$Trt
table(trtdata$Trt)
plot.default(trtdata$Trt,resid,las=1,main="Residuals vs Treatment Levels"
             ,ylab="Residuals",xlab="Treatments",axes=T,ann=par("ann"))
axis(1, at=1:16)
par(mfrow = c(1,1))
leveneTest(peelstrength, interaction(adhesive,primer,thickness))
leveneTest(peelstrength~factor(trtdata$Trt),data=maindata)
stargazer(leveneTest(peelstrength, interaction(adhesive,primer,thickness)),summary=F)

# Normality test

z<-qqnorm(resid)
cor(z$x,z$y) ###==0.99463
###both methods gave the same correlation
stargazer(shapiro.test(resid),type="text")
texreg(shapiro.test(resid))
f<-shapiro.test(resid)
stargazer(f)

n = 96
mse=12.3
ExpVals = sapply(1:n, function(k) 12.3 * qnorm((k-.375)/(n+.25)))
ExpVals
cor(ExpVals,sort(resid))
##0.9945999 Its almost the same as previous, so this method works too!!
hist(resid)

# Bonferroni OT

boxplot(rstudent(fit1),main="Rstudent box plot",ylab="Rstudent",las=1)
youtliers<-which(abs(rstudent(fit1))>qt(1-(0.05/(2*96)),96-16-1))
qt(1-(0.05/(2*75)),96-16-1)
rstudent(fit1)
outlierTest(fit1)


chart.Correlation(maindata, histogram=TRUE, pch=19)
vif(fit1)
stargazer(vif(fit1))

# Transformation

summary(peelstrength)
boxcox(peelstrength~factor(adhesive)*factor(primer)*factor(thickness),
       lambda=seq(-2,2,length = 13), plotit=T)
gmean<-exp(mean(log(peelstrength)))
sse<-c()
lambda<-c()
i<-1
for(lam in seq(-2,2,0.1)){
        if(lam !=0){
                tpeelstrength<-(peelstrength^lam-1)/(lam*gmean^(lam-1))
        }else{
                tpeelstrength<-log(peelstrength)*gmean
        }
        test<-anova(aov(tpeelstrength~factor(adhesive)*factor(primer)*factor(thickness)))
        sse[i]<-test['Residuals','Sum Sq']
        lambda[i]<-lam
        i<-i+1
}
cbind(lambda,sse)


# Trying log
logpeel<-log(peelstrength)
logfit<-fit2<-lm(logpeel~adhesive*primer*thickness,contrasts=c(adhesive=contr.sum,
                primer=contr.sum,thickness=contr.sum))
plot(logfit)
leveneTest(logpeel, interaction(adhesive,primer,thickness))
summary(logfit)
Anova(logfit,type=3) 

transpeel<-(peelstrength)^0.5

# transformed Fit
# Full model transformed
fit2<-lm(transpeel~adhesive*primer*thickness,contrasts=c(adhesive=contr.sum,
                                                  primer=contr.sum,thickness=contr.sum))
Anova(fit2,type=3) ##Type III SS
anova(fit2) ##Type I SS
summary(fit2)

# Residual analysis

par(mfrow = c(2,2))
transfits<-fitted.values(fit2)
transresid<-residuals(fit2)
plot(transfits,transresid,las=1,main="Residuals vs Fits (Transformed)", xlab="Fitted values"
     , ylab="Residuals")
plot(fit2)
plot.default(trtdata$Trt,transresid,las=1,main="Residuals vs Treatment Levels"
             ,ylab="Residuals",xlab="Treatments",axes=T,ann=par("ann"))
axis(1, at=1:16)
par(mfrow = c(1,1))

leveneTest(transpeel, interaction(adhesive,primer,thickness))
leveneTest(transpeel~factor(trtdata$Trt),data=maindata)
attach(transpeel)
transpeel<-data.frame(transpeel)
stargazer(leveneTest(transpeel, interaction(adhesive,primer,thickness), summary=F))
stargazer(leveneTest(peelstrength, interaction(adhesive,primer,thickness)),summary=F)
stargazer(leveneTest(fit2),summary=F)

data.class(peelstrength)
data.class(transpeel)
shapiro.test(transresid)

hist(resid)

# Bonferroni OT

boxplot(rstudent(fit2),main="Rstudent box plot",ylab="Rstudent",las=1)
youtliers<-which(abs(rstudent(fit2))>qt(1-(0.05/(2*96)),96-16-1))
qt(1-(0.05/(2*75)),96-16-1)
rstudent(fit1)
outlierTest(fit2)

bdata<-cbind(transpeel,adhesive,primer,thickness)
chart.Correlation(bdata, histogram=TRUE, pch=19)
vif(fit2)
stargazer(vif(fit2))

######################################################################################
########################### adehsive effect analysis #################################
######################################################################################

# Tukey value
T<-qtukey(0.95,16,80)/sqrt(2)

#Bonferroni contrasts
g=3
B<-qt(1-(0.05/(2*3)),96-16)
mse<-12.3
mu4..<-64.333
mu1..<-63.958
mu2..<-52.041
mu3..<-30.758

n4..<-24
n1..<-24
n2..<-24
n3..<-24

d1hat<-mu4..-mu1..
d2hat<-mu4..-mu2..
d3hat<-mu4..-mu3..
se_d1hat<-sqrt(mse*((1/24)+(1/24)))
se_d2hat<-sqrt(mse*((1/24)+(1/24)))
se_d3hat<-sqrt(mse*((1/24)+(1/24)))
ul<-d1hat+B*se_d1hat
ll<-d1hat-B*se_d1hat

ul<-d2hat+B*se_d2hat
ll<-d2hat-B*se_d2hat

ul<-d3hat+B*se_d3hat
ll<-d3hat-B*se_d3hat

L1hat<-((mu101+mu102+mu111+mu112)/4)-((mu201+mu202+mu211+mu212)/4)
L2hat<-((mu101+mu102+mu201+mu202)/4)-((mu111+mu112+mu211+mu212)/4)
se_Lhat<-sqrt(mse*(1/16)*(1/5+1/5+1/12+1/13+1/13+1/13+1/7+1/7))
ul<-L1hat+B*se_Lhat
ll<-L1hat-B*se_Lhat

ul<-L2hat+B*se_Lhat
ll<-L2hat-B*se_Lhat

######################################################################################
########################### treatment effect analysis ################################
######################################################################################
fit_aov<-aov(peelstrength~adhesive*primer*thickness,contrasts=c(adhesive=contr.sum,
      primer=contr.sum,thickness=contr.sum)) ##same as the (lm) model 

Anova(fit_aov,type=3) ##same as the (lm) modelfit2
par(mfrow = c(1,1))
fit_tukey<-TukeyHSD(fit_aov,conf.level=0.95)
plot(fit_tukey,sub="Tukey Honest Significant Differences",las=1)
summary(fit_aov)


# Line plot
summarySE(maindata)
tapply(trtdata$peelstrength,trtdata$Trt,mean)
meansbb<-tapply(trtdata$peelstrength,trtdata$Trt,mean)
sort(meansbb)


