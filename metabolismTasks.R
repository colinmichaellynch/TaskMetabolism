#clear workspace 
rm(list = ls())

#load required packages 
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(dunn.test)){install.packages("dunn.test")}
if(!require(dplyr)){install.packages("dplyr")}
if(!require(lme4)){install.packages("lme4")}
if(!require(MASS)){install.packages("MASS")}

setwd("~/Metabolic Scaling/Metabolism of tasks")

#analyze CO2 traces of tasks 
data = read.csv("CarbonDioxideofTasks.csv")

#count how often each task appears in the dataset 
table(data$Task.Manuel)

#here we calculate metabolic rate given the corrected CO2 value (mean CO2 - baseline mean CO2)
RQ = .86 
airFlow = 72 
Oxy = 20.1
data$MR = (((data$Corrected.CO2*airFlow)/RQ)*Oxy)/(60)

#Use mixed effects model to evaluate whether tasks have different metabolic rates. Task is the fixed effect, and colony ID is the random effect 
mdl1 = lmer(MR~Task.Manuel+(1|Colony.ID), data=data)
anova(mdl1) #check if task effect is significant 
plot(mdl1) #check for model adequacy. Residuals show data does not have constant variance 

#here we attempt to fix the data by performing a boxcox transformation 
boxCox = boxcox(lm(data$MR ~ 1))
lambda = boxCox$lambda[which.max(boxCox$objective)]
data$mrTransformed = (data$MR ^ lambda - 1) / lambda

mdl2 = lmer(mrTransformed~Task.Manuel+(1|Colony.ID), data=data)
plot(mdl2) #Data still does not have constant variance 

medians = c() #here we get the median value of each task, to be used later for a levene's test 
for(i in 1:nrow(data)){
  dataTemp = subset(data, Task.Manuel == data$Task.Manuel[i])
  medians[i] = median(dataTemp$MR)
}
data$medians = medians
data$Res = abs(data$MR-data$medians) #residual is the difference between the expected value of a task (median) and the actual value 

# Levene's test based on residuals from medians 
leveneTest = aov(Res~Task.Manuel,data)
summary(leveneTest)
TukeyHSD(leveneTest) #look at all pairwise comparisons 

kruskal.test(MR~Task.Manuel, data=data) #nonparametric test for task metabolsm 
dunn.test(data$MR, data$Task.Manuel) 

#visualize results (figure 2)
ggplot(data, aes(x = reorder(Task.Manuel, -MR), y = MR)) + geom_jitter(width=.2) + geom_boxplot(alpha=.5, coef = 10) + theme_classic() + theme(text = element_text(size = 14)) + xlab("Task") + ylab("Metabolic Rate (\u00b5W)") + scale_x_discrete(labels=c("B" = "Brood Care", "F" = "Foraging", "M" = "Maintenance", "R" = "Resting")) + geom_label(label="a",  x=1, y=410, label.size = NA, color = "black") + ylim(0, 410)+ geom_label(label="b",  x=2, y=185, label.size = NA, color = "black") +  geom_label(label="b",  x=3, y=185, label.size = NA, color = "black") + geom_label(label="c",  x=4, y=70, label.size = NA, color = "black") 

#get metabolic rate of resting to compare to other studies 
dataSub = subset(data, Task.Manuel == "R")
x = dataSub$MR
median(x) 

### toy model

#proportions for each task come from Holbrook, C. T., Eriksson, T. H., Overson, R. P., Gadau, J., & Fewell, J. H. (2013). Colony-size effects on task organization in the harvester ant Pogonomyrmex californicus. Insectes sociaux, 60, 191-201.

#get the median cost of each tasks
dataM = aggregate(data$MR, by = list(data$Task.Manuel), FUN = median)

#hypothetical colony sizes 
N1 = 50
N2 = 300

#get proportion of time spent on 4 tasks for both colony sizes 
pBroodCare1 = .122
pBroodCare2 = .1
mBroodCare = dataM$x[1]

pMain1 = .034 + .017
pMain2 = .032 + .036
mMain = dataM$x[3]

pFor1 = .008 + .23 + .184 
pFor2 = .007 + .231 + .232 
mFor = dataM$x[2]

#remaining time leftover is resting 
pRest1 = 1 - pBroodCare1 - pMain1 - pFor1
pRest2 = 1 - pBroodCare2 - pMain2 - pFor2
mRest = dataM$x[4]

M1 = N1*sum(pBroodCare1*mBroodCare+pMain1*mMain+pFor1*mFor+pRest1*mRest)
M2 = N2*sum(pBroodCare2*mBroodCare+pMain2*mMain+pFor2*mFor+pRest2*mRest)

s = (log(M2)-log(M1))/(log(N2)-log(N1))
s #the slope we are looking for 

### supplemental figure ###

rm(list = ls())

dataCO2 = read.csv("4PeriodCO2Trace.csv")

#figure S1
ggplot(dataCO2, aes(x = Time, y = CO2)) + geom_line() + geom_vline(xintercept = 2*60, color = "black", linetype = "dashed") + geom_vline(xintercept = 4*60, color = "black", linetype = "dashed") + geom_vline(xintercept = 6*60, color = "black", linetype = "dashed") + geom_vline(xintercept = (6+8)*60, color = "black", linetype = "dashed") + geom_label(label="Baseline",  x=40, y=2.985, label.size = NA, color = "black") + geom_label(label="Seeds",  x=181, y=2.985, label.size = NA, color = "black") + geom_label(label="Acclimation",  x=302, y=2.985, label.size = NA, color = "black") + geom_label(label="Observation Period",  x=600, y=2.985, label.size = NA, color = "black") + geom_label(label="Baseline",  x=930, y=2.985, label.size = NA, color = "black") + theme_classic() + theme(text = element_text(size = 14)) + xlab("Time (s)") + ylab(expression(CO[2]~(ppm))) 

#get the signal to noise ratio 
noise = max(dataCO2$CO2[c(1:76, 890:960)])-min(dataCO2$CO2[c(1:76, 890:960)])
signal = mean(dataCO2$CO2[360:840]) - mean(dataCO2$CO2[c(1:76, 890:960)])
signalNoiseRatio = signal/noise
signalNoiseRatio
