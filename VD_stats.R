#Problem 1

x <- c(15, 23, 12, 18, 9, 28, 11,10)
y <- c(25, 20, 35, 15, 40, 16, 10, 22, 18, 32)
hello <- c(x,y)
densityPlot(hello)
#part a
mean(x) #mean is 15.75

mean(y) #mean is 23.3

var(x) #sample variance is 46.21429

var(y) #sample variance is 92.67778

#estimator is var(x)/var(y)
var(x)/var(y) #0.4986555

#part b
#chi squared CI 95%
qf(1-.025,7,9) #4.197047
qf(1-.025,9,7) #4.823217

lcl <- var(x)/var(y) / qf(1-.025,7,9) #0.1188111
ucl <- var(x)/var(y) * qf(1-.025,9,7) #2.405124
c(lcl, ucl)
#part c
qqnorm(x)
qqline(x)
qqnorm(y)
qqline(y)
#assume a normal distribution



#problem 2
survey <- c(46.06, 64.41, 94.62, 63.41, 60.53,
            37.85, 110.40, 72.52, 60.68, 70.34,
            49.94, 37.47, 54.55, 82.33, 63.95,
            50.80, 46.34, 43.68, 73.86, 65.05,
            423.33, 478.83, 101.69, 227.06, 187.40)

densityPlot(survey) #since data is skewed use sign test
#part a
#Null hypothesis: the sample mean is 63.688
#ALt hypothesis: the sample mean is greater than 63.688

#part b
qqnorm(survey)
qqline(survey)
qqPlot(survey)
# does not follow a normal distribution
#t.test(survey, mu = 63.688,alternative = "greater") #parametric test
#wilcox.test(survey, mu = 63.688,alternative = "greater") #non parametric test
SIGN.test(survey, mu= 63.688, alternative = "greater")
conc_fried#part c
#p-value > 0.05, therefore fail to reject null hyp. 

#Problem 3
w_traf <- c(214, 159, 169, 202, 103, 119, 200, 109, 132, 142, 194, 104, 219, 119, 234)
wo_traf <- c(159, 135, 141, 101, 102, 168, 62, 167, 174, 159, 66, 118, 181, 171, 112)
traf <- data.frame(w_traf, wo_traf)
hello1 <- c(w_traf, wo_traf)
densityPlot(hello1)
#part a
##This is paired sample data
#part b
#H0: there is not a difference between pollution levels by closing streets to car traffic.
#Ha: there is a difference between pollution levels by closing streets to car traffic
qqplot(w_traf, wo_traf) #not normal 
qqPlot(traf)
qqnorm(wo_traf)
qqline(wo_traf)
qqnorm(w_traf)
qqline(w_traf)
#wilcox.test(traf$w_traf, traf$wo_traf, alternative = "two.sided", paired = TRUE)
t.test(traf$w_traf, traf$wo_traf, alternative = "two.sided", paired = TRUE)
mean(-wo_traf+  w_traf) #26.86667
sd(-wo_traf+  w_traf) #68.78213

#part c
abs(qnorm(0.1526/2)) #1.430408
1.430408 / sqrt(15) #0.3693298
#There is a small effect size from this test

#part d
lcl <- 26.86667 - qt(1-0.05/2,14)* 68.78213/ sqrt(15) #lcl
ucl <- 26.86667 + qt(1-0.05/2,14)* 68.78213/ sqrt(15) #ucl
c(lcl, ucl)
#since p-value >0.05, and CI includes 0 at 95% level, we fail to reject null hyp 
#and conclude that there is no difference between pollution levels. 


#Problem 4
conc <- data.frame("Dog" = as.factor(c(1,2,3,4,5,6,7,8,9,10)),
                   "Insoflurane" = c(0.28,0.51, 1.00, 0.39, 0.29, 0.36, 0.32, 0.69, 0.17, 0.33),
                   "Halothane" = c(0.30, 0.39, 0.63, 0.38, 0.21, 0.88, 0.39, 0.51, 0.32, 0.42),
                   "Cyclopropane" = c(1.07, 1.35, 0.69, 0.28, 1.24, 1.53, 0.49, 0.56, 1.02, 0.30))
#part a

#H0: There is no anesthetic effect on concentration
#Ha: There is anesthetic effect on concentration
qqnorm(conc$Insoflurane)
qqline(conc$Insoflurane)
qqnorm(conc$Halothane)
qqline(conc$Halothane)
qqnorm(conc$Cyclopropane)
qqline(conc$Cyclopropane)
boxplot(conc$Insoflurane, conc$Halothane, conc$Cyclopropane, names = c("A", "B", "C"))
#does not follow a normal distribution
densityPlot(conc_fried)
#cant do anova since not normal
#Friedman is the non parametric version of anova
#friedman.test(conc$Insoflurane, conc$Halothane, conc$Cyclopropane, groups = conc$Dog)
is.matrix(conc) #friedman only takes matrix
conc_fried <- as.matrix(conc[,-1]) #drop the dog factor
friedman.test(conc_fried, groups = conc$Dog)
#p-value is > 0.05, therefore there we fail to reject null hypothesis. There is 
#no anesthetic effect on concentration
kruskal.test(conc)

view(groups)
iso = c(0.28,0.51, 1.00, 0.39, 0.29, 0.36, 0.32, 0.69, 0.17, 0.33)
halo = c(0.30, 0.39, 0.63, 0.38, 0.21, 0.88, 0.39, 0.51, 0.32, 0.42)
cyclo = c(1.07, 1.35, 0.69, 0.28, 1.24, 1.53, 0.49, 0.56, 1.02, 0.30)
data1 = c(iso,halo,cyclo)
groups = factor(c(rep(1:10, each = 1),rep(1:10, each = 1), rep(1:10, each = 1) ))
#part b
anova1 <- aov(data1 ~ groups)
TukeyHSD(anova1)


#part b
#Since there is no anesthetic effect on concentration, none of the anesthetic treatments
#affects concentration differently. 


#Problem 5
dat <- c(12.0, 12.0, 15.0, 15.0, 10.0,
         15.0, 6.0, 12.0, 7.5, 7.0,
         15.0, 10.0, 15.0, 10.0, 15.0,
         25.0, 10.5, 8.0, 7.0, 12.0,
         7.0, 10.5)
#part a
library(boot)
median.fun <- function(dat, idx) median(dat[idx], na.rm = TRUE)
boot.out = boot(data = dat,statistic = median.fun,R = 1000)
boot.out
boot.ci(boot.out)
#Resample
B = 1000
M = NA
for(i in 1:B) {
  x = sample(dat,length(dat),replace=T)
  M[i] = median(x)
  }
median(M) #the median is 11.25
sd(M) #1.131318

#part b
#Confidence interval from Resamples
boot.ci(boot.out, type = "norm") #95%   ( 9.06, 13.43 ) 

#part c
hist(boot.out$t)
plot(density(boot.out$t))
#The distribution essentially only takes on a few values (majority of values are 10 and 12)

#part d
#Confidence Interval
boot.ci(boot.out, type = "perc") #95%   (10.0, 13.5 )

#part e
#Since there are no extreme scores in the data, mean is the better estimate for central 
#tendency in this case. Median would be better if there were major outliers, however, since the 
#data seems to be densely centered around 11, then mean would be a better estimate. 













