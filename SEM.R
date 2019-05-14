#SEM
install.packages("lavaan")
library("lavaan")
HS.model <- 'visual =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed =~ x7 + x8 + x9'
fit <- cfa(HS.model, data = HolzingerSwineford1939)
parTable(fit)
summary(fit, fit.measures=F)


model <- '
# measurement model
ind60 =~ x1 + x2 + x3
dem60 =~ y1 + y2 + y3 + y4
dem65 =~ y5 + y6 + y7 + y8
# regressions
dem60 ~ ind60
dem65 ~ ind60 + dem60
# residual correlations
y1 ~~ y5
y2 ~~ y4 + y6
y3 ~~ y7
y4 ~~ y8
y6 ~~ y8
'
fit <- sem(model, data=PoliticalDemocracy)
summary(fit, standardized=TRUE)



getwd()
dat = read.table("C:/Users/tienan/Documents/R/20190502SEM.txt",header = F,sep = "\t")
model <- '
# measurement model
demo =~ V1 + V2 + V3
livingHabit =~ V4 + V5 + V6 + V7 + V8 
MCI =~ V9 + V10 + V11 + V12 
# regressions
MCI ~ demo + livingHabit
'
fit <- sem(model, data=dat)
summary(fit, standardized=TRUE)
