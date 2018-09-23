library(MatchIt)
data(lalonde)
library(stddiff)
install.packages("stddiff")
treat=round(abs(rnorm(100)+1)*10,0) %% 2
numeric=round(abs(rnorm(100)+1)*10,0)
binary=round(abs(rnorm(100)+1)*10,0) %% 2
category=round(abs(rnorm(100)+1)*10,0) %% 3
data=data.frame(treat,numeric,binary,category)

?stddiff
stddiff.binary(data=data,gcol=1,vcol=c(3,4))

dat = read.table("TongrenSunpengRapidCoverage.txt",header = F,sep = "\t")

attach(dat)

dat_1 = dat[dat$V2==1|dat$V1==0,]

dat_1 = dat_1[dat_1$V2!=9,]

dat_1 = dat_1[dat_1$V18==0,]

dat_1 = na.omit(dat_1)

m.out = matchit(V1 ~ V3 + V4 + V5 + V6 + V7 + V8 +V9 + V10 + V11 + V12 + V13 + V14 + V15  , data = dat_1[,1:18], method ="exact")

m.out = matchit(V1 ~ V3 + V4 + V5 + V6 + V7 + V8 +V9 + V10 + V11 + V12 + V13 + V14 + V15  , data = dat_1[,1:18], method = "subclass")

m.out = matchit(V1 ~ V3 + V4 + V5 + V6 + V7 + V8 +V9 + V10 + V11 + V12 + V13 + V14 + V15  , data = dat_1[,1:18], method = "nearest")

install.packages("optmatch")
library(optmatch)

m.out = matchit(V1 ~ V3 + V4 + V5 + V6 + V7 + V8 +V9 + V10 + V11 + V12 + V13 + V14 + V15  , data = dat_1[,1:18], method = "optimal")

m.out = matchit(V1 ~ V3 + V4 + V5 + V6 + V7 + V8 +V9 + V10 + V11 + V12 + V13 + V14 + V15  , data = dat_1[,1:18], method = "full")

m.out = matchit(V1 ~ V3 + V4 + V5 + V6 + V7 + V8 +V9 + V10 + V11 + V12 + V13 + V14 + V15  , data = dat_1[,1:18], method = "optimal")




m.out = matchit(V1 ~ V3 + V4 + V5 + V6 + V7 + V8 +V9 + V10 + V11 + V12 + V13 + V14 + V15, data = dat_1[,1:15], method ="genetic")

m.out = matchit(V1 ~ V3 + V9 +  V11  + V14 + V15, data = dat_1[,1:15], method ="genetic")

m.out = matchit(V1 ~ V3 + V9 +  V11  + V14 + V15, data = dat_1[,1:15], method ="exact")


plot(m.out) 

summary(m.out)

m.data = match.data(m.out)


?matchit


matchit(treat ~ educ + black + hispan, data = lalonde,
        method = "exact")

matchit(treat ~ educ + black + hispan, data = lalonde,
        method = "exact")

dat_1$V1
dat_1$V2
dat_1$V3
dat_1$V4
dat_1$V5
dat_1$V6
dat_1$V7
dat_1$V8
dat_1$V9
dat_1$V10
dat_1$V11
dat_1$V12
dat_1$V13
dat_1$V14
dat_1$V15
dat_1$V18

na.omit(dat_1)






###########
install.packages("rgenoud")

m.out <- matchit(treat ~ age + educ + black + hispan + married + nodegree +
                  re74 + re75, data = lalonde, method = "genetic")


z.out <- zelig(Y ~ treat + x1 + x2, model = mymodel, data = m.data)

ks.test(dat[dat_1$V1==1,3],dat[dat_1$V1==0,3])

#=================================

tmpR = c()
j=1
for (i in c(3:15)){
  ks.result = ks.test(dat_1[dat_1$V1==1,i],dat_1[dat_1$V1==0,i])
  
  tmpR[j]=ks.result$p.value
  j=j+1
}

tmpR = c()
j=1
for (i in c(3:15)){
  ks.result = ks.test(m.data[m.data$V1==1,i],m.data[m.data$V1==0,i])
  
  tmpR[j]=ks.result$p.value
  j=j+1
}

tmpR


tmpR = c()
j=1
for (i in c(3,9,11,14,15)){
  ks.result = ks.test(dat_1[dat_1$V1==1,i],dat_1[dat_1$V1==0,i])
  
  tmpR[j]=ks.result$p.value
  j=j+1
}

tmpR = c()
j=1
for (i in c(3,9,11,14,15)){
  ks.result = ks.test(m.data[m.data$V1==1,i],m.data[m.data$V1==0,i])
  
  tmpR[j]=ks.result$p.value
  j=j+1
}

tmpR

m.out = matchit(V1 ~ V3 , data = dat_1[,1:15], method ="exact")
m.data = match.data(m.out)
m.out
ks.test(m.data[m.data$V1==1,3],m.data[m.data$V1==0,3])


m.out = matchit(V1 ~ V14 , data = dat_1[,1:15], method ="exact")
m.data = match.data(m.out)
m.out
ks.test(m.data[m.data$V1==1,14],m.data[m.data$V1==0,14])

dat_modified = dat_1[rownames(m.data),]

dat_modified = dat_modified[dat_modified$V3!=4,]

tmpP=c()
tmpF=c()
j=1
for (i in c(19:22,24:26)){
  chisqDataFrame = table(dat_modified$V1,dat_modified[,i])
#  res = chisq.test(chisqDataFrame)
  tmpR[j] = chisq.test(chisqDataFrame)$p.value
  tmpF[j] = fisher.test(chisqDataFrame)$p.value
  j=j+1
}


tmpR
tmpF

table(dat_modified$V1,dat_modified[,19])
table(dat_modified$V1,dat_modified[,20])
table(dat_modified$V1,dat_modified[,21])
table(dat_modified$V1,dat_modified[,22])
table(dat_modified$V1,dat_modified[,23])
table(dat_modified$V1,dat_modified[,24])
table(dat_modified$V1,dat_modified[,25])
table(dat_modified$V1,dat_modified[,26])


t.test(as.numeric(V16)~V1,data=dat_modified)
boxplot(as.numeric(V16)~V1,data=dat_modified)

t.test(as.numeric(V17)~V1,data=dat_modified)
boxplot(as.numeric(V17)~V1,data=dat_modified)

dat_modified$V1
dat_modified$V16


table(dat_modified$V1,dat_modified[,3]) # 
table(dat_modified$v1)



chisqDataFrame = table(dat_modified[dat_modified$V3!=4,1],dat_modified[dat_modified$V3!=4,3])
chisq.test(chisqDataFrame)$p.value

nrow(dat_modified)
