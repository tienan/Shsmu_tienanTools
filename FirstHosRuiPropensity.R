library(MatchIt)

dat = read.table("FirstHosRuiPropensity.txt",header = F,sep = "\t")

cor.test(dat$V1,dat$V5,method = "kendall")


dat_1 = dat[dat$V6==0,] # choose NOT TPO treatment

dat_1 = dat[dat$V7==1,] # choose < 100

ks.test(dat_1[dat_1$V5==1,1],dat_1[dat_1$V5==0,1])
?matchit


chisqDataFrame = table(dat_1[,3],dat_1[,5])
fisher.test(chisqDataFrame)$p.value
chisq.test(chisqDataFrame )


#"exact" (exact matching), "full" (full matching), "genetic" (genetic matching), "nearest" (nearest neighbor matching), "optimal" (optimal matching), and "subclass" (subclassification) are available. The default is "nearest". Note that within each of these matching methods, MatchIt offers a variety of options.
m.out = matchit(V5 ~ V1, data = dat_1, method ="optimal")  #

m.data = match.data(m.out)



ks.test(m.data[m.data$V5==1,1],m.data[m.data$V5==0,1])

chisqDataFrame = table(m.data[,3],m.data[,5])
fisher.test(chisqDataFrame)$p.value




dat_1 = dat[dat$V6==1,] # choose TPO treatment

m.out = matchit(V5 ~ V1   , data = dat_1, method ="exact")

m.data = match.data(m.out)

nrow(m.data)

ks.test(m.data[m.data$V3==1,1],m.data[m.data$V3==0,1])

chisqDataFrame = table(m.data[,3],m.data[,5])
fisher.test(chisqDataFrame)$p.value
chisqDataFrame


tmpR = c()
j=1
for (i in c(1:4)){
  ks.result = ks.test(m.data[m.data$V6==1,i],m.data[m.data$V6==0,i])
  
  tmpR[j]=ks.result$p.value
  j=j+1
}

tmpR


chisqDataFrame = table(m.data[,3],m.data[,5])
chisq.test(chisqDataFrame)$p.value

glm.fit=glm(V5~V1 + V3 +V2 + V4,data=m.data,family=binomial(link="logit"))



?predict.glm()

predict(glm.fit, m.data[,1:5],type ="response" )


summary(glm.fit)

install.packages("pROC")
install.packages("ggplot2") #下载ggplot2包
library(pROC) #调用pROC包
library(ggplot2) #调用ggplot2包以利用ggroc函数

roc(predict(glm.fit, m.data),m.data[,5])

library(ROCR)

data(ROCR.simple)

library(ROCR)
data(ROCR.simple)

## computing a simple ROC curve (x-axis: fpr, y-axis: tpr)
library(ROCR)
data(ROCR.simple)
pred <- prediction( ROCR.simple$predictions, ROCR.simple$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf)    

datROC = data.frame(as.numeric(predict(glm.fit, m.data[,1:5],type ="response" )),m.data[,5])
colnames(datROC) = c("value","label") 

roc(datROC$label,datROC$value,)  
  

pred <- prediction(as.numeric(predict(glm.fit, m.data[,1:5],type ="response" )), m.data[,5])
perf <- performance(pred,"tpr","fpr")
plot(perf,col=rainbow(10))


install.packages("pROC")
library(pROC)


