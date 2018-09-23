dat = read.table("Cancer_data.txt",header = T,sep="\t")
#计量变量的正态性判断
#BMI, HIP, WT
shapiro.test(dat[dat$Cancer==1,]$BMI)
shapiro.test(dat[dat$Cancer==0,]$BMI)

shapiro.test(dat[dat$Cancer==1,]$Hip)
shapiro.test(dat[dat$Cancer==0,]$Hip)

shapiro.test(dat[dat$Cancer==1,]$WT)
shapiro.test(dat[dat$Cancer==0,]$WT)

#样本间的比较

t.test(dat[dat$Cancer==1,]$WT,dat[dat$Cancer==0,]$WT)

ks.test(dat[dat$Cancer==1,]$BMI,dat[dat$Cancer==0,]$BMI,exact = F)
ks.test(dat[dat$Cancer==1,]$Hip,dat[dat$Cancer==0,]$Hip,exact = F)

boxplot(WT~Cancer,data=dat)
boxplot(Hip~Cancer,data=dat)
boxplot(BMI~Cancer,data=dat)

#计数资料的比较

tab = table(dat$Gender,dat$Cancer)

prop = prop.table(tab,1)

write.table(result,file="tmp.txt",quote = F,sep="\t",row.names = T,fileEncoding = "utf-8")

barplot(tab,xlab = "癌症",ylab="百分比",col=c("green","red"),legend=c("Male","Fema"))

chisq.test(tab,correct = F)

#其他比较
