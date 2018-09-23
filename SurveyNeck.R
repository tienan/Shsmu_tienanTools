neck = read.table("SurveyNeck.txt",header = T,sep = "\t")


tab = table(neck[,2],neck$outcomeNeck)

prop = prop.table(tab,1)

r = chisq.test(tab,correct = F)
s=c()
p=c()
s[1]=r$statistic
p[1]=r$p.value
j=2

result = cbind(tab[,1],prop[,1],tab[,2],prop[,2])

for (i in c(3:5,7:13)){

  tab = table(neck[,i],neck$outcomeNeck)
  
  prop = prop.table(tab,1)
  
  r = chisq.test(tab,correct = F)
  
  tmp = cbind(tab[,1],prop[,1],tab[,2],prop[,2])
  
  s[j]=r$statistic
  p[j]=r$p.value
  result = rbind(result,tmp)
  j=j+1
  
}
write.table(result,file="tmp.txt",quote = F,sep="\t",row.names = T,fileEncoding = "utf-8")


tab = table(neck$outcomeNeck,neck$Grade)
barplot(tab,xlab = "年级",ylab="百分比",col=c("green","red"),legend=c("颈部正常","颈部不适")) 

tab = table(neck$outcomeNeck,neck$school)
barplot(tab,xlab = "学校类型",ylab="百分比",col=c("green","red")) 



logit = glm(outcomeNeck~Grade+myopia+Burden+Reading,family=binomial(link='logit'),data=neck)
summary(logit)


logit = glm(outcomeNeck~Grade+myopia+Burden+Reading,family=binomial(link='logit'),data=neck[neck$Reading<4,])
summary(logit)

