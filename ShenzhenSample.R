#1：中山
#2：温一
#3：市二院

#####################zhongshan
remove(res)
b1="01-MT-T-";B1 = rep(b1,350)
b2="01-MT-C-";B2 = rep(b2,350)
b3="01-MT-S-";B3 = rep(b3,350)
set.seed(1)
a1 = sample(1:350,350)
set.seed(2)
a2 = sample(1:350,350)
set.seed(3)
a3 = sample(1:350,350)

b =cbind(a1,a2,a3)

res = cbind(c(1:350),paste(B1,b[,1],sep = ""),paste(B2,b[,2],sep = ""),paste(B3,b[,3],sep = ""))

source("Rgenetic/basicTools.R")
writeOutNoRowname(res,"zhongshan.txt")

#####################wenyi
remove(res)
b1="02-MT-T-";B1 = rep(b1,350)
b2="02-MT-C-";B2 = rep(b2,350)
b3="02-MT-S-";B3 = rep(b3,350)
set.seed(4)
a1 = sample(1:350,350)
set.seed(5)
a2 = sample(1:350,350)
set.seed(6)
a3 = sample(1:350,350)

b =cbind(a1,a2,a3)
res = cbind(c(1:350),paste(B1,b[,1],sep = ""),paste(B2,b[,2],sep = ""),paste(B3,b[,3],sep = ""))

source("Rgenetic/basicTools.R")
writeOutNoRowname(res,"wenyi.txt")


#####################shieryi
remove(res)
b1="03-MT-T-";B1 = rep(b1,350)
b2="03-MT-C-";B2 = rep(b2,350)
b3="03-MT-S-";B3 = rep(b3,350)
set.seed(7)
a1 = sample(1:350,350)
set.seed(8)
a2 = sample(1:350,350)
set.seed(9)
a3 = sample(1:350,350)

b =cbind(a1,a2,a3)

res = cbind(c(1:350),paste(B1,b[,1],sep = ""),paste(B2,b[,2],sep = ""),paste(B3,b[,3],sep = ""))

source("Rgenetic/basicTools.R")
writeOutNoRowname(res,"shieryi.txt")



#####################shanghai xinhua
#T1: 试验组；T2: 对照组
install.packages("randomizeR")
library(randomizeR)

install.packages("randomizr")
library("randomizr")

set.seed(343)
Group <- complete_ra(N=120,num_arms = 2)

Header = "Xihua"
#generate 5 Number

set.seed(343)
a3 = sample(1:10000)
no=paste(Header,"_",a3[1:120],"_",Z,sep = "")

write.csv(file = "tableRandom.csv",x = cbind(no,Group),sep = "/t")
