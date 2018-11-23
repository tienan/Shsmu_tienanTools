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

########################## First Hospital Liping re - seed 
library("randomizr")
setwd("D:/R/")
dat = read.table("randomCodeLiping.txt")
r1 = dat[c(1:41),1]
r2 = dat[,2]

r1p = paste(r1,collapse = "")
r2p = paste(r2,collapse = "")

t=1
for (i in 1:100000){
  set.seed(i)
  Group <-  simple_ra(N=41,num_arms = 2)
  rtp = paste(Group,collapse = "")
  if (rtp==r1p){
    t=i
    break()
  }
}

r1p==r1p
#################################频数一样
library("randomizr")
setwd("D:/R/")
dat = read.table("randomCodeLiping.txt")
r1 = dat[c(1:41),1]
r2 = dat[c(1:68),2]


r1p=paste(table(r1),collapse = "")
r2p=paste(table(r2),collapse = "")


r1p = 2219
r2p = 3236


#r1p = paste(r1,collapse = "")
#r2p = paste(r2,collapse = "")

t=1
for (i in 1:1000000){
  set.seed(i)
  Group <-  simple_ra(N=68,num_arms = 2)
  table(Group)
  rtp = paste(table(Group),collapse = "")
#  print(r2p)
#  print(rtp)
  if (rtp==r2p){
    t=i
    break()
  }
}



t=1
for (i in 1:10000){
  set.seed(i)
  Group <-  simple_ra(N=41,num_arms = 2)
  table(Group)
  rtp = paste(table(Group),collapse = "")
  #  print(r2p)
  #  print(rtp)
  if (rtp==r1p){
    t=i
    break()
  }
}



####41入组的样本数
set.seed(18)
Group_41 <-  simple_ra(N=182,num_arms = 2)
#table(Group_41)
table(Group_41[1:41])
####67入组的样本数
set.seed(2)
Group_67 <-  simple_ra(N=182,num_arms = 2)  
table(Group_67[1:68])
write.csv(file = "LipingRandom.csv",x = data.frame(Group_41,Group_67))


#甲癌测序的随机分配表
#影响因素：不同批次、试剂差异



library("randomizr")
setwd("D:/R/")
data(HairEyeColor)
HairEyeColor <- data.frame(HairEyeColor)
hec <- HairEyeColor[rep(1:nrow(HairEyeColor),
                        times = HairEyeColor$Freq), 1:3]
Z <- block_ra(blocks = hec$Hair)
table(Z, hec$Hair)
blocks <- with(hec, paste(Hair, Eye, Sex, sep = "_"))

dat = read.table("THCA_random.txt",header = T,sep = "\t")
blocks <- with(dat, paste(Batch,Kit,sep = "_"))

set.seed(343)
Z <- block_ra(blocks = blocks,num_arms = 3,prob_each = c(39/78,27/78,12/78))
table(blocks, Z)
Header = "THCA"
set.seed(343)
a3 = sample(1:10000,78)
no=paste(Header,"_",a3,sep = "")
write.csv(file = "THCA_tableRandom.csv",x = cbind(no,Z,dat$Batch,dat$Kit))

                                           
                                           
                                           