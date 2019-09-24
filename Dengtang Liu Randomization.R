d# Dengtang Liu Randomization

install.packages("randomizeR")
library("randomizeR")

params <- crPar(100)
rs <- genSeq(params)
rs
getRandList(rs)
plotSeq(rs, plotAllSeq = TRUE)
allSeqs <- getAllSeq(params)



install.packages("randomizr")
library("randomizr")
data(HairEyeColor)
HairEyeColor <- data.frame(HairEyeColor)
hec <- HairEyeColor[rep(1:nrow(HairEyeColor),
                        times = HairEyeColor$Freq), 1:3]

N <- nrow(hec)

# Fix the rownames
rownames(hec) <- NULL


# Set a seed for reproducability
set.seed(343)

# Create untreated and treated outcomes for all subjects
hec <- within(hec,{
  Y0 <- rnorm(n = N,mean = (2*as.numeric(Hair) + -4*as.numeric(Eye) + -6*as.numeric(Sex)), sd = 5)
  Y1 <- Y0 + 6*as.numeric(Hair) + 4*as.numeric(Eye) + 2*as.numeric(Sex)
})

# Calculate true ATE
with(hec, mean(Y1 - Y0))
#> [1] 25


Z <- simple_ra(N = N)
table(Z)

hec$Hair

hec <- within(hec,{
  Z_blocked <- complete_ra(N = N, m_each = c(100, 200, 292),
                           conditions = c("control", "placebo", "treatment"))
  id_var <- 1:nrow(hec)
})



#https://cran.r-project.org/web/packages/randomizr/vignettes/randomizr_vignette.html


#T1: 奥痰平；T2：利培酮；T3：氨磺必利；T4：阿立哌；T5：奋乃静
#7个中心

library("randomizr")
?seq
setwd("D:/R/")
# 1310 in all
hospital = c(rep(1,time=165),rep(2,time=165),rep(3,time=165),
             rep(4,time=165),rep(5,time=165),rep(6,time=165),
             rep(7,time=165),rep(8,time=165))


set.seed(343)
Z <- block_ra(blocks = hospital, num_arms = 5)
table(Z,hospital)
Header = "SMART"
#generate 5 Number

set.seed(343)
a3 = sample(1:10000,length(hospital))
no=paste(Header,"_","Phase1","_",a3,sep = "")
x = cbind(no,hospital,Z)
colnames(x)=c("No","Hospital","Group")
x = as.data.frame(x)
head(x)
table(x$Group)

write.table(file = "Dengtang_phase1_tableRandom.txt",x,sep = "\t")

#Phase_II - obey
#G =   floor(sum(as.numeric(table(x$Group)))/5*0.6+5)


#hospital = c(rep(1,time=25),rep(2,time=25),rep(3,time=25),rep(4,time=25),rep(5,time=25),rep(6,time=25),rep(7,time=25))

hospital =c(rep(1,time=550))

tmp = data.frame()
for (i in 1:4){
  hospital =c(rep(i,time=550/5))
  set.seed(i)
  Z <- block_ra(blocks = hospital, num_arms = 4,prob_each = c(1,1,1,1.6)/sum(c(1,1,1,1.6)))
  #sugroup order: 1. 奥氮平; 2. 利培酮; 3. 氨磺必利; 4. 阿立哌唑; 5. 氯氮平
  ?block_ra
  table(Z,hospital)
  set.seed(343+i)
  a3 = sample(1:10000,length(hospital))
  no=paste(Header,"_","Phase2","_",a3,sep = "")
  x = cbind(no,hospital,Z)
  colnames(x)=c("No","Group","subGroup")
  tmp = rbind(tmp,x)
}
hospital =c(rep(5,time=550/5))
set.seed(5)
Z <- block_ra(blocks = hospital, num_arms = 5,prob_each = c(1,1,1,1,1.6)/sum(c(1,1,1,1,1.6)))
#sugroup order: 1. 奥氮平; 2. 利培酮; 3. 氨磺必利; 4. 阿立哌唑; 5. 氯氮平
?block_ra
set.seed(343+5)
a3 = sample(1:10000,length(hospital))
table(Z,hospital)
no=paste(Header,"_","Phase2","_",a3,sep = "")
x = cbind(no,hospital,Z)
colnames(x)=c("No","Group","subGroup")
tmp = rbind(tmp,x)
#sugroup order: 奋乃静
tail(tmp)

write.table(file = "Dengtang_phase2_tableRandom_obey.txt",tmp,sep = "\t")

#Phase_II - Not obey
hospital =c(rep(1,time=130))

tmp = data.frame()
for (i in 1:4){
  hospital =c(rep(i,time=130/5))
  set.seed(i+10)
  Z <- block_ra(blocks = hospital, num_arms = 3,prob_each = c(1,1,1)/sum(c(1,1,1)))
  #sugroup order: 1. 奥氮平; 2. 利培酮; 3. 氨磺必利; 4. 阿立哌唑; 
  ?block_ra
  table(Z,hospital)
  set.seed(343+10+i)
  a3 = sample(1:10000,length(hospital))
  no=paste(Header,"_","Phase2","_",a3,sep = "")
  x = cbind(no,hospital,Z)
  colnames(x)=c("No","Group","subGroup")
  tmp = rbind(tmp,x)
}
hospital =c(rep(5,time=130/5))
set.seed(5)
Z <- block_ra(blocks = hospital, num_arms = 4,prob_each = c(1,1,1,1)/sum(c(1,1,1,1)))
#sugroup order: 1. 奥氮平; 2. 利培酮; 3. 氨磺必利; 4. 阿立哌唑; 5. 氯氮平
?block_ra
table(Z,hospital)
seed(343+10+5)
a3 = sample(1:10000,length(hospital))
no=paste(Header,"_","Phase2","_",a3,sep = "")
x = cbind(no,hospital,Z)
colnames(x)=c("No","Group","subGroup")
tmp = rbind(tmp,x)
#sugroup order: 奋乃静
write.table(file = "Dengtang_phase2_tableRandom_Not_obey.txt",tmp,sep = "\t")




#Phase_III 
hospital =c(rep(1,time=190))
set.seed(6)
Z <- block_ra(blocks = hospital, num_arms =2)
#Group order: 1. add; 2. Injection
?block_ra
table(Z,hospital)
seed(343+20+1)
a3 = sample(1:10000,length(hospital))
no=paste(Header,"_","Phase3","_",a3,sep = "")
x=data.frame()
x = cbind(no,Z)
colnames(x)=c("No","Group")
write.table(file = "Dengtang_phase3_tableRandom.txt",x,sep = "\t")

#import list 


#the 20190802 modify
#sugroup order: 1. 奥氮平; 2. 利培酮; 3. 氨磺必利; 4. 阿立哌唑; 5. 氯氮平
set.seed(343+10+5)
a3 = sample(1:10000,312)
no=paste("SmartCat","_","Phase2","_",a3,sep = "")
set.seed(3)
Z_blocked <- complete_ra(N = 312, m_each = c(52, 52, 52,52,104), conditions = c("T1", "T2", "T3","T4","T5"))
x = cbind(no,Z_blocked)
table(Z_blocked)
write.csv(file = "20190801Dengtang_phase2_tableRandom.csv",x)



###############
library("randomizr")
hospital =c(rep(1,time=550))

tmp = data.frame()
for (i in 1:4){
  hospital =c(rep(i,time=550/5))
  set.seed(i)
  Z <- block_ra(blocks = hospital, num_arms = 4,prob_each = c(1,1,1,2)/sum(c(1,1,1,2)))
  #sugroup order: 1. 奥氮平; 2. 利培酮; 3. 氨磺必利; 4. 阿立哌唑; 5. 氯氮平
  ?block_ra
  table(Z,hospital)
  set.seed(343+i)
  a3 = sample(1:10000,length(hospital))
  no=paste("SmartCat","_","Phase2","_",a3,sep = "")
  x = cbind(no,hospital,Z)
  colnames(x)=c("No","Group","subGroup")
  tmp = rbind(tmp,x)
}
#1.奋乃静
hospital =c(rep(5,time=550/5))
set.seed(5)
Z <- block_ra(blocks = hospital, num_arms = 5,prob_each = c(1,1,1,1,2)/sum(c(1,1,1,1,2)))
?block_ra
table(Z,hospital)
set.seed(343+5)
a3 = sample(1:10000,length(hospital))
no=paste("SmartCat","_","Phase2","_",a3,sep = "")
x = cbind(no,hospital,Z)
colnames(x)=c("No","Group","subGroup")
tmp = rbind(tmp,x)
write.csv(file = "20190801-1Dengtang_phase2_tableRandom.csv",tmp)

#############################################


