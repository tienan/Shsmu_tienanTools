# test
dat_1 = read.table("tmp",header = F, sep =  "\t")
head(dat_1)
datReso = dat_1

dat[-1,] # genes

# Main
source("Hosfunction.R")


#########test
for (k in 2:length(itemGap)){
print((itemGap[k-1]+1):itemGap[k])
}

############
#the format of  input file 
#rank
dat = read.table("principle_2018-rank.txt",header = F, sep =  "\t")
dat[is.na(dat)] <- 0
nrow(dat)
year = c("2011","2012","2013")

# 
dat = read.table("principle_rank_single.txt",header = F, sep =  "\t")
printRow(as.matrix(scoreHos(as.data.frame(dat[,1]))))
nrow(as.matrix(scoreHos(as.data.frame(dat[,1]))))
nrow(dat)

#Start 

source("Hosfunction.R")
dat = read.table("principle_2018.txt",header = F, sep =  "\t")
dat[is.na(dat)] <- 0
nrow(dat)


head(dat)
nrow(dat)
#year = c("2011","2012","2013")
year = c(2000)
class  = 1:7
#start = 5
itemGap = c(8,10,13,16,18,21,25,26,29,36,41,43,44)
itemGap = itemGap+5
itemGap = c(5,itemGap)
itemGap
resutTotal = data.frame(c(1:495))
l = 1 
m =  1
id = c()
datDealt=data.frame()
#
for (k in 2:length(itemGap)){
  resultsTmp = data.frame()
  for (j in 1:length(year)){
    for (i in 1:length(class)){
#      datTmp = dat[dat$V3==year[j]&dat$V5==class[i],]
      datTmp = dat[dat$V5==class[i],]
      datDealt = rbind(datDealt,datTmp)
#     id[m] = paste(year[j], class[i],collapse = "",sep = "_")
#     m=m+1
      #nrow(datTmp)
      #resource 5+
      datReso = as.data.frame(datTmp[,(itemGap[k-1]+1):itemGap[k]])
      #start=start+itemGap[k]
      result=as.data.frame(scoreHos(datReso))
      
      resultsTmp= rbind(resultsTmp,result)
    }
  }
  resutTotal=cbind(resutTotal,resultsTmp)
}


resutTotal_1=cbind(resutTotal,datDealt[c(1:495),])

?write.table
write.table(resutTotal_1[,-1],file = "resultItem.txt",
            quote = F,sep="\t",row.names = F)

resutTotal_2 = resutTotal_1[,-1]
#sd = apply(resutTotal,2,sd)


#dep = c(1,2)
gap = c(0,3,8,13)

dat_internal_medicine = resutTotal_2[resutTotal_2$V5==1,]
dat = dat_internal_medicine

dat_s = apply(dat[,c(1:3)],2,sd)
dat_p = apply(dat[,c(4:8)],2,sd)
dat_r = apply(dat[,c(9:13)],2,sd)

structureS=c()
pross=c()
res = c()


j=1
for (i in 1:nrow(dat_internal_medicine)){
  structureS[j] = sum(dat[j,c(1:3)] * dat_s)/sum(dat_s)
  j=j+1
}
j=1
for (i in 1:nrow(dat)){
  pross[j] = sum(dat[j,c(4:8)] * dat_p)/sum(dat_p)
  j=j+1
}
j=1
for (i in 1:nrow(dat)){
  res[j] = sum(dat[j,c(9:13)] * dat_r)/sum(dat_r)
  j=j+1
}

Level_1S = apply(cbind(structureS,pross,res),2,sd)


tscore = c()
j=1
for (i in 1:nrow(dat)){
  tscore[j] = sum(cbind(structureS,pross,res)[j,c(1:3)] * Level_1S)/sum(Level_1S)
  j=j+1
}

dat_internal_medicine=cbind(tscore,structureS,pross,res,dat)


#####


dat_surgery_medicine = resutTotal_2[resutTotal_2$V5==2,]
dat = dat_surgery_medicine

dat_s = apply(dat[,c(1:3)],2,sd)
dat_p = apply(dat[,c(4:8)],2,sd)
dat_r = apply(dat[,c(9:13)],2,sd)

structureS=c()
pross=c()
res = c()


j=1
for (i in 1:nrow(dat)){
  structureS[j] = sum(dat[j,c(1:3)] * dat_s)/sum(dat_s)
  j=j+1
}
j=1
for (i in 1:nrow(dat)){
  pross[j] = sum(dat[j,c(4:8)] * dat_p)/sum(dat_p)
  j=j+1
}
j=1
for (i in 1:nrow(dat)){
  res[j] = sum(dat [j,c(9:13)] * dat_r)/sum(dat_r)
  j=j+1
}

Level_1S = apply(cbind(structureS,pross,res),2,sd)


tscore = c()
j=1
for (i in 1:nrow(dat)){
  tscore[j] = sum(cbind(structureS,pross,res)[j,c(1:3)] * Level_1S)/sum(Level_1S)
  j=j+1
}

dat_surgery_medicine=cbind(tscore,structureS,pross,res,dat)

resultSurInter = rbind(dat_surgery_medicine,dat_internal_medicine)
ncol(dat_surgery_medicine)
ncol(dat_internal_medicine)

write.table(resultSurInter ,file = "resultItem.txt",
            quote = F,sep="\t",row.names = F)

nrow((resultSurInter))

########plot






#######calculate the score of each index Level 3
dat = read.table("principle_2018.txt",header = F, sep =  "\t")
datTmp = dat[dat$V5==1,]
head(datTmp)
source("Hosfunction.R")
score(dat[,c(6:ncol(dat))])
write.table(score(dat[,c(6:ncol(dat))]),file = "resultItemScore.txt",
            quote = F,sep="\t",row.names = F)
#######


# sum  calculation 

nrow(resutTotal)
ncol(resutTotal)
head(resutTotal)

#resutTotal = resutTotal[,-1]



#
dat = read.table("principle_20180_index.txt",header = F,sep = "\t")
#
resutTotal = dat

sdResult = apply(resutTotal, 2, sd)

weight = sdResult/sum(sdResult)

sumResult = (as.matrix(resutTotal) %*% as.matrix(weight))

for (i in 1:nrow(sumResult)){
  cat(sumResult[i])
  cat("\n")
}

# End
#weight calculation
ncol = ncol(datReso)
score = data.frame()
datTmp = data.frame()
datTmp = datReso

k = 1 # score count
for (i in 1:ncol){
  tmp = table(datTmp[,i])
  tmpScore=c()
  for (j in 1:length(tmp)){
    tmpScore[j] = (sum(tmp[1:j])/sum(tmp))
  }
  for (j in 1:length(datTmp[,i])){
    score[j,k] = tmpScore[names(tmp)==datTmp[j,i]]
  }
  k = k+1
  # tmpScore = tmpScore/max(tmpScore)
}
mean = apply(score, 2, mean)
weight = (max(mean)-mean)/(max(mean) - min(mean))
#apply(score, 2, sd)
result=c()
for (i in 1:nrow(score)){
  result[i] = sum(score[i,] * weight)/sum(weight)
}
for (i in 1:nrow(score)){
  cat(result[i])
  cat("\n")
}
plot(apply(score, 2, mean),  apply(score, 2, sd) )
cor.test(apply(score, 2, mean),  apply(score, 2, sd) )

length(result)



main = function(){
  dat = read.table("evaluatedDiscipline.txt",header = T,sep = "\t")
  score = data.frame() # Score File
  k = 1 # score count
  nrow = nrow(dat)
  nrow
  ncol = ncol(dat)
  ncol
  dat[,36]
  colnames(dat)
  # for (i in 9:ncol){
  #   table()
  # }
#  i = 9
  for (i in 9:ncol){
    tmp = table(dat[,i])
    tmpScore=c()
    for (j in 1:length(tmp)){
      tmpScore[j] = (sum(tmp[1:j])/sum(tmp))
    }
    for (j in 1:length(dat[,i])){
      score[j,k] = tmpScore[names(tmp)==dat[j,i]]
    }
    k = k+1
   # tmpScore = tmpScore/max(tmpScore)
  }
  mean = apply(score, 2, mean)
  weight = (max(mean)-mean)/(max(mean) - min(mean))
  #apply(score, 2, sd)
  result=c()
  for (i in 1:nrow(score)){
    result[i] = sum(score[i,] * weight)/sum(weight)
  }
  for (i in 1:nrow(score)){
    cat(result[i])
    cat("\n")
  }
  plot(apply(score, 2, mean),  apply(score, 2, sd) )
  cor.test(apply(score, 2, mean),  apply(score, 2, sd) )
  # i = 10
  # tmp = table(dat[,i])
  # tmpScore=c()
  # for (j in 1:length(tmp)){
  #   tmpScore[j] = (sum(tmp[1:j])/sum(tmp))
  # }
  # tmpScore = tmpScore/max(tmpScore)

  
   
}

scoringTable = function(datFram){
  ncol  = ncol(datFram)
}

scoring = function (){
  
  
}