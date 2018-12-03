#analysis model 

genesName = levels(as.factor(res$gene_id))

expValue = data.frame()

for (i in 1:length(geneName)){
  t = res[res$gene_id==geneName[i],]
  expValue[i,1]=geneName[i]
  expValue[i,2]=mean(t$normalized_count)
  expValue[i,3]=sd(t$normalized_count)
  expValue[i,4] = sd(t$normalized_count)/mean(t$normalized_count)
} 
#filter using CV coefficient
candidatGenes = expValue[expValue$V4<23.95 & expValue$V4>1,1] 


candidatGenesExp = res[res$gene_id==candidatGenes[1],]$normalized_count

for (i in 2:length(candidatGenes)){
  candidatGenesExp=cbind(candidatGenesExp,res[res$gene_id==candidatGenes[i],]$normalized_count)
}

candidatGenesExp=cbind(candidatGenesExp,res[res$gene_id==pdia3p[1,1],]$normalized_count)
candidatGenesExp=cbind(candidatGenesExp,res[res$gene_id==hmgb1[1,1],]$normalized_count)
hmgb1 = res[grepl("hmgb1",res$gene_id,ignore.case=T),c(3,4)]  
hmgb1 = res[grepl("pdia3p",res$gene_id,ignore.case=T),c(3,4)]  

na.omit(c(candidatGenes,hmgb1[1,1],pdia3p[1,1])) -> genens
colnames(candidatGenesExp)[c(8576,8577)] = c(pdia3p[1,1],hmgb1[1,1])


tmp = candidatGenesExp[,colnames(candidatGenesExp)%in%as.vector(genens)]





# using reg to find the gene

grepl("YY1",rownames(corMat))






write.table(res,file="LUAD_methylation.txt",quote = F,sep="\t",row.names = F,fileEncoding = "utf-8")

write.table(expValue,file="LUAD_expValue.txt",quote = F,sep="\t",row.names = F,fileEncoding = "utf-8")

install.packages("apcluster")

library(apclustera)





s = corSimMat(t(tmp), sel=NA, r=1, signed=TRUE, method="pearson")

apres1b = apcluster(s)


tiff(filename = "heatmap_luad.tif",width = 3000, height = 3000, units = "px", pointsize = 12,compression ="lzw",bg = "white", res = 300)
heatmap(apres1b,s)
dev.off()

# remove NA 
x[apply(x, 1, function(x) !all(is.na(x))),]
sign = ifelse(candidatGenesExp[,8577]>median(candidatGenesExp[,8576]),1,0)

sign = ifelse(dat[,8243]>median(dat[,8243]),1,0)


#
diff=data.frame()
for (i in  1:(8245-1)){
  r = ks.test(dat[dat$sign==1,i],dat[dat$sign==0,i])
  d = log(mean(dat[dat$sign==1,i])/mean(dat[dat$sign==0,i]))
  s =  r$statistic 
  p =  r$p.value
  diff[i,1] = d
  diff[i,2] = s
  diff[i,3] = p
}

names = colnames(dat)
rownames(diff) = names[1:length(names)-1]

diff[diff[,3]<0.05&abs(diff[,1])>0.5,]


##########################reshape dataframe 

colnum = nrow(candidatedGenesExp_s)/length(levels(as.factor(candidatedGenesExp_s[,1])))

rownum = length(levels(as.factor(candidatedGenesExp_s[,1])))

ind = seq(1,nrow(candidatedGenesExp_s),colnum)

SpecialGenesExp = candidatedGenesExp_s[seq(1,nrow(candidatedGenesExp_s),colnum),c(1,3)]

genesName = candidatedGenesExp_s[c(1:colnum),2]

for (i in 2:colnum){
  SpecialGenesExp = cbind(SpecialGenesExp,candidatedGenesExp_s[seq(i,nrow(candidatedGenesExp_s),colnum),3])
}

colnames(SpecialGenesExp) = c("sampleID",genesName)

# reshapge
resShape = res[seq(1,nrow(res)-1,nrow(res)/length(levels(as.factor(res[,1])))),c(1,3)]
len = nrow(res)/length(levels(as.factor(res[,1])))
for (i in 2:len){
  resShape = cbind(resShape,res[seq(i,nrow(res),nrow(res)/length(levels(as.factor(res[,1])))),3])
}
colnames(resShape) = c("sampleID",res[c(1:len),2]) 


resShape = resShape[order(resShape[,1]),]

cancerSign = clinic_s[,"cancer"]==1

datLUAD = annData.m[cancerSign,]

datLUAD = cbind(datLUAD,resShape[cancerSign,])

#continue 2 count
resShapeCount=c()
resShapeCount=resShape[,1]
for(i in 2:ncol(resShape)){
  median = median(resShape[,i])
  tmp = resShape[,i]
  tmp[resShape[,i]>median]=1
  tmp[resShape[,i]<=median]=0
  resShapeCount=cbind(resShapeCount,tmp)
}
colnames(resShapeCount) = colnames(resShape)


datLUAD$state=0
#datLUAD[datLUAD$vital_status=="Alive",]$state=0
datLUAD[datLUAD$vital_status=="Dead",]$state=1
######Survial
#remove(my.surv)
my.surv <- Surv(datLUAD$survialDat,datLUAD$state)



diffSurvRes = data.frame()
name = c()
for (i in 1:(ncol(datLUAD)-11-1)){
  diffSurvRes[i,1] = 0
  diffSurvRes[i,2] = 0
  if( length(levels(as.factor(resShapeCount[cancerSign,i+1])))>1){
  diffSurv = survdiff(my.surv ~ 
                        resShapeCount[cancerSign,i+1]
#                     datLUAD$pathologic_stage+
#                      datLUAD$tobacco_smoking_history+
#                     datLUAD$age_at_initial_pathologic_diagnosis+
#                      datLUAD[,i+11]
#  ,
#                    data=datLUAD
)
  diffSurvRes[i,1] = diffSurv$chisq
  diffSurvRes[i,2] = 1-pchisq(diffSurv$chisq,1)
  }
}

rowNames = colnames(datLUAD) 
colNames = c("Chisq","P","geneName")

diffSurvRes$name=rowNames[12:(length(colnames(datLUAD))-1)]

colnames(diffSurvRes)=colNames
writeOut(diffSurvRes,"LUADDiff")




orResult=data.frame()
names = c()
for (i in 1:(ncol(datLUAD)-11)){
coxResult = coxph(my.surv ~ 
                    datLUAD$gender+
                    datLU?AD$pathologic_stage+
                    datLUAD$tobacco_smoking_history+
                    datLUAD$age_at_initial_pathologic_diagnosis+
                    datLUAD[,i+11],
                  data=datLUAD)
tmp = summary(coxResult)
#ORvalue
orResult[i,1] = as.numeric(tmp$coefficients[nrow(tmp$coefficients),1])
#Z
orResult[i,2] = as.numeric(tmp$coefficients[nrow(tmp$coefficients),4])
#P
orResult[i,3] = as.numeric(tmp$coefficients[nrow(tmp$coefficients),5])
#names

}

rowNames = colnames(datLUAD) 
colNames = c("OR","Z","P","geneName")

orResult$name=rowNames[12:length(colnames(datLUAD))]

colnames(orResult)=colNames
writeOut(orResult,"LUADOR")
library(heatmap3)


###CancerSign
diff = data.frame()
for (i in  2:ncol(SpecialGenesExp)){
  r = ks.test(SpecialGenesExp[clinic_s$cancer==1,i],SpecialGenesExp[clinic_s$cancer==11,i])
  d = log(mean(SpecialGenesExp[clinic_s$cancer==1,i])/mean(SpecialGenesExp[clinic_s$cancer==11,i]))
  s =  r$statistic 
  p =  r$p.value
  diff[i,1] = d
  diff[i,2] = s
  diff[i,3] = p
}
names = colnames(SpecialGenesExp)
diff = diff[-1,]
rownames(diff) = names[2:length(names)]

sigDiff = diff[diff[,3]<0.05&abs(diff[,1])>0.5,]

write.table(sigDiff,quote = F,sep="\t",row.names = F,fileEncoding = "utf-8")



#key term
annData = data.frame(clinic_s[,c("gender","pathologic_T","pathologic_N","pathologic_M","pathologic_stage","tobacco_smoking_history","vital_status","days_to_last_followup","days_to_death","age_at_initial_pathologic_diagnosis","cancer")])
for(i in 1:7){
  annData[annData[,i] == "[Not_Available]",i] = 'Not_Sure'
  annData[annData[,i] == "[Unknown]",i] = 'Not_Sure'
  annData[annData[,i] == "[Discrepancy]",i] = 'Not_Sure'
  annData[,i]=as.factor(annData[,i]);
}

#annData[annData[,6]=="[Not_Available]",6] = ''
#annData[annData[,6]=="[Unknown]",6] = ''


annData[annData[,8]=="[Not_Available]",8] = 0
annData[,8] = as.numeric(annData[,8])

annData[annData[,9]=="[Not_Applicable]",9] = 0
annData[,9] = as.numeric(annData[,9])

annData[annData[,10]=="[Not_Available]",10] = ''
annData[,10] = as.numeric(annData[,10])

annData[,11] = as.factor(annData[,11])

annData$survialDat = annData[,8]+annData[,9]


annData.m = annData[,c(-8,-9)]


#heatmap3(t(log10(SpecialGenesExp.m1+0.01)))

tiff(filename = "heatmap3.tif",
     width = 10000, height = 20000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 300)
heatmap3(t(log10(SpecialGenesExp.m1+0.01)),ColSideCut=0.85,ColSideAnn=annData.m,ColSideFun=function(x)
showAnn(x),ColSideWidth=1.2,balanceColor=TRUE)

dev.off()


annData[annData=='[Not_Available]']="NA"
annData[annData$days_to_last_followup=='NA',]$days_to_last_followup = 0
annData[annData=='[Not_Applicable]']=0

ann<-data.frame(mtcars[,c("mpg","am","wt","gear")])
ann[,2]<-as.factor(ann[,2])
ann[,4]<-as.factor(ann[,4])



####################survial analysis


pdia3p = SpecialGenesExp.m1[,grepl("pdia3p",colnames(SpecialGenesExp.m1),ignore.case = T)]
pdia3pSign = ifelse(pdia3p>median(pdia3p),1,0)
pdia3pClinic = cbind(annData.m,pdia3p,pdia3pSign)
pdia3pClinic$status = ifelse(pdia3pClinic$vital_status=='Dead',1,0)
boxplot(data=pdia3pClinic,pdia3p~ pathologic_stage)
boxplot(data=pdia3pClinic,pdia3p~ pathologic_T)
boxplot(data=pdia3pClinic,pdia3p~ pathologic_N)
boxplot(data=pdia3pClinic,pdia3p~ pathologic_M)

c



heatmap3(t(log10(pdia3p)),ColSideCut=0.85,ColSideAnn=annData.m,ColSideFun=function(x)
  showAnn(x),ColSideWidth=1.2,balanceColor=TRUE)





