#source("C:/Users/Suvalue/Documents/LAPRealWorld/tools.R")

#dataOutDemoLab = dataIntergrateOut() 
#dataInDemoLab = dataIntergrateIn()
library(survival)

head(dataOutDemoLab)

dataDisOutSustain = read.csv("tating_disease_OUT(CHECK).csv")
dataDisOutNonSustain = read.csv("non_tating_disease_OUT(CHECK).csv")

ncol(dataDisOutSustain)

ncol(dataDisOutNonSustain)

dataDisOut =rbind(dataDisOutSustain[,-5],dataDisOutNonSustain)


intersect(dataDisOutSustain$MASTER_INDEX,dataDisOutNonSustain$MASTER_INDEX)

dataDisInSustain = read.csv("tating_disease_IN(CHECK).csv")
dataDisInNonSustain = read.csv("non_tating_disease_IN(CHECK).csv")

intersect(dataDisInSustain$MASTER_INDEX,dataDisInNonSustain$MASTER_INDEX)

ncol(dataDisInSustain)

ncol(dataDisInSustain)


dataDisIn =rbind(dataDisInSustain,dataDisInNonSustain)
nrow(dataDisIn)






#drugAnalysisRes = drugAnalysisRes()

names(dataOutDemoLab)

#dataIntegration_out
dataOutDemoLabIntegration = dataOutDemoLab$LPA
head(dataOutDemoLab$LDL_C)
LDL_CTemp = dataOutDemoLab$LDL_C[,c(1,2,3,4,10)]
head(LDL_CTemp)
dataOutDemoLabIntegration =  merge(dataOutDemoLabIntegration,LDL_CTemp,
                                   by.x = "MASTER_INDEX",by.y="MASTER_INDEX")
colnames(dataOutDemoLabIntegration)
#dataOutDemoLabIntegration = dataOutDemoLabIntegration[,-c(11:18)]
names(dataOutDemoLab)
#table(dataOutDemoLab$LPA$tating)
for ( i in 2:length(dataOutDemoLab)){
  dataOutDemoLabIntegration = merge(dataOutDemoLabIntegration,dataOutDemoLab[[i]][,c(1,2)],
                                    by.x = "MASTER_INDEX",by.y="MASTER_INDEX")
}
dataOutDemoLabIntegration = merge(dataOutDemoLabIntegration,dataDisOut,
                                  by.x = "MASTER_INDEX",by.y="MASTER_INDEX")

nrow(dataOutDemoLabIntegration)
head(dataOutDemoLabIntegration)
colnames(dataOutDemoLabIntegration)
colnames(dataOutDemoLabIntegration) [15:21]= c("LDL_C","HDL_C","APO_A","APO_B","TC" ,"TG", "CRP")
colnames(dataOutDemoLabIntegration)
dataOutDemoLabIntegration = dataOutDemoLabIntegration %>%  mutate(outcome.x=case_when(
  outcome.x == "1" ~ 1,
  TRUE ~ 0
))

dataOutDemoLabIntegration = dataOutDemoLabIntegration %>%  mutate(outcome.y=case_when(
  outcome.y == "1" ~ 1,
  TRUE ~ 0
))

head(dataOutDemoLabIntegration)

dataOutDemoLabIntegration$outcome.x 


dataInDemoLabIntegration = dataInDemoLab$LPA
LDL_CTemp = dataInDemoLab$LDL_C[,c(1,2,3,4,10)]

dataInDemoLabIntegration =  merge(dataInDemoLabIntegration ,LDL_CTemp,
                                  by.x = "MASTER_INDEX",by.y="MASTER_INDEX")
colnames(dataInDemoLabIntegration)
#dataInDemoLabIntegration = dataInDemoLabIntegration[,-c(11:18)]

for ( i in 2:length(dataInDemoLab)){
  dataInDemoLabIntegration = merge(dataInDemoLabIntegration,dataInDemoLab[[i]][,c(1,2)],
                                   by.x = "MASTER_INDEX",by.y="MASTER_INDEX")
}
dataInDemoLabIntegration = merge(dataInDemoLabIntegration,dataDisIn,
                                 by.x = "MASTER_INDEX",by.y="MASTER_INDEX")

nrow(dataInDemoLabIntegration)
colnames(dataOutDemoLabIntegration)

colnames(dataInDemoLabIntegration)

colnames(dataInDemoLabIntegration) [15:21]= c("LDL_C","HDL_C","APO_A","APO_B","TC" ,"TG", "CRP")

dataInDemoLabIntegration = dataInDemoLabIntegration %>%  mutate(outcome.x=case_when(
  TRUE ~ 0
))

dataInDemoLabIntegration = dataInDemoLabIntegration %>%  mutate(outcome.y=case_when(
  outcome.y == "1" ~ 1,
  TRUE ~ 0
))
head(dataInDemoLabIntegration)


sql = "
Select distinct a.MASTER_INDEX,datediff(MONTH,min(a.DIAG_TIME),max(a.DIAG_TIME)) month
from LPA_IN_0_5Years_2_Time a
join (SELECT * from   LDL_LAB_IN a where a.I_ITEM_ENAME like 'LPA'   ) b  on a.REG_CODE=b.REG_CODE 
group by a.MASTER_INDEX
order by a.MASTER_INDEX
"

followUPin= tbl_df(
  sqlQuery(cn,
           sql) )

followUP = rbind(followUPin,followUPout)
followUP
nrow(followUP)
tmpName = colnames(followUP)

by_master_index <- group_by(followUP, MASTER_INDEX)

followUP <- by_master_index %>% 
  summarise(month= max(month)
  )
followUP

##############################ALL No PSM dataset


##################################

nrow(dataOutDemoLabIntegration)
nrow(dataInDemoLabIntegration)
head(dataOutDemoLabIntegration)
head(dataInDemoLabIntegration)
7641+27944

nrow(dataOutDemoLabIntegrationPsm)
nrow(dataInDemoLabIntegrationPsm)
4888+ 16488





dataAllDemoLabIntegration = rbind(dataOutDemoLabIntegration,dataInDemoLabIntegration)
intersect(dataOutDemoLabIntegration$MASTER_INDEX,dataInDemoLabIntegration$MASTER_INDEX)

dataALLDemoLabIntegration_1 = merge(dataAllDemoLabIntegration,followUP,by.x = "MASTER_INDEX",by.y = "MASTER_INDEX")
nrow(dataALLDemoLabIntegration_1)
nrow(dataAllDemoLabIntegration)
dataOutDemoLabIntegration_1 = na.omit(dataALLDemoLabIntegration_1)
nrow(dataOutDemoLabIntegration_1)
head(dataOutDemoLabIntegration_1)

dataOutDemoLabIntegration_1 = 
  na.omit(dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$gender==3,])
breaks = c(0,45,65,150)
label = c(1,2,3)
dataOutDemoLabIntegration_1$ageGroup= cut(dataOutDemoLabIntegration_1$age,breaks = breaks,labels = label)
dataOutDemoLabIntegration_1 = na.omit(dataOutDemoLabIntegration_1)
nrow(dataOutDemoLabIntegration_1 )



#dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month>6,]
#dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<60,]

summary(dataOutDemoLabIntegration_1 )
breaks = c(0,36,60,200)
label = c(1,2,3)
dataOutDemoLabIntegration_1$monthGroup= cut(dataOutDemoLabIntegration_1$month,breaks = breaks,labels = label)
table(dataOutDemoLabIntegration_1$tating)
table(dataOutDemoLabIntegration_1$monthGroup,dataOutDemoLabIntegration_1$tating)
nrow(dataOutDemoLabIntegration_1)
#dataOutDemoLabIntegration_1= dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<36,]


####################disease modify
dataOutDemoLabIntegration_1 = 
  na.omit(dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$dis_group==5,])
nrow(dataOutDemoLabIntegration_1 )
dataOutDemoLabIntegration_1$dis_group.y = 0
dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$dis_group%in%c(1:4),]$dis_group.y = 1
dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$dis_group%in%c(6),]$dis_group.y = 2
table(dataOutDemoLabIntegration_1$dis_group)
table(dataOutDemoLabIntegration_1$dis_group.y)

#table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]

max(table(dataOutDemoLabIntegration_1$MASTER_INDEX))

#dataAllDemoLabIntegration = 
#  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]


dataOutDemoLabIntegration_1C = dataOutDemoLabIntegration_1
nrow(dataOutDemoLabIntegration_1C)
table(dataOutDemoLabIntegration_1C$tating)





dataOutDemoLabIntegrationPsm = dataOutDemoLabIntegration_1

colnames(dataOutDemoLabIntegrationPsm)[c(3,15:21)]
c(3,11:17)
rownames = c("Lp(a) at FT (Mean,ci95%)","LDL-C (Mean,ci95%)",
             "HDL-C (Mean,ci95%)","APO-A (Mean,ci95%)","APO-B (Mean,ci95%)",
             "TC (Mean,ci95%)","TG (Mean,ci95%)","CRP (Mean,ci95%)")
num = table(dataOutDemoLabIntegrationPsm$tating)

colnames = c(paste("Control","(","n=",num[1],")",sep = ""),
             paste("Treat","(","n=",num[2],")",sep = ""),
             "StaVal","PVal")

colnames(dataOutDemoLabIntegrationPsm)[c(3,15:21)]
tab_1 = data.frame()
for(i in c(3,15:21)){
  tmp = dataOutDemoLabIntegrationPsm[c(9,i)]
  tmp = as.data.frame(descri2(tmp)) # mean1+ci95%,mean2+ci95%,staticValue,Pvalue
  tab_1 = rbind(tab_1,tmp)
}

colnames(tab_1) = colnames
rownames(tab_1) = rownames
tab_1
colnames(dataOutDemoLabIntegrationPsm)[c(7,29,31,27,10,30)]
tab_2 = data.frame()
#gender, age, Disease in FE, Disease in FT,lp(a) outcome,month
for (i in c(7,29,31,27,10,30)){
  M = dataOutDemoLabIntegrationPsm[c(i,9)]
  tmp = decriCount(M)
  tab_2 = rbind(tab_2,tmp)
}
colnames(tab_2)=colnames(tab_1)
tab_2

rownames(tab_2) = c("Male","Female","<45","46-65",">65",
                    "Non-CVD-related disease(FE)","CVD-related disease(FE)",
                    "[0.5 - 3) years","[3 - 5) years","≥ 5 years",
                    "Non-CVD disease event (FU)","CVD disease event (FU)",
                    "Lp(a) decrease (FU)","Lp(a) increase (FU)"
                    # "LDL-C decrease (FU)","LDL-C increase (FU)",
                    
                    #"Modify Non-CVD disease event (FU)","Modify CVD disease event (FU)",
                    # "diabetes","dyslipidemia","hypertension","arteriosclerosis","Others"
                    #"Lp(a) decrease_1 (FU)","Lp(a) increase_1 (FU)",
                    #"LDL-C decrease_1 (FU)","LDL-C increase_1 (FU)"
)
tab_2


rbind(tab_1,tab_2)


write.csv(file ="ALLNoPSMPatientsDemo.csv" , rbind(tab_1,tab_2))



######################mainDataSet
#=================================outPatients

nrow(dataOutDemoLabIntegration)
nrow(dataInDemoLabIntegration)
head(dataOutDemoLabIntegration)
head(dataInDemoLabIntegration)
7641+27944

nrow(dataOutDemoLabIntegrationPsm)
nrow(dataInDemoLabIntegrationPsm)
4888+ 16488


intersect(dataOutDemoLabIntegration$MASTER_INDEX,dataInDemoLabIntegration$MASTER_INDEX)

dataOutDemoLabIntegrationRaw = dataOutDemoLabIntegration

dataOutDemoLabIntegration = dataOutDemoLabIntegration[!dataOutDemoLabIntegration$MASTER_INDEX%in%intersect(dataOutDemoLabIntegration$MASTER_INDEX,dataInDemoLabIntegration$MASTER_INDEX),]
nrow(dataOutDemoLabIntegration)


dataOutDemoLabIntegration# outpatients dataset 

###################outpatients

dataOutDemoLabIntegration_1 = na.omit(dataOutDemoLabIntegration[dataOutDemoLabIntegration$LDL_C>1.8,])
nrow(dataOutDemoLabIntegration_1 )
table(dataOutDemoLabIntegration_1$tating)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$gender==3,]
nrow(dataOutDemoLabIntegration_1 )

dataOutDemoLabIntegration_1 = 
  na.omit(dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$gender==3,])
breaks = c(0,45,65,150)
label = c(1,2,3)
dataOutDemoLabIntegration_1$ageGroup= cut(dataOutDemoLabIntegration_1$age,breaks = breaks,labels = label)
dataOutDemoLabIntegration_1 = na.omit(dataOutDemoLabIntegration_1)
nrow(dataOutDemoLabIntegration_1 )

dataOutDemoLabIntegration_1 = merge(dataOutDemoLabIntegration_1,followUP,by.x = "MASTER_INDEX",by.y = "MASTER_INDEX")
nrow(dataOutDemoLabIntegration_1 )

#dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month>6,]
#dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<60,]

summary(dataOutDemoLabIntegration_1 )
breaks = c(0,36,60,200)
label = c(1,2,3)
dataOutDemoLabIntegration_1$monthGroup= cut(dataOutDemoLabIntegration_1$month,breaks = breaks,labels = label)
table(dataOutDemoLabIntegration_1$tating)
table(dataOutDemoLabIntegration_1$monthGroup,dataOutDemoLabIntegration_1$tating)
nrow(dataOutDemoLabIntegration_1)
#dataOutDemoLabIntegration_1= dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<36,]


####################disease modify
dataOutDemoLabIntegration_1 = 
  na.omit(dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$dis_group==5,])
nrow(dataOutDemoLabIntegration_1 )
dataOutDemoLabIntegration_1$dis_group.y = 0
dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$dis_group%in%c(1:4),]$dis_group.y = 1
dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$dis_group%in%c(6),]$dis_group.y = 2
table(dataOutDemoLabIntegration_1$dis_group)
table(dataOutDemoLabIntegration_1$dis_group.y)

#table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]

max(table(dataOutDemoLabIntegration_1$MASTER_INDEX))

#dataAllDemoLabIntegration = 
#  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]


dataOutDemoLabIntegration_1C = dataOutDemoLabIntegration_1
nrow(dataOutDemoLabIntegration_1C)
table(dataOutDemoLabIntegration_1C$tating)



dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==1,]
dataOutDemoLabIntegrationPsm_1 = psmOut(dataOutDemoLabIntegration_1)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==2,]
dataOutDemoLabIntegrationPsm_2 = psmOut(dataOutDemoLabIntegration_1)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==3,]
dataOutDemoLabIntegrationPsm_3 = psmOut(dataOutDemoLabIntegration_1)

dataOutDemoLabIntegrationPsm= rbind(
  dataOutDemoLabIntegrationPsm_1,
  dataOutDemoLabIntegrationPsm_2,
  dataOutDemoLabIntegrationPsm_3
)



colnames(dataOutDemoLabIntegrationPsm)[c(3,11:17)]
c(3,11:17)
rownames = c("Lp(a) at FT (Mean,ci95%)","LDL-C (Mean,ci95%)",
             "HDL-C (Mean,ci95%)","APO-A (Mean,ci95%)","APO-B (Mean,ci95%)",
             "TC (Mean,ci95%)","TG (Mean,ci95%)","CRP (Mean,ci95%)")
num = table(dataOutDemoLabIntegrationPsm$tating)

colnames = c(paste("Control","(","n=",num[1],")",sep = ""),
             paste("Treat","(","n=",num[2],")",sep = ""),
             "StaVal","PVal")

colnames(dataOutDemoLabIntegrationPsm)[c(3,11:17)]
tab_1 = data.frame()
for(i in c(3,11:17)){
  tmp = dataOutDemoLabIntegrationPsm[c(9,i)]
  tmp = as.data.frame(descri2(tmp)) # mean1+ci95%,mean2+ci95%,staticValue,Pvalue
  tab_1 = rbind(tab_1,tmp)
}

colnames(tab_1) = colnames
rownames(tab_1) = rownames
tab_1
colnames(dataOutDemoLabIntegrationPsm)[c(7,28,31,27,10)]
tab_2 = data.frame()
#gender, age, Disease in FE, Disease in FT,lp(a) outcome,month
for (i in c(7,28,31,27,10)){
  M = dataOutDemoLabIntegrationPsm[c(i,9)]
  tmp = decriCount(M)
  tab_2 = rbind(tab_2,tmp)
}
colnames(tab_2)=colnames(tab_1)
tab_2
write.csv(file ="outPatientsDemo.csv" , rbind(tab_1,tab_2))

#function
####################################Inpatients


nrow((dataInDemoLabIntegration))
head(dataInDemoLabIntegration)
dataInDemoLabIntegration_1 = na.omit(dataInDemoLabIntegration[dataInDemoLabIntegration$LDL_C>1.8,])
nrow(dataInDemoLabIntegration_1 )
table(table(dataInDemoLabIntegration_1$tating))
dataInDemoLabIntegration_1 = 
  na.omit(dataInDemoLabIntegration_1[!dataInDemoLabIntegration_1$gender==3,])
nrow(dataInDemoLabIntegration_1 )
breaks = c(0,45,60,65,70,80,1000)
label = c(1,2,3,4,5,6)
dataInDemoLabIntegration_1$ageGroup= cut(dataInDemoLabIntegration_1$age,breaks = breaks,labels = label)
nrow(dataInDemoLabIntegration_1 )
table(dataInDemoLabIntegration_1$ageGroup)
sum(table(dataInDemoLabIntegration_1$ageGroup))
dataInDemoLabIntegration_1 = na.omit(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1 = dataInDemoLabIntegration_1[order(dataInDemoLabIntegration_1$tating,decreasing = T),]
nrow(dataInDemoLabIntegration_1 )

table(followUP$MASTER_INDEX)
followUP[followUP$MASTER_INDEX == 2644,]


dataInDemoLabIntegration_1 = merge(dataInDemoLabIntegration_1,followUP,by.x = "MASTER_INDEX",by.y = "MASTER_INDEX",all.x = T)
nrow(dataInDemoLabIntegration_1)
# dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month>6,]
# dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<60,]

table(dataInDemoLabIntegration$MASTER_INDEX)[table(dataInDemoLabIntegration$MASTER_INDEX)>=2]
#dataAllDemoLabIntegration = 
#  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]


breaks = c(0,36,60,200)
label = c(1,2,3)
dataInDemoLabIntegration_1$monthGroup= cut(dataInDemoLabIntegration_1$month,breaks = breaks,labels = label)
table(dataInDemoLabIntegration_1$tating)
table(dataInDemoLabIntegration_1$monthGroup,dataInDemoLabIntegration_1$tating)




####################disease modify
dataInDemoLabIntegration_1 = 
  na.omit(dataInDemoLabIntegration_1[!dataInDemoLabIntegration_1$dis_group==5,])
nrow(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1$dis_group.y=1
dataInDemoLabIntegration_1[dataInDemoLabIntegration_1$dis_group%in%c(1:4),]$dis_group.y = 1
dataInDemoLabIntegration_1[dataInDemoLabIntegration_1$dis_group%in%c(6),]$dis_group.y = 2
table(dataInDemoLabIntegration_1$dis_group)

dataInDemoLabIntegration_1C = dataInDemoLabIntegration_1
nrow(dataInDemoLabIntegration_1C )

table(dataInDemoLabIntegration_1C$ageGroup)
table(dataInDemoLabIntegration_1C$tating)

#write.csv(dataInDemoLabIntegration_1,file="dataInDemoLabIntegration_1.csv")
#table(dataInDemoLabIntegration_1C$tating)

dataInDemoLabIntegration_1C = (dataInDemoLabIntegration_1C[!dataInDemoLabIntegration_1C$MASTER_INDEX%in%intersect(dataOutDemoLabIntegrationPsm$master_id,dataInDemoLabIntegration_1C$master_id),])
nrow(dataInDemoLabIntegration_1C)
head(dataInDemoLabIntegration_1C)


dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==1,]
dataInDemoLabIntegrationPsm_1 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==2,]
dataInDemoLabIntegrationPsm_2 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==3,]
dataInDemoLabIntegrationPsm_3 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==4,]
dataInDemoLabIntegrationPsm_4 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==5,]
dataInDemoLabIntegrationPsm_5 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==6,]
dataInDemoLabIntegrationPsm_6 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegrationPsm = rbind(dataInDemoLabIntegrationPsm_1,
                                    dataInDemoLabIntegrationPsm_2,
                                    dataInDemoLabIntegrationPsm_3,
                                    dataInDemoLabIntegrationPsm_4,
                                    dataInDemoLabIntegrationPsm_5,
                                    dataInDemoLabIntegrationPsm_6
)
table(dataInDemoLabIntegrationPsm$tating)
table(dataOutDemoLabIntegrationPsm$tating)
breaks = c(0,45,65,150)
label = c(1,2,3)
dataInDemoLabIntegrationPsm$ageGroup= cut(dataInDemoLabIntegrationPsm$age,
                                          breaks = breaks,labels = label)
nrow(dataInDemoLabIntegrationPsm)
nrow(dataOutDemoLabIntegrationPsm)



dataAllDemoLabIntegrationPsm = rbind(dataOutDemoLabIntegrationPsm,dataInDemoLabIntegrationPsm)
nrow(dataAllDemoLabIntegrationPsm)
breaks = c(0,45,65,150)
label = c(1,2,3)
dataAllDemoLabIntegrationPsm$ageGroup= cut(dataAllDemoLabIntegrationPsm$age,
                                           breaks = breaks,labels = label)

table(dataOutDemoLabIntegrationPsm$master_id)[table(dataOutDemoLabIntegrationPsm$master_id)>=2]

table(dataInDemoLabIntegrationPsm$master_id)[table(dataInDemoLabIntegrationPsm$master_id)>=2]

intersect(dataInDemoLabIntegrationPsm$master_id,dataOutDemoLabIntegrationPsm$master_id)


table(dataAllDemoLabIntegrationPsm$master_id)[table(dataAllDemoLabIntegrationPsm$master_id)>=2]

#dataAllDemoLabIntegrationPsm = dataAllDemoLabIntegrationPsm[!dataAllDemoLabIntegrationPsm$MASTER_INDEX%in%names(table(CVDmodify$master_id)[table(CVDmodify$master_id)>=2]),]
nrow(dataAllDemoLabIntegrationPsm)
#CVDmodify = read.table("CVDmodified.txt",header = T,sep = "\t")


#dataAllDemoLabIntegrationPsm = merge(dataAllDemoLabIntegrationPsm,CVDmodify,
#                                     by.x = "MASTER_INDEX",by.y = "master_id",all.y  = F)
nrow(dataAllDemoLabIntegrationPsm)


dataAllDemoLabIntegration = rbind(dataOutDemoLabIntegration,dataInDemoLabIntegration)
nrow(dataAllDemoLabIntegration)

table(rbind(dataOutDemoLabIntegration_1C,dataInDemoLabIntegration_1C)$tating)

head(dataAllDemoLabIntegration)
colnames(dataAllDemoLabIntegration)
dataAllDemoLabIntegration = dataAllDemoLabIntegration[,c(1,24)]
head(dataAllDemoLabIntegration)
dataAllDemoLabIntegration = 
  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]
head(dataAllDemoLabIntegration)
head(dataAllDemoLabIntegrationPsm)

dataAllDemoLabIntegrationPsm_1=merge(dataAllDemoLabIntegrationPsm,dataAllDemoLabIntegration,
                                     by.x = "MASTER_INDEX",by.y = "MASTER_INDEX",all.x = F)
nrow(dataAllDemoLabIntegrationPsm_1)

dataAllDemoLabIntegrationPsm = dataAllDemoLabIntegrationPsm_1



rownames = c("Lp(a) at FT (Mean,CI95%)","LDL-C at FT (Mean,CI95%)","LDL-C (Mean,CI95%)",
             "HDL-C (Mean,CI95%)","APO-A (Mean,CI95%)","APO-B (Mean,CI95%)",
             "TC (Mean,CI95%)","TG (Mean,CI95%)","CRP (Mean,CI95%)")
num = table(dataAllDemoLabIntegrationPsm$tating)

colnames = c(paste("Control","(","n=",num[1],")",sep = ""),
             paste("Treat","(","n=",num[2],")",sep = ""),
             "StaVal","PVal")
colnames(dataAllDemoLabIntegrationPsm)
tab_1 = data.frame()
for(i in c(3,12,15:21)){
  tmp = dataAllDemoLabIntegrationPsm[c(9,i)]
  tmp = as.data.frame(descri2(tmp)) # mean1+ci95%,mean2+ci95%,staticValue,Pvalue
  tab_1 = rbind(tab_1,tmp)
}

colnames(tab_1) = colnames
rownames(tab_1) = rownames
tab_1
dataAllDemoLabIntegrationPsm$CRP
tab_2 = data.frame()
colnames(dataAllDemoLabIntegrationPsm)[ c(7,28,31,30,10,14,27,25)]
#gender,age,dis_history,follow_time,LPA_outcome,LDL_outcome,dis_follow,
for (i in c(7,28,31,30,10,14,27,25)){
  M = dataAllDemoLabIntegrationPsm[c(i,9)]
  tmp = decriCount(M)
  tab_2 = rbind(tab_2,tmp)
}
tab_2
colnames(tab_2)=colnames(tab_1)
nrow(tab_2)
rownames(tab_2) = c("Male","Female","<45","46-65",">65",
                    "Non-CVD-related disease(FE)","CVD-related disease(FE)",
                    "[0.5 - 3) years","[3 - 5) years","≥ 5 years",
                    "Lp(a) decrease (FU)","Lp(a) increase (FU)",
                    "LDL-C decrease (FU)","LDL-C increase (FU)",
                    "Non-CVD disease event (FU)","CVD disease event (FU)",
                    #"Modify Non-CVD disease event (FU)","Modify CVD disease event (FU)",
                    "diabetes","dyslipidemia","hypertension","arteriosclerosis","Others"
                    #"Lp(a) decrease_1 (FU)","Lp(a) increase_1 (FU)",
                    #"LDL-C decrease_1 (FU)","LDL-C increase_1 (FU)"
)
tab_2


rbind(tab_1,tab_2)


write.csv(file = "mainDataset.csv",rbind(tab_1,tab_2))

write.csv(file = "mainDatasetRaw.csv",dataAllDemoLabIntegrationPsm)

head(dataAllDemoLabIntegrationPsm)

dat = read.csv("mainDatasetRaw.csv")

head(dat)

dataAllDemoLabIntegrationPsm = dat[,-1]
nrow(dataAllDemoLabIntegrationPsm)

summary(dataAllDemoLabIntegrationPsm$month)
sort(dataAllDemoLabIntegrationPsm$month)
###############################

resOR=c()
resORCi=c()

colnames(dataAllDemoLabIntegrationPsm)

nrow( dataAllDemoLabIntegrationPsm)

# glmIn = glm(formula = outcome.x ~   tating, 
#             data = dataAllDemoLabIntegrationPsm, family = binomial)

glmIn = coxph(Surv(month,outcome.x) ~  tating,
            data = dataAllDemoLabIntegrationPsm)

?coxph

print(summary(glmIn))

tmp = summary(glmIn)

exp(tmp$coefficients)

resOR[1] = formatC(tmp$coefficients[1,2],format = "f",digits = 2)
resORCi[1] = paste("(",
                   formatC(exp(tmp$coefficients[1,1]-1.96*tmp$coefficients[1,3]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[1,1]+1.96*tmp$coefficients[1,3]),
                           format = "f",digits = 2),")",sep = ""
)

# glmIn = glm(formula = disFollow_group_new  ~   tating , 
#             data = dataAllDemoLabIntegrationPsm, family = binomial)

glmIn = coxph(Surv(month,disFollow_group)   ~   tating , 
            data = dataAllDemoLabIntegrationPsm)
print(summary(glmIn))
tmp = summary(glmIn)
resOR[2] = formatC(exp(tmp$coefficients[1,1]),format = "f",digits = 2)
resORCi[2] = paste("(",
                   formatC(exp(tmp$coefficients[1,1]-1.96*tmp$coefficients[1,2]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[1,1]+1.96*tmp$coefficients[1,2]),
                           format = "f",digits = 2),")",sep = ""
)
resOR
resORCi

# glmIn = glm(formula =  outcome.x ~ as.numeric(ageGroup) +  tating +CRP + as.numeric(monthGroup)+ LDL_C + dis_group + HDL_C + 
#               APO_A +  APO_B+as.numeric(monthGroup), 
#             data = dataAllDemoLabIntegrationPsm, family = binomial)

glmIn = coxph(Surv(month,outcome.x) ~  gender + as.numeric(ageGroup) +  tating +CRP + LDL_C + dis_group + HDL_C + 
                APO_A +  APO_B,
              data = dataAllDemoLabIntegrationPsm)

print(summary(glmIn))

tmp = summary(glmIn)
resOR[3] = formatC(exp(tmp$coefficients[2,1]),format = "f",digits = 2)
resORCi[3] = paste("(",
                   formatC(exp(tmp$coefficients[2,1]-1.96*tmp$coefficients[2,3]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[2,1]+1.96*tmp$coefficients[2,3]),
                           format = "f",digits = 2),")",sep = ""
)
resOR
resORCi
tmp = as.data.frame(tmp$coefficients)


tmp$OR=formatC(tmp$`exp(coef)`,format = "f",digits = 2)

tmp$`95%CI` = paste("(",formatC((tmp$`exp(coef)`-1.96*tmp$`se(coef)`),format = "f",digits = 2),","
                    ,formatC((tmp$`exp(coef)`+1.96*tmp$`se(coef)`),format = "f",digits = 2),")",sep="")
tmp
tmp[,5] = formatC(tmp[,5],format = "f",digits = 2)
tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)

rownames(tmp) = c("Sex","Age","statin use","CRP","LDL-C",
                  "CVD-related disease(FE)","HDL-C","APO-A","APO-B")
tmp = tmp[,c(6,7,4,5)]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")

write.csv(tmp,file = "statinMultiOutcome.csv")


# glmIn = glm(formula = disFollow_group_new ~ as.numeric(ageGroup) +  tating + CRP + as.numeric(monthGroup)
#             + LDL_C + dis_group.x + HDL_C + APO_A +  APO_B, 
#             data = dataAllDemoLabIntegrationPsm, family = binomial)

glmIn = coxph(Surv(month,disFollow_group)  ~ as.numeric(ageGroup) +  tating + CRP 
            + LDL_C + dis_group+ HDL_C + APO_A +  APO_B, 
            data = dataAllDemoLabIntegrationPsm)



print(summary(glmIn))
tmp = summary(glmIn)
resOR[4] = formatC(exp(tmp$coefficients[3,1]),format = "f",digits = 2)
resORCi[4] = paste("(",
                   formatC(exp(tmp$coefficients[3,1]-1.96*tmp$coefficients[3,2]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[3,1]+1.96*tmp$coefficients[3,2]),
                           format = "f",digits = 2),")",sep = ""
)

# tmp = as.data.frame(tmp$coefficients)
# 
# 
# tmp$OR=formatC(exp(tmp$Estimate),format = "f",digits = 2)
# 
# tmp$`95%ci` = paste("(",formatC(exp(tmp$Estimate-1.96*tmp$`Std. Error`),format = "f",digits = 2),","
#                     ,formatC(exp(tmp$Estimate+1.96*tmp$`Std. Error`),format = "f",digits = 2),")",sep="")
# 
# tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)
# tmp[,3] = formatC(tmp[,3],format = "f",digits = 2)
# 
# rownames(tmp) = c("Intercept","Age","statin use","CRP","Follow-up(months)","LDL-C",
#                   "CVD-related disease(FE)","HDL-C","APO-A","APO-B")
# tmp = tmp[,c(5,6,3,4)]
# 
# colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")
# 
# write.csv(tmp,file = "followMultiOutcome.csv")
# 
# write.csv(tmp,file = "lpaCVD.csv")

#nrow(dataAllDemoLabIntegrationPsm)

mitifactors =  data.frame(resOR,resORCi)

rownames(mitifactors) = c("Single statin","Single Dis","Muti statin","Muti Dis")

write.csv(mitifactors,file = "mitifactorsMainData.csv")


breaks = c(0,2.300, 3.540,13.230)
label = c(3,2,1)
dataAllDemoLabIntegrationPsm$LDL_CGroup= cut(dataAllDemoLabIntegrationPsm$LDL_C,breaks = breaks,labels = label)
dataAllDemoLabIntegrationPsm$LDL_CGroup



tabdrug = table(dataAllDemoLabIntegrationPsm$LDL_CGroup,dataAllDemoLabIntegrationPsm$disFollow_group)

tabdrug[,2]/tabdrug[,1]


fivenum(dataAllDemoLabIntegrationPsm$lastValue)

breaks = c(0, 179.00, 322.00,10000)
label = c(1,2,3)
dataAllDemoLabIntegrationPsm$LPAFlGroup= cut(as.numeric(dataAllDemoLabIntegrationPsm$median.x)+1,breaks = breaks,labels = label)
dataAllDemoLabIntegrationPsm$lastValue
dataAllDemoLabIntegrationPsm$LPAFlGroup

cbind(dataAllDemoLabIntegrationPsm$lastValue,dataAllDemoLabIntegrationPsm$LPAFlGroup)

breaks = c(0, 65, 150)
label = c(1,2)
dataAllDemoLabIntegrationPsm$AgeGroup2= cut(dataAllDemoLabIntegrationPsm$age,breaks = breaks,labels = label)



colnames( dataAllDemoLabIntegrationPsm)[c(9,10,27,15:21,7,14,29,35,31,33,34)]
# [1] "tating"          "outcome.x"       "disFollow_group" "LDL_C"           "HDL_C"          
# [6] "APO_A"           "APO_B"           "TC"              "TG"              "CRP"            
# [11] "gender"          "outcome.y"       "dis_group.y"     "disFollow_group" "monthGroup"     
# [16] "LPAFlGroup"      "AgeGroup2" 

subcovarience = dataAllDemoLabIntegrationPsm[c(9,10,27,15:21,7,14,29,35,31,33,34)]
subcovarienceMain = subcovarience
head(subcovarience )
subcovarience$dis_group.y
subcovarience [,12] = subcovarience [,12]+1

subcovarience [,14]
colnames(subcovarience)
colnames(subcovarience)

metaFor  = subgroup(subcovarience,12,"subgroupldlcDis",
                    c("LDL-C decrease","LDL-C increase"),"outcome.x")

tmp = subgroup(subcovarience,12,"subgroupldlcDis",
               c("LDL-C decrease","LDL-C increase"),"disFollow_group")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])

metaFor

tmp = subgroup(subcovarience,15,"subgrouplpaDis",
               c("Non-comorbidity history","Comorbidity history"),"outcome.x")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])



tmp = subgroup(subcovarience,15,"subgroupFolDisDis",
               c("Non-CVD disease history","CVD disease history"),"disFollow_group")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])
metaFor


tmp = subgroup(subcovarience,16,"subgrouplpaLDL_C",
               c("≥3.54mmol/L", "(2.30,3.54)mmol/L","≤2.30mmol/L"),"outcome.x")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])

tmp = subgroup(subcovarience,16,"subgroupFolDisLDL_C",
               c("≥3.54mmol/L", "(2.30,3.54)mmol/L","≤2.30mmol/L"),"disFollow_group")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])


tmp = subgroup(subcovarience,17,"subgrouplpaLLPAFL",
               c("≤179.00mg/L", "(179.00,322.00)mg/L","≥322.00mg/L"),"outcome.x")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])
tmp = subgroup(subcovarience,17,"subgroupFolDisLPAFL",
               c("≤179.00mg/L", "(179.00,322.00)mg/L","≥322.00mg/L"),"disFollow_group")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])


# tmp =subgroup(subcovarience,13,"subgrouplpaLMonth",
#               c("(0.5,3] years", "(3,5) year","≥5 year"),"outcome.x")
# metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
# metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])
# tmp =subgroup(subcovarience,13,"subgroupFolDisLMonth",
#               c("(0.5,3] years", "(3,5) year","≥5 year"),"disFollow_group")
# metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
# metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])

tmp = subgroup(subcovarience,11,"subgrouplpaSex",
               c("Male", "Female"),"outcome.x")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])
tmp = subgroup(subcovarience,11,"subgroupFolDisSex",
               c("Male", "Female"),"disFollow_group")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])


tmp = subgroup(subcovarience,14,"subgrouplpaSex",
               c("<65", "≥65"),"outcome.x")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])
tmp = subgroup(subcovarience,14,"subgroupFolDisSex",
               c("<65", "≥65"),"disFollow_group_new")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])

metaFor

#lpasubgroup = metaFor[[1]][c(1,2,3,4,7,8,11:13,17:19,23:25,29:31,35,36),]

#lpasubgroup = metaFor[[1]][c(1,2,3,4,7,8,11:13,17:19,23:25,29:30,33,34),]

nrow( metaFor[[1]])
lpasubgroup = metaFor[[1]][c(1,2,27,28,23,24,7,8,11:13,17:19,3,4),]


metaFor[[1]][,1]

Item = c("\n","\n",
         "Age","\n",
         "Sex","\n",
         "Comorbidity","\n",
         "LDL-C level","\n","\n",
         "Lp(a) level (FE)","\n","\n",
         "LDL-C outcome","\n"
)
length(Item)


# Item = c("\n","\n","LDL-C outcome","\n",
#         "CVD-related disease history","\n",
#          "Age","\n","\n",
#          "LDL-C level","\n","\n",
#          "Lp(a) level","\n","\n",
#          "Follow-up time","\n","\n",
#          "Sex","\n","\n")

lpasubgroup = cbind(Item,lpasubgroup)

lpasubgroup[2,3:6] = c("Lp(a) increase","Lp(a) decrease","Lp(a) increase","Lp(a) decrease")

metaFor[[1]][,1]
# #lpasubgroup = metaFor[[1]][c(1,2,33,34,29,30,7,8,23:25,11:13,17:19,3,4),]
# 
# CVDsubgroup = metaFor[[1]][c(1,2,35,36,31,32,9,10,26:28,14:16,20:22,5,6),]
# 
# CVDsubgroup  = cbind(Item,CVDsubgroup)
# 
# CVDsubgroup[2,3:6] = c("CVD event","Non-CVD event","CVD event","Non-CVD event")
# 

tiff(filename =  "subgroup-lpaInRish.tif",
     width = 9480, height = 3480, units = "px", pointsize = 12,
     compression ="lzw",
     bg = "white", res = 300)

forestplot(lpasubgroup, 
           txt_gp=fpTxtGp(
             label=gpar(cex=2),
             ticks=gpar(cex=2),
             xlab=gpar(cex = 1.7),
             title=gpar(cex = 1.8)),
           rbind( rep(NA,3),rep(NA,3),
                  metaFor[[2]][c(27,28,23,24,7,8,11:13,17:19,3,4)-2,]),
           col=clrs,
           colgap=unit(8,"mm"),
           lwd.ci=2, boxsize=0.5,
           zero=1,
           cex=2, lineheight = "auto",
           ci.vertices=TRUE, ci.vertices.height = 0.4
           #xlab="    <---Decrease--- lp(a) ---Increase--->"
)
dev.off()

tiff(filename =  "subgroup-CVDInRish.tif",
     width = 9480, height = 3480, units = "px", pointsize = 12,
     compression ="lzw",
     bg = "white", res = 300)

forestplot(CVDsubgroup, 
           txt_gp=fpTxtGp(
             label=gpar(cex=2),
             ticks=gpar(cex=2),
             xlab=gpar(cex = 1.7),
             title=gpar(cex = 1.8)),
           rbind( rep(NA,3),rep(NA,3),
                  metaFor[[2]][c(35,36,31,32,9,10,26:28,14:16,20:22,5,6)-2,]),
           col=clrs,
           colgap=unit(8,"mm"),
           lwd.ci=2, boxsize=0.5,
           zero=1,
           cex=2, lineheight = "auto",
           ci.vertices=TRUE, ci.vertices.height = 0.4
           #xlab="    <---Decrease--- lp(a) ---Increase--->"
)
dev.off()

dataAllDemoLabIntegrationPsm$monthGroup

dataAllDemoLabIntegrationPsm[dataAllDemoLabIntegrationPsm$monthGroup==1,]

table(dataAllDemoLabIntegrationPsm[dataAllDemoLabIntegrationPsm$monthGroup==1,]$tating,
      dataAllDemoLabIntegrationPsm[dataAllDemoLabIntegrationPsm$monthGroup==1,]$outcome.y)

table(dataAllDemoLabIntegrationPsm[dataAllDemoLabIntegrationPsm$monthGroup==2,]$tating,
      dataAllDemoLabIntegrationPsm[dataAllDemoLabIntegrationPsm$monthGroup==2,]$outcome.y)

table(dataAllDemoLabIntegrationPsm[dataAllDemoLabIntegrationPsm$monthGroup==3,]$tating,
      dataAllDemoLabIntegrationPsm[dataAllDemoLabIntegrationPsm$monthGroup==3,]$outcome.y)


breaks = c(0,2.300, 3.540,13.230)

colnames(dataAllDemoLabIntegrationPsm)

head(dataAllDemoLabIntegrationPsm)
#


#####################LPA assocaiton with CVD 

glmIn = coxph(Surv(month,disFollow_group)~ as.factor(LPAFlGroup)+age+lastValue.y +
              +dis_group+gender+HDL_C + APO_A + CRP +
                APO_B, data =dataAllDemoLabIntegrationPsm)

nrow(dataAllDemoLabIntegrationPsm)

tmp = summary(glmIn)

tmp

write()

cbind(dataAllDemoLabIntegrationPsm$firstValue.x,dataAllDemoLabIntegrationPsm$LPAFlGroup)


tmp = as.data.frame(tmp$coefficients)

tmp
tmp$OR=formatC(tmp$`exp(coef)`,format = "f",digits = 2)

tmp$`95%CI` = paste("(",formatC((tmp$`exp(coef)`-1.96*tmp$`se(coef)`),format = "f",digits = 2),","
                    ,formatC((tmp$`exp(coef)`+1.96*tmp$`se(coef)`),format = "f",digits = 2),")",sep="")
tmp
tmp[,5] = formatC(tmp[,5],format = "f",digits = 2)
tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)

rownames(tmp) = c("2","3","Age","LDL-C","CVD-related disease(FE)",
                  "Gender","HDL-C","APO-A","CRP","APO-B")
tmp = tmp[,c(6,7,4,5)]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")

write.csv(tmp,file = "LpaMultiOutcome.csv")





#############################################################################

#########drug analysis


sql = "
SELECT  [master_id]
,[reg_time]
,[UNITPRICE]
,[QUANTITY]
,[UNIT]
,[FEE]
,[I_ITEM_NAME]
FROM [DB_LPA_ZONG].[dbo].[TATING_ID_TT2_TYPE]
order by master_id,reg_time asc
"
drugInfo= tbl_df(
  sqlQuery(cn,
           sql) )
head(drugInfo)


drugInfo <- group_by(drugInfo, master_id)



drugInfoS  <-  drugInfo %>% 
  summarise(QUANTITY = sum(QUANTITY),time =(max(reg_time)-min(reg_time)), drugName = names(which.max(table(I_ITEM_NAME))),
            durg = length(table(I_ITEM_NAME)[table(I_ITEM_NAME)>0]),
            af = table(I_ITEM_NAME)[1],
            ff = table(I_ITEM_NAME)[2],
            lf = table(I_ITEM_NAME)[3],
            pf = table(I_ITEM_NAME)[4],
            rf = table(I_ITEM_NAME)[5],
            xf = table(I_ITEM_NAME)[6],
            index = max(table(I_ITEM_NAME))/sum(table(I_ITEM_NAME)))
nrow(drugInfoS)
head(drugInfoS)

drugInfoS$drugName


drugInfoTW <- group_by(drugInfo,EdrugName)




drugInfoTW_1  <-  drugInfoTW  %>% 
  summarise( StartTime = min(reg_time),EndTime = max(reg_time))

drugInfoTW_1

#write.csv("DrugTW.csv",x=drugInfoTW_1)


dataAllDemoLabIntegrationPsmdrug = merge(dataAllDemoLabIntegrationPsm,drugInfoS,
                                         by.x = "MASTER_INDEX",by.y = "master_id")

nrow(drugInfoS)
nrow(dataAllDemoLabIntegrationPsmdrug)
head(dataAllDemoLabIntegrationPsmdrug)


###lpa
dataAllDemoLabIntegrationPsmdrug$EdrugName = case_when(
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "阿托伐他汀" ~"Atorvastatin" ,
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "氟伐他汀" ~ "Fluvastatin",
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "洛伐他汀" ~ "Lovastatin",
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "普伐他汀" ~"Pravastatin",
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "瑞舒伐他汀" ~ "Rosuvastatin",
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "辛伐他汀" ~ "Simvastatin"
)

table(dataAllDemoLabIntegrationPsmdrug$EdrugName)
?pie
tiff(filename =  "drugDistribution.tif",
     width = 5480, height = 6480, units = "px", pointsize = 12,
     compression ="lzw",
     bg = "white", res = 300)
sales <- table(dataAllDemoLabIntegrationPsmdrug[dataAllDemoLabIntegrationPsmdrug$index==1,]$EdrugName)
names <- c("Atorvastatin","Fluvastatin","Lovastatin","Pravastatin", "Rosuvastatin", "Simvastatin")
per.sales <- paste(sales,"(",round(100 * sales / sum(sales)),"%",")",sep = "")
slice.col <- rainbow(6)
pie(sales,labels=per.sales,col= slice.col,cex=3)
legend("topright",names,cex=3, fill=slice.col)
dev.off()




# 
# drugInfoSM = drugInfoS[drugInfoS$index>0.8,]
# 
# table(drugInfoSM$drugName)
# 
# nrow(drugInfoS[drugInfoS$index>0.8,])
# 
# max
# 
# 
# table(drugInfoS$durg)
# 
# sum(table(drugInfoS$durg))
# 
# drugInfoAf = drugInfo[drugInfo$I_ITEM_NAME=="阿托伐他汀",]
# 
# barplot(table(abs(drugInfoAf$FEE)),xlab = "Count",ylab = "Fee")
# 
# drugInfoAf <- group_by(drugInfoAf, master_id)
# 
# head(drugInfoAf)
# 
# 
# 
# drugInfoAfS  <-  drugInfoAf %>% 
#   summarise(FEE = sum(FEE),durg = length(table(I_ITEM_NAME)[table(I_ITEM_NAME)>0]))
# 
# table(drugInfoAfS$durg)
# 
# tmp = drugInfo[drugInfo$master_id=="1449",]
# tmp = drugInfo[drugInfo$master_id=="1481",]
# 
# 
# names(which.max(table(tmp$I_ITEM_NAME)))
# 
# 
# 
# 
# table(drugInfoAf[drugInfoAf$master_id=="1449",]$I_ITEM_NAME)[table(drugInfoAf[drugInfoAf$master_id=="1449",]$I_ITEM_NAME)>0]
# 
# 
# head(drugInfoAfS)
# nrow(drugInfoAfS)
# barplot(table(abs(drugInfoAfS$FEE)),xlab = "Fee",ylab = "Count")
# 
# 
# table()
# 
# #阿托伐他汀   氟伐他汀   洛伐他汀   普伐他汀 瑞舒伐他汀   辛伐他汀 

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization
set.seed(123)
k.values <- 1:15

colnames(dataAllDemoLabIntegrationPsmdrug)
nrow(dataAllDemoLabIntegrationPsmdrug)
dataAllDemoLabIntegrationPsmdrug$time = as.numeric(dataAllDemoLabIntegrationPsmdrug$time)

wss <- function(k) {
  kmeans(na.omit(dataAllDemoLabIntegrationPsmdrug[,c(36,29)]), k, nstart = 10 )$tot.withinss
}
wss_values <- map_dbl(k.values, wss)

?tiff
tiff(filename = "kmeanGroupNum.tif",
     width = 1480, height = 1480, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res =300)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()
set.seed(123)
head( dataAllDemoLabIntegrationPsmdrug)
k4 <- kmeans( na.omit(dataAllDemoLabIntegrationPsmdrug[,c(36,29)]), centers = 4, nstart = 25)

k4$cluster = case_when(
  k4$cluster  ==  1 ~ 1,
  k4$cluster  ==  2 ~ 4,
  k4$cluster  ==  3 ~ 3,
  k4$cluster  ==  4 ~ 2
)

colnames(dataAllDemoLabIntegrationPsmdrug)[c(36,29)]=c("PCA-1(Dose)","PCA-2(Follow-up time)")

tiff(filename = "kmeanGroupDistribuiton.tif",
     width = 1480, height = 1480, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res =300)
fviz_cluster(k4, geom = "point",
             data =dataAllDemoLabIntegrationPsmdrug[,c(36,29)])+theme_classic()
dev.off()


group = k4$cluster
#varience
dataAllDemoLabIntegrationPsmdrug$group = group
colnames(dataAllDemoLabIntegrationPsmdrug)[c(36,29)]=c("Dose","Follow-up time")
dataAllDemoLabIntegrationPsmdrug$Dose
?par
tiff(filename = "kmeanGroupFeeTime.tif",
     width =3480, height = 1480, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res =300)
par(mfcol=c(1,2))
boxplot(Dose~group,data = dataAllDemoLabIntegrationPsmdrug,xlab = "Group",ylab = "Dose(Piece)")

boxplot(`Follow-up time`~group,data = dataAllDemoLabIntegrationPsmdrug,xlab = "Group",ylab = "Follow-up time(month)")

dataAllDemoLabIntegrationPsmDrug$`Follow-up time`

dev.off()
par(mfcol=c(1,1))

tmp = summary(lm(Dose~as.numeric(group), data =dataAllDemoLabIntegrationPsmdrug)) 
tmp
tmp$r.squared
tmp$coefficients[2,1]-1.96*tmp$coefficients[2,2]
tmp$coefficients[2,1]+1.96*tmp$coefficients[2,2]

head(dataAllDemoLabIntegrationPsmdrug)
tmp = summary(lm(`Follow-up time`~as.numeric(group), data =dataAllDemoLabIntegrationPsmdrug))   #  15.96
tmp
tmp$r.squared
tmp$coefficients[2,1]-1.96*tmp$coefficients[2,2]
tmp$coefficients[2,1]+1.96*tmp$coefficients[2,2]


#tmp = summary(lm(log(median.x+1)~group+age, data =dataAllDemoLabIntegrationPsmdrug))


tmp = summary(lm(log(median.x+1)~group+age+LDL_C + dis_group + HDL_C + APO_A +  APO_B+ 
                   as.numeric(monthGroup)+
                   as.numeric(gender), data =dataAllDemoLabIntegrationPsmdrug))  

tmp = summary(lm(log(median.x+1)~group, data =dataAllDemoLabIntegrationPsmdrug))  



tmp = as.data.frame(tmp$coefficients)



tmp$RS=formatC((tmp$Estimate),format = "f",digits = 2)


tmp$`95%CI` = paste(tmp$RS,"(",formatC((tmp$Estimate-1.96*tmp$`Std. Error`),format = "f",digits = 2),","
                    ,formatC((tmp$Estimate+1.96*tmp$`Std. Error`),format = "f",digits = 2),")",sep="")

tmp$low = formatC((tmp$Estimate-1.96*tmp$`Std. Error`),format = "f",digits = 2)

tmp$high = formatC((tmp$Estimate+1.96*tmp$`Std. Error`),format = "f",digits = 2)


tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)
tmp[,3] = formatC(tmp[,3],format = "f",digits = 2)

rownames(tmp) = c("Intercept","Group","Age","LDL-C","CVD-related disease","HDL-C",
                  "APO-A","APO-B","Follow-up time","Sex")



tmp = tmp[,c(5,6,3,4,7,8)]

durgIntension = tmp[2,]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")

write.csv(tmp,file = "groupLpa.csv")





tmp = summary(coxph(Surv(`Follow-up time`,disFollow_group)~group+age+LDL_C + dis_group + HDL_C + APO_A +  APO_B+ 
                         as.numeric(gender), data =dataAllDemoLabIntegrationPsmdrug))  


tmp = as.data.frame(tmp$coefficients)


tmp$HR=formatC((tmp$`exp(coef)`),format = "f",digits = 2)


tmp$`95%CI` = paste(tmp$HR,"(",formatC((tmp$`exp(coef)`-1.96*tmp$`se(coef)`),format = "f",digits = 2),","
                    ,formatC((tmp$`exp(coef)`+1.96*tmp$`se(coef)`),format = "f",digits = 2),")",sep="")


tmp[,5] = formatC(tmp[,5],format = "f",digits = 2)
tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)

rownames(tmp) = c("Group","Age","LDL-C","CVD-related disease","HDL-C",
                  "APO-A","APO-B","Sex")

durgIntension = tmp[2,]

tmp = tmp[,c(6,7,4,5)]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")


write.csv(tmp,file = "groupCVD.csv")

tmp = summary(lm(LDL_C~group + age + dis_group + HDL_C + APO_A +  APO_B+ as.numeric(monthGroup)+
                   as.numeric(gender), data =dataAllDemoLabIntegrationPsmdrug)) 

tmp = as.data.frame(tmp$coefficients)


tmp$RC=formatC((tmp$Estimate),format = "f",digits = 2)

tmp$`95%CI` = paste(tmp$RC,"(",formatC((tmp$Estimate-1.96*tmp$`Std. Error`),format = "f",digits = 2),","
                    ,formatC((tmp$Estimate+1.96*tmp$`Std. Error`),format = "f",digits = 2),")",sep="")

tmp$low = formatC((tmp$Estimate-1.96*tmp$`Std. Error`),format = "f",digits = 2)

tmp$high = formatC((tmp$Estimate+1.96*tmp$`Std. Error`),format = "f",digits = 2)



tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)
tmp[,3] = formatC(tmp[,3],format = "f",digits = 2)

rownames(tmp) = c("Intercept","Group","Age","CVD-related disease","HDL-C",
                  "APO-A","APO-B","Follow-up time","Sex")


tmp = tmp[,c(5,6,3,4,7,8)]

durgIntension = rbind(as.data.frame(durgIntension) ,as.character(tmp[2,]))



colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")

write.csv(tmp,file = "groupLDL-C.csv")


datFrame = durgIntension[,c(1,5,6,4)]
rownames(datFrame ) = c("Lp(a)","LDL-C")


colnames(datFrame) = c(colnames(HRQoL$Sweden),"Pvalue")
#coef       lower       upper

clrs <- fpColors(box="royalblue",line="darkblue", summary="royalblue")

tabletext <- cbind(c("Item",rownames(datFrame)),
                   c("RS(95%CI)",durgIntension[,2]),
                   c("Pvalue","<0.01","<0.01"))
tabletext 

datFrame$coef=as.numeric(datFrame$coef)
datFrame$lower=as.numeric(datFrame$lower)
datFrame$upper=as.numeric(datFrame$upper)
datFrame$Pvalue=as.numeric(datFrame$Pvalue)

tmp1 = rbind( rep(NA,3),as.data.frame(datFrame[,-4]))

tmp1$coef

tiff(filename = "LPALDL-CstatinDose.tiff",
     width = 5480, height = 3480, units = "px", pointsize = 12,
     compression ="lzw",
     bg = "white", res = 300)

forestplot(tabletext, 
            txt_gp=fpTxtGp(
              label=gpar(cex=2),
              ticks=gpar(cex=2),
             xlab=gpar(cex = 1.7),
              title=gpar(cex = 1.8)
           ),
           tmp1,
           col=clrs,
           #colgap=unit(4,"mm"),
          lwd.ci=0.1, 
          boxsize=0.1,
          # zero=1,
           #cex=2, lineheight = "auto",
           ci.vertices=TRUE, ci.vertices.height = 0.05
           # xlab="    <---Decrease--- CVD ---Increase--->"
)
dev.off()







nrow(dataAllDemoLabIntegrationPsmDrug)


groupCVD = 
  table(dataAllDemoLabIntegrationPsmdrug$group,dataAllDemoLabIntegrationPsmdrug$disFollow_group)
for(i in 1:4){
  groupCVD[i,] = groupCVD[i,]/sum(groupCVD[i,])
}

# table(dataAllDemoLabIntegrationPsmDrug$group)
# 
# nrow(dataAllDemoLabIntegrationPsmDrug)

# decriCount(dataAllDemoLabIntegrationPsmDrug[,c(23,28)])
# 
# write.csv(decriCount(dataAllDemoLabIntegrationPsmDrug[,c(23,28)]),file = "groupDisFol.csv")


groupCVD=
  as.data.frame(cbind(groupCVD[,1],groupCVD[,2]))

colnames(groupCVD)=c("No CVD event","CVD event")
groupCVD$group = c(1:4)
tiff(filename = "GroupCVD.tif",
     width = 1480, height = 1480, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res =300)
barplot(cbind(`No CVD event`, `CVD event`) ~ group, data = groupCVD[c(1:4),],
        col = c("green", "red"),
        legend = colnames(groupCVD)[c(1,2)], ylim = c(0, 1.5))
dev.off()


clusterF=data.frame()
clusterT = data.frame()
for(i in 1:4){
  tmp = dataAllDemoLabIntegrationPsmdrug[dataAllDemoLabIntegrationPsmdrug$group==i,]$Dose
  
  c1 = t(as.data.frame(formatC(c(mean(tmp),
                                 quantile(tmp,
                                          c(0.025, 0.975))),digits = 2,format = "f")))
  
  clusterF = rbind(clusterF,c1)
  tmp =   dataAllDemoLabIntegrationPsmdrug[dataAllDemoLabIntegrationPsmdrug$group==i,]$`Follow-up time`
  c1 = t(as.data.frame(formatC(c(mean(tmp),
                                 quantile(tmp,
                                          c(0.025, 0.975))),digits = 2,format = "f")))
  clusterT = rbind(clusterT,c1)
  
}
Dose = paste(clusterF[,1],"(",clusterF[,2],", ",clusterF[,3],")",sep = "")
`Follow-up time` = paste(clusterT[,1],"(",clusterT[,2],", ",clusterT[,3],")",sep = "")

groupTimeFee = data.frame(Dose,`Follow-up time`)
rownames(groupTimeFee) = c(paste("Group-1","(",nrow(dataAllDemoLabIntegrationPsmdrug[dataAllDemoLabIntegrationPsmdrug$group==1,]),")",sep = ""),
                           paste("Group-2","(",nrow(dataAllDemoLabIntegrationPsmdrug[dataAllDemoLabIntegrationPsmdrug$group==2,]),")",sep = ""),
                           paste("Group-3","(",nrow(dataAllDemoLabIntegrationPsmdrug[dataAllDemoLabIntegrationPsmdrug$group==3,]),")",sep = ""),
                           paste("Group-4","(",nrow(dataAllDemoLabIntegrationPsmdrug[dataAllDemoLabIntegrationPsmdrug$group==4,]),")",sep = "")
)

rownames(groupCVD)= rownames(groupTimeFee)

write.csv(groupTimeFee,file = "groupTimeFee.csv")

write.csv(groupCVD,file = "groupCVDe.csv")

###durg 

#dataAllDemoLabIntegrationPsmdrug08 = dataAllDemoLabIntegrationPsmdrug[dataAllDemoLabIntegrationPsmdrug$index>0.8,]

dataAllDemoLabIntegrationPsmdrug08 = dataAllDemoLabIntegrationPsmdrug[dataAllDemoLabIntegrationPsmdrug$index==1,]

nrow(dataAllDemoLabIntegrationPsmdrug08)


tmp = table(dataAllDemoLabIntegrationPsmdrug08$EdrugName,dataAllDemoLabIntegrationPsmdrug08$disFollow_group)

tabdrug[,2]/tabdrug[,1]

tabdrug = table(dataAllDemoLabIntegrationPsmdrug08$EdrugName,dataAllDemoLabIntegrationPsmdrug08$outcome.x)

tabdrug[,2]/tabdrug[,1]

tmp = cbind(tmp,tabdrug)

sum(tmp)

colnames(tmp) = c("CVD dec","CVD inc","Statin dec","Statin inc")

write.csv(tmp,"drugInt.csv")


dataAllDemoLabIntegrationPsmdrug$drugNmaeCodeDis = case_when(
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "阿托伐他汀" ~1 ,
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "氟伐他汀" ~ 2,
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "洛伐他汀" ~ 3,
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "普伐他汀" ~4,
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "瑞舒伐他汀" ~ 5,
  dataAllDemoLabIntegrationPsmdrug$drugName  ==  "辛伐他汀" ~ 6
)

tmp = summary(coxph(Surv(`Follow-up time`,outcome.x)~as.factor(drugNmaeCodeDis)+age+LDL_C+firstValue.x+dis_group+gender,
                data = dataAllDemoLabIntegrationPsmdrug))
tmp

tmp = summary(coxph(Surv(`Follow-up time`,disFollow_group)~as.factor(drugNmaeCodeDis)+age+LDL_C+firstValue.x+dis_group+gender,
                    data = dataAllDemoLabIntegrationPsmdrug))
tmp

tmp = as.data.frame(tmp$coefficients)


tmp$HR=formatC(tmp$`exp(coef)`,format = "f",digits = 2)

tmp$`95%CI` = paste("(",formatC((tmp$`exp(coef)`-1.96*tmp$`se(coef)`),format = "f",digits = 2),","
                    ,formatC((tmp$`exp(coef)`+1.96*tmp$`se(coef)`),format = "f",digits = 2),")",sep="")

tmp[,5] = formatC(tmp[,5],format = "f",digits = 2)
tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)

rownames(tmp) = c("Fluvastatin","Lovastatin","Pravastatin",
                  "Rosuvastatin","Sivastatin",
                  "Age","LDL-C","Lp(a) FE","CVD related disease history","Sex")


rownames(tmp) = c("Intercept",
                  "Age","LDL-C","Lp(a) Group",
                  "CVD related disease history","Sex","HDL-C",
                  "APO-A","APO-B")


tmp = tmp[,c(6,7,4,5)]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")

write.csv(tmp,file = "drugIntensityLpa.csv")

write.csv(tmp,file = "drugIntensityCVD.csv")

write.csv(tmp,file = "LPALevelCVD.csv")




tmp = summary(coxph(Surv(`Follow-up time`,outcome.x)~as.factor(drugNmaeCodeDis)+age+LDL_C+firstValue.x+dis_group+gender,
                    data = dataAllDemoLabIntegrationPsmdrug))
tmp = as.data.frame(tmp$coefficients)


tmp = summary(coxph(Surv(`Follow-up time`,disFollow_group)~as.factor(drugNmaeCodeDis)+age+LDL_C+firstValue.x+dis_group+gender,
                    data = dataAllDemoLabIntegrationPsmdrug))
tmp = as.data.frame(tmp$coefficients)


tmp$HR=formatC(tmp$`exp(coef)`,format = "f",digits = 2)

tmp$`95%CI` = paste("(",formatC((tmp$`exp(coef)`-1.96*tmp$`se(coef)`),format = "f",digits = 2),","
                    ,formatC((tmp$`exp(coef)`+1.96*tmp$`se(coef)`),format = "f",digits = 2),")",sep="")

tmp$`Pr(>|z|)`=formatC(tmp$`Pr(>|z|)`,format = "f",digits = 2)


tmp$low = (tmp$`exp(coef)`-1.96*tmp$`se(coef)`)
tmp$high = (tmp$`exp(coef)`+1.96*tmp$`se(coef)`)

drugLpa = table(dataAllDemoLabIntegrationPsmdrug08$drugName,dataAllDemoLabIntegrationPsmdrug08$outcome.x)

drugCVD = table(dataAllDemoLabIntegrationPsmdrug08$drugName,dataAllDemoLabIntegrationPsmdrug08$disFollow_group)

rownames(tmp) = c("Fluvastatin","Lovastatin","Pravastatin",
                  "Rosuvastatin","Sivastatin",
                  "Age","LDL-C","Lp(a) FE","CVD related disease history","Sex")
datFrame = tmp[c(1:5),c(6,8,9,5)]



colnames(datFrame) = c(colnames(HRQoL$Sweden),"Pvalue")
#coef       lower       upper

clrs <- fpColors(box="royalblue",line="darkblue", summary="royalblue")

tabletext <- cbind(c("Statin","\n","Atorvastatin (ref)",rownames(datFrame)),
                   c("Number","CVD event",drugCVD[,2]),
                   c("Number","Non-CVD event",drugCVD[,1]),
                   c("HR","\n", "1.00",datFrame[,"coef"]),
                   c("P","\n", "-",datFrame[,"Pvalue"]))
tabletext 

tiff(filename = "drugCVD.tiff",
     width = 5480, height = 3480, units = "px", pointsize = 12,
     compression ="lzw",
     bg = "white", res = 300)

forestplot(tabletext, 
           txt_gp=fpTxtGp(
             label=gpar(cex=2),
             ticks=gpar(cex=2),
             xlab=gpar(cex = 1.7),
             title=gpar(cex = 1.8)),
           rbind( rep(NA,3),rep(NA,3),rep(1,3),
                  datFrame[,-4]),
           col=clrs,
           colgap=unit(8,"mm"),
           lwd.ci=2, boxsize=0.5,
           zero=1,
           cex=2, lineheight = "auto",
           ci.vertices=TRUE, ci.vertices.height = 0.4,
           xlab="    <---Decrease--- CVD ---Increase--->"
)
dev.off()


tabletext <- cbind(c("Statin","\n","Atorvastatin (ref)",rownames(datFrame)),
                   c("Number","Lp(a) increase",drugLpa[,2]),
                   c("Number","Lp(a) decrease",drugLpa[,1]),
                   c("HR","\n", "1.00",datFrame[,"coef"]),
                   c("P","\n", "-",datFrame[,"Pvalue"]))
tabletext 

tiff(filename = "drugLpa.tiff",
     width = 5480, height = 3480, units = "px", pointsize = 12,
     compression ="lzw",
     bg = "white", res = 300)

forestplot(tabletext, 
           txt_gp=fpTxtGp(
             label=gpar(cex=2),
             ticks=gpar(cex=2),
             xlab=gpar(cex = 1.7),
             title=gpar(cex = 1.8)),
           rbind( rep(NA,3),rep(NA,3),rep(1,3),
                  datFrame[,-4]),
           col=clrs,
           colgap=unit(8,"mm"),
           lwd.ci=2, boxsize=0.5,
           zero=1,
           cex=2, lineheight = "auto",
           ci.vertices=TRUE, ci.vertices.height = 0.4,
           xlab="    <---Decrease--- Lp(a) ---Increase--->"
)
dev.off()










#################################################ALLDataSet
#=================================outPatients

nrow(dataOutDemoLabIntegration)
nrow(dataInDemoLabIntegration)
head(dataOutDemoLabIntegration)
head(dataInDemoLabIntegration)
7641+27944

nrow(dataOutDemoLabIntegrationPsm)
nrow(dataInDemoLabIntegrationPsm)
4888+ 16488


intersect(dataOutDemoLabIntegration$MASTER_INDEX,dataInDemoLabIntegration$MASTER_INDEX)

dataOutDemoLabIntegrationRaw = dataOutDemoLabIntegration

dataOutDemoLabIntegration = dataOutDemoLabIntegration[!dataOutDemoLabIntegration$MASTER_INDEX%in%intersect(dataOutDemoLabIntegration$MASTER_INDEX,dataInDemoLabIntegration$MASTER_INDEX),]
nrow(dataOutDemoLabIntegration)


dataOutDemoLabIntegration# outpatients dataset 

###################outpatients

dataOutDemoLabIntegration_1 = na.omit(dataOutDemoLabIntegration)
nrow(dataOutDemoLabIntegration_1 )
table(dataOutDemoLabIntegration_1$tating)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$gender==3,]
nrow(dataOutDemoLabIntegration_1 )

dataOutDemoLabIntegration_1 = 
  na.omit(dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$gender==3,])
breaks = c(0,45,65,150)
label = c(1,2,3)
dataOutDemoLabIntegration_1$ageGroup= cut(dataOutDemoLabIntegration_1$age,breaks = breaks,labels = label)
dataOutDemoLabIntegration_1 = na.omit(dataOutDemoLabIntegration_1)
nrow(dataOutDemoLabIntegration_1 )

dataOutDemoLabIntegration_1 = merge(dataOutDemoLabIntegration_1,followUP,by.x = "MASTER_INDEX",by.y = "MASTER_INDEX")
nrow(dataOutDemoLabIntegration_1 )

#dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month>6,]
#dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<60,]

summary(dataOutDemoLabIntegration_1 )
breaks = c(0,36,60,200)
label = c(1,2,3)
dataOutDemoLabIntegration_1$monthGroup= cut(dataOutDemoLabIntegration_1$month,breaks = breaks,labels = label)
table(dataOutDemoLabIntegration_1$tating)
table(dataOutDemoLabIntegration_1$monthGroup,dataOutDemoLabIntegration_1$tating)
nrow(dataOutDemoLabIntegration_1)
#dataOutDemoLabIntegration_1= dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<36,]


####################disease modify
dataOutDemoLabIntegration_1 = 
  na.omit(dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$dis_group==5,])
nrow(dataOutDemoLabIntegration_1 )
dataOutDemoLabIntegration_1$dis_group.y = 0
dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$dis_group%in%c(1:4),]$dis_group.y = 1
dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$dis_group%in%c(6),]$dis_group.y = 2
table(dataOutDemoLabIntegration_1$dis_group)
table(dataOutDemoLabIntegration_1$dis_group.y)

#table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]

max(table(dataOutDemoLabIntegration_1$MASTER_INDEX))

#dataAllDemoLabIntegration = 
#  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]


dataOutDemoLabIntegration_1C = dataOutDemoLabIntegration_1
nrow(dataOutDemoLabIntegration_1C)
table(dataOutDemoLabIntegration_1C$tating)



dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==1,]
dataOutDemoLabIntegrationPsm_1 = psmOut(dataOutDemoLabIntegration_1)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==2,]
dataOutDemoLabIntegrationPsm_2 = psmOut(dataOutDemoLabIntegration_1)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==3,]
dataOutDemoLabIntegrationPsm_3 = psmOut(dataOutDemoLabIntegration_1)

dataOutDemoLabIntegrationPsm= rbind(
  dataOutDemoLabIntegrationPsm_1,
  dataOutDemoLabIntegrationPsm_2,
  dataOutDemoLabIntegrationPsm_3
)



colnames(dataOutDemoLabIntegrationPsm)[c(3,11:17)]
c(3,11:17)
rownames = c("Lp(a) at FT (Mean,ci95%)","LDL-C (Mean,ci95%)",
             "HDL-C (Mean,ci95%)","APO-A (Mean,ci95%)","APO-B (Mean,ci95%)",
             "TC (Mean,ci95%)","TG (Mean,ci95%)","CRP (Mean,ci95%)")
num = table(dataOutDemoLabIntegrationPsm$tating)

colnames = c(paste("Control","(","n=",num[1],")",sep = ""),
             paste("Treat","(","n=",num[2],")",sep = ""),
             "StaVal","PVal")

colnames(dataOutDemoLabIntegrationPsm)[c(3,11:17)]
tab_1 = data.frame()
for(i in c(3,11:17)){
  tmp = dataOutDemoLabIntegrationPsm[c(9,i)]
  tmp = as.data.frame(descri2(tmp)) # mean1+ci95%,mean2+ci95%,staticValue,Pvalue
  tab_1 = rbind(tab_1,tmp)
}

colnames(tab_1) = colnames
rownames(tab_1) = rownames
tab_1
colnames(dataOutDemoLabIntegrationPsm)[c(7,28,31,27,10)]
tab_2 = data.frame()
#gender, age, Disease in FE, Disease in FT,lp(a) outcome,month
for (i in c(7,28,31,27,10)){
  M = dataOutDemoLabIntegrationPsm[c(i,9)]
  tmp = decriCount(M)
  tab_2 = rbind(tab_2,tmp)
}
colnames(tab_2)=colnames(tab_1)
tab_2
write.csv(file ="outPatientsDemo.csv" , rbind(tab_1,tab_2))


####################################Inpatients


nrow((dataInDemoLabIntegration))
head(dataInDemoLabIntegration)
dataInDemoLabIntegration_1 = na.omit(dataInDemoLabIntegration)
nrow(dataInDemoLabIntegration_1 )
table(table(dataInDemoLabIntegration_1$tating))
dataInDemoLabIntegration_1 = 
  na.omit(dataInDemoLabIntegration_1[!dataInDemoLabIntegration_1$gender==3,])
nrow(dataInDemoLabIntegration_1 )
breaks = c(0,45,60,65,70,80,1000)
label = c(1,2,3,4,5,6)
dataInDemoLabIntegration_1$ageGroup= cut(dataInDemoLabIntegration_1$age,breaks = breaks,labels = label)
nrow(dataInDemoLabIntegration_1 )
table(dataInDemoLabIntegration_1$ageGroup)
sum(table(dataInDemoLabIntegration_1$ageGroup))
dataInDemoLabIntegration_1 = na.omit(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1 = dataInDemoLabIntegration_1[order(dataInDemoLabIntegration_1$tating,decreasing = T),]
nrow(dataInDemoLabIntegration_1 )

table(followUP$MASTER_INDEX)
followUP[followUP$MASTER_INDEX == 2644,]


dataInDemoLabIntegration_1 = merge(dataInDemoLabIntegration_1,followUP,by.x = "MASTER_INDEX",by.y = "MASTER_INDEX",all.x = T)
nrow(dataInDemoLabIntegration_1)
# dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month>6,]
# dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<60,]

table(dataInDemoLabIntegration$MASTER_INDEX)[table(dataInDemoLabIntegration$MASTER_INDEX)>=2]
#dataAllDemoLabIntegration = 
#  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]


breaks = c(0,36,60,200)
label = c(1,2,3)
dataInDemoLabIntegration_1$monthGroup= cut(dataInDemoLabIntegration_1$month,breaks = breaks,labels = label)
table(dataInDemoLabIntegration_1$tating)
table(dataInDemoLabIntegration_1$monthGroup,dataInDemoLabIntegration_1$tating)




####################disease modify
dataInDemoLabIntegration_1 = 
  na.omit(dataInDemoLabIntegration_1[!dataInDemoLabIntegration_1$dis_group==5,])
nrow(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1$dis_group.y=1
dataInDemoLabIntegration_1[dataInDemoLabIntegration_1$dis_group%in%c(1:4),]$dis_group.y = 1
dataInDemoLabIntegration_1[dataInDemoLabIntegration_1$dis_group%in%c(6),]$dis_group.y = 2
table(dataInDemoLabIntegration_1$dis_group)

dataInDemoLabIntegration_1C = dataInDemoLabIntegration_1
nrow(dataInDemoLabIntegration_1C )

table(dataInDemoLabIntegration_1C$ageGroup)
table(dataInDemoLabIntegration_1C$tating)

#write.csv(dataInDemoLabIntegration_1,file="dataInDemoLabIntegration_1.csv")
#table(dataInDemoLabIntegration_1C$tating)

dataInDemoLabIntegration_1C = (dataInDemoLabIntegration_1C[!dataInDemoLabIntegration_1C$MASTER_INDEX%in%intersect(dataOutDemoLabIntegrationPsm$master_id,dataInDemoLabIntegration_1C$master_id),])
nrow(dataInDemoLabIntegration_1C)
head(dataInDemoLabIntegration_1C)


dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==1,]
dataInDemoLabIntegrationPsm_1 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==2,]
dataInDemoLabIntegrationPsm_2 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==3,]
dataInDemoLabIntegrationPsm_3 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==4,]
dataInDemoLabIntegrationPsm_4 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==5,]
dataInDemoLabIntegrationPsm_5 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==6,]
dataInDemoLabIntegrationPsm_6 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegrationPsm = rbind(dataInDemoLabIntegrationPsm_1,
                                    dataInDemoLabIntegrationPsm_2,
                                    dataInDemoLabIntegrationPsm_3,
                                    dataInDemoLabIntegrationPsm_4,
                                    dataInDemoLabIntegrationPsm_5,
                                    dataInDemoLabIntegrationPsm_6
)
table(dataInDemoLabIntegrationPsm$tating)
table(dataOutDemoLabIntegrationPsm$tating)
breaks = c(0,45,65,150)
label = c(1,2,3)
dataInDemoLabIntegrationPsm$ageGroup= cut(dataInDemoLabIntegrationPsm$age,
                                          breaks = breaks,labels = label)
nrow(dataInDemoLabIntegrationPsm)
nrow(dataOutDemoLabIntegrationPsm)



dataAllDemoLabIntegrationPsm = rbind(dataOutDemoLabIntegrationPsm,dataInDemoLabIntegrationPsm)
nrow(dataAllDemoLabIntegrationPsm)
breaks = c(0,45,65,150)
label = c(1,2,3)
dataAllDemoLabIntegrationPsm$ageGroup= cut(dataAllDemoLabIntegrationPsm$age,
                                           breaks = breaks,labels = label)

table(dataOutDemoLabIntegrationPsm$master_id)[table(dataOutDemoLabIntegrationPsm$master_id)>=2]

table(dataInDemoLabIntegrationPsm$master_id)[table(dataInDemoLabIntegrationPsm$master_id)>=2]

intersect(dataInDemoLabIntegrationPsm$master_id,dataOutDemoLabIntegrationPsm$master_id)


table(dataAllDemoLabIntegrationPsm$master_id)[table(dataAllDemoLabIntegrationPsm$master_id)>=2]

#dataAllDemoLabIntegrationPsm = dataAllDemoLabIntegrationPsm[!dataAllDemoLabIntegrationPsm$MASTER_INDEX%in%names(table(CVDmodify$master_id)[table(CVDmodify$master_id)>=2]),]
nrow(dataAllDemoLabIntegrationPsm)
#CVDmodify = read.table("CVDmodified.txt",header = T,sep = "\t")


#dataAllDemoLabIntegrationPsm = merge(dataAllDemoLabIntegrationPsm,CVDmodify,
#                                     by.x = "MASTER_INDEX",by.y = "master_id",all.y  = F)
nrow(dataAllDemoLabIntegrationPsm)


dataAllDemoLabIntegration = rbind(dataOutDemoLabIntegration,dataInDemoLabIntegration)
nrow(dataAllDemoLabIntegration)

table(rbind(dataOutDemoLabIntegration_1C,dataInDemoLabIntegration_1C)$tating)

head(dataAllDemoLabIntegration)
colnames(dataAllDemoLabIntegration)
dataAllDemoLabIntegration = dataAllDemoLabIntegration[,c(1,24)]
head(dataAllDemoLabIntegration)
dataAllDemoLabIntegration = 
  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]
head(dataAllDemoLabIntegration)
head(dataAllDemoLabIntegrationPsm)

dataAllDemoLabIntegrationPsm_1=merge(dataAllDemoLabIntegrationPsm,dataAllDemoLabIntegration,
                                     by.x = "MASTER_INDEX",by.y = "MASTER_INDEX",all.x = F)
nrow(dataAllDemoLabIntegrationPsm_1)

dataAllDemoLabIntegrationPsm = dataAllDemoLabIntegrationPsm_1



rownames = c("Lp(a) at FT (Mean,CI95%)","LDL-C at FT (Mean,CI95%)","LDL-C (Mean,CI95%)",
             "HDL-C (Mean,CI95%)","APO-A (Mean,CI95%)","APO-B (Mean,CI95%)",
             "TC (Mean,CI95%)","TG (Mean,CI95%)","CRP (Mean,CI95%)")
num = table(dataAllDemoLabIntegrationPsm$tating)

colnames = c(paste("Control","(","n=",num[1],")",sep = ""),
             paste("Treat","(","n=",num[2],")",sep = ""),
             "StaVal","PVal")
colnames(dataAllDemoLabIntegrationPsm)
tab_1 = data.frame()
for(i in c(3,12,15:21)){
  tmp = dataAllDemoLabIntegrationPsm[c(9,i)]
  tmp = as.data.frame(descri2(tmp)) # mean1+ci95%,mean2+ci95%,staticValue,Pvalue
  tab_1 = rbind(tab_1,tmp)
}

colnames(tab_1) = colnames
rownames(tab_1) = rownames
tab_1
dataAllDemoLabIntegrationPsm$CRP
tab_2 = data.frame()
colnames(dataAllDemoLabIntegrationPsm)[ c(7,28,31,30,10,14,27,25)]
#gender,age,dis_history,follow_time,LPA_outcome,LDL_outcome,dis_follow,
for (i in c(7,28,31,30,10,14,27,25)){
  M = dataAllDemoLabIntegrationPsm[c(i,9)]
  tmp = decriCount(M)
  tab_2 = rbind(tab_2,tmp)
}
tab_2
colnames(tab_2)=colnames(tab_1)
nrow(tab_2)
rownames(tab_2) = c("Male","Female","<45","46-65",">65",
                    "Non-CVD-related disease(FE)","CVD-related disease(FE)",
                    "[0.5 - 3) years","[3 - 5) years","≥ 5 years",
                    "Lp(a) decrease (FU)","Lp(a) increase (FU)",
                    "LDL-C decrease (FU)","LDL-C increase (FU)",
                    "Non-CVD disease event (FU)","CVD disease event (FU)",
                    #"Modify Non-CVD disease event (FU)","Modify CVD disease event (FU)",
                    "diabetes","dyslipidemia","hypertension","arteriosclerosis","Others"
                    #"Lp(a) decrease_1 (FU)","Lp(a) increase_1 (FU)",
                    #"LDL-C decrease_1 (FU)","LDL-C increase_1 (FU)"
)
tab_2


rbind(tab_1,tab_2)


write.csv(file = "ALLDataset.csv",rbind(tab_1,tab_2))

write.csv(file = "ALLDatasetRaw.csv",dataAllDemoLabIntegrationPsm)

dat = read.csv("ALLDatasetRaw.csv")

dataAllDemoLabIntegrationPsm=dat[,-1]

###############################

resOR=c()
resORCi=c()

colnames(dataAllDemoLabIntegrationPsm)

nrow( dataAllDemoLabIntegrationPsm)

glmIn = coxph(Surv(month,outcome.x) ~  tating,
              data = dataAllDemoLabIntegrationPsm)


print(summary(glmIn))

tmp = summary(glmIn)

resOR[1] = formatC(tmp$coefficients[1,2],format = "f",digits = 2)
resORCi[1] = paste("(",
                   formatC(exp(tmp$coefficients[1,1]-1.96*tmp$coefficients[1,3]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[1,1]+1.96*tmp$coefficients[1,3]),
                           format = "f",digits = 2),")",sep = ""
)
# glmIn = glm(formula = disFollow_group_new  ~   tating , 
#             data = dataAllDemoLabIntegrationPsm, family = binomial)

glmIn = glm(formula = disFollow_group  ~   tating , 
            data = dataAllDemoLabIntegrationPsm, family = binomial)
print(summary(glmIn))
tmp = summary(glmIn)
resOR[2] = formatC(exp(tmp$coefficients[2,1]),format = "f",digits = 2)
resORCi[2] = paste("(",
                   formatC(exp(tmp$coefficients[2,1]-1.96*tmp$coefficients[2,2]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[2,1]+1.96*tmp$coefficients[2,2]),
                           format = "f",digits = 2),")",sep = ""
)
resOR
resORCi

glmIn = coxph(Surv(month,outcome.x) ~  as.numeric(ageGroup) +  tating +CRP + LDL_C + dis_group + HDL_C + 
                APO_A +  APO_B,
              data = dataAllDemoLabIntegrationPsm)
print(summary(glmIn))

tmp = summary(glmIn)

resOR[3] = formatC(exp(tmp$coefficients[2,1]),format = "f",digits = 2)
resORCi[3] = paste("(",
                   formatC(exp(tmp$coefficients[2,1]-1.96*tmp$coefficients[2,3]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[2,1]+1.96*tmp$coefficients[2,3]),
                           format = "f",digits = 2),")",sep = ""
)
resOR
resORCi
tmp = as.data.frame(tmp$coefficients)


tmp$OR=formatC(tmp$`exp(coef)`,format = "f",digits = 2)

tmp$`95%CI` = paste("(",formatC((tmp$`exp(coef)`-1.96*tmp$`se(coef)`),format = "f",digits = 2),","
                    ,formatC((tmp$`exp(coef)`+1.96*tmp$`se(coef)`),format = "f",digits = 2),")",sep="")
tmp
tmp[,5] = formatC(tmp[,5],format = "f",digits = 2)
tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)

rownames(tmp) = c("Age","statin use","CRP","LDL-C",
                  "CVD-related disease(FE)","HDL-C","APO-A","APO-B")
tmp = tmp[,c(6,7,4,5)]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")

write.csv(tmp,file = "statinMultiOutcome.csv")


# glmIn = glm(formula = disFollow_group_new ~ as.numeric(ageGroup) +  tating + CRP + as.numeric(monthGroup)
#             + LDL_C + dis_group.x + HDL_C + APO_A +  APO_B, 
#             data = dataAllDemoLabIntegrationPsm, family = binomial)

glmIn = glm(formula = disFollow_group ~ as.numeric(ageGroup) +  tating + CRP + as.numeric(monthGroup)
            + LDL_C + dis_group+ HDL_C + APO_A +  APO_B, 
            data = dataAllDemoLabIntegrationPsm, family = binomial)



print(summary(glmIn))
tmp = summary(glmIn)
resOR[4] = formatC(exp(tmp$coefficients[3,1]),format = "f",digits = 2)
resORCi[4] = paste("(",
                   formatC(exp(tmp$coefficients[3,1]-1.96*tmp$coefficients[3,2]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[3,1]+1.96*tmp$coefficients[3,2]),
                           format = "f",digits = 2),")",sep = ""
)

tmp = as.data.frame(tmp$coefficients)


tmp$OR=formatC(exp(tmp$Estimate),format = "f",digits = 2)

tmp$`95%ci` = paste("(",formatC(exp(tmp$Estimate-1.96*tmp$`Std. Error`),format = "f",digits = 2),","
                    ,formatC(exp(tmp$Estimate+1.96*tmp$`Std. Error`),format = "f",digits = 2),")",sep="")

tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)
tmp[,3] = formatC(tmp[,3],format = "f",digits = 2)

rownames(tmp) = c("Intercept","Age","statin use","CRP","Follow-up(months)","LDL-C",
                  "CVD-related disease(FE)","HDL-C","APO-A","APO-B")
tmp = tmp[,c(5,6,3,4)]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")
# 
# write.csv(tmp,file = "followMultiOutcome.csv")
# 
# write.csv(tmp,file = "lpaCVD.csv")

nrow(dataAllDemoLabIntegrationPsm)

mitifactors =  data.frame(resOR,resORCi)

rownames(mitifactors) = c("Single statin","Single Dis","Muti statin","Muti Dis")

write.csv(mitifactors,file = "mitifactorsALLData.csv")

#################################################NormalDataSet
#=================================outPatients

nrow(dataOutDemoLabIntegration)
nrow(dataInDemoLabIntegration)
head(dataOutDemoLabIntegration)
head(dataInDemoLabIntegration)
7641+27944

nrow(dataOutDemoLabIntegrationPsm)
nrow(dataInDemoLabIntegrationPsm)
4888+ 16488


intersect(dataOutDemoLabIntegration$MASTER_INDEX,dataInDemoLabIntegration$MASTER_INDEX)

dataOutDemoLabIntegrationRaw = dataOutDemoLabIntegration

dataOutDemoLabIntegration = dataOutDemoLabIntegration[!dataOutDemoLabIntegration$MASTER_INDEX%in%intersect(dataOutDemoLabIntegration$MASTER_INDEX,dataInDemoLabIntegration$MASTER_INDEX),]
nrow(dataOutDemoLabIntegration)


dataOutDemoLabIntegration# outpatients dataset 

###################outpatients

dataOutDemoLabIntegration_1 = na.omit(dataOutDemoLabIntegration[dataOutDemoLabIntegration$LDL_C<=1.8,])
nrow(dataOutDemoLabIntegration_1 )
table(dataOutDemoLabIntegration_1$tating)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$gender==3,]
nrow(dataOutDemoLabIntegration_1 )

dataOutDemoLabIntegration_1 = 
  na.omit(dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$gender==3,])
breaks = c(0,45,65,150)
label = c(1,2,3)
dataOutDemoLabIntegration_1$ageGroup= cut(dataOutDemoLabIntegration_1$age,breaks = breaks,labels = label)
dataOutDemoLabIntegration_1 = na.omit(dataOutDemoLabIntegration_1)
nrow(dataOutDemoLabIntegration_1 )

dataOutDemoLabIntegration_1 = merge(dataOutDemoLabIntegration_1,followUP,by.x = "MASTER_INDEX",by.y = "MASTER_INDEX")
nrow(dataOutDemoLabIntegration_1 )

#dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month>6,]
#dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<60,]

summary(dataOutDemoLabIntegration_1 )
breaks = c(0,36,60,200)
label = c(1,2,3)
dataOutDemoLabIntegration_1$monthGroup= cut(dataOutDemoLabIntegration_1$month,breaks = breaks,labels = label)
table(dataOutDemoLabIntegration_1$tating)
table(dataOutDemoLabIntegration_1$monthGroup,dataOutDemoLabIntegration_1$tating)
nrow(dataOutDemoLabIntegration_1)
#dataOutDemoLabIntegration_1= dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<36,]


####################disease modify
dataOutDemoLabIntegration_1 = 
  na.omit(dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$dis_group==5,])
nrow(dataOutDemoLabIntegration_1 )
dataOutDemoLabIntegration_1$dis_group.y = 0
dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$dis_group%in%c(1:4),]$dis_group.y = 1
dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$dis_group%in%c(6),]$dis_group.y = 2
table(dataOutDemoLabIntegration_1$dis_group)
table(dataOutDemoLabIntegration_1$dis_group.y)

#table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]

max(table(dataOutDemoLabIntegration_1$MASTER_INDEX))

#dataAllDemoLabIntegration = 
#  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]


dataOutDemoLabIntegration_1C = dataOutDemoLabIntegration_1
nrow(dataOutDemoLabIntegration_1C)
table(dataOutDemoLabIntegration_1C$tating)



dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==1,]
dataOutDemoLabIntegrationPsm_1 = psmOut(dataOutDemoLabIntegration_1)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==2,]
dataOutDemoLabIntegrationPsm_2 = psmOut(dataOutDemoLabIntegration_1)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==3,]
dataOutDemoLabIntegrationPsm_3 = psmOut(dataOutDemoLabIntegration_1)

dataOutDemoLabIntegrationPsm= rbind(
  dataOutDemoLabIntegrationPsm_1,
  dataOutDemoLabIntegrationPsm_2,
  dataOutDemoLabIntegrationPsm_3
)



colnames(dataOutDemoLabIntegrationPsm)[c(3,11:17)]
c(3,11:17)
rownames = c("Lp(a) at FT (Mean,ci95%)","LDL-C (Mean,ci95%)",
             "HDL-C (Mean,ci95%)","APO-A (Mean,ci95%)","APO-B (Mean,ci95%)",
             "TC (Mean,ci95%)","TG (Mean,ci95%)","CRP (Mean,ci95%)")
num = table(dataOutDemoLabIntegrationPsm$tating)

colnames = c(paste("Control","(","n=",num[1],")",sep = ""),
             paste("Treat","(","n=",num[2],")",sep = ""),
             "StaVal","PVal")

colnames(dataOutDemoLabIntegrationPsm)[c(3,11:17)]
tab_1 = data.frame()
for(i in c(3,11:17)){
  tmp = dataOutDemoLabIntegrationPsm[c(9,i)]
  tmp = as.data.frame(descri2(tmp)) # mean1+ci95%,mean2+ci95%,staticValue,Pvalue
  tab_1 = rbind(tab_1,tmp)
}

colnames(tab_1) = colnames
rownames(tab_1) = rownames
tab_1
colnames(dataOutDemoLabIntegrationPsm)[c(7,28,31,27,10)]
tab_2 = data.frame()
#gender, age, Disease in FE, Disease in FT,lp(a) outcome,month
for (i in c(7,28,31,27,10)){
  M = dataOutDemoLabIntegrationPsm[c(i,9)]
  tmp = decriCount(M)
  tab_2 = rbind(tab_2,tmp)
}
colnames(tab_2)=colnames(tab_1)
tab_2
write.csv(file ="outPatientsDemo.csv" , rbind(tab_1,tab_2))


####################################Inpatients


nrow((dataInDemoLabIntegration))
head(dataInDemoLabIntegration)
dataInDemoLabIntegration_1 = na.omit(dataInDemoLabIntegration[dataInDemoLabIntegration$LDL_C<=1.8,])
nrow(dataInDemoLabIntegration_1 )
table(table(dataInDemoLabIntegration_1$tating))
dataInDemoLabIntegration_1 = 
  na.omit(dataInDemoLabIntegration_1[!dataInDemoLabIntegration_1$gender==3,])
nrow(dataInDemoLabIntegration_1 )
breaks = c(0,45,60,65,70,80,1000)
label = c(1,2,3,4,5,6)
dataInDemoLabIntegration_1$ageGroup= cut(dataInDemoLabIntegration_1$age,breaks = breaks,labels = label)
nrow(dataInDemoLabIntegration_1 )
table(dataInDemoLabIntegration_1$ageGroup)
sum(table(dataInDemoLabIntegration_1$ageGroup))
dataInDemoLabIntegration_1 = na.omit(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1 = dataInDemoLabIntegration_1[order(dataInDemoLabIntegration_1$tating,decreasing = T),]
nrow(dataInDemoLabIntegration_1 )

table(followUP$MASTER_INDEX)
followUP[followUP$MASTER_INDEX == 2644,]


dataInDemoLabIntegration_1 = merge(dataInDemoLabIntegration_1,followUP,by.x = "MASTER_INDEX",by.y = "MASTER_INDEX",all.x = T)
nrow(dataInDemoLabIntegration_1)
# dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month>6,]
# dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<60,]

table(dataInDemoLabIntegration$MASTER_INDEX)[table(dataInDemoLabIntegration$MASTER_INDEX)>=2]
#dataAllDemoLabIntegration = 
#  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]


breaks = c(0,36,60,200)
label = c(1,2,3)
dataInDemoLabIntegration_1$monthGroup= cut(dataInDemoLabIntegration_1$month,breaks = breaks,labels = label)
table(dataInDemoLabIntegration_1$tating)
table(dataInDemoLabIntegration_1$monthGroup,dataInDemoLabIntegration_1$tating)




####################disease modify
dataInDemoLabIntegration_1 = 
  na.omit(dataInDemoLabIntegration_1[!dataInDemoLabIntegration_1$dis_group==5,])
nrow(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1$dis_group.y=1
dataInDemoLabIntegration_1[dataInDemoLabIntegration_1$dis_group%in%c(1:4),]$dis_group.y = 1
dataInDemoLabIntegration_1[dataInDemoLabIntegration_1$dis_group%in%c(6),]$dis_group.y = 2
table(dataInDemoLabIntegration_1$dis_group)

dataInDemoLabIntegration_1C = dataInDemoLabIntegration_1
nrow(dataInDemoLabIntegration_1C )

table(dataInDemoLabIntegration_1C$ageGroup)
table(dataInDemoLabIntegration_1C$tating)

#write.csv(dataInDemoLabIntegration_1,file="dataInDemoLabIntegration_1.csv")
#table(dataInDemoLabIntegration_1C$tating)

dataInDemoLabIntegration_1C = (dataInDemoLabIntegration_1C[!dataInDemoLabIntegration_1C$MASTER_INDEX%in%intersect(dataOutDemoLabIntegrationPsm$master_id,dataInDemoLabIntegration_1C$master_id),])
nrow(dataInDemoLabIntegration_1C)
head(dataInDemoLabIntegration_1C)


dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==1,]
dataInDemoLabIntegrationPsm_1 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==2,]
dataInDemoLabIntegrationPsm_2 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==3,]
dataInDemoLabIntegrationPsm_3 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==4,]
dataInDemoLabIntegrationPsm_4 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==5,]
dataInDemoLabIntegrationPsm_5 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==6,]
dataInDemoLabIntegrationPsm_6 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegrationPsm = rbind(dataInDemoLabIntegrationPsm_1,
                                    dataInDemoLabIntegrationPsm_2,
                                    dataInDemoLabIntegrationPsm_3,
                                    dataInDemoLabIntegrationPsm_4,
                                    dataInDemoLabIntegrationPsm_5,
                                    dataInDemoLabIntegrationPsm_6
)
table(dataInDemoLabIntegrationPsm$tating)
table(dataOutDemoLabIntegrationPsm$tating)
breaks = c(0,45,65,150)
label = c(1,2,3)
dataInDemoLabIntegrationPsm$ageGroup= cut(dataInDemoLabIntegrationPsm$age,
                                          breaks = breaks,labels = label)
nrow(dataInDemoLabIntegrationPsm)
nrow(dataOutDemoLabIntegrationPsm)



dataAllDemoLabIntegrationPsm = rbind(dataOutDemoLabIntegrationPsm,dataInDemoLabIntegrationPsm)
nrow(dataAllDemoLabIntegrationPsm)
breaks = c(0,45,65,150)
label = c(1,2,3)
dataAllDemoLabIntegrationPsm$ageGroup= cut(dataAllDemoLabIntegrationPsm$age,
                                           breaks = breaks,labels = label)

table(dataOutDemoLabIntegrationPsm$master_id)[table(dataOutDemoLabIntegrationPsm$master_id)>=2]

table(dataInDemoLabIntegrationPsm$master_id)[table(dataInDemoLabIntegrationPsm$master_id)>=2]

intersect(dataInDemoLabIntegrationPsm$master_id,dataOutDemoLabIntegrationPsm$master_id)


table(dataAllDemoLabIntegrationPsm$master_id)[table(dataAllDemoLabIntegrationPsm$master_id)>=2]

#dataAllDemoLabIntegrationPsm = dataAllDemoLabIntegrationPsm[!dataAllDemoLabIntegrationPsm$MASTER_INDEX%in%names(table(CVDmodify$master_id)[table(CVDmodify$master_id)>=2]),]
nrow(dataAllDemoLabIntegrationPsm)
#CVDmodify = read.table("CVDmodified.txt",header = T,sep = "\t")


#dataAllDemoLabIntegrationPsm = merge(dataAllDemoLabIntegrationPsm,CVDmodify,
#                                     by.x = "MASTER_INDEX",by.y = "master_id",all.y  = F)
nrow(dataAllDemoLabIntegrationPsm)


dataAllDemoLabIntegration = rbind(dataOutDemoLabIntegration,dataInDemoLabIntegration)
nrow(dataAllDemoLabIntegration)

table(rbind(dataOutDemoLabIntegration_1C,dataInDemoLabIntegration_1C)$tating)

head(dataAllDemoLabIntegration)
colnames(dataAllDemoLabIntegration)
dataAllDemoLabIntegration = dataAllDemoLabIntegration[,c(1,24)]
head(dataAllDemoLabIntegration)
dataAllDemoLabIntegration = 
  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]
head(dataAllDemoLabIntegration)
head(dataAllDemoLabIntegrationPsm)

dataAllDemoLabIntegrationPsm_1=merge(dataAllDemoLabIntegrationPsm,dataAllDemoLabIntegration,
                                     by.x = "MASTER_INDEX",by.y = "MASTER_INDEX",all.x = F)
nrow(dataAllDemoLabIntegrationPsm_1)

dataAllDemoLabIntegrationPsm = dataAllDemoLabIntegrationPsm_1



rownames = c("Lp(a) at FT (Mean,CI95%)","LDL-C at FT (Mean,CI95%)","LDL-C (Mean,CI95%)",
             "HDL-C (Mean,CI95%)","APO-A (Mean,CI95%)","APO-B (Mean,CI95%)",
             "TC (Mean,CI95%)","TG (Mean,CI95%)","CRP (Mean,CI95%)")
num = table(dataAllDemoLabIntegrationPsm$tating)

colnames = c(paste("Control","(","n=",num[1],")",sep = ""),
             paste("Treat","(","n=",num[2],")",sep = ""),
             "StaVal","PVal")
colnames(dataAllDemoLabIntegrationPsm)
tab_1 = data.frame()
for(i in c(3,12,15:21)){
  tmp = dataAllDemoLabIntegrationPsm[c(9,i)]
  tmp = as.data.frame(descri2(tmp)) # mean1+ci95%,mean2+ci95%,staticValue,Pvalue
  tab_1 = rbind(tab_1,tmp)
}

colnames(tab_1) = colnames
rownames(tab_1) = rownames
tab_1
dataAllDemoLabIntegrationPsm$CRP
tab_2 = data.frame()
colnames(dataAllDemoLabIntegrationPsm)[ c(7,28,31,30,10,14,27,25)]
#gender,age,dis_history,follow_time,LPA_outcome,LDL_outcome,dis_follow,
for (i in c(7,28,31,30,10,14,27,25)){
  M = dataAllDemoLabIntegrationPsm[c(i,9)]
  tmp = decriCount(M)
  tab_2 = rbind(tab_2,tmp)
}
tab_2
colnames(tab_2)=colnames(tab_1)
nrow(tab_2)
rownames(tab_2) = c("Male","Female","<45","46-65",">65",
                    "Non-CVD-related disease(FE)","CVD-related disease(FE)",
                    "[0.5 - 3) years","[3 - 5) years","≥ 5 years",
                    "Lp(a) decrease (FU)","Lp(a) increase (FU)",
                    "LDL-C decrease (FU)","LDL-C increase (FU)",
                    "Non-CVD disease event (FU)","CVD disease event (FU)",
                    #"Modify Non-CVD disease event (FU)","Modify CVD disease event (FU)",
                    "diabetes","dyslipidemia","hypertension","arteriosclerosis","Others"
                    #"Lp(a) decrease_1 (FU)","Lp(a) increase_1 (FU)",
                    #"LDL-C decrease_1 (FU)","LDL-C increase_1 (FU)"
)
tab_2


rbind(tab_1,tab_2)


write.csv(file = "NormalDataset.csv",rbind(tab_1,tab_2))

write.csv(file = "NormalDatasetRaw.csv",dataAllDemoLabIntegrationPsm)

dat =read.csv("NormalDatasetRaw.csv")

dataAllDemoLabIntegrationPsm = dat[,-1]

###############################

resOR=c()
resORCi=c()

colnames(dataAllDemoLabIntegrationPsm)

nrow( dataAllDemoLabIntegrationPsm)

glmIn = coxph(Surv(month,outcome.x) ~  tating,
              data = dataAllDemoLabIntegrationPsm)

print(summary(glmIn))

tmp = summary(glmIn)

exp(tmp$coefficients)

resOR[1] = formatC(tmp$coefficients[1,2],format = "f",digits = 2)
resORCi[1] = paste("(",
                   formatC(exp(tmp$coefficients[1,1]-1.96*tmp$coefficients[1,3]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[1,1]+1.96*tmp$coefficients[1,3]),
                           format = "f",digits = 2),")",sep = ""
)

glmIn = glm(formula = disFollow_group  ~   tating , 
            data = dataAllDemoLabIntegrationPsm, family = binomial)
print(summary(glmIn))
tmp = summary(glmIn)
resOR[2] = formatC(exp(tmp$coefficients[2,1]),format = "f",digits = 2)
resORCi[2] = paste("(",
                   formatC(exp(tmp$coefficients[2,1]-1.96*tmp$coefficients[2,2]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[2,1]+1.96*tmp$coefficients[2,2]),
                           format = "f",digits = 2),")",sep = ""
)
resOR
resORCi
glmIn = coxph(Surv(month,outcome.x) ~  as.numeric(ageGroup) +  tating +CRP + LDL_C + dis_group + HDL_C + 
                APO_A +  APO_B,
              data = dataAllDemoLabIntegrationPsm)

print(summary(glmIn))

tmp = summary(glmIn)
resOR[3] = formatC(exp(tmp$coefficients[2,1]),format = "f",digits = 2)
resORCi[3] = paste("(",
                   formatC(exp(tmp$coefficients[2,1]-1.96*tmp$coefficients[2,3]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[2,1]+1.96*tmp$coefficients[2,3]),
                           format = "f",digits = 2),")",sep = ""
)
resOR
resORCi
tmp = as.data.frame(tmp$coefficients)


tmp$OR=formatC(tmp$`exp(coef)`,format = "f",digits = 2)

tmp$`95%CI` = paste("(",formatC((tmp$`exp(coef)`-1.96*tmp$`se(coef)`),format = "f",digits = 2),","
                    ,formatC((tmp$`exp(coef)`+1.96*tmp$`se(coef)`),format = "f",digits = 2),")",sep="")
tmp
tmp[,5] = formatC(tmp[,5],format = "f",digits = 2)
tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)

rownames(tmp) = c("Age","statin use","CRP","LDL-C",
                  "CVD-related disease(FE)","HDL-C","APO-A","APO-B")
tmp = tmp[,c(6,7,4,5)]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")

write.csv(tmp,file = "statinMultiOutcome.csv")

# glmIn = glm(formula = disFollow_group_new ~ as.numeric(ageGroup) +  tating + CRP + as.numeric(monthGroup)
#             + LDL_C + dis_group.x + HDL_C + APO_A +  APO_B, 
#             data = dataAllDemoLabIntegrationPsm, family = binomial)

glmIn = glm(formula = disFollow_group ~ as.numeric(ageGroup) +  tating + CRP + as.numeric(monthGroup)
            + LDL_C + dis_group+ HDL_C + APO_A +  APO_B, 
            data = dataAllDemoLabIntegrationPsm, family = binomial)



print(summary(glmIn))
tmp = summary(glmIn)
resOR[4] = formatC(exp(tmp$coefficients[3,1]),format = "f",digits = 2)
resORCi[4] = paste("(",
                   formatC(exp(tmp$coefficients[3,1]-1.96*tmp$coefficients[3,2]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[3,1]+1.96*tmp$coefficients[3,2]),
                           format = "f",digits = 2),")",sep = ""
)

tmp = as.data.frame(tmp$coefficients)


tmp$OR=formatC(exp(tmp$Estimate),format = "f",digits = 2)

tmp$`95%ci` = paste("(",formatC(exp(tmp$Estimate-1.96*tmp$`Std. Error`),format = "f",digits = 2),","
                    ,formatC(exp(tmp$Estimate+1.96*tmp$`Std. Error`),format = "f",digits = 2),")",sep="")

tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)
tmp[,3] = formatC(tmp[,3],format = "f",digits = 2)

rownames(tmp) = c("Intercept","Age","statin use","CRP","Follow-up(months)","LDL-C",
                  "CVD-related disease(FE)","HDL-C","APO-A","APO-B")
tmp = tmp[,c(5,6,3,4)]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")
# 
# write.csv(tmp,file = "followMultiOutcome.csv")
# 
# write.csv(tmp,file = "lpaCVD.csv")

nrow(dataAllDemoLabIntegrationPsm)

mitifactors =  data.frame(resOR,resORCi)

rownames(mitifactors) = c("Single statin","Single Dis","Muti statin","Muti Dis")

write.csv(mitifactors,file = "mitifactorsNormalData.csv")



#################################################WellComplianceDataSet
#=================================outPatients

nrow(dataOutDemoLabIntegration)
nrow(dataInDemoLabIntegration)
head(dataOutDemoLabIntegration)
head(dataInDemoLabIntegration)
7641+27944

nrow(dataOutDemoLabIntegrationPsm)
nrow(dataInDemoLabIntegrationPsm)
4888+ 16488


intersect(dataOutDemoLabIntegration$MASTER_INDEX,dataInDemoLabIntegration$MASTER_INDEX)

dataOutDemoLabIntegrationRaw = dataOutDemoLabIntegration

dataOutDemoLabIntegration = dataOutDemoLabIntegration[!dataOutDemoLabIntegration$MASTER_INDEX%in%intersect(dataOutDemoLabIntegration$MASTER_INDEX,dataInDemoLabIntegration$MASTER_INDEX),]
nrow(dataOutDemoLabIntegration)


dataOutDemoLabIntegration# outpatients dataset 

###################outpatients

dataOutDemoLabIntegration_1 = na.omit(dataOutDemoLabIntegration[dataOutDemoLabIntegration$LDL_C>1.8,])
nrow(dataOutDemoLabIntegration_1 )

dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1[!(dataOutDemoLabIntegration_1$tating==1&dataOutDemoLabIntegration_1$outcome.y==1),]

nrow(dataOutDemoLabIntegration_1 )
table(dataOutDemoLabIntegration_1$tating)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$gender==3,]
nrow(dataOutDemoLabIntegration_1 )

dataOutDemoLabIntegration_1 = 
  na.omit(dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$gender==3,])
breaks = c(0,45,65,150)
label = c(1,2,3)
dataOutDemoLabIntegration_1$ageGroup= cut(dataOutDemoLabIntegration_1$age,breaks = breaks,labels = label)
dataOutDemoLabIntegration_1 = na.omit(dataOutDemoLabIntegration_1)
nrow(dataOutDemoLabIntegration_1 )

dataOutDemoLabIntegration_1 = merge(dataOutDemoLabIntegration_1,followUP,by.x = "MASTER_INDEX",by.y = "MASTER_INDEX")
nrow(dataOutDemoLabIntegration_1 )

#dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month>6,]
#dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<60,]

summary(dataOutDemoLabIntegration_1 )
breaks = c(0,36,60,200)
label = c(1,2,3)
dataOutDemoLabIntegration_1$monthGroup= cut(dataOutDemoLabIntegration_1$month,breaks = breaks,labels = label)
table(dataOutDemoLabIntegration_1$tating)
table(dataOutDemoLabIntegration_1$monthGroup,dataOutDemoLabIntegration_1$tating)
nrow(dataOutDemoLabIntegration_1)
#dataOutDemoLabIntegration_1= dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<36,]


####################disease modify
dataOutDemoLabIntegration_1 = 
  na.omit(dataOutDemoLabIntegration_1[!dataOutDemoLabIntegration_1$dis_group==5,])
nrow(dataOutDemoLabIntegration_1 )
dataOutDemoLabIntegration_1$dis_group.y = 0
dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$dis_group%in%c(1:4),]$dis_group.y = 1
dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$dis_group%in%c(6),]$dis_group.y = 2
table(dataOutDemoLabIntegration_1$dis_group)
table(dataOutDemoLabIntegration_1$dis_group.y)

#table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]

max(table(dataOutDemoLabIntegration_1$MASTER_INDEX))

#dataAllDemoLabIntegration = 
#  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]


dataOutDemoLabIntegration_1C = dataOutDemoLabIntegration_1
nrow(dataOutDemoLabIntegration_1C)
table(dataOutDemoLabIntegration_1C$tating)



dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==1,]
dataOutDemoLabIntegrationPsm_1 = psmOut(dataOutDemoLabIntegration_1)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==2,]
dataOutDemoLabIntegrationPsm_2 = psmOut(dataOutDemoLabIntegration_1)
dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1C[dataOutDemoLabIntegration_1C$ageGroup==3,]
dataOutDemoLabIntegrationPsm_3 = psmOut(dataOutDemoLabIntegration_1)

dataOutDemoLabIntegrationPsm= rbind(
  dataOutDemoLabIntegrationPsm_1,
  dataOutDemoLabIntegrationPsm_2,
  dataOutDemoLabIntegrationPsm_3
)



colnames(dataOutDemoLabIntegrationPsm)[c(3,11:17)]
c(3,11:17)
rownames = c("Lp(a) at FT (Mean,ci95%)","LDL-C (Mean,ci95%)",
             "HDL-C (Mean,ci95%)","APO-A (Mean,ci95%)","APO-B (Mean,ci95%)",
             "TC (Mean,ci95%)","TG (Mean,ci95%)","CRP (Mean,ci95%)")
num = table(dataOutDemoLabIntegrationPsm$tating)

colnames = c(paste("Control","(","n=",num[1],")",sep = ""),
             paste("Treat","(","n=",num[2],")",sep = ""),
             "StaVal","PVal")

colnames(dataOutDemoLabIntegrationPsm)[c(3,11:17)]
tab_1 = data.frame()
for(i in c(3,11:17)){
  tmp = dataOutDemoLabIntegrationPsm[c(9,i)]
  tmp = as.data.frame(descri2(tmp)) # mean1+ci95%,mean2+ci95%,staticValue,Pvalue
  tab_1 = rbind(tab_1,tmp)
}

colnames(tab_1) = colnames
rownames(tab_1) = rownames
tab_1
colnames(dataOutDemoLabIntegrationPsm)[c(7,28,31,27,10)]
tab_2 = data.frame()
#gender, age, Disease in FE, Disease in FT,lp(a) outcome,month
for (i in c(7,28,31,27,10)){
  M = dataOutDemoLabIntegrationPsm[c(i,9)]
  tmp = decriCount(M)
  tab_2 = rbind(tab_2,tmp)
}
colnames(tab_2)=colnames(tab_1)
tab_2
write.csv(file ="outPatientsDemo.csv" , rbind(tab_1,tab_2))


####################################Inpatients


nrow((dataInDemoLabIntegration))
head(dataInDemoLabIntegration)
dataInDemoLabIntegration_1 = na.omit(dataInDemoLabIntegration[dataInDemoLabIntegration$LDL_C>1.8,])
nrow(dataInDemoLabIntegration_1 )
dataInDemoLabIntegration_1 = dataInDemoLabIntegration_1[!(dataInDemoLabIntegration_1$tating==1&dataInDemoLabIntegration_1$outcome.y==1),]
nrow(dataInDemoLabIntegration_1 )

table(table(dataInDemoLabIntegration_1$tating))
dataInDemoLabIntegration_1 = 
  na.omit(dataInDemoLabIntegration_1[!dataInDemoLabIntegration_1$gender==3,])
nrow(dataInDemoLabIntegration_1 )
breaks = c(0,45,60,65,70,80,1000)
label = c(1,2,3,4,5,6)
dataInDemoLabIntegration_1$ageGroup= cut(dataInDemoLabIntegration_1$age,breaks = breaks,labels = label)
nrow(dataInDemoLabIntegration_1 )
table(dataInDemoLabIntegration_1$ageGroup)
sum(table(dataInDemoLabIntegration_1$ageGroup))
dataInDemoLabIntegration_1 = na.omit(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1 = dataInDemoLabIntegration_1[order(dataInDemoLabIntegration_1$tating,decreasing = T),]
nrow(dataInDemoLabIntegration_1 )

table(followUP$MASTER_INDEX)
followUP[followUP$MASTER_INDEX == 2644,]


dataInDemoLabIntegration_1 = merge(dataInDemoLabIntegration_1,followUP,by.x = "MASTER_INDEX",by.y = "MASTER_INDEX",all.x = T)
nrow(dataInDemoLabIntegration_1)
# dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month>6,]
# dataOutDemoLabIntegration_1 =dataOutDemoLabIntegration_1[dataOutDemoLabIntegration_1$month<60,]

table(dataInDemoLabIntegration$MASTER_INDEX)[table(dataInDemoLabIntegration$MASTER_INDEX)>=2]
#dataAllDemoLabIntegration = 
#  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]


breaks = c(0,36,60,200)
label = c(1,2,3)
dataInDemoLabIntegration_1$monthGroup= cut(dataInDemoLabIntegration_1$month,breaks = breaks,labels = label)
table(dataInDemoLabIntegration_1$tating)
table(dataInDemoLabIntegration_1$monthGroup,dataInDemoLabIntegration_1$tating)




####################disease modify
dataInDemoLabIntegration_1 = 
  na.omit(dataInDemoLabIntegration_1[!dataInDemoLabIntegration_1$dis_group==5,])
nrow(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1$dis_group.y=1
dataInDemoLabIntegration_1[dataInDemoLabIntegration_1$dis_group%in%c(1:4),]$dis_group.y = 1
dataInDemoLabIntegration_1[dataInDemoLabIntegration_1$dis_group%in%c(6),]$dis_group.y = 2
table(dataInDemoLabIntegration_1$dis_group)

dataInDemoLabIntegration_1C = dataInDemoLabIntegration_1
nrow(dataInDemoLabIntegration_1C )

table(dataInDemoLabIntegration_1C$ageGroup)
table(dataInDemoLabIntegration_1C$tating)

#write.csv(dataInDemoLabIntegration_1,file="dataInDemoLabIntegration_1.csv")
#table(dataInDemoLabIntegration_1C$tating)

dataInDemoLabIntegration_1C = (dataInDemoLabIntegration_1C[!dataInDemoLabIntegration_1C$MASTER_INDEX%in%intersect(dataOutDemoLabIntegrationPsm$master_id,dataInDemoLabIntegration_1C$master_id),])
nrow(dataInDemoLabIntegration_1C)
head(dataInDemoLabIntegration_1C)


dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==1,]
dataInDemoLabIntegrationPsm_1 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==2,]
dataInDemoLabIntegrationPsm_2 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==3,]
dataInDemoLabIntegrationPsm_3 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==4,]
dataInDemoLabIntegrationPsm_4 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==5,]
dataInDemoLabIntegrationPsm_5 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegration_1  = dataInDemoLabIntegration_1C[dataInDemoLabIntegration_1C$ageGroup==6,]
dataInDemoLabIntegrationPsm_6 = PsmIn(dataInDemoLabIntegration_1)
dataInDemoLabIntegrationPsm = rbind(dataInDemoLabIntegrationPsm_1,
                                    dataInDemoLabIntegrationPsm_2,
                                    dataInDemoLabIntegrationPsm_3,
                                    dataInDemoLabIntegrationPsm_4,
                                    dataInDemoLabIntegrationPsm_5,
                                    dataInDemoLabIntegrationPsm_6
)
table(dataInDemoLabIntegrationPsm$tating)
table(dataOutDemoLabIntegrationPsm$tating)
breaks = c(0,45,65,150)
label = c(1,2,3)
dataInDemoLabIntegrationPsm$ageGroup= cut(dataInDemoLabIntegrationPsm$age,
                                          breaks = breaks,labels = label)
nrow(dataInDemoLabIntegrationPsm)
nrow(dataOutDemoLabIntegrationPsm)



dataAllDemoLabIntegrationPsm = rbind(dataOutDemoLabIntegrationPsm,dataInDemoLabIntegrationPsm)
nrow(dataAllDemoLabIntegrationPsm)
breaks = c(0,45,65,150)
label = c(1,2,3)
dataAllDemoLabIntegrationPsm$ageGroup= cut(dataAllDemoLabIntegrationPsm$age,
                                           breaks = breaks,labels = label)

table(dataOutDemoLabIntegrationPsm$master_id)[table(dataOutDemoLabIntegrationPsm$master_id)>=2]

table(dataInDemoLabIntegrationPsm$master_id)[table(dataInDemoLabIntegrationPsm$master_id)>=2]

intersect(dataInDemoLabIntegrationPsm$master_id,dataOutDemoLabIntegrationPsm$master_id)


table(dataAllDemoLabIntegrationPsm$master_id)[table(dataAllDemoLabIntegrationPsm$master_id)>=2]

#dataAllDemoLabIntegrationPsm = dataAllDemoLabIntegrationPsm[!dataAllDemoLabIntegrationPsm$MASTER_INDEX%in%names(table(CVDmodify$master_id)[table(CVDmodify$master_id)>=2]),]
nrow(dataAllDemoLabIntegrationPsm)
#CVDmodify = read.table("CVDmodified.txt",header = T,sep = "\t")


#dataAllDemoLabIntegrationPsm = merge(dataAllDemoLabIntegrationPsm,CVDmodify,
#                                     by.x = "MASTER_INDEX",by.y = "master_id",all.y  = F)
nrow(dataAllDemoLabIntegrationPsm)


dataAllDemoLabIntegration = rbind(dataOutDemoLabIntegration,dataInDemoLabIntegration)
nrow(dataAllDemoLabIntegration)

table(rbind(dataOutDemoLabIntegration_1C,dataInDemoLabIntegration_1C)$tating)

head(dataAllDemoLabIntegration)
colnames(dataAllDemoLabIntegration)
dataAllDemoLabIntegration = dataAllDemoLabIntegration[,c(1,24)]
head(dataAllDemoLabIntegration)
dataAllDemoLabIntegration = 
  dataAllDemoLabIntegration[!dataAllDemoLabIntegration$MASTER_INDEX%in%names(table(dataAllDemoLabIntegration$MASTER_INDEX)[table(dataAllDemoLabIntegration$MASTER_INDEX)>=2]),]
head(dataAllDemoLabIntegration)
head(dataAllDemoLabIntegrationPsm)

dataAllDemoLabIntegrationPsm_1=merge(dataAllDemoLabIntegrationPsm,dataAllDemoLabIntegration,
                                     by.x = "MASTER_INDEX",by.y = "MASTER_INDEX",all.x = F)
nrow(dataAllDemoLabIntegrationPsm_1)

dataAllDemoLabIntegrationPsm = dataAllDemoLabIntegrationPsm_1



rownames = c("Lp(a) at FT (Mean,CI95%)","LDL-C at FT (Mean,CI95%)","LDL-C (Mean,CI95%)",
             "HDL-C (Mean,CI95%)","APO-A (Mean,CI95%)","APO-B (Mean,CI95%)",
             "TC (Mean,CI95%)","TG (Mean,CI95%)","CRP (Mean,CI95%)")
num = table(dataAllDemoLabIntegrationPsm$tating)

colnames = c(paste("Control","(","n=",num[1],")",sep = ""),
             paste("Treat","(","n=",num[2],")",sep = ""),
             "StaVal","PVal")
colnames(dataAllDemoLabIntegrationPsm)
tab_1 = data.frame()
for(i in c(3,12,15:21)){
  tmp = dataAllDemoLabIntegrationPsm[c(9,i)]
  tmp = as.data.frame(descri2(tmp)) # mean1+ci95%,mean2+ci95%,staticValue,Pvalue
  tab_1 = rbind(tab_1,tmp)
}

colnames(tab_1) = colnames
rownames(tab_1) = rownames
tab_1
dataAllDemoLabIntegrationPsm$CRP
tab_2 = data.frame()
colnames(dataAllDemoLabIntegrationPsm)[ c(7,28,31,30,10,14,27,25)]
#gender,age,dis_history,follow_time,LPA_outcome,LDL_outcome,dis_follow,
for (i in c(7,28,31,30,10,14,27,25)){
  M = dataAllDemoLabIntegrationPsm[c(i,9)]
  tmp = decriCount(M)
  tab_2 = rbind(tab_2,tmp)
}
tab_2
colnames(tab_2)=colnames(tab_1)
nrow(tab_2)
rownames(tab_2) = c("Male","Female","<45","46-65",">65",
                    "Non-CVD-related disease(FE)","CVD-related disease(FE)",
                    "[0.5 - 3) years","[3 - 5) years","≥ 5 years",
                    "Lp(a) decrease (FU)","Lp(a) increase (FU)",
                    "LDL-C decrease (FU)","LDL-C increase (FU)",
                    "Non-CVD disease event (FU)","CVD disease event (FU)",
                    #"Modify Non-CVD disease event (FU)","Modify CVD disease event (FU)",
                    "diabetes","dyslipidemia","hypertension","arteriosclerosis","Others"
                    #"Lp(a) decrease_1 (FU)","Lp(a) increase_1 (FU)",
                    #"LDL-C decrease_1 (FU)","LDL-C increase_1 (FU)"
)
tab_2


rbind(tab_1,tab_2)


write.csv(file = "WellComDataset.csv",rbind(tab_1,tab_2))

write.csv(file = "WellComDatasetRaw.csv",dataAllDemoLabIntegrationPsm)

dat = read.csv("WellComDatasetRaw.csv")

dataAllDemoLabIntegrationPsm = dat[,-1]

###############################

resOR=c()
resORCi=c()

colnames(dataAllDemoLabIntegrationPsm)

nrow( dataAllDemoLabIntegrationPsm)

glmIn = coxph(Surv(month,outcome.x) ~  tating,
              data = dataAllDemoLabIntegrationPsm)
print(summary(glmIn))

tmp = summary(glmIn)

resOR[1] = formatC(tmp$coefficients[1,2],format = "f",digits = 2)
resORCi[1] = paste("(",
                   formatC(exp(tmp$coefficients[1,1]-1.96*tmp$coefficients[1,3]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[1,1]+1.96*tmp$coefficients[1,3]),
                           format = "f",digits = 2),")",sep = ""
)

# glmIn = glm(formula = disFollow_group_new  ~   tating , 
#             data = dataAllDemoLabIntegrationPsm, family = binomial)

glmIn = glm(formula = disFollow_group  ~   tating , 
            data = dataAllDemoLabIntegrationPsm, family = binomial)
print(summary(glmIn))
tmp = summary(glmIn)
resOR[2] = formatC(exp(tmp$coefficients[2,1]),format = "f",digits = 2)
resORCi[2] = paste("(",
                   formatC(exp(tmp$coefficients[2,1]-1.96*tmp$coefficients[2,2]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[2,1]+1.96*tmp$coefficients[2,2]),
                           format = "f",digits = 2),")",sep = ""
)
resOR
resORCi

glmIn = coxph(Surv(month,outcome.x) ~  as.numeric(ageGroup) +  tating +CRP + LDL_C + dis_group + HDL_C + 
                APO_A +  APO_B,
              data = dataAllDemoLabIntegrationPsm)

print(summary(glmIn))

tmp = summary(glmIn)
resOR[3] = formatC(exp(tmp$coefficients[2,1]),format = "f",digits = 2)
resORCi[3] = paste("(",
                   formatC(exp(tmp$coefficients[2,1]-1.96*tmp$coefficients[2,3]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[2,1]+1.96*tmp$coefficients[2,3]),
                           format = "f",digits = 2),")",sep = ""
)
resOR
resORCi
tmp = as.data.frame(tmp$coefficients)


tmp$OR=formatC(tmp$`exp(coef)`,format = "f",digits = 2)

tmp$`95%CI` = paste("(",formatC((tmp$`exp(coef)`-1.96*tmp$`se(coef)`),format = "f",digits = 2),","
                    ,formatC((tmp$`exp(coef)`+1.96*tmp$`se(coef)`),format = "f",digits = 2),")",sep="")
tmp
tmp[,5] = formatC(tmp[,5],format = "f",digits = 2)
tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)

rownames(tmp) = c("Age","statin use","CRP","LDL-C",
                  "CVD-related disease(FE)","HDL-C","APO-A","APO-B")
tmp = tmp[,c(6,7,4,5)]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")

write.csv(tmp,file = "statinMultiOutcome.csv")


# glmIn = glm(formula = disFollow_group_new ~ as.numeric(ageGroup) +  tating + CRP + as.numeric(monthGroup)
#             + LDL_C + dis_group.x + HDL_C + APO_A +  APO_B, 
#             data = dataAllDemoLabIntegrationPsm, family = binomial)

glmIn = glm(formula = disFollow_group ~ as.numeric(ageGroup) +  tating + CRP + as.numeric(monthGroup)
            + LDL_C + dis_group+ HDL_C + APO_A +  APO_B, 
            data = dataAllDemoLabIntegrationPsm, family = binomial)



print(summary(glmIn))
tmp = summary(glmIn)
resOR[4] = formatC(exp(tmp$coefficients[3,1]),format = "f",digits = 2)
resORCi[4] = paste("(",
                   formatC(exp(tmp$coefficients[3,1]-1.96*tmp$coefficients[3,2]),
                           format = "f",digits = 2),",",
                   formatC(exp(tmp$coefficients[3,1]+1.96*tmp$coefficients[3,2]),
                           format = "f",digits = 2),")",sep = ""
)

tmp = as.data.frame(tmp$coefficients)


tmp$OR=formatC(exp(tmp$Estimate),format = "f",digits = 2)

tmp$`95%ci` = paste("(",formatC(exp(tmp$Estimate-1.96*tmp$`Std. Error`),format = "f",digits = 2),","
                    ,formatC(exp(tmp$Estimate+1.96*tmp$`Std. Error`),format = "f",digits = 2),")",sep="")

tmp[,4] = formatC(tmp[,4],format = "f",digits = 2)
tmp[,3] = formatC(tmp[,3],format = "f",digits = 2)

rownames(tmp) = c("Intercept","Age","statin use","CRP","Follow-up(months)","LDL-C",
                  "CVD-related disease(FE)","HDL-C","APO-A","APO-B")
tmp = tmp[,c(5,6,3,4)]

colnames(tmp)[c(3,4)] = c("Zvalue","Pvalue")
# 
# write.csv(tmp,file = "followMultiOutcome.csv")
# 
# write.csv(tmp,file = "lpaCVD.csv")

nrow(dataAllDemoLabIntegrationPsm)

mitifactors =  data.frame(resOR,resORCi)

rownames(mitifactors) = c("Single statin","Single Dis","Muti statin","Muti Dis")

write.csv(mitifactors,file = "mitifactorsWellComData.csv")



breaks = c(0,2.300, 3.540,13.230)
label = c(3,2,1)
dataAllDemoLabIntegrationPsm$LDL_CGroup= cut(dataAllDemoLabIntegrationPsm$LDL_C,breaks = breaks,labels = label)
dataAllDemoLabIntegrationPsm$LDL_CGroup



tabdrug = table(dataAllDemoLabIntegrationPsm$LDL_CGroup,dataAllDemoLabIntegrationPsm$disFollow_group)

tabdrug[,2]/tabdrug[,1]


fivenum(dataAllDemoLabIntegrationPsm$lastValue)

breaks = c(0, 179.00, 322.00,10000)
label = c(1,2,3)
dataAllDemoLabIntegrationPsm$LPAFlGroup= cut(as.numeric(dataAllDemoLabIntegrationPsm$median.x)+1,breaks = breaks,labels = label)
dataAllDemoLabIntegrationPsm$lastValue
dataAllDemoLabIntegrationPsm$LPAFlGroup

cbind(dataAllDemoLabIntegrationPsm$lastValue,dataAllDemoLabIntegrationPsm$LPAFlGroup)

breaks = c(0, 65, 150)
label = c(1,2)
dataAllDemoLabIntegrationPsm$AgeGroup2= cut(dataAllDemoLabIntegrationPsm$age,breaks = breaks,labels = label)



colnames( dataAllDemoLabIntegrationPsm)[c(9,10,27,15:21,7,14,29,35,31,33,34)]
# [1] "tating"          "outcome.x"       "disFollow_group" "LDL_C"           "HDL_C"          
# [6] "APO_A"           "APO_B"           "TC"              "TG"              "CRP"            
# [11] "gender"          "outcome.y"       "dis_group.y"     "disFollow_group" "monthGroup"     
# [16] "LPAFlGroup"      "AgeGroup2" 

subcovarience = dataAllDemoLabIntegrationPsm[c(9,10,27,15:21,7,14,29,35,31,33,34)]
subcovarienceMain = subcovarience
head(subcovarience )
subcovarience$dis_group.y
subcovarience [,12] = subcovarience [,12]+1

subcovarience [,14]

colnames(subcovarience)

metaFor  


metaFor = subgroup(subcovarience,15,"subgrouplpaDis",
               c("Non-CVD disease history","CVD disease history"),"outcome.x")


tmp = subgroup(subcovarience,15,"subgroupFolDisDis",
               c("Non-CVD disease history","CVD disease history"),"disFollow_group")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])
metaFor


tmp = subgroup(subcovarience,16,"subgrouplpaLDL_C",
               c("≥3.54mmol/L", "(2.30,3.54)mmol/L","≤2.30mmol/L"),"outcome.x")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])

tmp = subgroup(subcovarience,16,"subgroupFolDisLDL_C",
               c("≥3.54mmol/L", "(2.30,3.54)mmol/L","≤2.30mmol/L"),"disFollow_group")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])


tmp = subgroup(subcovarience,17,"subgrouplpaLLPAFL",
               c("≤179.00mg/L", "(179.00,322.00)mg/L","≥322.00mg/L"),"outcome.x")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])
tmp = subgroup(subcovarience,17,"subgroupFolDisLPAFL",
               c("≤179.00mg/L", "(179.00,322.00)mg/L","≥322.00mg/L"),"disFollow_group")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])


# tmp =subgroup(subcovarience,13,"subgrouplpaLMonth",
#               c("(0.5,3] years", "(3,5) year","≥5 year"),"outcome.x")
# metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
# metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])
# tmp =subgroup(subcovarience,13,"subgroupFolDisLMonth",
#               c("(0.5,3] years", "(3,5) year","≥5 year"),"disFollow_group")
# metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
# metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])

tmp = subgroup(subcovarience,11,"subgrouplpaSex",
               c("Male", "Female"),"outcome.x")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])
tmp = subgroup(subcovarience,11,"subgroupFolDisSex",
               c("Male", "Female"),"disFollow_group")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])


tmp = subgroup(subcovarience,14,"subgrouplpaSex",
               c("<65", "≥65"),"outcome.x")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])
tmp = subgroup(subcovarience,14,"subgroupFolDisSex",
               c("<65", "≥65"),"disFollow_group_new")
metaFor[[1]]=rbind(metaFor[[1]],tmp[[1]][c(-1,-2),])
metaFor[[2]]=rbind(metaFor[[2]],tmp[[2]])

metaFor

#lpasubgroup = metaFor[[1]][c(1,2,3,4,7,8,11:13,17:19,23:25,29:31,35,36),]

#lpasubgroup = metaFor[[1]][c(1,2,3,4,7,8,11:13,17:19,23:25,29:30,33,34),]

nrow( metaFor[[1]])
lpasubgroup = metaFor[[1]][c(1,2,23,24,19,20,3,4,7:9,13:15),]


metaFor[[1]][,1]

Item = c("\n","\n",
         "Age","\n",
         "Sex","\n",
         "CVD-related disease history","\n",
         "LDL-C level","\n","\n",
         "Lp(a) level (FE)","\n","\n"
)
length(Item)


# Item = c("\n","\n","LDL-C outcome","\n",
#          "CVD-related disease history","\n",
#          "Age","\n","\n",
#          "LDL-C level","\n","\n",
#          "Lp(a) level","\n","\n",
#          "Follow-up time","\n","\n",
#          "Sex","\n","\n")

lpasubgroup = cbind(Item,lpasubgroup)

lpasubgroup[2,3:6] = c("Lp(a) increase","Lp(a) decrease","Lp(a) increase","Lp(a) decrease")

metaFor[[1]][,1]
# #lpasubgroup = metaFor[[1]][c(1,2,33,34,29,30,7,8,23:25,11:13,17:19,3,4),]
# 
# CVDsubgroup = metaFor[[1]][c(1,2,35,36,31,32,9,10,26:28,14:16,20:22,5,6),]
# 
# CVDsubgroup  = cbind(Item,CVDsubgroup)
# 
# CVDsubgroup[2,3:6] = c("CVD event","Non-CVD event","CVD event","Non-CVD event")
# 

tiff(filename =  "subgroup-lpaInRishWell.tif",
     width = 9480, height = 3480, units = "px", pointsize = 12,
     compression ="lzw",
     bg = "white", res = 300)

forestplot(lpasubgroup, 
           txt_gp=fpTxtGp(
             label=gpar(cex=2),
             ticks=gpar(cex=2),
             xlab=gpar(cex = 1.7),
             title=gpar(cex = 1.8)),
           rbind( rep(NA,3),rep(NA,3),
                  metaFor[[2]][c(23,24,19,20,3,4,7:9,13:15)-2,]),
           col=clrs,
           colgap=unit(8,"mm"),
           lwd.ci=2, boxsize=0.5,
           zero=0.5,
           cex=2, lineheight = "auto",
           ci.vertices=TRUE, ci.vertices.height = 0.4
           #xlab="    <---Decrease--- lp(a) ---Increase--->"
)
dev.off()





#Function 
psmOut = function(dataOutDemoLabIntegration_1){
  dataOutDemoLabIntegration_1 = dataOutDemoLabIntegration_1[order(dataOutDemoLabIntegration_1$tating,decreasing = T),]
  #attach(dataOutDemoLabIntegration_1)
  dataOutDemoLabIntegration_1
  t_ind = dataOutDemoLabIntegration_1$tating
  dist_mat = NULL
  subset_weight = 1
  mom_covs = cbind(
    dataOutDemoLabIntegration_1$firstValue.x,
    dataOutDemoLabIntegration_1$LDL_C
    # dataOutDemoLabIntegration_1$HDL_C,
    # dataOutDemoLabIntegration_1$APO_A,
    # dataOutDemoLabIntegration_1$APO_B
  ) 
  mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
  mom = list(covs = mom_covs, tols = mom_tols)
  
  fine_covs = cbind(
    dataOutDemoLabIntegration_1$gender,
    dataOutDemoLabIntegration_1$dis_group.y,
    dataOutDemoLabIntegration_1$monthGroup
  )
  fine = list(covs = fine_covs)
  
  t_max = 60*5 
  solver = "glpk" 
  approximate = 1 
  solver = list(name = solver, t_max = t_max, approximate = approximate, round_cplex = 0, trace = 0)
  out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight, 
               mom = mom, fine = fine,solver = solver)
  
  t_id = out$t_id
  length(t_id)
  c_id = out$c_id
  length(c_id)
  intersect(t_id,c_id)
  length(out$group_id)
  meantab(mom_covs, t_ind, t_id, c_id)
  nrow(dataOutDemoLabIntegration_1)
  for (i in 1:ncol(fine_covs)) {
    print(finetab(fine_covs[, i], t_id, c_id))
  }
  nrow(mom_covs)
  intersect(t_id,c_id)
  #dataOutDemoLabIntegrationPsm = dataOutDemoLabIntegration_1[c(t_id, c_id[t_ind[c_id]!=1]),]
  dataOutDemoLabIntegrationPsm = dataOutDemoLabIntegration_1[c(t_id, c_id),]
  return(dataOutDemoLabIntegrationPsm)
}

PsmIn  <- function(dataInDemoLabIntegration_1) {
  dataInDemoLabIntegration_1 = dataInDemoLabIntegration_1[order(dataInDemoLabIntegration_1$tating,decreasing = T),]
  
  t_ind = dataInDemoLabIntegration_1$tating
  #write.csv(file = "id.csv",t_ind)
  dist_mat = NULL
  subset_weight = 1
  mom_covs = cbind(
    dataInDemoLabIntegration_1$firstValue.x,
    dataInDemoLabIntegration_1$LDL_C
  ) 
  #write.csv(mom_covs,file = "psmIn.csv")
  mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
  mom = list(covs = mom_covs, tols = mom_tols)
  
  fine_covs = cbind(
    dataInDemoLabIntegration_1$gender,
    dataInDemoLabIntegration_1$dis_group,
    dataInDemoLabIntegration_1$monthGroup
  )
  fine = list(covs = fine_covs)
  
  t_max = 60*5 
  solver = "glpk" 
  approximate = 1 
  solver = list(name = solver, t_max = t_max, approximate = approximate, round_cplex = 0, trace = 0)
  out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight, 
               mom = mom, fine = fine,solver = solver)
  
  t_id = out$t_id
  print(length(t_id))
  c_id = out$c_id
  print(length(c_id))
  print(intersect(t_id,c_id))
  print(meantab(mom_covs, t_ind, t_id, c_id))
  for (i in 1:ncol(fine_covs)) {
    print(finetab(fine_covs[, i], t_id, c_id))
  }
  dataInDemoLabIntegrationPsm = dataInDemoLabIntegration_1[c(t_id, c_id),]
  return(dataInDemoLabIntegrationPsm)
}


decriCount = function(M){
  #M = dataAllDemoLabIntegrationPsm[c(25,9)]
  TM = table(M)
  TM
  Mpertage=TM
  Mpertage
  for (i in 1:ncol(TM)) {
    Mpertage[,i] = c(TM[,i]/sum(TM[,i]))
  }
  M_1 = ( matrix( paste(as.matrix(TM),"(",formatC(as.matrix(Mpertage)*100,format = "f",digits = 2),"%",")",sep = ""),
                  nrow = nrow(TM),dimnames = dimnames(TM)))
  Q = chisq.test(TM)
  res=cbind(M_1,Q$statistic,Q$p.value)
  return(res)
}


descri2 = function(tmp){
  #tmp = tmp 
  t = t.test(tmp[,2]~tmp[,1],data=tmp)
  m_ctrl = formatC(c(mean(tmp[tmp[,1]==0,2]),
                     quantile(tmp[tmp[,1]==0,2], c(0.025, 0.975))),
                   format = "f",digits = 2)
  r_1 =as.character(paste( m_ctrl [1],"(", m_ctrl [2],", ", m_ctrl[3],")",sep = ""))
  m_treat = formatC(c(mean(tmp[tmp[,1]==1,2]),
                      quantile(tmp[tmp[,1]==1,2], c(0.025, 0.975))),
                    format = "f",digits = 2)
  r_2 =as.character(paste( m_treat[1],"(",m_treat [2],", ",m_treat[3],")",sep = ""))
  staVal = as.character(formatC(abs(t$statistic),format = "f",digits = 2))
  pVal = as.character(formatC(abs(t$p.value),format = "f",digits = 2))
  res = data.frame(r_1,r_2,staVal,pVal)
  return(res)
}


subgroup = function(subcovarience,index,fileName, name,outcomeIndex){
  #index=12 ; index=16； index=14
  #name = c("Non-CVD disease history","CVD disease history")
  # name =c("(0.5,3] years", "(3,5) year","≥5 year")
  # fileName = "subgroupDis"
  #outcomeIndex ="disFollow_group"
  #outcomeIndex = "outcome.x"
  #name= c("LDL-C decrease","LDL-C increase")
  head(subcovarience)
  {
    if(outcomeIndex=="outcome.x") 
      subcovarience_1 = subcovarience[,c(-index,-3)] else
        subcovarience_1 = subcovarience[,c(-index,-2)]
  }
  
  head(subcovarience_1)
  mean = c()
  lower = c()
  upper=c()
  statinIncP= c()
  statinInc=c()
  statinDec=c()
  nonStatinIncP = c()
  nonStatinInc = c()
  nonStatinDec=c()
  for (i in as.numeric(levels(as.factor(subcovarience[,index])))){
    
    tmpData =subcovarience_1[subcovarience[,index]==i,]
    head(tmpData)
    tab = table(tmpData[,2],tmpData$tating) # the first = row, the second = colomn
    statinIncP[i]= tab[2,2]/sum(tab[,2])
    nonStatinIncP[i] = tab[2,1]/sum(tab[,1])
    statinInc[i] = tab[2,2]
    nonStatinInc[i] = tab[2,1]
    statinDec[i] = tab[1,2]
    nonStatinDec[i] =tab[1,1]
    glmOut = coxph(Surv(month,tmpData[,2]) ~ ., 
                 data =tmpData[,-2])
    tmp = summary(glmOut)
    mean[i]=tmp$coefficients["tating",][2]
    lower[i] = tmp$coefficients["tating",][2]-1.96*tmp$coefficients["tating",][3]
    upper[i] = tmp$coefficients["tating",][2]+1.96*tmp$coefficients["tating",][3]
  }
  
  
  datFrame = as.data.frame(cbind(mean,lower,upper))
  datFrame
  rownames(datFrame) =name
  
  colnames(datFrame) = c("coef","lower","upper")
  
  clrs <- fpColors(box="royalblue",line="darkblue", summary="royalblue")
  
  tabletext <- cbind(c("Subgroup","\n",rownames(datFrame)),
                     c("No statin-used","Case increase number",paste(nonStatinInc,"(",sprintf("%.2f",nonStatinIncP*100),"%",")",sep = "")),
                     c("\n","Case decrease number",paste(nonStatinDec,"(",sprintf("%.2f",(1-nonStatinIncP)*100),"%",")",sep = "")),
                     c("Statin-used","Case increase number",paste(statinInc,"(",sprintf("%.2f",statinIncP*100),"%",")",sep = "")),
                     c("\n","Case decrease number",paste(statinDec,"(",sprintf("%.2f",(1-statinIncP)*100),"%",")",sep = "")),
                     c("HR","\n",sprintf("%.2f", datFrame[,"coef"])))
  filename = paste(fileName, gsub(pattern = ":|-| ",replacement = "",Sys.time()),"lp(a)",".tif",sep = "_")
  # tiff(filename =  filename,
  #      width = 5480, height = 3480, units = "px", pointsize = 12,
  #      compression ="lzw",
  #      bg = "white", res = 300)
  # # pdf(file  =  filename,
  # #      width = 5480, height = 3480, pointsize = 12) 
  # #xlab = 
  # forestplot(tabletext, 
  #            txt_gp=fpTxtGp(
  #              label=gpar(cex=2),
  #              ticks=gpar(cex=2),
  #              xlab=gpar(cex = 1.7),
  #              title=gpar(cex = 1.8)),
  #            rbind( rep(NA,3),rep(NA,3),
  #                   datFrame),
  #            col=clrs,
  #            colgap=unit(8,"mm"),
  #            lwd.ci=2, boxsize=0.5,
  #            zero=1,
  #            cex=2, lineheight = "auto",
  #            ci.vertices=TRUE, ci.vertices.height = 0.4
  #            #xlab="    <---Decrease--- lp(a) ---Increase--->"
  # )
  # dev.off()
  return(list(tabletext,datFrame))
  #return(head( head(subcovarience_1)))
}

#


