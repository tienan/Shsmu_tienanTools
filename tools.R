#tools
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization
library(designmatch)
library(dplyr)
library(metafor)
library(forestplot)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
#install.packages("factoextra")
library(factoextra) # clustering algorithms & visualization
library(designmatch)
library(metafor)
library(forestplot)
trendAnalysis = function(rawData,sign){
#rawData = df_lpa_tating
  df_lpa = rawData
  df_lpa$sex = as.character(df_lpa$sex)
  #男1女2其他
  df_lpa[grepl(pattern = "男",x=df_lpa$sex),]$sex = "1"
  df_lpa[grepl(pattern = "女",x=df_lpa$sex),]$sex = "2"
  df_lpa[!df_lpa$sex%in%c(1,2),]$sex = "3"
  df_lpa$RESULTS=as.numeric(as.character(df_lpa$RESULTS))
  df_lpa = na.omit(df_lpa)
  df_lpa$order =c(1:nrow(df_lpa))
  by_master_index <- group_by(df_lpa, MASTER_INDEX)
  nrow(by_master_index)
  models <- by_master_index %>% do(mod = lm( RESULTS~order+age, data = .))
  tatingSign = sign(as.data.frame(summarise(models, rsq = as.data.frame(summary(mod)$coefficients)[2,1])))
  table(tatingSign)
  #tatingSign
  #-1    1 
  #955 2166 
  #1306 1815 年龄调整后
  #1851/1306
  basicData <- by_master_index %>% 
    summarise(median= median(RESULTS),firstValue=first(RESULTS),lastValue= last(RESULTS),
              age=min(age),ORG_CODE=first(ORG_CODE),gender=first(sex),PROVINCE=last(PROVINCE))
  res = basicData
  res$tating = sign
  res$outcome = tatingSign
  res =as.data.frame(res)
  return(res)
}

diseaseAnalysis = function(){
##out patients  sustain 
  data = tbl_df(
    sqlQuery(cn,"
Select distinct a.MASTER_INDEX,a.ICD,a.DISEASE_NAME 
from Demo_Drug_ICD_Disease_OUT a  join
             (Select distinct  MASTER_INDEX From TATING_OUT ) b on a.MASTER_INDEX=b.MASTER_INDEX
             order by  a.MASTER_INDEX
             ") )
  OUT_tating = data
head(data)


##out patients  non sustain   
  data = tbl_df(
    sqlQuery(cn,"
Select distinct a.MASTER_INDEX,a.ICD,a.DISEASE_NAME 
from Demo_Drug_ICD_Disease_OUT a  join
             (Select distinct  MASTER_INDEX From NON_TATING_OUT ) b on a.MASTER_INDEX=b.MASTER_INDEX
             order by  a.MASTER_INDEX
             ") )
head(data)
OUT_non_tating = data
 
##IN patients  sustain 
  data = tbl_df(
    sqlQuery(cn,"
             Select distinct a.MASTER_INDEX,a.ICD,a.DISEASE_NAME 
             from Demo_Drug_ICD_Disease_IN a  join
             (Select distinct  MASTER_INDEX From TATING_IN ) b on a.MASTER_INDEX=b.MASTER_INDEX
             order by  a.MASTER_INDEX
             ") )
  head(data)
  IN_tating = data
 
  
 
  data = tbl_df(
    sqlQuery(cn,"
  Select distinct a.MASTER_INDEX,a.ICD,a.DISEASE_NAME 
  from Demo_Drug_ICD_Disease_IN a  join
  (Select distinct  MASTER_INDEX From NON_TATING_IN ) b on a.MASTER_INDEX=b.MASTER_INDEX
  order by  a.MASTER_INDEX
            ") )
  head(data)
  IN_non_tating = data
  
  intersect(IN_tating$MASTER_INDEX,IN_non_tating$MASTER_INDEX)
  intersect(OUT_tating$MASTER_INDEX,OUT_non_tating$MASTER_INDEX)
  
  
  write.csv(x=DiseaseSep(OUT_tating),  file = paste("tating_disease_OUT",
                                                gsub(pattern = ":|-| ",replacement = "",Sys.time()),"csv",sep = "-"))
  write.csv(x=DiseaseSep(OUT_non_tating),  file = paste("non_tating_disease_OUT",
                                                gsub(pattern = ":|-| ",replacement = "",Sys.time()),"csv",sep = "-"))
  write.csv(x=DiseaseSep(IN_tating),  file = paste("tating_disease_IN",
                                                gsub(pattern = ":|-| ",replacement = "",Sys.time()),"csv",sep = "-"))
  write.csv(x=DiseaseSep(IN_non_tating),  file = paste("non_tating_disease_IN",
                                                gsub(pattern = ":|-| ",replacement = "",Sys.time()),"csv",sep = "-"))
  
  tating_disease_OUT = read.csv("tating_disease_OUT(CHECK).csv")
  non_tating_disease_OUT = read.csv("non_tating_disease_OUT(CHECK).csv")
  tating_disease_IN = read.csv("tating_disease_IN(CHECK).csv")
  non_tating_disease_IN = read.csv("non_tating_disease_IN(CHECK).csv")
  
  #tmp = DiseaseSep(data[c(1:10000),])

}

disChisquare = function(){
  #dat_1 = tmp$dis_group
  #dat_2 = tmp$dis_group
  tating_disease_OUT = read.csv("tating_disease_OUT(CHECK).csv")
  non_tating_disease_OUT = read.csv("non_tating_disease_OUT(CHECK).csv")
  tating_disease_IN = read.csv("tating_disease_IN(CHECK).csv")
  non_tating_disease_IN = read.csv("non_tating_disease_IN(CHECK).csv")
  dat_1 = tating_disease_OUT$dis_group
  dat_2 = non_tating_disease_OUT$dis_group
  tmp_1 = table( dat_1)
  tmp_2 = table( dat_2)
  T <- as.table( rbind(tmp_1,tmp_2))
  dimnames(T) <- list(Sustain = c("Used", "Non-used"),
                      Dis = c("Diabetes","Dyslipidemia","Hypertension ","Pre-ASCVD","CVD","Other"))
  Tpertage = t(data.frame(T[1,]/sum( T[1,])))
  Tpertage = as.data.frame(rbind(Tpertage,t(data.frame(T[2,]/sum( T[2,])))))
  rownames(Tpertage) =  c("Used", "Non-used")
  Tpertage
  Tsq <- chisq.test(T)
  T_1 = t( matrix( paste(as.matrix(T),"(",formatC(as.matrix(Tpertage)*100,format = "f",digits = 2),"%",")",sep = ""),
                   nrow = nrow(T),dimnames = dimnames(T)))
  write.csv(x=T_1,file = "dis_OUT.csv")
  
  dat_1 = tating_disease_OUT$disFollow_group
  dat_2 = non_tating_disease_OUT$disFollow_group
  tmp_1 = table( dat_1)
  tmp_2 = table( dat_2)
  T <- as.table( rbind(tmp_1,tmp_2))
  dimnames(T) <- list(Sustain = c("Used", "Non-used"),
                      Dis = c("Non-ASCVD","ASCVD"))
  Tpertage = t(data.frame(T[1,]/sum( T[1,])))
  Tpertage = as.data.frame(rbind(Tpertage,t(data.frame(T[2,]/sum( T[2,])))))
  rownames(Tpertage) =  c("Used", "Non-used")
  Tpertage
  Tsq <- chisq.test(T)
  T_1_1 = t( matrix( paste(as.matrix(T),"(",formatC(as.matrix(Tpertage)*100,format = "f",digits = 2),"%",")",sep = ""),
                   nrow = nrow(T),dimnames = dimnames(T)))
  write.csv(x=T_1_1,file = "Outcome_OUT.csv")
  
  dat_1 = tating_disease_IN$dis_group
  dat_2 = non_tating_disease_IN$dis_group
  tmp_1 = table( dat_1)
  tmp_2 = table( dat_2)
  T <- as.table( rbind(tmp_1,tmp_2))
  dimnames(T) <- list(Sustain = c("Used", "Non-used"),
                      Dis = c("Diabetes","Dyslipidemia","Hypertension ","Pre-ASCVD","CVD","Other"))
  Tpertage = t(data.frame(T[1,]/sum( T[1,])))
  Tpertage = as.data.frame(rbind(Tpertage,t(data.frame(T[2,]/sum( T[2,])))))
  rownames(Tpertage) =  c("Used", "Non-used")
  Tpertage
  Tsq <- chisq.test(T)
  T_2 = t( matrix( paste(as.matrix(T),"(",formatC(as.matrix(Tpertage)*100,format = "f",digits = 2),"%",")",sep = ""),
                   nrow = nrow(T),dimnames = dimnames(T)))
  write.csv(x=T_2,file = "dis_IN.csv")
  
  dat_1 = tating_disease_IN$disFollow_group
  dat_2 = non_tating_disease_IN$disFollow_group
  tmp_1 = table( dat_1)
  tmp_2 = table( dat_2)
  T <- as.table( rbind(tmp_1,tmp_2))
  dimnames(T) <- list(Sustain = c("Used", "Non-used"),
                      Dis = c("Non-ASCVD","ASCVD"))
  Tpertage = t(data.frame(T[1,]/sum( T[1,])))
  Tpertage = as.data.frame(rbind(Tpertage,t(data.frame(T[2,]/sum( T[2,])))))
  rownames(Tpertage) =  c("Used", "Non-used")
  Tpertage
  Tsq <- chisq.test(T)
  T_2_1 = t( matrix( paste(as.matrix(T),"(",formatC(as.matrix(Tpertage)*100,format = "f",digits = 2),"%",")",sep = ""),
                     nrow = nrow(T),dimnames = dimnames(T)))
  write.csv(x=T_2_1,file = "Outcome_IN.csv")
  
  
}

DiseaseSep = function(data){
#  library(dplyr)
  dis_data = data
  
  by_master_index <- group_by(dis_data, MASTER_INDEX)
 
  Dis <- by_master_index %>% 
    summarise(master_id = first(MASTER_INDEX),history= first(DISEASE_NAME)
    )
  Dis = Dis %>%  mutate(dis_group=case_when(
    grepl(pattern =  "糖尿病|消渴",Dis$history) ~ 1,
    grepl(pattern = "血脂|脂血",Dis$history) ~  2,
    grepl(pattern ="高血压",Dis$history) ~ 3,
    grepl(pattern = "((血管)|(动脉)&(硬化))",Dis$history) ~ 4 ,
    grepl(pattern = "
((脑)&(出血))|((脑)&(梗死))|((脑)&(梗塞))|((心)&(梗))|((脑)&(供血))|冠心病|((动脉粥样)&(心脏病))|(脑血管)|((心)&(梗死))|((心)&(梗塞))|(脑血管)
",Dis$history) ~ 5,
    TRUE ~ 6
  ))
  
   table(Dis$dis_group)
  
  disFollow=c()
  for (i in 1:nrow(Dis)){
    tmp = as.character(dis_data[dis_data$MASTER_INDEX==Dis$MASTER_INDEX[i],]$DISEASE_NAME)
    tmp = paste0(tmp[-1],collapse = ",")#去除第一的诊断信息
    tmp
    disFollow[i]=tmp
  }
  Dis$disFollow = disFollow
  Dis = Dis %>%  mutate(disFollow_group=case_when(
    grepl(pattern = "
((脑)&(出血))|((脑)&(梗死))|((脑)&(梗塞))|((心)&(梗))|((脑)&(供血))|冠心病|((动脉粥样)&(心脏病))|(脑血管)|((心)&(梗死))|((心)&(梗塞))|(脑血管)
          ",Dis$history) ~ 1,
    TRUE ~ 0
  ))
  return(Dis)
}


drugAnalysisRes = function(){
  drug_out = tbl_df(
    sqlQuery(cn,"
             Select distinct a.MASTER_INDEX,a.ICD,a.DISEASE_NAME,a.I_ITEM_NAME,a.FEE
             from Demo_Drug_ICD_Disease_OUT a  join
             (Select distinct  MASTER_INDEX From TATING_OUT ) b on a.MASTER_INDEX=b.MASTER_INDEX
             where a.I_ITEM_NAME like '%他汀%'
             order by  a.MASTER_INDEX
             "))
  
  drug_time_out =  tbl_df(
    sqlQuery(cn,"
             select master_index,datediff(month,min(diag_time),max(diag_time)) time  from Demo_Drug_ICD_Disease_OUT 
             group by master_index
             "))
  
  
  
  
  drug_in = tbl_df(
    sqlQuery(cn,"
             Select distinct a.MASTER_INDEX,a.ICD,a.DISEASE_NAME,a.I_ITEM_NAME,a.FEE
             from Demo_Drug_ICD_Disease_IN a  join
             (Select distinct  MASTER_INDEX From TATING_IN ) b on a.MASTER_INDEX=b.MASTER_INDEX
             where a.I_ITEM_NAME like '%他汀%'
             order by  a.MASTER_INDEX
             "))
  drug_time_in =  tbl_df(
    sqlQuery(cn,"
             select master_index,datediff(month,min(diag_time),max(diag_time)) time  from Demo_Drug_ICD_Disease_IN 
             group by master_index
             "))
  
  intersect(drug_in$MASTER_INDEX,drug_out$MASTER_INDEX)
  
  by_master_index_out = 
      by_master_index_out <- group_by(drug_time_out, master_index)
  
  
  basicData_out <-  by_master_index_out %>% 
    summarise(time = sum(time))
  
  
  by_master_index_in <- group_by(drug_time_in, master_index)
  
  
  basicData_in <-  by_master_index_in %>% 
    summarise(time = sum(time))
  
  
  by_master_index_drug_out <- group_by(drug_out, MASTER_INDEX)
  
  
  basicData_drug_out <-  by_master_index_drug_out %>% 
    summarise( FEE= sum(FEE),drugNmae = first(I_ITEM_NAME) )
  
  by_master_index_drug_in <- group_by(drug_in, MASTER_INDEX)
  
  
  basicData_drug_in <-  by_master_index_drug_in %>% 
    summarise( FEE= sum(FEE),drugNmae = first(I_ITEM_NAME) )
  
  drug_out_fee_time = merge( basicData_drug_out,basicData_out,by.x = "MASTER_INDEX",by.y = "master_index")
  drug_in_fee_time = merge( basicData_drug_in,basicData_in,by.x = "MASTER_INDEX",by.y = "master_index")
  
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  #install.packages("factoextra")
  library(factoextra) # clustering algorithms & visualization
  set.seed(123)
  k.values <- 1:15
  wss <- function(k) {
    kmeans(na.omit(drug_out_fee_time[,c(2,4)]), k, nstart = 10 )$tot.withinss
  }
  wss_values <- map_dbl(k.values, wss)
  plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  
  
  k4 <- kmeans( drug_out_fee_time[,c(2,4)], centers = 6, nstart = 25)
  fviz_cluster(k4, geom = "point",data = drug_out_fee_time[,c(2,4)])+theme_classic()
  
  group = k4$cluster
  #varience
  drug_out_fee_time = cbind(drug_out_fee_time,group)
  
  
  set.seed(123)
  k.values <- 1:15
  drug_in_fee_time = na.omit( drug_in_fee_time)
  wss <- function(k) {
    kmeans(na.omit(drug_in_fee_time[,c(2,4)]), k, nstart = 10 )$tot.withinss
  }
  wss_values <- map_dbl(k.values, wss)
  plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  k4 <- kmeans(drug_in_fee_time[,c(2,4)], centers = 6, nstart = 25)
  fviz_cluster(k4, geom = "point",data = drug_in_fee_time[,c(2,4)])+theme_classic()
  group = k4$cluster
  drug_in_fee_time = cbind(drug_in_fee_time,group)
  drug_in_fee_time = drug_in_fee_time[drug_in_fee_time$time<50,]
  drugAnalysisRes = list(drug_out_fee_time,drug_in_fee_time)
  names(drugAnalysisRes) = c("drug_out_fee_time","drug_in_fee_time")
  
 return(drugAnalysisRes)
  
  
  }


drugAnalysis = function(){
  drug_out = tbl_df(
    sqlQuery(cn,"
  Select distinct a.MASTER_INDEX,a.ICD,a.DISEASE_NAME,a.I_ITEM_NAME,a.FEE
  from Demo_Drug_ICD_Disease_OUT a  join
  (Select distinct  MASTER_INDEX From TATING_OUT ) b on a.MASTER_INDEX=b.MASTER_INDEX
  where a.I_ITEM_NAME like '%他汀%'
  order by  a.MASTER_INDEX
             "))
  
  drug_time_out =  tbl_df(
    sqlQuery(cn,"
select master_index,datediff(month,min(diag_time),max(diag_time)) time  from Demo_Drug_ICD_Disease_OUT 
group by master_index
 "))
  
  
  
  
  drug_in = tbl_df(
    sqlQuery(cn,"
  Select distinct a.MASTER_INDEX,a.ICD,a.DISEASE_NAME,a.I_ITEM_NAME,a.FEE
  from Demo_Drug_ICD_Disease_IN a  join
  (Select distinct  MASTER_INDEX From TATING_IN ) b on a.MASTER_INDEX=b.MASTER_INDEX
  where a.I_ITEM_NAME like '%他汀%'
  order by  a.MASTER_INDEX
             "))
  drug_time_in =  tbl_df(
    sqlQuery(cn,"
select master_index,datediff(month,min(diag_time),max(diag_time)) time  from Demo_Drug_ICD_Disease_IN 
group by master_index
 "))
  

  tmp_1 = table( drug_out$I_ITEM_NAME)
  tmp_2 = table( drug_in$I_ITEM_NAME)
  T <- as.table( rbind(tmp_1,tmp_2))
  dimnames(T) <- list(Sustain = c("Outpatient", "Inpatient"),
                      Dis = c("Atorvastatin","Fluvastatin","Lovastatin","Pivastatin","Pravastatin",
                              "Rosuvastatin","Simvastatin"))
  Tpertage = t(data.frame(T[1,]/sum( T[1,])))
  Tpertage = as.data.frame(rbind(Tpertage,t(data.frame(T[2,]/sum( T[2,])))))
  rownames(Tpertage) =  c("Outpatient", "Inpatient")
  Tpertage
  Tsq <- chisq.test(T)
  Drug = t( matrix( paste(as.matrix(T),"(",formatC(as.matrix(Tpertage)*100,format = "f",digits = 2),"%",")",sep = ""),
                     nrow = nrow(T),dimnames = dimnames(T)))
  write.csv(x= Drug,file = "Drug.csv")
  
  by_master_index_out <- group_by(drug_time_out, master_index)
  
  basicData_out <-  by_master_index_out %>% 
    summarise(time = sum(time))
  
  by_master_index_in <- group_by(drug_time_in, master_index)
  
  basicData_in <-  by_master_index_in %>% 
    summarise(time = sum(time))
  
  by_master_index_drug_out <- group_by(drug_out, MASTER_INDEX)
  
  basicData_drug_out <- by_master_index_drug_out  %>% 
    summarise( FEE= sum(FEE),drugNmae = first(I_ITEM_NAME) )
  
  by_master_index_drug_in <- group_by(drug_in, MASTER_INDEX)
  
  basicData_drug_in <-  by_master_index_drug_in %>% 
    summarise( FEE= sum(FEE),drugNmae = first(I_ITEM_NAME) )
  
  
  
  
  drug_out_fee_time = merge( basicData_drug_out,basicData_out,by.x = "MASTER_INDEX",by.y = "master_index")
  nrow(drug_out_fee_time)
  drug_in_fee_time = merge( basicData_drug_in,basicData_in,by.x = "MASTER_INDEX",by.y = "master_index")
  nrow(drug_in_fee_time)
  #install.packages("tidyverse")
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  #install.packages("factoextra")
  library(factoextra) # clustering algorithms & visualization
  set.seed(123)
  k.values <- 1:15
  wss <- function(k) {
    kmeans(na.omit(drug_out_fee_time[,c(2,4)]), k, nstart = 10 )$tot.withinss
  }
  wss_values <- map_dbl(k.values, wss)
  plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  
  
  k4 <- kmeans( drug_out_fee_time[,c(2,4)], centers = 6, nstart = 25)
  fviz_cluster(k4, geom = "point",data = drug_out_fee_time[,c(2,4)])+theme_classic()
  
  group = k4$cluster
  #varience
  drug_out_fee_time = cbind(drug_out_fee_time,group)
  
  
  set.seed(123)
  k.values <- 1:15
  drug_in_fee_time = na.omit( drug_in_fee_time)
  wss <- function(k) {
    kmeans(na.omit(drug_in_fee_time[,c(2,4)]), k, nstart = 10 )$tot.withinss
  }
  wss_values <- map_dbl(k.values, wss)
  plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  k4 <- kmeans(drug_in_fee_time[,c(2,4)], centers = 6, nstart = 25)
  fviz_cluster(k4, geom = "point",data = drug_in_fee_time[,c(2,4)])+theme_classic()
  group = k4$cluster
  drug_in_fee_time = cbind(drug_in_fee_time,group)
  drug_in_fee_time = drug_in_fee_time[drug_in_fee_time$time<30,]
  table(drug_in_fee_time$group)
  result = list( drug_in_fee_time,drug_out_fee_time)
  names(result) = c("drug_in","drug_out")
  drug_level(drug_in_fee_time)
  drug_level(drug_out_fee_time)
  
  write.csv(x= drug_level(drug_in_fee_time),file = "drug_in_fee_time_sum.csv")
  write.csv(x= drug_level(drug_out_fee_time),file = "drug_out_fee_time_sum.csv")
  return(result)
  


}

drug_level =function(data){
  #data = drug_out_fee_time
  r_1=c()
  r_2=c()
  cout=c()

  for (i in 1:6){
    tmp_1 =data[data$group==i,]$FEE
    tmp_2 = data[data$group==i,]$time
    t=t.test(tmp_1)
    r = formatC(c(t$estimate,t$conf.int[c(1,2)]),format = "f",digits = 2)
    r = paste(r[1],"(",r[2],",",r[3],")",sep = "")
    r_1=c(r_1,r)
    t=t.test(tmp_2)
    r = formatC(c(t$estimate,t$conf.int[c(1,2)]),format = "f",digits = 2)
    r = paste(r[1],"(",r[2],",",r[3],")",sep = "")
    r_2=c(r_2,r)
    cout=c(cout,length(tmp_1))
  }
  drug_level = as.data.frame(cbind(cout,cbind(r_1,r_2)))
  dimnames(drug_level) = list(Level = c(2,6,3,5,1,4),Item = c("Number","Fee","Time Window"))
  return(drug_level)
  
} 




colculateIndex = function(input,areaCateg){
  #input=result
  # result = fivenumLAP_OUT("LPA")
  input[[1]]$PROVINCE = substr(input[[1]]$PROVINCE,1,2)
  input[[2]]$PROVINCE = substr(input[[2]]$PROVINCE,1,2)
  data_1 = merge(input[[1]],areaCateg,by.x = "PROVINCE",by.y = "ares",all.x=T)
  data_2 = merge(input[[2]],areaCateg,by.x = "PROVINCE",by.y = "ares",all.x=T)
  breaks = c(0,45,70,150)
  label = c(1,2,3)
  data_1$age_group = cut(data_1$age,breaks = breaks,labels = label)
  data_2$age_group = cut(data_2$age,breaks = breaks,labels = label)
  age_1 = table( data_1$age_group)
  age_2 = table( data_2$age_group)
  M <- as.table( rbind(age_1, age_2))
  M
  dimnames(M) <- list(Sustain = c("Used", "Non-used"),
                      age = c("≤45","46-70", "≥71"))
  Mpertage = t(data.frame(M[1,]/sum( M[1,])))
  Mpertage = as.data.frame(rbind(Mpertage,t(data.frame(M[2,]/sum( M[2,])))))
  rownames(Mpertage) =  c("Used", "Non-used")
  Mpertage
  Xsq <- chisq.test(M)
  Xsq
  Mpertage 
  
  gender_1 = table( data_1$gender)[1:2]
  gender_2 = table( data_2$gender)[1:2]
  G <- as.table( rbind(gender_1,gender_2))
  dimnames(G) <- list(Sustain = c("Used", "Non-used"),
                      age = c("Male","Female"))
  Gpertage = t(data.frame(G[1,]/sum( G[1,])))
  Gpertage = as.data.frame(rbind(Gpertage,t(data.frame(G[2,]/sum( G[2,])))))
  rownames(Gpertage) =  c("Used", "Non-used")
  Gpertage
  Gsq <- chisq.test(G)
  Gsq
  index_1 = c(fivenum(data_1$firstValue)[c(3,2,4)],fivenum(data_1$lastValue)[c(3,2,4)],
              fivenum(data_1$median)[c(3,2,4)]) 
  index_2 =  c(fivenum(data_2$firstValue)[c(3,2,4)],fivenum(data_2$lastValue)[c(3,2,4)],
               fivenum(data_2$median)[c(3,2,4)]) 
  
  t_1 = t.test(data_1$firstValue,data_2$firstValue)
  t_2 = t.test(data_1$lastValue,data_2$lastValue)
  t_3 = t.test(data_1$median,data_2$median)
  t_pair_1 = t.test(data_1$firstValue,data_1$lastValue,paired = T)
  t_pair_2 = t.test(data_2$firstValue,data_2$lastValue,paired = T)
  
  area_1 = table(data_1$code)
  area_2 = table(data_2$code)
  A <- as.table( rbind(area_1,area_2))
  dimnames(A) <- list(Sustain = c("Used", "Non-used"),
                      Area = c("North","South","Northwest","Xizang-Qinghai"))
  Apertage = t(data.frame(A[1,]/sum( A[1,])))
  Apertage = as.data.frame(rbind(Apertage,t(data.frame(A[2,]/sum( A[2,])))))
  rownames(Apertage) =  c("Used", "Non-used")
  Apertage
  Asq <- fisher.test(A)
  Asq  
  
  change_1 = table(data_1$outcome)
  change_2 = table(data_2$outcome)
  C <- as.table( rbind(change_1,change_2))
  dimnames(C) <- list(Sustain = c("Used", "Non-used"),
                      Outcome = c("decrease","unchange","increase"))
  Cpertage = t(data.frame(C[1,]/sum( C[1,])))
  Cpertage = as.data.frame(rbind(Cpertage,t(data.frame(C[2,]/sum( C[2,])))))
  rownames(Cpertage) =  c("Used", "Non-used")
  Cpertage
  Csq <-  chisq.test(C)
  Csq  
  
  C_1= t( matrix( paste(as.matrix(C),"(",formatC(as.matrix(Cpertage)*100,format = "f",digits = 2),"%",")",sep = ""),
                  nrow = nrow(C),dimnames = dimnames(C)))
  
  M_1 =t( matrix( paste(as.matrix(M),"(",formatC(as.matrix(Mpertage)*100,format = "f",digits = 2),"%",")",sep = ""),
                  nrow = nrow(M),dimnames = dimnames(M)))
  
  G_1 = t( matrix( paste(as.matrix(G),"(",formatC(as.matrix(Gpertage)*100,format = "f",digits = 2),"%",")",sep = ""),
                   nrow = nrow(G),dimnames = dimnames(G)))
  A_1 = t( matrix( paste(as.matrix(A),"(",formatC(as.matrix(Apertage)*100,format = "f",digits = 2),"%",")",sep = ""),
                   nrow = nrow(A),dimnames = dimnames(A)))
  
  #summary(data_1$firstValue)
 

 
  #t=t.test(data_1$firstValue)
  r_1_f = formatC( c(mean(data_1$firstValue),quantile(data_1$firstValue, c(0.025, 0.975))),
                   format = "f",digits = 2)
  r_1 = paste(r_1_f[1],"(",r_1_f[2],",",r_1_f[3],")",sep = "")
  
  #t=t.test(data_1$lastValue)
  r_1_l = formatC(c(mean(data_1$lastValue),quantile(data_1$lastValue, c(0.025, 0.975))),
                  format = "f",digits = 2)
  r_1 =rbind(r_1,paste(r_1_l[1],"(",r_1_l[2],", ",r_1_l[3],")",sep = ""))
  
  #t=t.test(data_1$median)
  r_1_m =  formatC(c(mean(data_1$median),quantile(data_1$median, c(0.025, 0.975))),
                   format = "f",digits = 2)
  r_1 =as.data.frame(rbind(r_1,paste(r_1_m[1],"(",r_1_m[2],", ",r_1_m[3],")",sep = "")))
 
  #t=t.test(data_2$firstValue)
  r_2_f = formatC(c(mean(data_2$firstValue),quantile(data_2$firstValue, c(0.025, 0.975)))
                  ,format = "f",digits = 2)
  r_2 = paste(r_2_f[1],"(",r_2_f[2],", ",r_2_f[3],")",sep = "")
  
  #t=t.test(data_2$lastValue)
  r_2_l = formatC(c(mean(data_2$lastValue),quantile(data_2$lastValue, c(0.025, 0.975))),
                  format = "f",digits = 2)
  r_2 =rbind(r_2,paste(r_2_l[1],"(",r_2_l[2],", ",r_2_l[3],")",sep = ""))
  
  #t=t.test(data_2$median)
  r_2_m =  formatC(c(mean(data_2$median),quantile(data_2$median, c(0.025, 0.975))),
                   format = "f",digits = 2)
  r_2 =as.data.frame(rbind(r_2,paste(r_2_m[1],"(",r_2_m[2],", ",r_2_m[3],")",sep = "")))
  
  r_value = cbind(r_1,r_2)  
  
  rownames(r_value) = c("First_value","Last_value","Median_value")
  colnames(r_value) = c("Sustain-used","Non-sustain-used")
  
  
  
  res = list(Age=M_1,Age_chi=Xsq,
             Gender=G_1,Gender_chi = Gsq,
             Area = A_1,Area_chi = Asq,
             Change = C_1,Change_chi = Csq,
             Result_first = t_1,
             Result_last = t_2,
             Result_medien = t_3,
             Result_pair_sustain_used = t_pair_1,
             Result_pair_non_sustain_used = t_pair_2,
             Demo = r_value
  )
  return(res)
}


fivenumLAP_OUT = function(index){
  #index="LPA"
  sqlQuery(cn,"drop table tmp") 
  ## sql code to select the index
  sql = "
  Select distinct a.MASTER_INDEX,
  a.BIRTHDAY,
  a.DIAG_TIME,
  a.SEX,
  a.ORG_CODE,
  b.RESULTS,
  b.ITEM_ENAME,
  b.REFRANGE,
  b.I_ITEM_ENAME
  into tmp
  from  LPA_CHECK0_5YEARS_2times_OUT a
  join (SELECT * from   LDL_LAB_OUT a where a.I_ITEM_ENAME like '%$index%'   ) b  on a.REG_CODE=b.REG_CODE 
  "
  sql <- fn$identity(sql)
  sqlQuery(cn,sql)
  #tating
  sqlQuery(cn,"drop table result") 
  sql = "
  Select distinct a.MASTER_INDEX,a.sex, a.RESULTS, ORG_CODE, datediff(year,(a.birthday),(a.DIAG_TIME))age  
  into result
  from  tmp a join 
  (Select distinct  MASTER_INDEX From TATING_OUT ) b on a.MASTER_INDEX=b.MASTER_INDEX
  order by  a.MASTER_INDEX
  "
  sqlQuery(cn,sql)
  
  
  data = tbl_df(
    sqlQuery(cn,"
             Select  distinct a.*,b.PROVINCE from result a,pub_Patient b where a.MASTER_INDEX = b.MASTER_INDEX  order by a.MASTER_INDEX
             ") )
  data_trend_tating = trendAnalysis(data,1)
  #head(data_trend)
  
 
  
  #non_tating
  sqlQuery(cn,"drop table result") 
  sql = "
Select distinct a.MASTER_INDEX,a.sex, a.RESULTS, ORG_CODE, datediff(year,(a.birthday),(a.DIAG_TIME))age  
into result
from  tmp a join 
(Select distinct  MASTER_INDEX From NON_TATING_OUT ) b on a.MASTER_INDEX=b.MASTER_INDEX
order by  a.MASTER_INDEX
"
  sqlQuery(cn,sql)
  data = tbl_df(
    sqlQuery(cn,"
Select  distinct a.*,b.PROVINCE from result a,pub_Patient b where a.MASTER_INDEX = b.MASTER_INDEX  order by a.MASTER_INDEX
") )
  data_trend_non_tating = trendAnalysis(data,0)
  #head(data_trend)
  

  list_data <- list( data_trend_tating,data_trend_non_tating)
  
  return(list_data)
  
} 


fivenumLAP_IN = function(index){
  #index="LPA"
  sqlQuery(cn,"drop table tmp") 
  ## sql code to select the index
  sql = "
  Select distinct a.MASTER_INDEX,
  a.BIRTHDAY,
  a.DIAG_TIME,
  a.SEX,
  a.ORG_CODE,
  b.RESULTS,
  b.ITEM_ENAME,
  b.REFRANGE,
  b.I_ITEM_ENAME
  into tmp
  from  LPA_IN_0_5Years_2_Time a
  join (SELECT * from LDL_LAB_IN a where a.I_ITEM_ENAME like '%$index%'   ) b  on a.REG_CODE=b.REG_CODE
  "
  sql <- fn$identity(sql)
  sqlQuery(cn,sql)
  #tating
  sqlQuery(cn,"drop table result") 
  sql = "
 Select distinct a.MASTER_INDEX,a.sex, a.RESULTS, ORG_CODE, datediff(year,(a.birthday),(a.DIAG_TIME))age  
  into result
  from  tmp a join 
  (Select distinct  MASTER_INDEX From TATING_IN ) b on a.MASTER_INDEX=b.MASTER_INDEX
  order by  a.MASTER_INDEX
  "
  sqlQuery(cn,sql)
  
  
  data = tbl_df(
    sqlQuery(cn,"
             Select  distinct a.*,b.PROVINCE from result a,pub_Patient b where a.MASTER_INDEX = b.MASTER_INDEX  order by a.MASTER_INDEX
             ") )
  data_trend_tating = trendAnalysis(data,1)
  #head(data_trend)
  
  
  
  #non_tating
  sqlQuery(cn,"drop table result") 
  sql = "
  Select distinct a.MASTER_INDEX,a.sex, a.RESULTS, ORG_CODE, datediff(year,(a.birthday),(a.DIAG_TIME))age  
  into result
  from  tmp a join 
  (Select distinct  MASTER_INDEX From NON_TATING_IN ) b on a.MASTER_INDEX=b.MASTER_INDEX
  order by  a.MASTER_INDEX
  "
  sqlQuery(cn,sql)
  data = tbl_df(
    sqlQuery(cn,"
             Select  distinct a.*,b.PROVINCE from result a,pub_Patient b where a.MASTER_INDEX = b.MASTER_INDEX  order by a.MASTER_INDEX
             ") )
  data_trend_non_tating = trendAnalysis(data,0)
  #head(data_trend)
  

  list_data <- list( data_trend_tating,data_trend_non_tating)
  
  return(list_data)
  
} 


dataIntergrateOut = function(){
  result_out=list()
  index = c("LPA","LDL_C","HDL_C","APO-A","APO-B","TC","TG","CRP")
  #outpatients
  for (i in 1:length(index)){
    tmp = fivenumLAP_OUT(index[i]) # if we input the wrong name of index, the terminal happened
    #tmp_1 = as.matrix(tmp[[1]])
    #tmp_2 = as.matrix(tmp[[2]])
    
    tmp_1 = tmp[[1]][,c(1:9)]
    tmp_2 = tmp[[2]][,c(1:9)]
    tmp_12 = rbind(tmp_1,tmp_2)
    outcome = c(as.character(tmp[[1]][,10]$rsq),as.character(tmp[[2]][,10]$rsq))
    tmp_12$outcome= outcome
    
    result_out[i] = list(tmp_12)
  }
  names(result_out) = index
  return(result_out)
}


dataIntergrateIn = function(){
  result_in=list()
  index = c("LPA","LDL_C","HDL_C","APO-A","APO-B","TC","TG","CRP")
  #outpatients
  for (i in 1:length(index)){
    tmp = fivenumLAP_IN(index[i])
    tmp_1 = tmp[[1]][,c(1:9)]
    tmp_2 = tmp[[2]][,c(1:9)]
    tmp_12 = rbind(tmp_1,tmp_2)
    outcome = c(as.character(tmp[[1]][,10]$rsq),as.character(tmp[[2]][,10]$rsq))
    tmp_12$outcome= outcome
    
    result_in[i] = list(tmp_12)
  }
  names(result_in) = index
  return(result_in)
}

library(dplyr)
library("RODBC")
library(Hmisc)
#install.packages("dplyr")
library(dplyr)
#library(dbplyr)
library(data.table)
library(odbc)
#devtools::install_github("dgrtwo/broom")
library(broom)
#if(!require(broom)){install.packages("broom")}
library(gsubfn) 

#library(RMSSQL)
#library(dplyr.mssql)
#library(Rmisc)
cn <- odbcDriverConnect(
  connection=
    "Driver={SQL Server Native Client 11.0};server=localhost; 
    database=DB_LPA_ZONG;
    trusted_connection=yes;")
#source("C:\\Users\\Suvalue\\Documents\\LAPRealWorld\\tools.R",encoding="utf-8")
areaCateg = read.table("C:\\Users\\Suvalue\\Documents\\LAPRealWorld\\areaCode.txt",header = T,sep=",")
#head(areaCateg)
#colnames(areaCateg)




# tab_1 = function(data){
#   
# }


