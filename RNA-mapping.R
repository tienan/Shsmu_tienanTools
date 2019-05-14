#########################SHSMU


#########################




#index build 
setwd("/LVM1/LVM/hg38/")
system("salmon index -t hg38.transcript.fa  -i  salmon_transcript_index  --type quasi -k 31", intern = TRUE)
?system2
#getwd()
#system2("salmon index -t hg38.transcript.fa  -i  salmon_transcript_index  --type quasi -k 31",stdout = T)
setwd("/LVM1/LVM1/2016-10-2-cancer-cell-fq-2/")
?grepl
transcripts_index = "/LVM1/LVM/hg38/salmon_transcript_index"
Libtype = "--libType"
read_1 = "-1"
read_2 = "-2"
list = dir()
list_1 = list[grepl("*val*",list)]
#salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq -o transcripts_quant
source("https://bioconductor.org/biocLite.R")
biocLite("update")
biocLite("limma")
biocLite("edgeR")


tmp_name=c()
tmp=c()
for (i in 1:length(list_1)){
  tmp = unlist(strsplit(list_1[i],"_R"))
  # tmp = unlist(strsplit(tmp[1],"_"))
  # tmp= paste(tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],sep="_")}
  tmp_name[i] = tmp[1]
  system(paste("salmon","quant -i",transcripts_index,"--libType U",
                   read_1 ,list_1[i],
                   read_2,list_1[i+1],"-o",paste("salom",tmp_name[i],sep = "_"),sep=" "), intern = TRUE)
  i=i+1
}

rm(quant_list ,quant_list_1,quant_list_2)

system("find -name 'quant.sf' > quant_list.txt")
quant_list = (read.table("quant_list.txt"))
quant_list_1 = as.character(quant_list[order(quant_list),])
quant_list_2 = quant_list_1[grepl(quant_list_1,pattern = "sal")] 



# quant_list_2[82] = "./salom_A549_siRNA_rcc2_as_2_S44/quant.sf"
# quant_list_2[79] = "./salom_A549_siRNA_con_3_S36/quant.sf"      
# quant_list_2[80] = "./salom_A549_siRNA_hcg18_3_S42/quant.sf"    
# quant_list_2[81] = "./salom_A549_siRNA_loc146800_3_S39/quant.sf"
rm(salmon_exp)
salmon_exp = read.table(quant_list_2[1],sep = "\t",header = T)[,c(1,4)]
for (i in 2: length(quant_list_2)){
  tmp = read.table(quant_list_2[i],sep = "\t",header = T)
  salmon_exp = cbind(salmon_exp,tmp[4])
}

rm(tmp)
name = NULL
for (i in 1: length(quant_list_2)){
  tmp = unlist(strsplit(quant_list_2[i],split = "_"))
  name[i] = paste(tmp[2],tmp[3],tmp[4],sep = "_")
}
colnames(salmon_exp) = c("GeneName",name)
head(salmon_exp)
write.csv(salmon_exp,file = "salmon_exp_1st_seq.csv")
colnames(salmon_exp)

#1
name = "1_p10vsp0.csv"
tmp= salmon_exp[,c(1:8,65:72)+3]
group = NULL
group[1:8]="P10"
group[9:16]="P0"
t1 = edgeR::DGEList(tmp,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = name)






##################################second Seq 

setwd("/S2600T-1/20170930/")
list = dir()
list_1 = list[grepl("*val*",list)]
list_1 = list_1[grepl("fq",list_1)]
tmp_name=c()
tmp=c()
for (i in 1:length(list_1)){
  tmp = unlist(strsplit(list_1[i],"_R"))
  # tmp = unlist(strsplit(tmp[1],"_"))
  # tmp= paste(tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],sep="_")}
  tmp_name[i] = tmp[1]
  system(paste("salmon","quant -i",transcripts_index,"--libType U",
               read_1 ,list_1[i],
               read_2,list_1[i+1],"-o",paste("salom",tmp_name[i],sep = "_"),sep=" "), intern = TRUE)
  i=i+1
}

system("find -name 'quant.sf' > quant_list.txt")
quant_list = (read.table("quant_list.txt"))
quant_list_1 = as.character(quant_list[order(quant_list),])
quant_list_2 = quant_list_1[grepl(quant_list_1,pattern = "sal")] 

rm(salmon_exp)
salmon_exp = read.table(quant_list_2[1],sep = "\t",header = T)[,c(1:4)]
for (i in 2: length(quant_list_2)){
  tmp = read.table(quant_list_2[i],sep = "\t",header = T)
  salmon_exp = cbind(salmon_exp,tmp[4])
}
library("limma")

name_1 = strsplit2(quant_list_2,split = "./")[,2]
colnames(salmon_exp) = c(colnames(salmon_exp)[1:3],name_1)
write.csv(salmon_exp,file = "base2bexp.csv")
library(edgeR)


p10_con_vs_dl = salmon_exp[,c(36:43,4:11)]
group = NULL
group[1:8]="Con"
group[9:16]="DL"
t1 = edgeR::DGEList(p10_con_vs_dl,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = "p10_con_vs_dl.csv")

p30_con_vs_dl = salmon_exp[,c(44:51,12:19)]
group = NULL
group[1:8]="Con"
group[9:16]="DL"
t1 = edgeR::DGEList(p30_con_vs_dl,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = "p30_con_vs_dl.csv")

p10_CT_con_vs_dl = salmon_exp[,c(116:123,76:83)]
group = NULL
group[1:8]="Con"
group[9:16]="DL"
t1 = edgeR::DGEList(p10_CT_con_vs_dl,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = "p10_CT_con_vs_dl.csv")

p30_CT_con_vs_dl = salmon_exp[,c(124:131,84:91)]
group = NULL
group[1:8]="Con"
group[9:16]="DL"
t1 = edgeR::DGEList(p30_CT_con_vs_dl,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = "p30_CT_con_vs_dl.csv")

############################################

#1
name = "CC_p10vsp0.csv"
tmp= salmon_exp[,c(1:8,65:72)+3]
group = NULL
group[1:8]="P10"
group[9:16]="P0"
t1 = edgeR::DGEList(tmp,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = name)

#2
name = "CC_p30vsp0.csv"
tmp= salmon_exp[,c(9:16,65:72)+3]
group = NULL
group[1:8]="P30"
group[9:16]="P0"
t1 = edgeR::DGEList(tmp,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = name)

#3
name = "CC_p10vsp30.csv"
tmp= salmon_exp[,c(9:16,1:8)+3]
group = NULL
group[1:8]="P30"
group[9:16]="P10"
t1 = edgeR::DGEList(tmp,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = name)


#4
name = "CT_p10vsp0.csv"
tmp= salmon_exp[,c(73:80,105:112)+3]
group = NULL
group[1:8]="P10"
group[9:16]="P0"
t1 = edgeR::DGEList(tmp,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = name)

#5
name = "CT_p30vsp0.csv"
tmp= salmon_exp[,c(81:88,105:112)+3+72]
group = NULL
group[1:8]="P30"
group[9:16]="P0"
t1 = edgeR::DGEList(tmp,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = name)

#6
name = "CT_p10vsp30.csv"
tmp= salmon_exp[,c(81:88,73:80)+3]
group = NULL
group[1:8]="P30"
group[9:16]="P10"
t1 = edgeR::DGEList(tmp,group = as.factor(group),genes = salmon_exp$Name)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
dat = as.data.frame(cbind(t3$genes,t3$table))
write.csv(x = dat,file = name)

######################################################
setwd("/LVM1/LVM1/2016-10-2-cancer-cell-fq-2/")
?grepl
transcripts_index = "/LVM1/LVM/hg38/salmon_transcript_index"
Libtype = "--libType"
read_1 = "-1"
read_2 = "-2"
list = dir()
list_1 = list[grepl("*val*",list)]



######################################################