#Transcript_of_RNA_Signature_2

######################Setting the work direction
getwd()
setwd("../DL/")
######################Reading FPKM value
geneNmae = read.table("uniqGeneName",header = F)
file_list = dir(pattern = "*.fpkm*")
dat_tmp = read.table(file_list[1],header = T,sep = "\t")
tmp=c()
tmp_name=dat_tmp$tracking_id
for (i in 1:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp = dat_tmp[,c(1,10)]
  dat_tmp = dat_tmp[dat_tmp $tracking_id%in%geneNmae$V1,]
  dat_tmp = dat_tmp[order(dat_tmp$tracking_id),]
  dat_tmp = summarise(group_by(dat_tmp,tracking_id),maxFPKM=max(FPKM))
  tmp_name = intersect(tmp_name,dat_tmp$tracking_id)
  tmp[i] = nrow(dat_tmp)
  #head(dat_tmp)
  #tmp_file = cbind(tmp_file,dat_tmp[,2])
}
length(tmp_name)

tmp_file = read.table(file_list[1],header = T,sep = "\t")
#head(tmp_file)
#nrow(tmp_file)
#colnames(tmp_file)
tmp_file = tmp_file[,c(1,10)]
tmp_file = tmp_file[tmp_file$tracking_id%in%tmp_name,]
tmp_file = tmp_file[order(tmp_file$tracking_id),]
head(tmp_file)
library(dplyr)
tmp_file = summarise(group_by(tmp_file,tracking_id),maxFPKM=max(FPKM))
nrow(tmp_file)

for (i in 2:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp = dat_tmp[,c(1,10)]
  dat_tmp = dat_tmp[dat_tmp $tracking_id%in%tmp_name,]
  dat_tmp = dat_tmp[order(dat_tmp$tracking_id),]
  dat_tmp = summarise(group_by(dat_tmp,tracking_id),maxFPKM=max(FPKM))
  #nrow(dat_tmp)
  #head(dat_tmp)
  tmp_file = cbind(tmp_file,dat_tmp[,2])
}
colname = gsub(pattern = "_genes.fpkm_tracking",replacement = "",file_list) 
colname = gsub(pattern = "1975",replacement = "H1975",colname) 
tmp = strsplit2(x = colname,split = "_")
group = as.factor(paste(tmp[,1],tmp[,2],sep="_"))
colnames(tmp_file) = c("Gene",colname)
condition = ifelse(grepl(pattern = "con",colname),"Con","DL")
library(edgeR)
dat_1 = tmp_file[,-1]
rownames(dat_1) = tmp_file[,1]

library(edgeR)

y = DGEList(dat_1,group = group,genes = dat_1)
library(org.Hs.eg.db)
rownames(y)
y$genes$Symbol<-mapIds(org.Hs.eg.db,rownames(y),keytype="SYMBOL",column="ENTREZID")


head(y$genes)
y<-y[!is.na(y$genes$Symbol), ]
dim(y)
rownames(y)=y$genes$Symbol
y$genes$Symbol<-mapIds(org.Hs.eg.db,y$genes$Symbol,keytype="ENTREZID",column="SYMBOL")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")
#library(org.Mm.eg.db)

#y$genes$Symbol = mapIds(org.Hs.eg.db,rownames(y),keytype = "ENTREZID",column = "SYMBOL")
#group=factor(name)
#table(group)
#y = y[!is.na(y$genes$Symbol ),]
#dim(y)
#keep 
y=calcNormFactors(y)
y$samples

pch<- c(0,1,2,15)
colors<- rep(c("blue","red"),2)
plotMDS(y,col=colors[group],pch=pch[group])
legend("topright",legend=levels(group),pch=pch,col=colors,ncol=2)


design<-model.matrix(~0+group)
colnames(design)<-levels(group)
design
y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
group
dlcon<-makeContrasts(H1975_dl-H1975_con,levels=design)
res<-glmQLFTest(fit,contrast=dlcon)
res$table[res$table$PValue<0.05,]

tr<-glmTreat(fit,contrast=dlcon,lfc=log2(1.5))
tr$table[tr$table$PValue<0.05,]
is.de<-decideTestsDGE(tr)
summary(is.de)
is.de<-decideTestsDGE(res,lfc = 0.8)
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")

# library(Hmisc)
# capitalize(tr$genes$Symbol)
# tr$genes$Symbol = capitalize(tr$genes$Symbol)
# tr$genes$Symbol = toupper(tr$genes$Symbol)
# 
# mapIds(org.Hs.eg.db,rownames(y),keytype="SYMBOL",column="ENTREZID")
# 
#

go<-goana(res,species="Hs",FDR = 0.1)
summary(res)
topTags(res)
res$table$PValue = res$table$PValue/7
keg<-kegga(res,species="Hs")
topKEGG(keg,n=15,truncate=34)
library(GO.db)
Rkeys(org.Hs.egGO2ALLEGS)

term = select(GO.db,keys=Rkeys(org.Hs.egGO2ALLEGS),columns="TERM")
head(term)  
term[grepl(pattern = "tnf",term$TERM,ignore.case = T),]


load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
idx<-ids2indices(Hs.c2,id=rownames(y))
dlvscom<-makeContrasts(H1975_dl-H1975_con,levels=design)
cam<-camera(y, idx, design,contrast=dlvscom,inter.gene.cor=0.01)
head(cam,14)

?decideTestsDGE
is.de<-
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")

topTags(res)
res$genes$Symbol = mapIds(org.Hs.eg.db,rownames(res),keytype = "ENTREZID",column = "SYMBOL")
is.de<-decideTestsDGE(res)
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")
logCPM<-cpm(y,prior.count=2,log=TRUE)
rownames(logCPM)<-y$genes$Symbol
colnames(logCPM)<-paste(y$samples$group,1:3,sep="-")
tr<-glmTreat(fit,contrast=ASAL,lfc=log2(1.5))
o<-order(tr$table$PValue)
nrow(tr)
is.de<-decideTestsDGE(tr)
summary(is.de)
tr
logCPM<-logCPM[is.de!=0,]
tr[is.de!=0,]$table




