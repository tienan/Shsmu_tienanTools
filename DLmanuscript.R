getwd()
setwd("~/DL/")
################################
geneNmae = read.table("uniqGeneName",header = F)
file_list = dir(pattern = "*.fpkm*")
dat_tmp = read.table(file_list[1],header = T,sep = "\t")
tmp=c()
tmp_name=dat_tmp$tracking_id
library(dplyr)
for (i in 1:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp = dat_tmp[,c(1,10)]
  dat_tmp = dat_tmp[dat_tmp $tracking_id%in%geneNmae$V1,]
  dat_tmp = dat_tmp[order(dat_tmp$tracking_id),]
  dat_tmp = summarise(group_by(dat_tmp,tracking_id),maxFPKM=max(FPKM))
  tmp_name = intersect(tmp_name,dat_tmp$tracking_id)
  tmp[i] = nrow(dat_tmp)
}
tmp_file = read.table(file_list[1],header = T,sep = "\t")
tmp_file = tmp_file[,c(1,10)]
tmp_file = tmp_file[tmp_file$tracking_id%in%tmp_name,]
tmp_file = tmp_file[order(tmp_file$tracking_id),]
tmp_file = summarise(group_by(tmp_file,tracking_id),maxFPKM=max(FPKM))
for (i in 2:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp = dat_tmp[,c(1,10)]
  dat_tmp = dat_tmp[dat_tmp $tracking_id%in%tmp_name,]
  dat_tmp = dat_tmp[order(dat_tmp$tracking_id),]
  dat_tmp = summarise(group_by(dat_tmp,tracking_id),maxFPKM=max(FPKM))
  tmp_file = cbind(tmp_file,dat_tmp[,2])
}
colname = gsub(pattern = "_genes.fpkm_tracking",replacement = "",file_list) 
colname = gsub(pattern = "1975",replacement = "H1975",colname) 
library(limma)
tmp = strsplit2(x = colname,split = "_")
group = as.factor(paste(tmp[,1],tmp[,2],sep="_"))
colnames(tmp_file) = c("Gene",colname)
condition = ifelse(grepl(pattern = "con",colname),"Con","DL")
library(edgeR)
dat_1 = tmp_file[,-1]
rownames(dat_1) = tmp_file[,1]
y = DGEList(dat_1,group = group,genes = dat_1)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
pch<- c(0,1,2,15)
colors<- rep(c("blue","red"),2)
plotMDS(y,col=colors[group],pch=pch[group],labels=rownames(y$samples))
legend("topright",legend=levels(group),pch=pch,col=colors,ncol=2)

dat_2 = dat_1[,c(-3,-9)]
####################H1975F
H1975F=dat_2[,1:5]
y = DGEList(H1975F,group = group[c(1:2,4:6)],genes = H1975F)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
t2 = edgeR::estimateCommonDisp(y)
t3 = edgeR::exactTest(t2,pair=c("H1975_con","H1975_dl"))
tableDEA_H1975F <- edgeR::topTags(t3, n = nrow(t3$table))$table
DegH1975F <- tableDEA_H1975F[tableDEA_H1975F$FDR <= 0.05,c(6:9)]

####################A549F
A549F=dat_2[,6:10]
y = DGEList(A549F,group = group[c(7:8,10:12)],genes = A549F)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
t2 = edgeR::estimateCommonDisp(y)
t3 = edgeR::exactTest(t2,pair=c("A549_con","A549_dl"))
tableDEA_A549F <- edgeR::topTags(t3, n = nrow(t3$table))$table
DegA549F <- tableDEA_A549F[tableDEA_A549F$FDR <= 0.05,c(6:9)]


############################
getwd()
setwd("~/DL/")
dat_Third = read.csv("ExpressDLThird.csv")
dat_3 = dat_Third[,-1]
rownames(dat_3) = dat_Third[,1]
##############A549S
A549S=dat_3[,c(4:6,10:12)]
group = as.factor(gsub(pattern = "[0-9]",replacement = "",colnames(A549S) ))
library(org.Hs.eg.db)
A549S$Symbol=mapIds(org.Hs.eg.db,rownames(A549S),keytype="ENTREZID",column="SYMBOL")
A549S<-A549S[!is.na(A549S$Symbol), ]
library(dplyr)
A549S<- A549S%>%
  group_by(Symbol) %>% 
  summarise_all(max)
A549S=data.frame(A549S)
rownames(A549S)=A549S[,1]
library(edgeR)
y = DGEList(A549S[,-1],group = group,genes = A549S[,-1])
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
t2 = edgeR::estimateCommonDisp(y)
t3 = edgeR::exactTest(t2,pair=c("AS","AL"))
tableDEA_A549S <- edgeR::topTags(t3, n = nrow(t3$table))$table
DegA549S <- tableDEA_A549S[tableDEA_A549S$FDR <= 0.05,c(7:10)]

###################################### H1975S
H1975S=dat_3[,c(16:18,22:24)]
group = as.factor(gsub(pattern = "[0-9]",replacement = "",colnames(H1975S) ))
H1975S$Symbol=mapIds(org.Hs.eg.db,rownames(H1975S),keytype="ENTREZID",column="SYMBOL")
H1975S<-H1975S[!is.na(H1975S$Symbol), ]
H1975S<- H1975S%>%
  group_by(Symbol) %>% 
  summarise_all(max)
H1975S=data.frame(H1975S)
rownames(H1975S)=H1975S[,1]
y = DGEList(H1975S[,-1],group = group,genes = H1975S[,-1])
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
t2 = edgeR::estimateCommonDisp(y)
t3 = edgeR::exactTest(t2,pair=c("HS","HL"))
tableDEA_H1975S <- edgeR::topTags(t3, n = nrow(t3$table))$table
DegH1975S <- tableDEA_H1975S[tableDEA_H1975S$FDR <= 0.05,c(7:10)]


##################################
DegF=union(rownames(tableDEA_A549F[tableDEA_A549F$FDR<= 0.05,]),rownames(tableDEA_H1975F[tableDEA_H1975F$FDR<= 0.05,]))
DegS=union(rownames(tableDEA_A549S[tableDEA_A549S$FDR<= 0.05,]),rownames(tableDEA_H1975S[tableDEA_H1975S$FDR<= 0.05,]))
Deg=union(DegF,DegS)
#################### ͬ????
Deg_A549F=tableDEA_A549F[rownames(tableDEA_A549F)%in%Deg,]
Deg_A549S=tableDEA_A549S[rownames(tableDEA_A549S)%in%Deg,]
Deg_H1975F=tableDEA_H1975F[rownames(tableDEA_H1975F)%in%Deg,]
Deg_H1975S=tableDEA_H1975S[rownames(tableDEA_H1975S)%in%Deg,]

IntDegA549Up=intersect(rownames(Deg_A549F[Deg_A549F$logFC>0,]),rownames(Deg_A549S[Deg_A549S$logFC>0,]))
IntDegA549Down=intersect(rownames(Deg_A549F[Deg_A549F$logFC<0,]),rownames(Deg_A549S[Deg_A549S$logFC<0,]))
IntDegH1975Up=intersect(rownames(Deg_H1975F[Deg_H1975F$logFC>0,]),rownames(Deg_H1975S[Deg_H1975S$logFC>0,]))
IntDegH1975Down=intersect(rownames(Deg_H1975F[Deg_H1975F$logFC<0,]),rownames(Deg_H1975S[Deg_H1975S$logFC<0,]))
DegUp=intersect(IntDegA549Up,IntDegH1975Up)
DegDown=intersect(IntDegA549Down,IntDegH1975Down)
Deg=union(DegUp,DegDown)
Deg_A549F=tableDEA_A549F[rownames(tableDEA_A549F)%in%Deg,]
Deg_A549S=tableDEA_A549S[rownames(tableDEA_A549S)%in%Deg,]
Deg_H1975F=tableDEA_H1975F[rownames(tableDEA_H1975F)%in%Deg,]
Deg_H1975S=tableDEA_H1975S[rownames(tableDEA_H1975S)%in%Deg,]
write.csv(Deg_A549F,"E:/??Ŀ/DL/Deg_A549F.csv")
write.csv(Deg_A549S,"E:/??Ŀ/DL/Deg_A549S.csv")
write.csv(Deg_H1975F,"E:/??Ŀ/DL/Deg_H1975F.csv")
write.csv(Deg_H1975S,"E:/??Ŀ/DL/Deg_H1975S.csv")

############??ɽͼ######
library(ggplot2)
library(RColorBrewer)
tableDEA_A549F$threshold[tableDEA_A549F$FDR<0.05&tableDEA_A549F$logFC>0]="up"
tableDEA_A549F$threshold[tableDEA_A549F$FDR<0.05&tableDEA_A549F$logFC<(-0)]="down"
tableDEA_A549F$threshold[(tableDEA_A549F$FDR>0.05)|(tableDEA_A549F$logFC<=0)&tableDEA_A549F$logFC>=(-0)]="normal"
tableDEA_A549F$threshold=factor(tableDEA_A549F$threshold,levels = c("up","down","normal"),ordered = T)
p <- ggplot(data = tableDEA_A549F, 
            aes(x = logFC, 
                y = -log10(FDR))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=threshold)) +
  scale_color_manual(values=c("red","blue", "grey"))+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

tableDEA_A549S$threshold[tableDEA_A549S$FDR<0.05&tableDEA_A549S$logFC>0]="up"
tableDEA_A549S$threshold[tableDEA_A549S$FDR<0.05&tableDEA_A549S$logFC<(-0)]="down"
tableDEA_A549S$threshold[(tableDEA_A549S$FDR>0.05)|(tableDEA_A549S$logFC<=0)&tableDEA_A549S$logFC>=(-0)]="normal"
tableDEA_A549S$threshold=factor(tableDEA_A549S$threshold,levels = c("up","down","normal"),ordered = T)
p <- ggplot(data = tableDEA_A549S, 
            aes(x = logFC, 
                y = -log10(FDR))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=threshold)) +
  scale_color_manual(values=c("red","blue", "grey"))+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

tableDEA_H1975F$threshold[tableDEA_H1975F$FDR<0.05&tableDEA_H1975F$logFC>0]="up"
tableDEA_H1975F$threshold[tableDEA_H1975F$FDR<0.05&tableDEA_H1975F$logFC<(-0)]="down"
tableDEA_H1975F$threshold[(tableDEA_H1975F$FDR>0.05)|(tableDEA_H1975F$logFC<=0)&tableDEA_H1975F$logFC>=(-0)]="normal"
tableDEA_H1975F$threshold=factor(tableDEA_H1975F$threshold,levels = c("up","down","normal"),ordered = T)
p <- ggplot(data = tableDEA_H1975F, 
            aes(x = logFC, 
                y = -log10(FDR))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=threshold)) +
  scale_color_manual(values=c("red","blue", "grey"))+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

tableDEA_H1975S$threshold[tableDEA_H1975S$FDR<0.05&tableDEA_H1975S$logFC>0]="up"
tableDEA_H1975S$threshold[tableDEA_H1975S$FDR<0.05&tableDEA_H1975S$logFC<(-0)]="down"
tableDEA_H1975S$threshold[(tableDEA_H1975S$FDR>0.05)|(tableDEA_H1975S$logFC<=0)&tableDEA_H1975S$logFC>=(-0)]="normal"
tableDEA_H1975S$threshold=factor(tableDEA_H1975S$threshold,levels = c("up","down","normal"),ordered = T)
p <- ggplot(data = tableDEA_H1975S, 
            aes(x = logFC, 
                y = -log10(FDR))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=threshold)) +
  scale_color_manual(values=c("red","blue", "grey"))+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p 
##############??ͼ###
library(pheatmap)
heatdata <- Deg_A549F[,c(1:5)]
annotation_col = data.frame(Type = factor(rep(c("Control","DL"), c(2,3))))
rownames(annotation_col) = colnames(heatdata)
ann_colors = list(Type = c(Control="blue",DL= "firebrick"))
pheatmap(scale(heatdata), 
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col,
         show_rownames = F,
         show_colnames = F,
         scale = "row",
         color =colorRampPalette(c("blue", "white","red"))(100),
         cellwidth = 40, cellheight = 0.1,
         fontsize = 10)

heatdata <- Deg_H1975F[,c(1:5)]
annotation_col = data.frame(Type = factor(rep(c("Control","DL"), c(2,3))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(scale(heatdata), 
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col,
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = F,
         show_colnames = F,
         scale = "row",
         cellwidth = 40, cellheight = 0.1,
         fontsize = 10)


heatdata <- Deg_A549S[,c(1:6)]
colnames(heatdata)=c("DL-1","DL-2","DL-3","Control-1","Control-2","Control-3")
heatdata=heatdata[,order(colnames(heatdata))]
annotation_col = data.frame(Type = factor(rep(c("Control","DL"), c(3,3))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(scale(heatdata), 
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col, 
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = F,
         show_colnames = F,
         scale = "row",
         cellwidth = 40, cellheight = 0.1,
         fontsize = 10)

heatdata <- Deg_H1975S[,c(1:6)]
colnames(heatdata)=c("DL-1","DL-2","DL-3","Control-1","Control-2","Control-3")
heatdata=heatdata[,order(colnames(heatdata))]
annotation_col = data.frame(Type = factor(rep(c("Control","DL"), c(3,3))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(scale(heatdata), 
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col, 
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = F,
         show_colnames = F,
         scale = "row",
         cellwidth = 40, cellheight = 0.1,
         fontsize = 10)
dev.off()
###############YuGene#######
Deg_A549F=Deg_A549F[order(rownames(Deg_A549F)),]
Deg_A549S=Deg_A549S[order(rownames(Deg_A549S)),]
Deg_H1975F=Deg_H1975F[order(rownames(Deg_H1975F)),]
Deg_H1975S=Deg_H1975S[order(rownames(Deg_H1975S)),]
alltable=cbind(Deg_A549F[,c(1:5)],Deg_H1975F[,c(1:5)],Deg_A549S[,c(1:6)],Deg_H1975S[,c(1:6)])
library(YuGene)
YuGene.transformed <- YuGene(alltable, progressBar = TRUE)
write.csv(YuGene.transformed,"E:/??Ŀ/DL/YuGene.transformed.csv")
YuGene.transformed_1 = read.csv("E:/??Ŀ/DL/YuGene.transformed.csv",header = T)
rownames(YuGene.transformed_1)=YuGene.transformed_1[,1]
heatdata=YuGene.transformed_1[,-1]
annotation_col = data.frame(Type = factor(rep(c("Control","DL"), c(10,12))),
                            Trial=factor(c(rep(c("A549#1"),2),
                                           rep(c("H1975#1"),2),
                                           rep(c("A549#2"),3),
                                           rep(c("H1975#2"),3),
                                           rep(c("A549#1"),3),
                                           rep(c("H1975#1"),3),
                                           rep(c("A549#2"),3),
                                           rep(c("H1975#2"),3))))
rownames(annotation_col) = colnames(heatdata)
ann_colors = list(Type = c(Control="blue",DL= "firebrick"),
                  Trial = c("A549#1" = "#1B9E77", "H1975#1" = "#D95F02","A549#2" = "#7570B3", "H1975#2" = "#E7298A"))
pheatmap(scale(heatdata), 
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col,
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = F,
         show_colnames = F,
         scale = "row",
         cellwidth = 40, cellheight = 0.2,
         fontsize = 10)

################################ ????????#########
library(clusterProfiler)
genelist=mapIds(org.Hs.eg.db,Deg,keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T)
write.csv(go,"~/R/????/ʵ??/Go.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"~/R/????/ʵ??/Kegg.csv")
library(enrichplot)
dotplot(go, showCategory=30)+ ggtitle("dotplot for GO")
dotplot(kegg, showCategory=30)+ ggtitle("dotplot for KEGG")

#################### TCGA########
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
LUADRnaseqSE <- GDCdownload(query.exp.hg38)
LUADRnaseqSE <- GDCprepare(query.exp.hg38)
rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
exp.hg38.values <- assay(LUADRnaseqSE)
Exp = exp.hg38.values[,substring(colnames(exp.hg38.values),14,16)=="01A"]
Exp=Exp[which(rowSums(Exp)>0),]
cli.query <- GDCquery(project =c("TCGA-LUAD"),file.type = "xml",data.category = "Clinical") 
GDCdownload(cli.query)
cli <- GDCprepare_clinic(cli.query, clinical.info = "patient")
samples=cli[,c("bcr_patient_barcode","gender","vital_status","days_to_death","days_to_last_followup","race_list","age_at_initial_pathologic_diagnosis","tobacco_smoking_history","stage_event_pathologic_stage")]
samples$survival_time<-ifelse(samples$vital_status=='Alive',
                              samples$days_to_last_followup,samples$days_to_death)
samples=samples[!(is.na(samples$survival_time)),]
samples$status=ifelse(samples$vital_status=="Dead",1,0)
samples=samples[!(substring(samples$stage_event_pathologic_stage,7,9)==""),]
samples$stage=ifelse(substring(samples$stage_event_pathologic_stage,7,9)=="IV"|
                       substring(samples$stage_event_pathologic_stage,7,9)=="III",2,1)
samples=samples[!duplicated(samples),]
samples=samples[!(samples$race_list==""),]
samples=samples[!(is.na(samples$tobacco_smoking_history)),]
Deg=read.csv("~/Downloads/项目/DL/Deg_A549F.csv")
DegExp=as.data.frame(Exp[rownames(Exp)%in%Deg$X,])
ExpCan=as.data.frame(t(DegExp))
ExpCan=log2(ExpCan+1)
ExpCan$SAMPID=substring(rownames(ExpCan),1,12)
ExpCanCli=merge(ExpCan,samples,by.x="SAMPID",by.y="bcr_patient_barcode")
ExpCanCli=log2(ExpCanCli+1)
library(tableone)
stable <- CreateTableOne(vars=c("gender","race_list","age_at_initial_pathologic_diagnosis","tobacco_smoking_history","status"),
                         strata="stage",data=ExpCanCli,
                         factorVars=c("gender","race_list","tobacco_smoking_history","status"))
stable=print(stable,showAllLevels = TRUE)
write.csv(stable,"E:/??Ŀ/DL/table1.csv")

####################################################################### cox
#############stage1###########
library("survival")
library("survminer")
stage1=ExpCanCli[ExpCanCli$stage==1,]
for (i in 1:length(colnames(stage1[,2:1866]))){
  group=ifelse(as.numeric(as.character(stage1[,i+1]))>median(as.numeric(as.character(stage1[,i+1]))),1,0)
  stage1[,i+1]=group
}
P=list()
HR=list()
CI=matrix(ncol = 2,nrow = 1865)
colnames(CI) <- c("Lower", "Higher")
for (i in 1:length(colnames(stage1[,2:1866]))){
  fit =  coxph(Surv(survival_time, status)~as.numeric(as.character(stage1[,i+1])),data=stage1)
  HR[i] <- round(exp(coef(fit)), 2)
  CI[i,] <- round(exp(confint(fit)), 2)
  P[i] <- round(coef(summary(fit))[,5], 3)
}
table1=as.data.frame(cbind(HR, CI, P))
rownames(table1)=colnames(stage1[,2:1866])
hubGeneUp1=table1[intersect(rownames(table1[table1$P<0.05&table1$HR<1,]),DegUp),]
hubGeneDown1=table1[intersect(rownames(table1[table1$P<0.05&table1$HR>1,]),DegDown),]
############stage2#########
stage2=ExpCanCli[ExpCanCli$stage==2,]
for (i in 1:length(colnames(stage2[,2:1866]))){
  group=ifelse(as.numeric(as.character(stage2[,i+1]))>median(as.numeric(as.character(stage2[,i+1]))),1,0)
  stage2[,i+1]=group
}
P=list()
HR=list()
CI=matrix(ncol = 2,nrow = 1865)
colnames(CI) <- c("Lower", "Higher")
for (i in 1:length(colnames(stage2[,2:1866]))){
  fit =  coxph(Surv(survival_time, status)~as.numeric(as.character(stage2[,i+1])),data=stage2)
  HR[i] <- round(exp(coef(fit)), 2)
  CI[i,] <- round(exp(confint(fit)), 2)
  P[i] <- round(coef(summary(fit))[,5], 3)
}
table2=as.data.frame(cbind(HR, CI, P))
rownames(table2)=colnames(stage2[,2:1866])
hubGeneUp2=table2[intersect(rownames(table2[table2$P<0.05&table2$HR<1,]),DegUp),]
hubGeneDown2=table2[intersect(rownames(table2[table2$P<0.05&table2$HR>1,]),DegDown),]
###########
hubGeneUp=intersect(rownames(hubGeneUp2),rownames(hubGeneUp1))
hubGeneDown=intersect(rownames(hubGeneDown1),rownames(hubGeneDown2))
fit1 <- coxph(Surv(survival_time,status)~FZD3, data=stage1)
summary(fit1)
ggsurvplot(survfit(Surv(survival_time,status)~FZD3,data=stage1),data=stage1,pval = T,linetype = c('solid', 'dashed'),palette=c("red","black"))
ggsurvplot(survfit(Surv(survival_time,status)~FZD3,data=stage2),data=stage2,pval = T,linetype = c('solid', 'dashed'),palette=c("red","black"))
ggsurvplot(survfit(Surv(survival_time,status)~SCAND2P,data=stage1),data=stage1,pval = T,linetype = c('solid', 'dashed'),palette=c("red","black"))
ggsurvplot(survfit(Surv(survival_time,status)~SCAND2P,data=stage2),data=stage2,pval = T,linetype = c('solid', 'dashed'),palette=c("red","black"))
ggsurvplot(survfit(Surv(survival_time,status)~MTURN,data=stage1),data=stage1,pval = T,linetype = c('solid', 'dashed'),palette=c("red","black"))
ggsurvplot(survfit(Surv(survival_time,status)~MTURN,data=stage2),data=stage2,pval = T,linetype = c('solid', 'dashed'),palette=c("red","black"))
ggsurvplot(survfit(Surv(survival_time,status)~PRC1,data=stage1),data=stage1,pval = T,linetype = c('solid', 'dashed'),palette=c("red","black"))
ggsurvplot(survfit(Surv(survival_time,status)~PRC1,data=stage2),data=stage2,pval = T,linetype = c('solid', 'dashed'),palette=c("red","black"))
library(forestplot)
stage1$race_list=ifelse(stage1$race_list=="WHITE","White","Others")
stage1 <- within(stage1, {
  gender <- factor(gender, labels = c('Female', 'Male'))
  tobacco_smoking_history <- factor(tobacco_smoking_history , 
                                    labels = c('Lifelong Non-smoker', 'Current smoker',
                                               "Current reformed smoker for > 15 years",
                                               "Current reformed smoker for < 15 years",
                                               "Current reformed smoker, duration not specified"))
})
colnames(stage1)[c(1872,1873)]=c("age","smoking_status")
f=coxph(Surv(survival_time,status)~FZD3+gender+age+race_list+smoking_status,data = stage1)
ggforest(f,data = stage1,cpositions = c(0.05, 0.15, 0.35),fontsize = 1.2, 
         refLabel = 'reference',noDigits = 2)
f=coxph(Surv(survival_time,status)~SCAND2P+gender+age+race_list+smoking_status,data = stage1)
ggforest(f,data = stage1,cpositions = c(0.05, 0.15, 0.35),fontsize = 1.2, 
         refLabel = 'reference',noDigits = 2)
f=coxph(Surv(survival_time,status)~MTURN+gender+age+race_list+smoking_status,data = stage1)
ggforest(f,data = stage1,cpositions = c(0.05, 0.15, 0.35),fontsize = 1.2, 
         refLabel = 'reference',noDigits = 2)
f=coxph(Surv(survival_time,status)~PRC1+gender+age+race_list+smoking_status,data = stage1)
ggforest(f,data = stage1,cpositions = c(0.05, 0.15, 0.35),fontsize = 1.2, 
         refLabel = 'reference',noDigits = 2)
stage2$race_list=ifelse(stage2$race_list=="WHITE","White","Others")
stage2 <- within(stage2, {
  gender <- factor(gender, labels = c('Female', 'Male'))
  tobacco_smoking_history <- factor(tobacco_smoking_history , 
                                    labels = c('Lifelong Non-smoker', 'Current smoker',
                                               "Current reformed smoker for > 15 years",
                                               "Current reformed smoker for < 15 years"))
})
colnames(stage2)[c(1872,1873)]=c("age","smoking_status")
f=coxph(Surv(survival_time,status)~FZD3+gender+age+race_list+smoking_status,data = stage2)
ggforest(f,data = stage2,cpositions = c(0.05, 0.15, 0.35),fontsize = 1.2, 
         refLabel = 'reference',noDigits = 2)
f=coxph(Surv(survival_time,status)~SCAND2P+gender+age+race_list+smoking_status,data = stage2)
ggforest(f,data = stage2,cpositions = c(0.05, 0.15, 0.35),fontsize = 1.2, 
         refLabel = 'reference',noDigits = 2)
f=coxph(Surv(survival_time,status)~MTURN+gender+age+race_list+smoking_status,data = stage2)
ggforest(f,data = stage2,cpositions = c(0.05, 0.15, 0.35),fontsize = 1.2, 
         refLabel = 'reference',noDigits = 2)
f=coxph(Surv(survival_time,status)~PRC1+gender+age+race_list+smoking_status,data = stage2)
ggforest(f,data = stage2,cpositions = c(0.05, 0.15, 0.35),fontsize = 1.2, 
         refLabel = 'reference',noDigits = 2)

forest1=rbind(hubGeneUp1[hubGeneUp,],hubGeneDown1[hubGeneDown,])
forest2=rbind(hubGeneUp2[hubGeneUp,],hubGeneDown2[hubGeneDown,])
############### survival########
library("rms") 
library(Hmisc)
library(grid)
library(lattice)
library(Formula)
library(ggplot2)
stage1=ExpCanCli[ExpCanCli$stage==1,]
stage1$race_list=ifelse(stage1$race_list=="WHITE","White","Others")
stage2=ExpCanCli[ExpCanCli$stage==2,]
stage2$race_list=ifelse(stage2$race_list=="WHITE","White","Others")
dd<-datadist(stage1)
options(datadist="dd")
f=cph(Surv(survival_time,status)~FZD3+PRC1+MTURN,data = stage1,
      x=T,y=T,surv=T)
survival=Survival(f)
nom=nomogram(f,fun=list(function(x) survival(1095, x),
                        function(x) survival(1825, x)),
             funlabel=c("3-year Overall Survival","5-year Overall Survival"),lp=T,
             maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
plot(nom,xfrac=.8,lplabel="Linear Predictor",varname.label=T,varname.label.sep="=",ia.space=.2,
     tck=NA,tcl=-0.20,lmap=0.3,points.label="Points",total.points.label="Total Points",
     total.sep.page=F,cap.labels=F,cex.var=1.6,cex.axis=1.05,lwd=5,label.every=1,
     col.grid=gray(c(0.8,0.95)))
fp<-predict(f)
cindex1=1-rcorr.cens(fp,Surv(stage1$survival_time,stage1$status))
library(nomogramEx)
riskFormul = nomogramEx(nomo=nom,np=2,digit=9)
riskScore=c()
for(i in 1:nrow(stage1)){
  s1 = 0 * stage1[i,"FZD3"] ^2 + -2.289307739 * stage1[i,"FZD3"] + 14.880500305
  s2 = 0 * stage1[i,"PRC1"] ^2 + 1.219272413 * stage1[i,"PRC1"] + 0 
  s3 = -3.846153846 * stage1[i,"MTURN"] + 100 
  riskScore[i] = s1+s2+s3
}
stage1$riskScoreGroup = ifelse(riskScore>median(riskScore),1,0)
stage1$survival_time=stage1$survival_time/365
fit1 <- coxph(Surv(survival_time,status)~riskScoreGroup, data=stage1)
summary(fit1)
ggsurvplot(survfit(Surv(stage1$survival_time,stage1$status)~stage1$riskScoreGroup),data=stage1,pval = T,linetype = c('solid', 'dashed'),palette=c("black","red"))
stage1 <- within(stage1, {
  gender <- factor(gender, labels = c('Female', 'Male'))
  riskScoreGroup <- factor(riskScoreGroup , labels = c('Low-risk score', 'High-risk score'))
  tobacco_smoking_history <- factor(tobacco_smoking_history , 
                                    labels = c('Lifelong Non-smoker', 'Current smoker',
                                               "Current reformed smoker for > 15 years",
                                               "Current reformed smoker for < 15 years",
                                               "Current reformed smoker, duration not specified"))
})
colnames(stage1)[c(1872,1873)]=c("age","smoking_status")
model <- coxph( Surv(survival_time,status) ~  riskScoreGroup+ gender+
                  age+smoking_status+race_list , data =  stage1 )
ggforest(model,data = stage1,cpositions = c(0.05, 0.15, 0.43),fontsize = 1.2, 
         refLabel = 'reference',noDigits = 2)

dd<-datadist(stage2)
options(datadist="dd")
f=cph(Surv(survival_time,status)~FZD3+PRC1+MTURN,data = stage2,
      x=T,y=T,surv=T)
survival=Survival(f)
nom=nomogram(f,fun=list(function(x) survival(1095, x),
                        function(x) survival(1825, x)),
             funlabel=c("3-year Overall Survival","5-year Overall Survival"),lp=T,
             maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
plot(nom,xfrac=.8,lplabel="Linear Predictor",varname.label=T,varname.label.sep="=",ia.space=.2,
     tck=NA,tcl=-0.20,lmap=0.3,points.label="Points",total.points.label="Total Points",
     total.sep.page=F,cap.labels=F,cex.var=1.6,cex.axis=1.05,lwd=5,label.every=1,
     col.grid=gray(c(0.8,0.95)))
fp<-predict(f)
cindex1=1-rcorr.cens(fp,Surv(stage2$survival_time,stage2$status))
riskFormul = nomogramEx(nomo=nom,np=2,digit=9)
riskScore=c()
for(i in 1:nrow(stage2)){
  s1 = -22.222222222 * stage2[i,"FZD3"] + 100
  s2 = 0 * stage2[i,"PRC1"] ^2 + 1.298157901 * stage2[i,"PRC1"] + 0 
  s3 = 0 * stage2[i,"MTURN"] ^3 + 0 * stage2[i,"MTURN"] ^2 + -0.194141945 * stage2[i,"MTURN"] + 11.648516696 
  riskScore[i] = s1+s2+s3
}
stage2$riskScoreGroup = ifelse(riskScore>median(riskScore),1,0)
stage2$survival_time=stage2$survival_time/365
fit1 <- coxph(Surv(survival_time,status)~riskScoreGroup, data=stage2)
summary(fit1)
ggsurvplot(survfit(Surv(stage2$survival_time,stage2$status)~stage2$riskScoreGroup),data=stage2,pval = T,linetype = c('solid', 'dashed'),palette=c("black","red"))
stage2 <- within(stage2, {
  gender <- factor(gender, labels = c('Female', 'Male'))
  riskScoreGroup <- factor(riskScoreGroup , labels = c('Low-risk score', 'High-risk score'))
  tobacco_smoking_history <- factor(tobacco_smoking_history , 
                                    labels = c('Lifelong Non-smoker', 'Current smoker',
                                               "Current reformed smoker for > 15 years",
                                               "Current reformed smoker for < 15 years"))
})
colnames(stage2)[c(1872,1873)]=c("age","smoking_status")
model <- coxph( Surv(survival_time,status) ~  riskScoreGroup+ gender+
                  age+smoking_status+race_list , data =  stage2 )
ggforest(model,data = stage2,cpositions = c(0.05, 0.15, 0.43),fontsize = 1.2, 
         refLabel = 'reference',noDigits = 2)

############### heatmap########
library(pheatmap)
heatmap=read.csv("heatmap.csv",header = T)
rownames(heatmap)=heatmap[,1]
heatmap=heatmap[,-1]
annotation_col = data.frame(Type = factor(rep(c("HR","logFC"), c(2,4))))
colnames(heatmap)=c("Early stage","Advanced stage",'A549#1',"A549#2","H1975#1","H1975#2")
rownames(annotation_col) = colnames(heatmap)
pheatmap(scale(heatmap), 
         cluster_rows = F,
         cluster_cols = T,
         annotation_col = annotation_col,
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = T,
         show_colnames = T,
         scale = "row")

##########
library(GEOquery)
gset=getGEO("GSE13213",GSEMatrix = T,AnnotGPL = T)
expSet=exprs(gset[[1]])
pdata=pData(gset[[1]])
fdata=fData(gset[[1]])
symbol=fdata[,c(1,3)]
expSet=data.frame(expSet)
expSet$ID=rownames(expSet)
expSet=merge(expSet,symbol,by.x = "ID",by.y = "ID")
boxplot(expSet[,c(2:118)],col="blue")
expSet[,c(2:118)]=normalizeBetweenArrays(expSet[,c(2:118)])
boxplot(expSet[,c(2:118)],col="blue")
expSet<-expSet[!(expSet$`Gene symbol`==""), ]
library(dplyr)
expSet<- aggregate(expSet[,c(2:118)],by=list(expSet$`Gene symbol`),FUN=max)
rownames(expSet)=expSet$Group.1
expSet=expSet[,-1]
hubGEO=expSet[c("FZD3","MTURN","PRC1"),]
hubGEO=hubGEO+5
samples=pdata[,c("geo_accession","Survival (days):ch1","Status:ch1","Stage (Pathological ):ch1","Smoking (BI):ch1","Sex:ch1","Age:ch1","Histology:ch1")]
samples=samples[!(is.na(samples$`Survival (days):ch1`)),]
hubgene=as.data.frame(t(hubGEO))
hubgene$patient=rownames(hubgene)
hubCli=merge(hubgene,samples,by.x="patient",by.y="geo_accession")
hubCli$`Survival (days):ch1`=as.numeric(as.character(hubCli$`Survival (days):ch1`))
hubCli$status=ifelse(hubCli$`Status:ch1`=="Alive",0,1)
hubCli[is.na(hubCli)]=0.000001
hubCli$stage=ifelse(substring(hubCli$`Stage (Pathological ):ch1`,1,3)=="III",2,1)
stage1=hubCli[hubCli$stage==1,]
stage2=hubCli[hubCli$stage==2,]
riskScore=c()
for(i in 1:nrow(stage1)){
  s1 = 0 * stage1[i,"FZD3"] ^2 + -2.289307739 * stage1[i,"FZD3"] + 14.880500305
  s2 = 0 * stage1[i,"PRC1"] ^2 + 1.219272413 * stage1[i,"PRC1"] + 0 
  s3 = -3.846153846 * stage1[i,"MTURN"] + 100 
  riskScore[i] = s1+s2+s3
}
stage1$riskScoreGroup = ifelse(riskScore>median(riskScore),1,0)
stage1$`Survival (days):ch1`=stage1$`Survival (days):ch1`/365
fit1 <- coxph(Surv(`Survival (days):ch1`,status)~riskScoreGroup, data=stage1)
summary(fit1)
ggsurvplot(survfit(Surv(`Survival (days):ch1`,status)~riskScoreGroup,data=stage1),data=stage1,pval = T,linetype = c('solid', 'dashed'),palette=c("black","red"))
library("survRM2")
obj=rmst2(stage1$`Survival (days):ch1`, stage1$status, stage1$riskScoreGroup)
plot(obj, xlab="Years", ylab="Probability")
riskScore=c()
for(i in 1:nrow(stage2)){
  s1 = -22.222222222 * stage2[i,"FZD3"] + 100
  s2 = 1.298157901 * stage2[i,"PRC1"] + 0 
  s3 = -0.194141945 * stage2[i,"MTURN"] + 11.648516696 
  riskScore[i] = s1+s2+s3
}
stage2$riskScoreGroup = ifelse(riskScore>median(riskScore),1,0)
stage2$`Status:ch1`[which(stage2$`Survival (days):ch1`>1000)] <- 0
stage2$`Survival (days):ch1`[which(stage2$`Survival (days):ch1`>1000)] <- 1000
fit1 <- coxph(Surv(`Survival (days):ch1`,status)~riskScoreGroup, data=stage2)
summary(fit1)
ggsurvplot(survfit(Surv(`Survival (days):ch1`,status)~riskScoreGroup,data=stage2),data=stage2,pval = T,linetype = c('solid', 'dashed'),palette=c("black","red"))
##########
library(clusterProfiler)
SYMBOL=c("FZD3","MTURN","PRC1")
genelist=mapIds(org.Hs.eg.db,SYMBOL,keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T)
write.csv(go,"E:/??Ŀ/DL/Go.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"E:/??Ŀ/DL/Kegg.csv")
library(enrichplot)
dotplot(go, showCategory=10)
dotplot(kegg, showCategory=10)
cnetplot(go, showCategory = 20)
cnetplot(kegg, showCategory = 20)
