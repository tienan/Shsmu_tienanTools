library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
query.exp.hg38 <- GDCquery(project = "TCGA-STAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts")
STADRnaseqSE <- GDCdownload(query.exp.hg38)
STADRnaseqSE <- GDCprepare(query.exp.hg38)
rownames(STADRnaseqSE) <- values(STADRnaseqSE)$external_gene_name
exp.hg38.values <- assay(STADRnaseqSE)
Exp = exp.hg38.values[,substring(colnames(exp.hg38.values),14,16)=="01A"|
                        substring(colnames(exp.hg38.values),14,16)=="11A"]
sign =c()
patient_id = colnames(Exp)
for (i in 1:length(colnames(Exp))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "-")
  patient_id[i] = tmp
}
colnames(Exp)=patient_id
ExpCan = Exp[,sign == "01A"]
ExpNor = Exp[,sign == "11A"]
SAMPID=colnames(ExpCan)
ExpCan=rbind(ExpCan,SAMPID)
ExpCan=data.frame(t(ExpCan))
SAMPID=colnames(ExpNor)
ExpNor=rbind(ExpNor,SAMPID)
ExpNor=data.frame(t(ExpNor))
cli.query <- GDCquery(project =c("TCGA-STAD"),file.type = "xml",
                      data.category = "Clinical") 
GDCdownload(cli.query)
cli <- GDCprepare_clinic(cli.query, clinical.info = "patient")
samples=cli[,c("bcr_patient_barcode","gender","vital_status","days_to_death",
               "days_to_last_followup","race_list","age_at_initial_pathologic_diagnosis",
               "stage_event_pathologic_stage")]
samples=samples[!(is.na(samples$days_to_last_followup)&is.na(samples$days_to_death)),]
samples[samples$vital_status=="Dead","survival_time"]=samples[samples$vital_status
                                                              =="Dead","days_to_death"]
samples[samples$vital_status=="Alive","survival_time"]=samples[samples$vital_status
                                                               =="Alive","days_to_last_followup"]
samples$status=ifelse(samples$vital_status=="Dead",1,0)
samples=samples[!(substring(samples$stage_event_pathologic_stage,7,9)==""),]
samples$stage=ifelse(substring(samples$stage_event_pathologic_stage,7,9)=="IV"|
                       substring(samples$stage_event_pathologic_stage,7,9)=="III",2,1)
samples=samples[!duplicated(samples),]
samples=samples[!(samples$race_list==""),]
ExpCanCli=merge(ExpCan,samples,by.x="SAMPID",by.y="bcr_patient_barcode")
ExpCanCli$SAMPID=paste(ExpCanCli$SAMPID,"01A",sep = "-")
ExpNorCli=merge(ExpNor,samples,by.x="SAMPID",by.y="bcr_patient_barcode")
ExpNorCli$SAMPID=paste(ExpNorCli$SAMPID,"11A",sep = "-")
ExpCli=rbind(ExpCanCli,ExpNorCli)
ExpCli$group=ifelse(substring(ExpCli$SAMPID,14,16)=="01A",1,0)
library(tableone)
stable <- CreateTableOne(vars=c("gender","race_list","age_at_initial_pathologic_diagnosis"),
                         strata="group",data=ExpCli,
                         factorVars=c("gender","race_list"))
stable=print(stable,showAllLevels = TRUE)
write.csv(stable,"E:/项目/STAD/table1.csv")
rownames(ExpCli)=ExpCli[,1]
ExpCli=t(ExpCli[,-1])
ExpM=matrix(as.numeric(unlist(ExpCli[1:56493,])),nrow=nrow(ExpCli[1:56493,]), dimnames = list(rownames(ExpCli[1:56493,]),colnames(ExpCli[1:56493,])))
library(DESeq2)
qualifiedExp=c()
for (genes_in_sheet in rownames(ExpM)) {
  qualification = ExpM[genes_in_sheet,] <= 10
  if (sum(qualification) < 0.8*length(ExpM)) {
    qualifiedExp = append(qualifiedExp,genes_in_sheet)
  }
}##去除在80%样本以上counts<10的基因
ExpM=ExpM[qualifiedExp,]
condition <- factor(ExpCli[56504,])
colData <- data.frame(row.names =colnames(ExpM),condition)
colData$condition=relevel(colData$condition,ref = "0")
dds <-DESeqDataSetFromMatrix(ExpM,colData,design=~condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, "condition")
expVst <- as.data.frame(assay(vsd))
par(cex=0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(expVst)
cols=rainbow(n.sample*1.2)
par(mfrow=c(2,1))
boxplot(expVst,col=cols,las=2)
res <- results(dds, tidy=TRUE) 
res=na.omit(res)
Deg=res[res$padj<0.05&abs(res$log2FoldChange)>1,]
write.csv(Deg,"E:/项目/STAD/Deg.csv")

library(ggplot2)
library(RColorBrewer)
res$threshold[res$padj<0.05&res$log2FoldChange>1]="up"
res$threshold[res$padj<0.05&res$log2FoldChange<(-1)]="down"
res$threshold[(res$padj>0.05)|(res$log2FoldChange<=1)&res$log2FoldChange>=(-1)]="normal"
res$threshold=factor(res$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(res$log2FoldChange,-res$log2FoldChange)
theme_set(theme_bw())
p <- ggplot(res,aes(log2FoldChange,-1*log10(padj),
                             color =threshold ))+geom_point()+
  xlim(-10,10) +  labs(x="log2(FoldChange)",y="-log10(padj)")
p <- p + scale_color_manual(values =c('up'="red","normal"="grey","down"="blue"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-log2(2),log2(2)),linetype=4)
p <- p +theme(plot.title = element_text(size = 25,face = 'bold', vjust = 0.5, hjust = 0.5),
              legend.title = element_blank(),
              legend.text = element_text(size = 18, face = 'bold'),
              legend.position = 'right',
              legend.key.size=unit(0.8,'cm'),
              axis.ticks.x=element_blank(),
              axis.text.x=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.text.y=element_text(size = 15,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.title.x = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
              axis.title.y = element_text(size = 20,face = 'bold', vjust = 0.5, hjust = 0.5),
              panel.background = element_rect(fill = 'transparent',colour = 'black'),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_rect(fill = 'transparent',colour = 'black'))
p

library(pheatmap)
heatdata <- expVst[Deg$row,]
group <- ifelse(substring(colnames(heatdata),14,16)=="01A","cancer","normal")
annotation_col = data.frame(Type = factor(rep(c("cancer","normal"), c(297,24))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(t(scale(t(heatdata))),cluster_rows = T,show_rownames = F,show_colnames = F,
         cluster_cols = T,annotation_col = annotation_col,
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         legend = T,legend_breaks = c(-2,0,2),breaks = unique(seq(-2,2,length=100)),
         cutree_rows = 2,cutree_cols = 2,
         lagend_labels = c("≤2","0","≥2"))
library("survival")
library("survminer")
DegExp=expVst[rownames(expVst)%in%Deg$row,which(substring(colnames(expVst),14,16)=="01A")]
surExp=matrix(nrow = 10971,ncol = 297)
for (i in 1:length(rownames(DegExp))){
  group=ifelse(DegExp[i,]>median(as.numeric(DegExp[i,])),1,0)
  surExp[i,]=group
}
colnames(surExp)=colnames(expVst[,1:297])
rownames(surExp)=rownames(DegExp)
surExp=rbind(surExp,ExpCli[56494:56504,1:297])
surExp=as.data.frame(t(surExp))
P=list()
HR=list()
CI=matrix(ncol = 2,nrow = 10971)
colnames(CI) <- c("Lower", "Higher")
for (i in 1:10971){
  fit =  coxph(Surv(as.numeric(survival_time), as.numeric(status))~surExp[,i],data=surExp)
  HR[i] <- round(exp(coef(fit)), 2)
  CI[i,] <- round(exp(confint(fit)), 2)
  P[i] <- round(coef(summary(fit))[,5], 3)
}
surTab=as.data.frame(cbind(HR, CI, P))
rownames(surTab)=colnames(surExp[,1:10971])
surTab$HR=unlist(surTab$HR)
surTab$Lower=unlist(surTab$Lower)
surTab$Higher=unlist(surTab$Higher)
surTab$P=unlist(surTab$P)
surGene=surTab[surTab$P<0.05,]
write.csv(surGene,"E:/项目/STAD/surGene.csv")

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
genelist=mapIds(org.Hs.eg.db,rownames(surGene),keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T,ont="ALL")
write.csv(go,"E:/项目/STAD/GoSur.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"E:/项目/STAD/KeggSur.csv")
library(enrichplot)
dotplot(go, showCategory=30)+ ggtitle("dotplot for GO")
dotplot(kegg, showCategory=30)+ ggtitle("dotplot for KEGG")

query <- GDCquery(project = "TCGA-STAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "miRNA Expression Quantification", 
                  workflow.type = "BCGSC miRNA Profiling")
GDCdownload(query, 
            method = "api", 
            files.per.chunk = 20)
mir_exp<- GDCprepare(query = query,
                     summarizedExperiment=F)
rownames(mir_exp)=mir_exp[,1]
mirExp=mir_exp[,substring(colnames(mir_exp),6,10)=="count"]
sign =c()
patient_id = colnames(mirExp)
for (i in 1:length(colnames(mirExp))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "_"))
  patient_id[i] = tmp[3]
}
for (i in 1:length(patient_id)){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  tmp = paste(tmp[1],tmp[2],tmp[3],tmp[4],sep = "-")
  patient_id[i] = tmp
}
colnames(mirExp)=patient_id
mirCan = mirExp[,sign == "01A"]
mirNor = mirExp[,sign == "11A"]
SAMPID=colnames(mirCan)
mirCan=rbind(mirCan,SAMPID)
mirCan=data.frame(t(mirCan))
SAMPID=colnames(mirNor)
mirNor=rbind(mirNor,SAMPID)
mirNor=data.frame(t(mirNor))
mirCanCli=merge(mirCan,samples,by.x="X1882",by.y="bcr_patient_barcode")
mirCanCli$X1882=paste(mirCanCli$X1882,"01A",sep = "-")
mirNorCli=merge(mirNor,samples,by.x="X1882",by.y="bcr_patient_barcode")
mirNorCli$X1882=paste(mirNorCli$X1882,"11A",sep = "-")
mirCli=rbind(mirCanCli,mirNorCli)
mirCli$group=ifelse(substring(mirCli$X1882,14,16)=="01A",1,0)
stable1 <- CreateTableOne(vars=c("gender","race_list","age_at_initial_pathologic_diagnosis"),
                         strata="group",data=mirCli,
                         factorVars=c("gender","race_list"))
stable1=print(stable1,showAllLevels = TRUE)
write.csv(stable1,"E:/项目/STAD/table1.csv")
rownames(mirCli)=mirCli[,1]
mirCli=mirCli[,-1]
mirM=t(mirCli[,1:1881])
mirM=matrix(as.numeric(unlist(mirM)),nrow=nrow(mirM), dimnames = list(rownames(mirM),colnames(mirM)))
qualifiedExp=c()
for (genes_in_sheet in rownames(mirM)) {
  qualification = mirM[genes_in_sheet,] <= 10
  if (sum(qualification) < 0.8*length(mirM)) {
    qualifiedExp = append(qualifiedExp,genes_in_sheet)
  }
}
mirM=mirM[qualifiedExp,]
condition <- factor(mirCli[,1892])
colData <- data.frame(row.names =colnames(ExpM),condition)
colData$condition=relevel(colData$condition,ref = "0")
dds <-DESeqDataSetFromMatrix(mirM,DataFrame(condition),design=~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsd, "condition")
mirVst <- as.data.frame(assay(vsd))
res <- results(dds, tidy=TRUE) 
res=na.omit(res)
mirDeg=subset(res,padj<0.05&abs(res$log2FoldChange)>1)
write.csv(mirDeg,"E:/项目/STAD/mirDeg.csv")
heatdata=mirVst[mirDeg$row,]
group <- ifelse(substring(colnames(heatdata),14,16)=="01A","cancer","normal")
annotation_col = data.frame(Type = factor(rep(c("cancer","normal"), c(343,32))))
rownames(annotation_col) = colnames(heatdata)
pheatmap(t(scale(t(heatdata))),cluster_rows = T,show_rownames = F,show_colnames = F,
         cluster_cols = T,annotation_col = annotation_col,
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         legend = T,legend_breaks = c(-2,0,2),breaks = unique(seq(-2,2,length=100)),
         cutree_rows = 2,cutree_cols = 2,
         lagend_labels = c("≤2","0","≥2"))
