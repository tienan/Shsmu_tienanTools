getwd()
setwd("../R/")
DL=read.csv("DL.csv")
CAN=read.csv("CAN.csv")
library(biomaRt)


library(org.Hs.eg.db)

citation(org.Hs.eg.db)

DL_gene_id = mapIds(org.Hs.eg.db,keys =as.character(DL$gene),keytype="SYMBOL",column="ENTREZID")
DL$gene_id = DL_gene_id
DL$threshold[DL$q_value<0.1&(DL$log2.fold_change.)>=0.5]="up"
DL$threshold[DL$q_value<0.1&DL$log2.fold_change.<=-0.5]="down"
DL$threshold[(DL$q_value>0.1)|(DL$log2.fold_change.<0.5)&DL$log2.fold_change.>-0.5]="normal"
DL$threshold=factor(DL$threshold,levels = c("up","down","normal"),ordered = T)
DL_Dif=DL[DL$q_value<0.1,c(2,3,10,12,14)]
DL_Dif
head(DL_Dif)
nrow(DL_Dif[DL_Dif$threshold=="up",])
nrow(DL_Dif[DL_Dif$threshold=="down",])

library()
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T)
write.csv(go,"~/Downloads/DL_CAN/Go.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"~/R/TC/HUB GENE/Kegg.csv")
library(enrichplot)
dotplot(go, showCategory=30)+ ggtitle("dotplot for GO")
dotplot(kegg, showCategory=30)+ ggtitle("dotplot for KEGG")


CAN$threshold[CAN$q_value<0.05&CAN$log2.fold_change.>=0.5]="up"
CAN$threshold[CAN$q_value<0.05&CAN$log2.fold_change.<=-0.5]="down"
CAN$threshold[(CAN$q_value>0.05)|(CAN$log2.fold_change.<0.5)&CAN$log2.fold_change.>-0.5]="normal"
CAN$threshold=factor(CAN$threshold,levels = c("up","down","normal"),ordered = T)

CAN_Dif=CAN[CAN$q_value<0.05,c(2,3,10,12,14)]

CAN_Dif$gene_id = mapIds(org.Hs.eg.db,keys =as.character(CAN_Dif$gene),keytype="REFSEQ",column="ENTREZID")


DL_U_Can_D=merge(DL_Dif[DL_Dif$log2.fold_change.>0,],CAN_Dif[CAN_Dif$log2.fold_change.<0,],by.x = "gene_id",by.y = "gene_id")
DL_D_Can_U=merge(DL_Dif[DL_Dif$log2.fold_change.<0,],CAN_Dif[CAN_Dif$log2.fold_change.>0,],by.x = "gene_id",by.y = "gene_id")
DL_Can=rbind(DL_U_Can_D,DL_D_Can_U)
DL_Can_1 = na.omit(DL_Can)
write.csv(DL_Can_1 ,"DL_Can_1.csv")
library("pheatmap")
rownames(DL_Can_1) = DL_Can_1[,2]
pheatmap(DL_Can_1[,c(3,7)], 
         cluster_rows = F,
         cluster_cols = F,
         annotation_legend=TRUE, 
         color =colorRampPalette(c("blue", "white","red"))(100),
         scale = "column",
         show_rownames = T,
         show_colnames = F,
         #display_numbers = T,
         number_color = "black",
        #cellwidth = 70, cellheight = 15,
         fontsize = 15)


library(org.Hs.eg.db)
library(dplyr)
CAN_gene_id = mapIds(org.Hs.eg.db,keys =as.character(CAN$gene),keytype="REFSEQ",column="ENTREZID")
# keytypes()
CAN$gene = CAN_gene_id 
mapped_seqs <- mappedkeys(x)
head(mapped_seqs)

# ensembl = useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
# ipro1 = getBM(attributes=c("refseq_mrna","hgnc_symbol"), 
#              filters="refseq_mrna",
#              values=CAN$gene, 
#              mart=ensembl)
# ipro2 = getBM(attributes=c("refseq_ncrna","hgnc_symbol"), 
#              filters="refseq_ncrna",
#              values=CAN$gene, 
#              mart=ensembl)
# colnames(ipro1)=c("ref","symbol")
# colnames(ipro2)=c("ref","symbol")
# ipro=rbind(ipro1,ipro2)
# CAN=merge(CAN,ipro,by.x = "gene",by.y = "ref")
# CAN=CAN[!CAN$symbol=="",]


library(ggplot2)
library(RColorBrewer)
DL$threshold[DL$q_value<0.1&(DL$log2.fold_change.)>=0.5]="up"
DL$threshold[DL$q_value<0.1&DL$log2.fold_change.<=-0.5]="down"
DL$threshold[(DL$q_value>0.1)|(DL$log2.fold_change.<0.5)&DL$log2.fold_change.>-0.5]="normal"
DL$threshold=factor(DL$threshold,levels = c("up","down","normal"),ordered = T)
DL_Dif=DL[DL$q_value<0.1,c(3,10,12)]
write.csv(DL_Dif,"DL_Dif.csv")
BiocManager::install("clusterProfiler")
BiocManager::install("DO.db")
library(clusterProfiler)
library(org.Hs.eg.db)
#genelist=mapIds(org.Hs.eg.db,as.character(DL_Can$gene),keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=DL_Dif$gene_id,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T)

go@result

dotplot(go)+ ggtitle("Dotplot for GO")

kegg@result
kegg=enrichKEGG(DL_Dif$gene_id,organism = "hsa",pvalueCutoff = 0.3)
dotplot(kegg, showCategory=30)+ ggtitle("Dotplot for KEGG")



write.csv(go,"~/Downloads/DL_CAN/Go.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"~/R/TC/HUB GENE/Kegg.csv")
BiocManager::install("enrichplot")
library(enrichplot)
dotplot(go, showCategory=30)+ ggtitle("dotplot for GO")
dotplot(kegg, showCategory=30)+ ggtitle("dotplot for KEGG")




x_lim=max(DL$log2.fold_change.,-DL$log2.fold_change.)
theme_set(theme_bw())
p <- ggplot(DL,aes(log2.fold_change.,-1*log10(q_value),
                        color =threshold ))+geom_point()+
  xlim(-8,8) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
p <- p + scale_color_manual(values =c('up'="red","normal"="grey","down"="blue"))+
  geom_hline(yintercept=-log10(0.1),linetype=4)+
  geom_vline(xintercept=c(-0.5,0.5),linetype=4)
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
CAN$threshold[CAN$q_value<0.05&CAN$log2.fold_change.>=0.5]="up"
CAN$threshold[CAN$q_value<0.05&CAN$log2.fold_change.<=-0.5]="down"
CAN$threshold[(CAN$q_value>0.05)|(CAN$log2.fold_change.<0.5)&CAN$log2.fold_change.>-0.5]="normal"
CAN$threshold=factor(CAN$threshold,levels = c("up","down","normal"),ordered = T)
x_lim=max(CAN$log2.fold_change.,-CAN$log2.fold_change.)
theme_set(theme_bw())
p <- ggplot(CAN,aes(CAN$log2.fold_change.,-1*log10(CAN$q_value),
                   color =threshold ))+geom_point()+
  xlim(-16,16) +  labs(x="log2(FoldChange)",y="-log10(FDR)")
p <- p + scale_color_manual(values =c('up'="red","normal"="grey","down"="blue"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-0.5,0.5),linetype=4)
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


#DL_Dif=DL[DL$q_value<0.1,c(3,10,12)]
CAN_Dif=CAN[CAN$q_value<0.05,c(14,10,12)]
DL_U_Can_D=merge(DL_Dif[DL_Dif$log2.fold_change.>0,],CAN_Dif[CAN_Dif$log2.fold_change.<0,],by.x = "gene",by.y = "symbol")
DL_D_Can_U=merge(DL_Dif[DL_Dif$log2.fold_change.<0,],CAN_Dif[CAN_Dif$log2.fold_change.>0,],by.x = "gene",by.y = "symbol")
DL_Can=rbind(DL_U_Can_D,DL_D_Can_U)
write.csv(DL_Can,"~/Downloads/DL_CAN/DL_Can.csv")
library(clusterProfiler)
library(org.Hs.eg.db)
genelist=mapIds(org.Hs.eg.db,as.character(DL_Can$gene),keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05,readable = T)
write.csv(go,"~/Downloads/DL_CAN/Go.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
kegg=setReadable(kegg,org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegg,"~/R/TC/HUB GENE/Kegg.csv")
library(enrichplot)
dotplot(go, showCategory=30)+ ggtitle("dotplot for GO")
dotplot(kegg, showCategory=30)+ ggtitle("dotplot for KEGG")

rownames(DL_Can)=DL_Can[,1]
DL_Can=DL_Can[,-1]
colnames(DL_Can)=c("DL_A549","q1","Cancer_Noncancer","q2")
library(pheatmap)
pheatmap(DL_Can[,c(1,3)], 
         cluster_rows = F,
         cluster_cols = F,
         annotation_legend=TRUE, 
         color =colorRampPalette(c("blue", "white","red"))(100),
         show_rownames = T,
         show_colnames = T,
         display_numbers = T,
         number_color = "black",
         cellwidth = 70, cellheight = 15,
         fontsize = 15)
