t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
write.csv(x=tmp,file =  paste(fileName,".csv",sep = "",collapse = ""))
write.csv(x=fgseaRes[fgseaRes$pval<0.1,c(1:5)],file = paste(fileName,"Gsea.csv",sep = "",collapse = ""))
write.csv(x=pathway,file = paste(fileName,"enrichKEGG.csv",sep = "",collapse = ""))
tmp
fileName = "A549DConvsPM2.5"
t1 = edgeR::DGEList(  datDLExp [,c(1:3,10:12)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
write.csv(x=tmp,file =  paste(fileName,".csv",sep = "",collapse = ""))
write.csv(x=fgseaRes[fgseaRes$pval<0.1,c(1:5)],file = paste(fileName,"Gsea.csv",sep = "",collapse = ""))
write.csv(x=pathway,file = paste(fileName,"enrichKEGG.csv",sep = "",collapse = ""))
tmp
#A549DConvsPM2.5
fileName = "A549DConvsPM2.5"
t1 = edgeR::DGEList(  datDLExp [,c(1:3,10:12)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
tableDEA$geneName
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
mapping
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp
tmp
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes
fgseaRes
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
ranks
plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
names(ranks)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
fileName = "A549dmsovsPM2.5"
t1 = edgeR::DGEList(  datDLExp [,c(10:12,7:9)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
#attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
#plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
write.csv(x=tmp,file =  paste(fileName,".csv",sep = "",collapse = ""))
write.csv(x=fgseaRes[fgseaRes$pval<0.1,c(1:5)],file = paste(fileName,"Gsea.csv",sep = "",collapse = ""))
write.csv(x=pathway,file = paste(fileName,"enrichKEGG.csv",sep = "",collapse = ""))
tmp
pathway
#1975msovsPM2.5
fileName = "1975msovsPM2.5"
t1 = edgeR::DGEList(  datDLExp [,c(22:24,19:21)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
#attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
#plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
write.csv(x=tmp,file =  paste(fileName,".csv",sep = "",collapse = ""))
write.csv(x=fgseaRes[fgseaRes$pval<0.1,c(1:5)],file = paste(fileName,"Gsea.csv",sep = "",collapse = ""))
write.csv(x=pathway,file = paste(fileName,"enrichKEGG.csv",sep = "",collapse = ""))
tmp
#1975ConvsPM2.5
fileName = "1975ConvsPM2.5"
t1 = edgeR::DGEList(  datDLExp [,c(13:15,19:21)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
#attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
#plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
write.csv(x=tmp,file =  paste(fileName,".csv",sep = "",collapse = ""))
write.csv(x=fgseaRes[fgseaRes$pval<0.1,c(1:5)],file = paste(fileName,"Gsea.csv",sep = "",collapse = ""))
write.csv(x=pathway,file = paste(fileName,"enrichKEGG.csv",sep = "",collapse = ""))
#A549DConvsPM2.5
fileName = "A549DConvsPM2.5"
t1 = edgeR::DGEList(  datDLExp [,c(1:3,7:9)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
#attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
#plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
write.csv(x=tmp,file =  paste(fileName,".csv",sep = "",collapse = ""))
write.csv(x=fgseaRes[fgseaRes$pval<0.1,c(1:5)],file = paste(fileName,"Gsea.csv",sep = "",collapse = ""))
write.csv(x=pathway,file = paste(fileName,"enrichKEGG.csv",sep = "",collapse = ""))
tmp
#1975DmsovsDL
fileName = "1975DmsovsDL"
t1 = edgeR::DGEList(  datDLExp [,c(22:24,16:18)],group = as.factor(group))
#1975DmsovsDL
fileName = "1975DmsovsDL"
t1 = edgeR::DGEList(  datDLExp [,c(22:24,16:18)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
#attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
#plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
write.csv(x=tmp,file =  paste(fileName,".csv",sep = "",collapse = ""))
write.csv(x=fgseaRes[fgseaRes$pval<0.1,c(1:5)],file = paste(fileName,"Gsea.csv",sep = "",collapse = ""))
write.csv(x=pathway,file = paste(fileName,"enrichKEGG.csv",sep = "",collapse = ""))
tmp
fileName = "A549DmsovsDL"
t1 = edgeR::DGEList(  datDLExp [,c(10:12,4:6)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
#attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
#plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
write.csv(x=tmp,file =  paste(fileName,".csv",sep = "",collapse = ""))
write.csv(x=fgseaRes[fgseaRes$pval<0.1,c(1:5)],file = paste(fileName,"Gsea.csv",sep = "",collapse = ""))
write.csv(x=pathway,file = paste(fileName,"enrichKEGG.csv",sep = "",collapse = ""))
tmp
pathway
#1975ConvsDL
fileName = "1975ConvsDL"
t1 = edgeR::DGEList(  datDLExp [,c(16:19,19:21)],group = as.factor(group))
#1975ConvsDL
fileName = "1975ConvsDL"
t1 = edgeR::DGEList(  datDLExp [,c(16:18,19:21)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
#attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
#plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
write.csv(x=tmp,file =  paste(fileName,".csv",sep = "",collapse = ""))
write.csv(x=fgseaRes[fgseaRes$pval<0.1,c(1:5)],file = paste(fileName,"Gsea.csv",sep = "",collapse = ""))
write.csv(x=pathway,file = paste(fileName,"enrichKEGG.csv",sep = "",collapse = ""))
tmp
fileName = "A549ConvsDL"
t1 = edgeR::DGEList(  datDLExp [,c(1:3,4:6)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =rownames(tableDEA)
mygenes = strsplit2(x = rownames(tableDEA),split = "[.]")[,1]
mygenes
#attributes( hsmart )
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'refseq_mrna',
values =  mygenes,
mart =ensembl
)
mapping = as.data.frame(mapping)
tmp = merge( tableDEA  ,mapping,by.x = "geneName",by.y= "refseq_mrna" )
tmp = na.omit(tmp)
ranks = tmp$logFC
names(ranks)=tmp$entrezgene_id
ranks =sort(ranks,decreasing = T)
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes[fgseaRes$pval<0.2,c(1:5,8)]
ncol(fgseaRes)
#plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)
library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.2,]
write.csv(x=tmp,file =  paste(fileName,".csv",sep = "",collapse = ""))
write.csv(x=fgseaRes[fgseaRes$pval<0.1,c(1:5)],file = paste(fileName,"Gsea.csv",sep = "",collapse = ""))
write.csv(x=pathway,file = paste(fileName,"enrichKEGG.csv",sep = "",collapse = ""))
}
tmp
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'hgnc_symbol',
values =  "PDIA3P1",
mart =ensembl
)
mapping
mapping <- getBM(
attributes = c('entrezgene_id', 'hgnc_symbol',"refseq_mrna"),
filters = 'hgnc_symbol',
values =  "PDIA3P",
mart =ensembl
)
mapping
datDLExp
rownames(datDLExp)
rownames(datDLExp)="NR_002305.1"
datDLExp[rownames(datDLExp)=="NR_002305.1",]
datDl
head( datDl)
datDLExp[rownames(datDLExp)=="NR_002305.1",]
datDl
head( datDl)
datDLExp[rownames(datDLExp)=="NR_002305",]
datDLExp[rownames(datDLExp)=="NR_002305.1",]
dat=read.csv("ExpressDLThird.csv")
head(dat)



##############################################################
getwd()
dat=read.csv("/media/tienan/0000678400004823/R/DL/ExpressDLThird.csv")
head(dat)
library(edgeR)
dat_1 = dat[,-1]
rownames(dat_1) = dat[,1]
name = strsplit2(colnames(dat)[-1],split = "[0-9]")
y = DGEList(dat_1,group = name,genes = dat_1)
y
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")
library(org.Mm.eg.db)
library(org.Hs.eg.db)
y$genes$Symbol = mapIds(org.Hs.eg.db,rownames(y),keytype = "ENTREZID",column = "SYMBOL")
group=factor(name)
table(group)
y = y[!is.na(y$genes$Symbol ),]
dim(y)
keep 

y=calcNormFactors(y)
y$samples
group
pch<- c(0,1,3,2,6,15,16,17)
colors<- rep(c("darkgreen","red","blue"),2)
plotMDS(y,col=colors[group],pch=pch[group])
legend("topright",legend=levels(group),pch=pch,col=colors,ncol=2)
design<-model.matrix(~0+group)
colnames(design)<-levels(group)
design
y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
head(design)
ASAL<-makeContrasts(AL-AS,levels=design)
res<-glmQLFTest(fit,contrast=ASAL)
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

install.packages("gplots")
library(gplots)
col.pan<-colorpanel(100,"blue","white","red")
heatmap.2(logCPM,col=col.pan,Rowv=TRUE,scale="none",
          trace="none",dendrogram="both",cexRow=1,cexCol=1.4,density.info="none",margin=c(10,9),lhei=c(2,10),lwid=c(2,6))

go<-goana(tr,species="Hs")
topGO(go)

keg<-kegga(tr,species="Hs")
topKEGG(keg,n=15,truncate=34)

library(GO.db)
cyt.go<-c("GO:0032465","GO:0000281","GO:0000920")
term<-select(GO.db,keys=cyt.go,columns="TERM")


Rkeys(org.Hs.egGO2ALLEGS)[Rkeys(org.Hs.egGO2ALLEGS)%in%cyt.go]
cyt.go.genes<-as.list(org.Hs.egGO2ALLEGS[Rkeys(org.Hs.egGO2ALLEGS)[Rkeys(org.Hs.egGO2ALLEGS)%in%cyt.go]])
fry(y,index=cyt.go.genes,design=design,contrast=ASAL)
##########################
setwd("../THCA/wes/")
file_list = dir(pattern = "*.vcf")
for (i in 1:length(file_list)){
  file_in = file_list[i]
  file_out = paste(file_list[i],"maf",sep = ".") 
  cmd = paste("vcf2maf.pl --input-vcf ",file_in,"--output-maf ",file_out,
              "--ref-fasta /home/tienan/.vep/homo_sapiens/87_GRCh38/hg38.fa --filter-vcf \"\" --vep-path ~/Downloads/ensembl-vep --vep-data /home/tienan/.vep",seq=" ")
  system(cmd)
}
paste(file_list[1],"maf",sep = ".") 

#system("salmon index -t hg38.transcript.fa  -i  salmon_transcript_index  --type quasi -k 31", intern = TRUE)



#vcf2maf.pl --input-vcf  homo_sapiens_GRCh38.vcf --output-maf homo_sapiens_GRCh38.maf --ref-fasta /home/tienan/.vep/homo_sapiens/87_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --filter-vcf "" --vep-path ~/Downloads/ensembl-vep --vep-data /home/tienan/.vep




