# DL manuscript code
# Figure 1. Significant genetic profile by DL 


library("edgeR")
library("gdata")
library("heatmaply")
library("ggplot2")
library("genefilter")
library("methylumi")


############################################################
#Generating the profile of DEGs 
#############################################################
getwd()
setwd("/home/tienan/R/Shsmu_tienanTools/")
setwd("../DL/")
# read different exp gene data 
genes_1975_dl_con = read.table("1975_con_dl/gene_exp.diff",header = T,sep = "\t")
genes_A549_dl_con = read.table("A549_con_dl/gene_exp.diff",header = T,sep = "\t")

#Union#################################
diff_gene = union(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]$gene,
                  genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]$gene)
######################intersect
diff_gene = intersect(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]$gene,
                 genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]$gene)

###############################
# extract the foldchange of diff exp gene
diff_gene_fold = cbind(as.character(genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$gene),
                       as.character(genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$gene),
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$q_value,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$q_value
)

# extract the data with the same direction 
diff_gene_filer_1 = diff_gene_fold[as.numeric(diff_gene_fold[,3])*as.numeric(diff_gene_fold[,5])>0,]


gene_name = as.data.frame(sort(diff_gene_filer_1[,1]))
colnames(gene_name)="gene_name"


file_list = dir(pattern = "*.fpkm*")

#extract the data of target gene
?data.frame
dat_tmp = read.table(file_list[i],header = T,sep = "\t")
merge_tmp = merge(dat_tmp,gene_name,by.x = "gene_id",by.y  = "gene_name")
tmp_file = merge_tmp 
for (i in 2:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  merge_tmp = merge(dat_tmp,gene_name,by.x = "gene_id",by.y  = "gene_name")
  tmp_file =cbind(tmp_file,as.data.frame(merge_tmp$FPKM))
}


colnames(tmp_file[,c(1,10,14:24)]) = c("Gene",gsub("_genes.fpkm_tracking", "", file_list))

tmp_file_dat = tmp_file[,c(1,10,14:24)]
colnames(tmp_file_dat) = c("Gene",gsub("_genes.fpkm_tracking", "", file_list))

library(pheatmap)
rownames(tmp_file_dat) = tmp_file_dat[,1]
pheatmap(tmp_file_dat[,c(2:7)],
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = T,
         cluster_cols = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)

pheatmap(tmp_file_dat[,c(8:13)],
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = T,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)

pheatmap(tmp_file_dat[,c(2:13)],
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = T,
         cluster_cols = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)

write.csv(file = "../Figures/UnionDLGenesName.csv",x = tmp_file_dat[,1])

write.csv(file = "../Figures/IntersectDLGenesName.csv",x = tmp_file_dat[,1])
##########################################Non-samoke Female data###############################
#Based on the sequencing of A549 cell lines, PDIA3P1 is the target of research
###########################################TCGA set#############################################


#####################################################################
#TCGA dataset of LUAD to generate the KEGG and GSEA pthway (geneset)
#####################################################################

BiocManager::install("TCGAbiolinks")
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)

query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
LUADRnaseqSE <- GDCprepare(query.exp.hg38)

rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
LUADRnaseqSEvalues <- assay(LUADRnaseqSE)
ncol(LUADRnaseqSEvalues)
head(exp.hg38.values)


query.exp.hg38 <- GDCquery(project = "TCGA-THCA", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
query.exp.hg38 <- GDCprepare(query.exp.hg38)
THCARnaseqSE = assay(query.exp.hg38)
rownames(THCARnaseqSE) <- values(query.exp.hg38 )$external_gene_name
THCARnaseqSEValue <- assay(THCARnaseqSE)
head(THCARnaseqSE)
ncol(THCARnaseqSE)
write.csv(x = THCARnaseqSE,file = "THCARnaseqV")


query.exp.hg19 <- GDCquery(project = "TCGA-THCA", 
                           data.category = "Gene expression", 
                           data.type = "Gene expression quantification",
                           legacy = TRUE)
GDCdownload(query.exp.hg38)
exp.hg19 <- GDCprepare(query = query.exp.hg19)
exp.values <- assay(query.exp.hg19)


sign = ifelse(exp.hg38.values_caner[rownames(exp.hg38.values_caner) == targetGene,]
              >median(exp.hg38.values_caner[rownames(exp.hg38.values_caner) == targetGene,]),1,0)




#pdia3p1 Sign#########################################################################

targetGene="PDIA3P1"

##selecting the cancer sample

exp.hg38.values_caner = exp.hg38.values[,grep("01A",colnames(exp.hg38.values))]

?ifelse
# high:1 low:0
sign = ifelse(exp.hg38.values_caner[rownames(exp.hg38.values_caner) == targetGene,]
              >median(exp.hg38.values_caner[rownames(exp.hg38.values_caner) == targetGene,]),1,0)

length(sign)#519

#explanation: the first sign is the ref, namely the denominator####################################
t1 = edgeR::DGEList(exp.hg38.values_caner,group = as.factor(sign))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
write.csv(tableDEA,file = "PDIA3P1TCGA.csv")
getwd()

####################################CHECK
?boxplot
dataCheck = exp.hg38.values_caner[rownames(exp.hg38.values_caner) == "NPY",]
boxplot(log(dataCheck) ~ sign,col = "lightgray")


###########################################GSEA results
BiocManager::install("fgsea")
library(fgsea)
#load("../human_H_v5p2.rdata")
load("../human_c2_v5p2.rdata") # relative conprehensive

#input logFC vaule and entrezgene
rank = tableDEA$logFC
row.names(tableDEA)

mygenes = row.names(tableDEA)
mapping <- getBM(
  attributes = c('entrezgene', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = mygenes,
  mart = hsmart
)

names(rank) = mapping[,1]

ranks = rank[na.omit(names(rank))]

barplot(sort(ranks, decreasing = T))
pathwaysH <- Hs.c2
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)

################## result 
TCGAGsea = fgseaRes[order(padj, -abs(NES)), ] ################## result

topUp <- TCGAGsea %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- TCGAGsea %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway], 
              ranks, 
              TCGAGsea, 
              gseaParam = 0.5)
################## result 




head(fgseaRes[order(padj, -abs(NES)), ], n=10)
plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)

kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')

TCGApathway  =  kk@result[kk@result$pvalue<0.1,2]
#####################################################################
#non-female dataset generate the KEGG and GSEA pthway (geneset)
#####################################################################

getwd()
setwd("../non_female_smoker_pm2.5/")
dir()

file_list = dir(pattern = "^[1-9]")
i=1
tmp = read.table(file_list[i],header = T,sep = "\t")
dat = tmp[,c(1,10)]
head(dat)
for (i in 2:length(file_list)){
  tmp = read.table(file_list[i],header = T,sep = "\t")
  dat =cbind(dat,tmp[,10])
}
colnames(dat) = c("Genes",gsub("_genes.fpkm_tracking", "", file_list))


targetGene="PDIA3P1"

##selecting the cancer sample

dat_caner = dat[,grep("A1",colnames(dat))]

dat_caner = as.matrix(dat_caner)

head(dat_caner)

sign = ifelse(dat_caner[dat[,1] == targetGene,]
              >median(dat_caner[dat[,1] == targetGene,]),1,0)

#explanation: the "second" sign is the ref, namely the denominator
t1 = edgeR::DGEList(dat_caner,group = as.factor(sign))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)

?boxplot
dataCheck = dat_caner[28868,]
boxplot(log(dataCheck) ~ sign,col = "lightgray")

tableDEA$geneName = dat[rownames(tableDEA),1]

write.csv(tableDEA,file = "PDIA3P1non_female_smoker_pm2.5.csv")




load("../human_c2_v5p2.rdata") # relative conprehensive

#input logFC vaule and entrezgene
rank = tableDEA$logFC
row.names(tableDEA)

mygenes = dat[row.names(tableDEA),1]
mapping <- getBM(
  attributes = c('entrezgene', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = mygenes,
  mart = hsmart
)

names(rank) = mapping[,1]

ranks = rank[na.omit(names(rank))]

barplot(sort(ranks, decreasing = T))
pathwaysH <- Hs.c2
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)

################## result 
NonFemaleNonSmakerGsea = fgseaRes[order(padj, -abs(NES)), ] ################## result
topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)

kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')

NonFemaleNonSmakerpathway  =  kk@result[kk@result$pvalue<0.05,2]
################## result 

#####################################################

library(GEOquery)

gset <- getGEO("GSE86958", GSEMatrix =TRUE, AnnotGPL=TRUE )

getwd()
setwd("../GSE86958//")
dir()




file_list = dir(pattern = "_Tumor_fpkm.txt")
i=1
tmp = read.table(file_list[i],header = T,sep = "\t")
dat = tmp[,c(2,6)]
head(dat)
for (i in 1:length(file_list)){
  tmp = read.table(file_list[i],header = T,sep = "\t")
#  dat =cbind(dat,tmp[,6])
  print(nrow(tmp))  
}

for (i in c(2:4)){
  tmp = read.table(file_list[i],header = T,sep = "\t")
  dat =cbind(dat,tmp[,6])
}

colnames(dat) = c("Genes",gsub("_fpkm.txt", "", file_list[c(1:4)]))



tmp = read.table(file_list[5],header = T,sep = "\t")
dat_1 = tmp[,c(2,6)]
tmp = read.table(file_list[6],header = T,sep = "\t")
dat_1 =cbind(dat_1,tmp[,6])
colnames(dat_1) = c("Genes",gsub("_fpkm.txt", "", file_list[c(5,6)]))

targetGene="PDIA3P"

dat_caner = as.matrix(dat[,c(2:5)])

sign = ifelse(dat_caner[dat[,1] == targetGene,]
              >median(dat_caner[dat[,1] == targetGene,]),1,0)

#explanation: the first sign is the ref, namely the denominator
t1 = edgeR::DGEList(dat_caner,group = as.factor(sign))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$PValue <= 0.01, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)


?boxplot
dataCheck = dat_caner[12005,]
boxplot(log(dataCheck) ~ sign,col = "lightgray")

tableDEA$geneName = dat[rownames(tableDEA),1]

write.csv(tableDEA,file = "PDIA3P1GSE86958.csv")


load("../human_c2_v5p2.rdata") # relative conprehensive

#input logFC vaule and entrezgene
rank = tableDEA$logFC
row.names(tableDEA)

mygenes = dat[row.names(tableDEA),1]
mapping <- getBM(
  attributes = c('entrezgene', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = mygenes,
  mart = hsmart
)

names(rank) = mapping[,1]

ranks = rank[na.omit(names(rank))]

barplot(sort(ranks, decreasing = T))
pathwaysH <- Hs.c2
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)

################## result 
GSE86958Gsea = fgseaRes[order(padj, -abs(NES)), ] ################## result
topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)

kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')

GSE86958pathway  =  kk@result[kk@result$pvalue<0.05,2]


################## result 


intersect(intersect(NonFemaleNonSmakerpathway,GSE86958pathway),TCGApathway)

intersect(GSE86958pathway,TCGApathway)

intersect(intersect(NonFemaleNonSmakerGsea[,1],TCGAGsea[,1]),GSE86958Gsea[,1])

pathwaysH 

############################################################KEGG


BiocManager::install("clusterProfiler")
library(clusterProfiler)
search_kegg_organism('hsa', by='kegg_code')
?enrichKEGG


head(kk, n=10)
data(geneList, package='DOSE')
browseKEGG(kk, 'hsa05034')
BiocManager::install("pathview")

library(pathview)

ranks


pathview(gene.data = ranks, 
         pathway.id = "hsa05034", 
         species = "hsa", 
         limit = list(gene=5, cpd=1))

data(examplePathways)
data(exampleRanks)
#######################################
library(biomaRt)

hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")

hsmart

# Object of class 'Mart':
#   Using the ENSEMBL_MART_ENSEMBL BioMart database
#   Using the hsapiens_gene_ensembl dataset
mygenes <- c("TNF", "IL6", "IL1B", "IL10", "CRP", "TGFB1", "CXCL8")
mapping <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = mygenes,
  mart = hsmart
)

library(GEOquery)

GSM118720 <- getGEO('GSM118720')
GSM118721 <- getGEO('GSM118721')

GSM118720 <- getGEO(filename=system.file("extdata/GSM118720.soft",package="GeneExpressionSignature"))
#control gene-expression profiles
GSM118721 <- getGEO(filename=system.file("extdata/GSM118721.soft",package="GeneExpressionSignature"))

control <- as.matrix(as.numeric(Table(GSM118721)[,2]))
head(control)
treatment <- as.matrix(as.numeric(Table(GSM118720)[,2]))
ranked_list <-getRLs(control,treatment)
data(exampleSet)
show(exampleSet)






library(dplyr)


topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway], 
              rank, 
              fgseaRes, 
              gseaParam = 0.5)

BiocManager::install("goseq")
library(goseq)
supportedOrganisms() %>% filter(str_detect(Genome, "mm"))

genes =as.integer(tableDEA$FDR < 0.01 & !is.na(tableDEA$FDR))
tableDEA[tableDEA$FDR < 0.01 & !is.na(tableDEA$FDR),]
names(genes)=rownames(tableDEA[tableDEA$FDR < 0.01 & !is.na(tableDEA$FDR),])


mapping <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values =names(genes),
  mart = hsmart
)
names(genes) = mapping[,1]

genes = genes[na.omit(names(genes))]



pwf <- nullp(genes, "hg19", "ensGene")
?nullp

getlength(names(genes),'hg19','ensGene')

BiocManager::install("org.Hs.eg.db")

goResults <- goseq(pwf, "hg19","ensGene", test.cats=c("GO:CC"))
library(ggplot2)
goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")


