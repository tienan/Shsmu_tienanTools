# DL diff genes#############################################################################################
setwd("D:/R/DL/")
genes_1975_dl_con = read.table("1975_con_dl/gene_exp.diff",header = T,sep = "\t")
genes_A549_dl_con = read.table("A549_con_dl/gene_exp.diff",header = T,sep = "\t")
diff_gene = 
  union(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1&abs(genes_1975_dl_con$log2.fold_change.)>0.6,]$gene,
            genes_A549_dl_con[genes_A549_dl_con$q_value<0.1&abs(genes_A549_dl_con$log2.fold_change.)>0.6,]$gene
            )
diff_gene[order(diff_gene)]

diff_gene_fold = cbind(as.character(genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$gene),
                       as.character(genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$gene),
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$q_value,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$q_value
)

# extract the data with the same direction 
diff_gene_filer_1 = diff_gene_fold[as.numeric(diff_gene_fold[,3])*as.numeric(diff_gene_fold[,5])>0,]
diff_gene_filer_1[order(diff_gene_filer_1[,1]),1]
diff_gene_filer_1  = as.data.frame(diff_gene_filer_1 )
colnames(diff_gene_filer_1)=c("GeneName_1","GeneName_2","log2FoldChange_1","qval_1","log2FoldChange_2","qval_2")
# Cancer vs benign diff genes LUAD##########################################################################
source("https://bioconductor.org/biocLite.R")
BiocManager::install("SummarizedExperiment")
BiocManager::install("TCGAbiolinks")
BiocManager::install("zoo")
BiocManager::install("limma")
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
getwd()
getwd()
query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  legacy = TRUE)
# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
LUADRnaseqSE <- GDCprepare(query)

LUADMatrix <- assay(LUADRnaseqSE ,"raw_count") # 
###############
# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
#LUADRnaseq_CorOutliers <- TCGAanalyze_Preprocessing(LUADRnaseqSE)
#################################################################### one way to detect DEGs
library(edgeR)
library(limma)
patient_id = colnames(LUADMatrix )
for (i in 1:length(patient_id)){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "_")
  patient_id[i] = tmp
}


gene_name_exp_dif = LUADMatrix[,sign == "01A"|sign == "11A"]

t1 = edgeR::DGEList(gene_name_exp_dif,group = as.factor(sign[sign == "01A"|sign == "11A"]))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)

?exactTest

tmp_id = rownames(t3) 
tmp=NULL
tmp_name = NULL
for (i in 1:length(tmp_id)){
  tmp = strsplit2(tmp_id[i],split = "\\|")
  tmp_name[i] = tmp[1] 
}

rownames(t3) = tmp_name 
t4 = cbind(tmp_name,t3$table)
t4[t4$tmp_name=="PDIA3P",]


#############################################################

dataNorm <- TCGAanalyze_Normalization(tabDF = LUADMatrix, geneInfo =  geneInfo)
# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.05 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])

dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA=="PDIA3P",]$mRNA="PDIA3P1"

# diff gene analysis

taget_gene = intersect(dataDEGsFiltLevel$mRNA,diff_gene_filer_1[,1])

taget_gene[order(taget_gene)]

?merge
merge_data = merge(dataDEGsFiltLevel,diff_gene_filer_1,by.x="mRNA",by.y="GeneName_1")

dl_cancer_gene = merge_data[merge_data$logFC*
  as.numeric(as.character(merge_data$log2FoldChange_2))<0,]

nrow(dl_cancer_gene)

dl_cancer_gene$mRNA


#A549 LOC 有差异  Pdia3p1
#H1975  LOC441666 diff  LOC low exp & non-invasion dataset no sig diff
#loc 46, 52 
###################################fig.1 
source("https://bioconductor.org/biocLite.R")
BiocManager::install("ggplot2")
library(ggplot2)
###
set.seed(1234)
dat <- data.frame(cond = factor(rep(c("A","B"), each=200)), 
                  rating = c(rnorm(200),rnorm(200, mean=.8)))
# View first few rows
ggplot(dat, aes(x=rating, fill=cond)) + geom_density(alpha=.3)

###################### map diff IDs 
BiocManager::install("gage")
BiocManager::install("pathview")
library(gage)
library(pathview)
data(egSymb)
data(bods) 
head(egSymb)
refseq.data <- sim.mol.data(mol.type = "gene", id.type = "REFSEQ",
                            nexp = 2, nmol = 1000)
id.map.refseq <- id2eg(ids = as.character(dl_cancer_gene$mRNA),ategory =
                          "REFSEQ", org = "Hsa")
?eg2id()
data(gene.idtype.list)
head(gene.idtype.list)

id2eg(ids = as.character(dl_cancer_gene$mRNA), 
      category = gene.idtype.list[1], org = "Hs",geneannot.map(in.))

geneannot.map(in.ids = as.character(dl_cancer_gene$mRNA) , 
              in.type = gene.idtype.list[1],out.type = gene.idtype.list[9],
              org="Hs", pkg.name=NULL,
              unique.map=TRUE, na.rm=TRUE, keep.order=TRUE)

###################################tpm DL 
setwd("D:/R/DL/")
tpm_data = read.csv("salmon_exp_1st_seq.csv")
head(tpm_data)
diff_gene_filer_1_intersection#0.1,intersect, 

DL_genes = tpm_data[tpm_data$GeneName%in%diff_gene_filer_1_intersection[,1],c(1:8)]
rownames(DL_genes) = DL_genes[,2]


library(pheatmap)

pheatmap(DL_genes[,c(3:8)],
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)
DL_genes[,c(3:8)]
