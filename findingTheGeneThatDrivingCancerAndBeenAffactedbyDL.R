#finding the Gene that driving cancer and been affacted by DL
getwd()
setwd("../DL/")
# read different exp gene data 
genes_1975_dl_con = read.table("1975_con_dl/gene_exp.diff",header = T,sep = "\t")
genes_A549_dl_con = read.table("A549_con_dl/gene_exp.diff",header = T,sep = "\t")
#cutoff < 0.1 to get diff exp gene data
diff_H1975 = genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]
diff_A549 = genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]
# intersect of gene name 
diff_gene = intersect(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]$gene,
                      genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]$gene)
# extract the foldchange of diff exp gene
diff_gene_fold = cbind(as.character(genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$gene),
                       as.character(genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$gene),
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$log2.fold_change.)
# extract the data with the same direction 
diff_gene_filer_1 = diff_gene_fold[as.numeric(diff_gene_fold[,3])*as.numeric(diff_gene_fold[,4])>0,]
gene_name = as.data.frame(sort((diff_gene_filer_1[,1])))

DLGene = diff_gene_filer_1[,c(1,3)]
colnames(DLGene)=c("gene_name","FoldChange")

file_list = dir(pattern = "*.fpkm*")
tmp_file=c()
tmp_file = data.frame(1:23)
for (i in 1:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp$gene_id=dat_tmp$gene_id
  merge_tmp = merge(dat_tmp,gene_name,by.x = "gene_id",by.y  = "gene_name")
  tmp_file =cbind(tmp_file,as.data.frame(merge_tmp$FPKM))
}
tmp_file$gene_id = merge_tmp$gene_id
colnames(tmp_file)=c("id",file_list,"gene_id")


profileDL = cbind(tmp_file[,grepl(pattern = "dl",colnames(tmp_file))],tmp_file[,grepl(pattern = "con",colnames(tmp_file))])
rownames(profileDL)=toupper(tmp_file$gene_id)
#colnames(profileDL) = ()

#install.packages("pheatmap")
library(pheatmap)

pheatmap(profileDL,
         clustering_distance_cols  = "euclidean",
         show_colnames =   T,
         scale = "row",
         cluster_cols = F,
         cluster_rows = F
)





################################################hg38
#hg 38 RNA se
library("TCGAbiolinks")
library(SummarizedExperiment)




query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM")


LUADRnaseqSE <- GDCprepare(query)

LUADMatrix <- assay(LUADRnaseqSE) 

rownames(LUADMatrix) <- values(LUADRnaseqSE)$external_gene_name

# 首先，将刚才GDCprepare好的数据进行normalization，用normalization()
# 这里注意geneInfo=geneInfoHT，default其实是geneInfo，但由于我们前面选择的是HTseq，所以要选择geneInfoHT

dataNorm <- TCGAanalyze_Normalization(tabDF = LUADMatrix, geneInfo =  geneInfoHT)
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

##############################
rownames(LUADMatrix) <- values(LUADRnaseqSE)$external_gene_name

patient_id = colnames(LUADMatrix )
for (i in 1:length(patient_id)){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "_")
  patient_id[i] = tmp
}


gene_name_exp_dif = LUADMatrix[,sign == "01A"|sign == "11A"]




LUADProfile = LUADMatrix[rownames(LUADMatrix)%in%toupper(tmp_file$gene_id),order(sign)]

pheatmap(LUADProfile,
         clustering_distance_cols  = "euclidean",
         show_colnames =   F,
         scale = "row",
         cluster_cols = F,
         cluster_rows = F
)


t1 = edgeR::DGEList(gene_name_exp_dif,group = as.factor(sign[sign == "01A"|sign == "11A"]))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)

logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]

tableDEA <- tableDEA[tableDEA$PValue <= 0.1, ]
tableDEA$external_gene_name=rownames(tableDEA)
CanerRelativeGene=tableDEA[,c(5,1)]




##############################

CanerRelativeGene = dataDEGsFiltLevel[,c(8,2)]



DLCancerGene = merge(CanerRelativeGene,DLGene,by.x ="external_gene_name" ,by.y="gene_name")

DLCancerGene[DLCancerGene$logFC*as.numeric(as.character(DLCancerGene$FoldChange)) >0,]






library("biomaRt")
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")
listAttributes(snpmart )

affyids=DLCancerGene[DLCancerGene$logFC*as.numeric(as.character(DLCancerGene$FoldChange)) >0,1]



genelist = getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                     'start_position', 'end_position', 'band'),
      filters = 'hgnc_symbol', 
      values = affyids, 
      mart = ensembl)

snpRisk = getBM(attributes = c('refsnp_id','allele','chrom_start','ensembl_gene_stable_id','chr_name'), 
                        filters = c('chr_name','start','end'), 
                        values = list(genelist[1,2],genelist[1,3],genelist[1,4]), #VGF 
                        mart = snpmart)
for (i in 3:nrow(genelist)){
  tmp = getBM(attributes = c('refsnp_id','allele','chrom_start','ensembl_gene_stable_id','chr_name'), 
              filters = c('chr_name','start','end'), 
              #     values = list(11,313506,315272), #IFITM1 
              #     values = list(19,55364382,55370463), #IL11
              #    values = list( 1,1001138,1014540 ), #ISG15
              values = list(genelist[i,2],genelist[i,3],genelist[i,4]), #VGF 
              mart = snpmart)
  snpRisk=rbind(snpRisk,tmp)
}

write.csv(snpRisk,file = "SNPRISK.csv")

genelist[1,2]

# hgnc_symbol chromosome_name start_position end_position   band
# 1      IFITM1              11         313506       315272  p15.5
# 2        IL11              19       55364382     55370463 q13.42
# 3       ISG15               1        1001138      1014540 p36.33
# 4         VGF               7      101162509    101165569  q22.1




VGF = getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_strand'), 
     filters = c('chr_name','start','end'), 
#     values = list(11,313506,315272), #IFITM1 
#     values = list(19,55364382,55370463), #IL11
#    values = list( 1,1001138,1014540 ), #ISG15
     values = list(7,101162509,101165569), #VGF 
      mart = snpmart)






