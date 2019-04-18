# DL manuscript code
# Figure 1. Significant genetic profile by DL 
library("edgeR")
library("gdata")
library("heatmaply")
library("ggplot2")
library("genefilter")
library("methylumi")
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



