#PLOT
setwd("D:/r/TongrenZhang/")
getwd()
gene_list = read.table("geneList.txt",header = T,sep="\t")
head(gene_list)
ref_seq_list = read.table("refseqlist.txt")
ref_2_geneSymbol = read.table("conv_3F7AE0CAF5791540436862887.txt",header = T,sep = "\t")
head(ref_2_geneSymbol )
dat_control = read.table("SSC_neo.fpkm.genes.refseq_id.no_NR.txt",header = T,sep="\t")
head(dat_control)
dat_exp = read.table("SSC_VG_UV.txt",header = T, sep="\t")
head(dat_exp)
dat_protein = read.table("SSC-proteinome-exp49_prot sorted by gmean C4.xlsx")

dat_exp_merge_geneSymbol=merge(ref_2_geneSymbol,dat_exp,by.x = "From",by.y = "refseq_id_1")
nrow(dat_exp_merge_geneSymbol)
dat_exp_contrl_merge_geneSymbol = merge(dat_exp_merge_geneSymbol,dat_control,by.x = "To",by.y="Tracking_id")

dat_exp_contrl_merge_geneSymbol_target = merge(dat_exp_contrl_merge_geneSymbol,gene_list,by.x = "To",by.y="Genes")


dat_exp_contrl_merge = cbind(dat_exp,dat_control)

dat_exp_contrl_merge_ontarget = merge(dat_exp_contrl_merge,gene_list,by.x = "Gene_id",by.y = "Genes")
nrow(dat_exp_contrl_merge_geneSymbol_target)


dat_exp_contrl_merge[dat_exp_contrl_merge$Gene_id == "Lig4",1]

?merge

dat_exp_merge_geneSymbol$To
dat_control$Tracking_id

nrow(dat_control)
nrow(dat_exp)

dat_rna = cbind(dat_control$Gene_name,dat_control$SSC_neo_1,dat_control$SSC_neo_2,dat_exp[,c(9:14)])
head(dat_rna)
#dat_exp_contrl_merge_ontarget = merge(dat_rna,gene_list,by.x = "Gene_id",by.y = "Genes")

colnames(dat_rna)=c("geneName","Ctr1","Ctr2","Exp1_24","Exp2_24","Exp3_24","Exp1_48","Exp2_48","Exp3_48")

dat_rna$geneName = tolower(dat_rna$geneName)

gene_list$Genes = tolower(gene_list$Genes )

dat_rna_ontarget = merge(dat_rna,gene_list,by.x = "geneName",by.y = "Genes")

rownames(dat_rna_ontarget) = paste(dat_rna_ontarget$geneName,dat_rna_ontarget$Function,c(1:length(dat_rna_ontarget$Function)),sep = "_")
  
rownames(dat_rna_ontarget) <- gsub("_gene", "", rownames(dat_rna_ontarget))
rownames(dat_rna_ontarget) <- gsub("self_renewal", "", rownames(dat_rna_ontarget))
rownames(dat_rna_ontarget) <- gsub("__", "_SR_", rownames(dat_rna_ontarget))


library(pheatmap)
tiff(filename = "Figure-1.tif",
     width = 1800, height = 3000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
)
pheatmap(dat_rna_ontarget[,c(2:9)],clustering_distance_rows = "correlation",scale="column")
dev.off()

getwd()

library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listFilters(mart)[grep("hgnc",as.character(listFilters(mart)$description)),]
listDatasets(mart)[grep("mouse",as.character(listDatasets(mart)$description)),]
head(listDatasets(mart))


?useMart

mart = useDataset()
?useDataset

mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset="",host="www.ensembl.org")
attr <- c("refseq_mrna",'hgnc_symbol',
          'chromosome_name','start_position','end_position')
head(listAttributes(mart))
listAttributes(mart)[grep("mouse",as.character(listAttributes(mart)$description)),]

ids <-ref_seq_list 
ids = c("NM_001270493", "NM_001185076","NM_001185075","NM_001185082","NM_021446
")
mart<- useDataset("mcaroli_gene_ensembl", useMart("ensembl"))
xx <- getBM(attributes=attr,filters = "refseq_mrna", values = ids,
            mart =mart)


dat_protein = read.table("protein.txt",header = T,sep = "\t")
library(ggplot2)
ggplot(dat_protein,aes(pValue,abs(logGmean))) + geom_point()
+ xlab("P_Value") + ylab("logGmean")



plot(dat_protein$p._val,dat_protein$pValue)

install.packages("gcookbook")
library(gcookbook)


install.packages("ggthemes")
install.packages("Cairo")

library(ggthemes)
library(Cairo)

dat_protein = read.table("protein.txt",header = T,sep = "\t")
dat_protein$threshold <- as.factor(ifelse(dat_protein$pValue< 0.05 & abs(dat_protein$logGmean) >=1.5,ifelse(dat_protein$logGmean> 1.5 ,'Up','Down'),'Not'))
tiff(filename = "Figure-2.tif",
     width = 2000, height = 2000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
)
ggplot(data=dat_protein, 
aes(x=logGmean, y =-log10(pValue), 
    colour=threshold,fill=threshold)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_point(alpha=1, size=1.8) +
  xlim(c(-4, 4)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=15),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=15),
        axis.text.y = element_text(face="bold",  color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=15),
        axis.title.y = element_text(face="bold",color="black", size=15))+
  labs(x="log2 (fold change)",y="-log10 (p-value)",title="Volcano picture of proteomics")

dev.off()





#### ################################################guofuyin - zhangjian
library(ggplot2)

setwd("D:/R/guofuyin-jianzhang/")
dir()
dat = read.table("dat_1.txt",header = T,sep="\t")
head(dat)
dat_for_use = dat[,c(4:7,12:16)]
rownames(dat_for_use) = dat[,1]
install.packages("pheatmap")
library(pheatmap)
tiff(filename = "Figure-1.tiff",
     width = 2000, height = 2000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
     )
pheatmap(dat_for_use,clustering_distance_rows = "correlation",scale="column")
dev.off()


data <- read.csv("volcano_2.txt",header = T)
data
# 设置颜色域

tiff(filename = "Figure-1B.tiff",
     width = 2000, height = 2000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
)
ggplot(data=data, 
       aes(x=log2FC, y =log10pvalue, 
           colour=significant,fill=significant))+
  geom_point(alpha=1, size=3)+
  geom_text(aes(y=log10pvalue+0.02,label=id))+
  scale_color_manual(values=c("blue","red"))+
  xlim(c(-3, 3)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-0.58,0.58),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
    theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=15),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=15),
        axis.text.y = element_text(face="bold",  color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=15),
        axis.title.y = element_text(face="bold",color="black", size=15))+
  labs(x="log2 (fold change)",y="-log10 (p-value)",title="Volcano picture of proteomics")
dev.off() 


ah




ggplot(data=data, 
       aes(x=log2FC, y =log10pvalue, 
           colour=significant,fill=significant))+
  geom_point(alpha=1, size=3)+
  geom_text(aes(y=log10pvalue+0.2,label=id))


