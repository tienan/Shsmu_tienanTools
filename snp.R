source("https://bioconductor.org/biocLite.R")
#install.packages("BiocManager")
#BiocManager::install("signeR")
library(signeR)
library(VariantAnnotation)
library(Rsamtools)
#biocLite("BSgenome.Hsapiens.UCSC.hg19")

#iocManager::install("maftools")
require(maftools)
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')

laml = read.maf(maf = laml.maf, clinicalData = laml.clin)


plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plot.new()
oncoplot(maf = laml, top = 10, fontSize = 12)
oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1'))
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
lollipopPlot(maf = laml, gene = 'DNMT3A', AACol = 'Protein_Change', showMutationRate = TRUE)
mafSurvival(maf = laml, genes = 'DNMT3A', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)

library(plyr)





install.packages("ukbtools")
install.packages("devtools")
library("devtools")
# library(ukbtools)
my_ukb_data = ukb_df("ukbxxxx", path = "D:/R/")
devtools::install_github("kenhanscombe/ukbtools", dependencies = TRUE)
options(unzip = "internal")
install.packages("plotr")
install.packages("reprex")
devtools::install_github("tienan/lxctk",args ="--no-test-load")
install.packages("Rcpp", dependencies=TRUE, INSTALL_opts = c('--no-lock'))

devtools::install_github("lixiangchun/lxctk", dependencies = TRUE)

library(lxctk)
source("lxctk-master/R/sortDataFrame.R")
source("lxctk-master/R/oncoprinter.R")

xx <- sort.data.frame.by.index(x, nmfsubtypes)
oncoprinter(xx, cols, legend.panel)

legend.panel <- data.frame(V1=c(7,6,5,4,3,2),
                           V2=c("Nonsense","Splice site","Frame shift",
                                "Inframe indel","Missense","Syn."))
# The molecular subtype
nmfsubtypes <- RegularMutatedGC$nmf_clustid
# Mutation matrix of SMGs across samples.
x <- RegularMutatedGC[, 15:45] # Columns are SMGs.
# Sort x column-by-column for better visualization
x <- sortDataFrame(x, decreasing = TRUE)
# Color codes will be used in visualizing matrix 'x'.
cols <- c('white','grey88','#644B39','forestgreen',
          '#FF8B00','#9867CC','#DB1C00')
oncoprinter(x, cols, legend.panel)


##bigmemory
devtools::install_github("kaneplusplus/bigmemory")
devtools::install_github("kaneplusplus/bigtabulate")
devtools::install_github("kaneplusplus/biganalytics")
devtools::install_github("kaneplusplus/bigalgebra")
devtools::install_github("kaneplusplus/synchronicity")



x <- read.big.matrix("ALLtraining.txt", sep = "\t", type = "integer",
                     shared = TRUE, col.names = c("movie", "customer", "rating",
                                                    "year", "month"))


library(lxctk)
library(plotr)


install.packages('deconstructSigs',dependencies = T)

install.packages('deconstructSigs')
# dependencies 'BSgenome', 'BSgenome.Hsapiens.UCSC.hg19' 
BiocInstaller::biocLite('BSgenome')
BiocInstaller::biocLite('BSgenome.Hsapiens.UCSC.hg19')
## https://github.com/raerose01/deconstructSigsfile1='TCGA.STAD.muse..somatic.maf.gz'
TCGA.STAD.muse=read.table(file1,sep = '\t',quote="",header = T)
TCGA.STAD.muse[1:5,1:15]
## data frame including 5 columns: sample.ID,chr,pos,ref,alt sample.mut.ref <- data.fram(Sample='TCGA.STAD.muse',                                 chr = TCGA.STAD.muse[,5],                                 pos = TCGA.STAD.muse[,6],                                 ref = TCGA.STAD.muse[,11],                                 alt = TCGA.STAD.muse[,13])sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref,                                 sample.id = "Sample",                                 chr = "chr",                                 pos = "pos",                                 ref = "ref",                                 alt = "alt")class(sigs.input)sample_1 = whichSignatures(tumor.ref = sigs.input,                            signatures.ref = signatures.nature2013,                            sample.id = 'TCGA.STAD.muse',                            contexts.needed = TRUE,                           tri.counts.method = 'exome')# Plot example                  plot_example <- whichSignatures(tumor.ref = sigs.input  ,                         signatures.ref = signatures.nature2013,                        sample.id = 'TCGA.STAD.muse' )# Plot outputplotSignatures(plot_example, sub = 'example')
#mutect,muse,vanscan,somaticsniper这4款软件call 到的somatic mutation文件。

BiocInstaller::biocLite('BSgenome')
BiocInstaller::biocLite('BSgenome.Hsapiens.UCSC.hg19')
BiocInstaller::biocLite('BSgenome.Hsapiens.UCSC.hg38')

# 注意 ... 位置可能不同
file1='TCGA.STAD.muse..somatic.maf.gz'
laml = read.maf(maf ="TCGA.LUAD.varscan.acb6852e-dd48-4ca5-80f2-3d1a2c7d7ceb.DR-10.0.somatic.maf")

TCGA.STAD.muse=read.table("TCGA.LUAD.varscan.acb6852e-dd48-4ca5-80f2-3d1a2c7d7ceb.DR-10.0.somatic.maf",sep = '\t',quote="",header = T)
#TCGA.STAD.muse=read.table(file1,sep = '\t',quote="",header = T)
TCGA.STAD.muse[1:5,1:15]
## data frame including 5 columns: sample.ID,chr,pos,ref,alt 
library(BSgenome.Hsapiens.UCSC.hg38)
sample.mut.ref <- data.frame(Sample='TCGA.STAD.muse', 
                            chr = TCGA.STAD.muse[,5], 
                            pos = TCGA.STAD.muse[,6], 
                            ref = TCGA.STAD.muse[,11], 
                            alt = TCGA.STAD.muse[,13])
sigs.input <- mut.to.sigs.input(mut.ref  = sample.mut.ref, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
class(sigs.input)

sample_1 = whichSignatures(tumor.ref = sigs.input, 
                           signatures.ref = signatures.nature2013, 
                          sample.id = "TCGA.STAD.muse", 
                           contexts.needed = TRUE,
                           tri.counts.method = 'exome')

# Plot example                  
plot_example <- whichSignatures(tumor.ref = sigs.input  ,  
                                signatures.ref = signatures.nature2013, 
                                sample.id = 'TCGA.STAD.muse' )

# Plot output
plotSignatures(sample_1 , sub = 'example')

makePie(sample_1, sub = 'example')
library(deconstructSigs)

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = "TCGA.LUAD.varscan.acb6852e-dd48-4ca5-80f2-3d1a2c7d7ceb.DR-10.0.somatic.maf",clinicalData = clinical_LUAD)
mut =prepareMutSig(maf = laml)




source("./TCGA2STAT/R/TCGA2STAT.R")
rnaseq.ov <- getTCGA(disease="OV", data.type="RNASeq", type="RPKM")

rnaseq.thca <- getTCGA(disease="THCA", data.type="RNASeq", type="RPKM")

?getTaskCallbackNames

Sys.which("tar")
Sys.setenv(R_GZIPCMD ="C:\\Rtools\\bin\\gzip.exe")
Sys.setenv(TAR ="C:\\Rtools\\bin\\tar.exe")


devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
library("TCGAbiolinks")
clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
colnames(clinical_LUAD)

clinical_names=c("Tumor_Sample_Barcode","FAB_classification","days_to_last_followup","Overall_Survival_Status")

clinical_LUAD$bcr_patient_barcode,
clinical_LUAD$tumor_stage,
clinical_LUAD$days_to_last_follow_up,
clinical_LUAD$days_to_death


#https://rdrr.io/bioc/maftools/man/prepareMutSig.html
#http://www.bio-info-trainee.com/2518.html 教程
#https://github.com/raerose01/deconstructSigs
#https://gdc.cancer.gov/access-data/gdc-data-transfer-tool  数据下载
# https://gdc.cancer.gov/access-data/obtaining-access-controlled-data 数据获取教程
