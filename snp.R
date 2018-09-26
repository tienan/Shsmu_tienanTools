source("https://bioconductor.org/biocLite.R")
install.packages("BiocManager")
BiocManager::install("signeR")
library(signeR)
library(VariantAnnotation)
library(Rsamtools)

BiocManager::install("maftools")
require(maftools)
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')




plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
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
devtools::install_github("lixiangchun/lxctk",args ="--no-test-load")
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