#multiMiR requires the XML and RCurl packages from R. The packages can be installed from R by typing:
install.packages("XML")
install.packages("RCurl")

# To install multiMiR, first install suggested package BiocStyle
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("BiocStyle")

a# Now install devtools (for installing from GitHub repositories)
install.packages("devtools")
library(devtools)

# Now install the development version of the multiMiR package
devtools::install_github("kechrislab/multimir", build_vignettes = TRUE)

source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("multiMiR")
library(multiMiR)
library("biomaRt")

BiocManager::install(version = "3.8")
BiocManager::install("biomaRt", version = "3.8")

ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)



ipro = getBM(attributes=c("refseq_mrna","interpro","interpro_description"), 
             filters="refseq_mrna",
             values=refseqids, 
             mart=ensembl)

genes = list.multimir("gene", limit = 10)

refseqids = c("ENO1")
ipro = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), 
             filters="hgnc_symbol",
             values=refseqids, 
             mart=ensembl)
ipro

ipro = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), 
             filters="ensembl_gene_id",
             values=refseqids, 
             mart=ensembl)

?getBM


example1 <- get_multimir(mirna = 'hsa-miR-539-5p', summary = TRUE)
example1@summary
