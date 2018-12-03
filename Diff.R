dat = read.table("LUSC_eset(1).txt",header = F,sep = "\t")

dat = dat[,-463]

ncol(dat)

group = dat[1,]
group[group=="11"]=0
group=group[-1]

geneData = dat[-1,-1]
rownames(geneData) = dat[-1,1]
head(geneData)

head(dat)

library(limma)
?limma
limmaUsersGuide(view=T)
source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")

targets <- readTargets("targets.txt")
?eBayes

#example
set.seed(2016)
sigma2 <- 0.05 / rchisq(100, df=10) * 10
y <- matrix(rnorm(100*6,sd=sqrt(sigma2)),100,6)

Group=c(0,0,0,1,1,1)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))


View(dat[-1,])

design = cbind(Intercept=1,Group=as.numeric(group))
fit <- lmFit(geneData,design)




?lmFit
?eBayes
?topTable
?treat
#  Moderated t-statistic
fit <- eBayes(fit)
res = topTable(fit,number = 20600,sort.by = "p")
?topTable
source("Rgenetic/basicTools.R")
writeOut(res,"LUSCDiff_2.txt")

nrow(geneData)


head(geneData[,c(1:10)])

figureName = "MAGEC2"

tiff(filename = figureName,
     width = 1000, height = 1000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 300
     )

boxplot(genesExp~group,
        data = data.frame(genesExp = 
                            as.numeric(geneData[rownames(geneData)==figureName,-463]),
                                                group = as.character(group[-463])))

dev.off()






head(data.frame(as.numeric(geneData[1,]),as.character(group)))


data = data.frame(genesExp = 
                    as.numeric(geneData[rownames(geneData)=="CENPA",-463]),
                  group = as.character(group[-463]))



?manova()




######code of result 
dat = read.table("LUSC_eset(1).txt",header = F,sep = "\t")
design = cbind(Intercept=1,Group=as.numeric(group))
fit <- lmFit(geneData,design)
design = cbind(Intercept=1,Group=as.numeric(group))
fit <- lmFit(geneData,design)
design = cbind(Intercept=1,Group=as.numeric(group))
fit <- lmFit(geneData,design)
dat = dat[,-463]
ncol(dat)
group = dat[1,]
group[group=="11"]=0
group=group[-1]
geneData = dat[-1,-1]
rownames(geneData) = dat[-1,1]
head(geneData)
head(dat)
library(limma)
design = cbind(Intercept=1,Group=as.numeric(group))
fit <- lmFit(geneData,design)
tfit <- treat(fit,lfc=log2(1.1))
res = topTreat(tfit,number = 20600,coef=2)
writeOut(res,"LUSCDiff_2.txt")


writeOut = function(dat,name){
  write.table(dat,file=name,quote = F,sep="\t",row.names = T,fileEncoding = "utf-8")
}