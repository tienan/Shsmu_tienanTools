install.packages("rms")  #这一步一劳永逸。安装好后，再次启动R时，无需再次输入代码，下同
install.packages("Hmisc")
install.packages("lattice")
install.packages("survival")
install.packages("Formula")
install.packages("ggplot2")
install.packages("foreign")

library("rms")  #这一步一劳永逸。安装好后，再次启动R时，无需再次输入代码，下同
library("Hmisc")
library("lattice")
library("survival")
library("Formula")
library("ggplot2")
library("foreign")
library(pROC)

F1<-read.table("ruijinNongram_1.txt",header = T,sep = "\t")
attach(F1)
ddist<-datadist(Size,pT,Location,Ulcer,NerveInvasion,VascularInvasion)
options(datadist='ddist')
logi<-lrm(LNM~Size+pT+Location+Ulcer+NerveInvasion+VascularInvasion, x=TRUE, y=TRUE)
nomo<-nomogram(logi, fun=plogis, fun.at=c(.001, .01, .05, seq(.1,.9,by=.1), .95, .99, .999),
               lp=F, funlabel="LNM")
plot(nomo)



cal<-calibrate(logi, method="boot", B=1000, bw=FALSE, rule="p",
               type="individual", sls=.05, aics=0, force=NULL, estimates=T, pr=FALSE,
               smoother="lowess", digits=NULL)
plot(cal,scat1d.opts=list(nhistSpike=500))



Predi<-predict(logi,F1,type="lp")
roc1<-roc(F1$LNM,Predi,legacy.axes=TRUE)
plot(roc1,font=2,legacy.axes=TRUE)
ci.auc(roc1)
auc(roc1)



dat = read.table("ruijinNongram_1.txt",header = T,sep = "\t")


#logi<-lrm(Outcome~Age+pT+pN+Grade+ molecularSubtype+geneRs, data=g1,x=TRUE, y=TRUE)

logi<-lrm(LNM ~ Size+pT+Location+Ulcer+NerveInvasion+LymphInvasion+VascularInvasion, data=dat,x=TRUE, y=TRUE)

nomo<-nomogram(logi,  fun=plogis,  fun.at=c(.001,  .01,  .05,  seq(.1,.9,  by=.1), .95, .99, .999),lp=F, funlabel="Outcome")


tiff(filename = "Fig-1.tif",
     width =800, height = 800, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 300)
plot(nomo)
dev.off()

Pred = predict(logi,g1,type="lp")
roc1=roc(g1$Outcome,Predi,legend.axes=T)

plot(roc1,font=2,legacy,axes=T)


g2=F1[id[(1*126):(2*126)],] 

g3=F1[id[(2*126):(3*126)],]

g4=F1[id[(3*126):(4*126)],]  

