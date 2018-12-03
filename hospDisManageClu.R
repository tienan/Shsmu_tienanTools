install.packages("tm")
install.packages("pheatmap")
install.packages("gplots")
install.packages("devtools")
install.packages("heatmap3")
install.packages("wavelets")
                 
library(heatmap3)
library(tm)
library(cluster)
library(pheatmap)
library("gplots")
library("devtools")
library(ggplot2)
#source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
library(wavelets) 

#Data import
dat = read.table("hospDisManage.txt",header = T,sep = "\t")
#Data check
head(dat)
dat = dat[,-13] # remove year
datClusterC50_900 = dat[dat$icd=="C50.900",] 
#Data grouped by age 
datClusterC50_900$ageGroup =  cut(datClusterC50_900$age,breaks=c(0,40,50,100),labels=c('1','2','3'))

#The cost and los distribution of different age 
summary(aov(datCluster$cost~ageGroup,data=datCluster))
boxplot(datCluster$cost~ageGroup,data=datCluster,ylim=c(-1,100000))
summary(aov(datCluster$los~ageGroup,data=datCluster))
boxplot(datCluster$los~ageGroup,data=datCluster,ylim=c(-1,50))

#Hierarchical Clustering) 
hc <- hclust(dist(datClusterC50_900[,c(1,2,3,4,5,6,7,9,10,12,14,16:23)]), "ave")
plot(hc)
# 
wss <- 2:29
for (i in 2:29) wss[i] <- sum(kmeans(datClusterC50_900[,c(1,2,3,4,5,6,7,9,10,12,14,16:23)],centers=i,nstart=25)$withinss)
print(wss)
#Wavelet to find "singular point" 
wt <- dwt(wss, filter="Daubechies", boundary="periodic")
#Kmean group 
kfit <- kmeans(datClusterC50_900[,c(1,2,3,4,5,6,7,9,10,12,14,16:23)],6, nstart=100)
datClusterC50_900$class =  kfit$cluster
#Mean and std calculation
for (i in 1:6){
  mean(datClusterC50_900Oder[datClusterC50_900Oder$class==i,]$cost)
  sd(datClusterC50_900Oder[datClusterC50_900Oder$class==i,]$cost)
  
  mean(datClass[datClass$cla==i,]$los)
  sd(datClass[datClass$cla==i,]$los)
}
#Boxplot
boxplot(cost ~ class, data = datClusterC50_900, col = "lightgray")
boxplot(los ~ class, data = datClusterC50_900, col = "lightgray")
#data clean 
datClusterC50_900Oder = data.frame()
datClusterC50_900Oder=datClusterC50_900[datClusterC50_900$class==1,]
tmp = datClusterC50_900[datClusterC50_900$class==3,]
tmp$class=2
datClusterC50_900Oder=rbind(datClusterC50_900Oder,tmp)
tmp = datClusterC50_900[datClusterC50_900$class==6,]
tmp$class=3
datClusterC50_900Oder=rbind(datClusterC50_900Oder,tmp)
tmp = datClusterC50_900[datClusterC50_900$class==2,]
tmp$class=4
datClusterC50_900Oder=rbind(datClusterC50_900Oder,tmp)
tmp = datClusterC50_900[datClusterC50_900$class==4,]
tmp$class=5
datClusterC50_900Oder=rbind(datClusterC50_900Oder,tmp)
boxplot(cost ~ class, data = datClusterC50_900Oder, col = "lightgray")
boxplot(los ~ class, data = datClusterC50_900Oder, col = "lightgray",ylim(-1,50))
#PCA analysis
ClusterC50_900 = datClusterC50_900Oder[,c(1,2,3,4,5,6,7,9,10,12,14,16:23)]
C50_900.pr <- princomp(ClusterC50_900, cor = T)
summary(C50_900.pr, loadings = T)
pca_data <- predict(C50_900.pr)
datClusterC50_900Oder$C1 = pca_data[,1]
datClusterC50_900Oder$C2 = pca_data[,2]
qplot(C1,C2, data = datClusterC50_900Oder, colour = class)
#heatmap
result<-heatmap3(t(scale(datCluster[,-10])),ColSideAnn=Ann,ColSideFun=function(x)
  showAnn(x),ColSideWidth=0.8,col=colorRampPalette(c("green",
                                                     "black", "red"))(1024),RowAxisColors=1,legendfun=function() showLegend(legend=c("Low","High"),col=c("chartreuse4","firebrick")),verbose=TRUE)




#x=sample(seq(1921,1997),1000,replace = T)
#y=cut(x,breaks=c(1920,1940,1970,1997),labels=c('young','mid-life','old'))
#table(y)
#hist(x, breaks = c(1920,1940,1970,1997))


datCluster$ageGroup =  cut(datCluster$age,breaks=c(0,40,50,100),labels=c('1','2','3'))


?a()
summary(aov(datCluster$cost~ageGroup,data=datCluster))
boxplot(datCluster$cost~ageGroup,data=datCluster,ylim=c(-1,100000))

summary(aov(datCluster$los~ageGroup,data=datCluster))
boxplot(datCluster$los~ageGroup,data=datCluster,ylim=c(-1,50))

?plot

# 
?hclust
?cmdscale
require(graphics)

hc <- hclust(dist(datCluster), "ave")


 
loc <- cmdscale(eurodist)
x <- loc[, 1]
y <- -loc[, 2] # reflect so North is at the top




head(datClusterC50_900)
wss <- 2:29
for (i in 2:29) wss[i] <- sum(kmeans(datClusterC50_900[,c(1,2,3,4,5,6,7,9,10,12,14,16:23)],centers=i,nstart=25)$withinss)
wss

plot(2:29,wss[2:29],type="b",xaxt="n",xlab="Number of Clusters",ylab="Within groups sum of squares")
axis(1,at=seq(1,30,1))

kfit <- kmeans(datClusterC50_900[,c(1,2,3,4,5,6,7,9,10,12,14,16:23)],6, nstart=100)

datClusterC50_900$class =  kfit$cluster

i=1
print(table(datClusterC50_900[datClusterC50_900$class==i,]$operation))
head(datClusterC50_900)

sink("result.txt")
sink()
kfit$size

mean(datClass[datClass$cla==i,]$cost)
sd(datClass[datClass$cla==i,]$cost)

#kfit$size

i=1
mean(datClusterC50_900Oder[datClusterC50_900Oder$class==i,]$cost)
sd(datClusterC50_900Oder[datClusterC50_900Oder$class==i,]$cost)

i=6
mean(datClass[datClass$cla==i,]$los)
sd(datClass[datClass$cla==i,]$los)



datClusterC50_900[order(datClusterC50_900$class),]$class

boxplot(cost ~ class, data = datClusterC50_900, col = "lightgray")
boxplot(los ~ class, data = datClusterC50_900, col = "lightgray")

datClusterC50_900Oder = data.frame()
datClusterC50_900Oder=datClusterC50_900[datClusterC50_900$class==1,]
tmp = datClusterC50_900[datClusterC50_900$class==3,]
tmp$class=2
datClusterC50_900Oder=rbind(datClusterC50_900Oder,tmp)
tmp = datClusterC50_900[datClusterC50_900$class==6,]
tmp$class=3
datClusterC50_900Oder=rbind(datClusterC50_900Oder,tmp)
tmp = datClusterC50_900[datClusterC50_900$class==2,]
tmp$class=4
datClusterC50_900Oder=rbind(datClusterC50_900Oder,tmp)
tmp = datClusterC50_900[datClusterC50_900$class==4,]
tmp$class=5
datClusterC50_900Oder=rbind(datClusterC50_900Oder,tmp)
boxplot(cost ~ class, data = datClusterC50_900Oder, col = "lightgray")
boxplot(los ~ class, data = datClusterC50_900Oder, col = "lightgray",ylim(-1,50))

print(table(datClusterC50_900Oder[datClusterC50_900Oder$class==5,]$operation))


cor.test(as.numeric(datClusterC50_900Oder[datClusterC50_900Oder$class==5,]$los),
    as.numeric(datClusterC50_900Oder[datClusterC50_900Oder$class==5,]$cost))



ClusterC50_900 = datClusterC50_900Oder[,c(1,2,3,4,5,6,7,9,10,12,14,16:23)]

C50_900.pr <- princomp(ClusterC50_900, cor = T)
summary(C50_900.pr, loadings = T)
pca_data <- predict(C50_900.pr)

datClusterC50_900Oder$C1 = pca_data[,1]
datClusterC50_900Oder$C2 = pca_data[,2]

qplot(C1,C2, data = datClusterC50_900Oder, colour = class)


plot(datClusterC50_900Oder$C1,datClusterC50_900Oder$C2)

#datClusterC50_900Oder=rbind(datClusterC50_900Oder,datClusterC50_900[datClusterC50_900$class==2,])
#datClusterC50_900Oder=rbind(datClusterC50_900Oder,datClusterC50_900[datClusterC50_900$class==6,])
#datClusterC50_900Oder[datClusterC50_900Oder$class==6,]$class=2


boxplot(cost ~ class, data = datClusterC50_900Oder, col = "lightgray")


?pie

for (i in 1:5)
  
  print(table(datClass[datClass$cla==6,]$icd))

print(table(datClass[datClass$cla==6,]$operation))

print(table(datClass$operation))


for (i in 1:6) print(table(datClass[datClass$cla==i,]$icd))

i=6
mean(datClass[datClass$cla==i,]$cost)
sd(datClass[datClass$cla==i,]$cost)

i=6
mean(datClass[datClass$cla==i,]$los)
sd(datClass[datClass$cla==i,]$los)





boxplot(cost ~ cla, data = datClassOrd, col = "lightgray")
boxplot(los ~ cla, data = datClassOrd, col = "lightgray")


datClassOrd = datClass[datClass$cla==1,]
for (i in c(5,6,2,3) ) datClassOrd = rbind(datClassOrd,datClass[datClass$cla==i,])

datClassOrd[datClassOrd$cla==5,]$cla = 15
datClassOrd[datClassOrd$cla==6,]$cla = 16
datClassOrd[datClassOrd$cla==2,]$cla = 4
datClassOrd[datClassOrd$cla==3,]$cla = 5
datClassOrd[datClassOrd$cla==15,]$cla = 2
datClassOrd[datClassOrd$cla==16,]$cla = 3


head(datClassOrd)
pheatmap( scale(datClassOrd[,c(1:7,9:22)])) # scale: stardlizd by 

pheatmap( scale(datClassOrd[,c(7,12,9,22)]),cluster_rows = F, cluster_cols = F) # scale: stardlizd by

install.packages("CHAID", repos="http://R-Forge.R-project.org")
install.packages("partykit")
library("CHAID")
head(datClassOrdSam[,c(7,12,15,16)])

set.seed(1)

datClassOrdSam <- datClassOrd[sample(1:nrow(datClassOrd), 1000),]

table(datClassOrdSam$cla)


table(datClassOrd$cla)

RowSideColors<-colorRampPalette(c("chartreuse4", "white", "firebrick"))(40)

result<-heatmap3(t(scale(datClassOrdSam[,c(7,12)])),ColSideCut=1.2,ColSideAnn=Ann,ColSideFun=function(x)
  showAnn(x),ColSideWidth=0.8,col=colorRampPalette(c("green",
"black", "red"))(1024),RowAxisColors=1,legendfun=function() showLegend(legend=c("Low","High"),col=c("chartreuse4","firebrick")),verbose=TRUE)



result<-heatmap3(t(scale(datCluster[,-10])),ColSideAnn=Ann,ColSideFun=function(x)
  showAnn(x),ColSideWidth=0.8,col=colorRampPalette(c("green",
"black", "red"))(1024),RowAxisColors=1,legendfun=function() showLegend(legend=c("Low","High"),col=c("chartreuse4","firebrick")),verbose=TRUE)





result$cutTable


Ann = data.frame(datClass$cla)
result<-heatmap3(t(scale(datCluster[,-10])),ColSideCut=1.8,ColSideAnn=Ann,ColSideFun=function(x)
  showAnn(x),ColSideWidth=0.8,col=colorRampPalette(c("green",
"black", "red"))(1024),RowAxisColors=1,legendfun=function() showLegend(legend=c("Low","High"),col=c("chartreuse4","firebrick")),verbose=TRUE)

Ann = data.frame(datClass$icd)
result<-heatmap3(t(scale(datCluster[,-10])),ColSideCut=1.8,ColSideAnn=Ann,ColSideFun=function(x)
  showAnn(x),ColSideWidth=0.8,col=colorRampPalette(c("green",
                                                     "black", "red"))(1024),RowAxisColors=1,legendfun=function() showLegend(legend=c("Low","High"),col=c("chartreuse4","firebrick")),verbose=TRUE)

summary(aov(datClass$cost~datClass$cla))

summary(aov(datClass$los~datClass$cla))
