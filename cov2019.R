setwd("../R")
remotes::install_github("GuangchuangYu/nCov2019")
library(nCov2019)
x=load_nCov2019()
head(x[])
#province summary
data=summary(x)[,1:5]
# write.csv(data,"~/R/nCov2019/data.csv")
# getwd()
# dat = read.table("../R/Cov2019Data.txt",header = T,sep = "\t")
# head(x)
# head(dat)
# row.names(dat)=dat[,1]

install.packages("pheatmap")
library(pheatmap)

pheatmap(dat[,c(2)],
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = F,
         scale="column",
         fontsize_row = 0.1, 
         fontsize_col = 15)



# install.packages("cluster")
# 
# head(x)
# x = as.data.frame(x[,c(1:6)])
# 
# md = melt(x,id=c("province","city","time"),measure=c("cum_confirm"))
# head(md)
# 
# 
# md_1 = dcast(md,province+time~variable,sum,margins = c("province","variable"))
# md_1 = md_1[,c(1:3)]
# head(md_1)

#code 

library(cluster)
library(reshape2)
library(reshape2)
library(sjPlot)
#library(effects)
library(sjPlot)
library(ggplot2)
#library(comato)
library(NbClust)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization


md_2 = dcast(data[,c(1:3)],province~time)
md_2[is.na(md_2)] <- 0
head(md_2)

rownames(md_2) = md_2[,1]
md_3 = as.data.frame(md_2[,c(2:ncol(md_2))])
head(as.data.frame(md_3))
md_4=scale(md_3)

distance <- get_dist(md_4)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(md_4, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


md_5 = md_4[rownames(md_4)!="湖北",]

wss <- function(k) {
  kmeans(md_5, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


k4 <- kmeans(md_3, centers = 4, nstart = 25)
fviz_cluster(k4, data = md_3)
sed.seed(123)
?fviz_cluster
k4 <- kmeans(md_4, centers = 4, nstart = 25)
fviz_cluster(k4, data = md_3,ggtheme = theme_classic(),repel = F,geom = c("point"))
# 
# install.packages("sjPlot")
# install.packages("effects")
# install.packages("comato")
# install.packages("reshape2")
# 
# BiocManager::install("sjPlot")
# devtools::install_github('cdeterman/gpuR', ref = 'develop')

remove.packages("sjPlot")
sj(data, steps = 15, show.diff = FALSE)
install.packages("Rtools")
BiocManager::install("Rtools")
install.packages("sjPlot")
library(devtools)
devtools::install_github("strengejacke/sjPlot")

install.packages("NbClust")
install.packages("factoextra")

cluster =  as.data.frame(k4$cluster)
cluster = cbind(row.names(cluster),cluster)

colnames(cluster)= c("province","cluster")

data

data_last_day = data[data$time=="2020-02-13",]
data_last_day = data_last_day[order(data_last_day$province),]
data_last_day_cluster = as.data.frame(cbind(data_last_day,as.data.frame(k4$cluster)))


library("RODBC")
library(Hmisc)
library(dplyr)
library(dbplyr)
library(data.table)
#library(odbc)
data_last_day_cluster$cum_dead = as.numeric(data_last_day_cluster$cum_dead)
by_cluster <- group_by(data_last_day_cluster, k4$cluster)

basicData <- by_cluster %>% 
  summarise(
    median_confirm= median(cum_confirm),P1_confirm=fivenum(cum_confirm)[2],P3_confirm=fivenum(cum_confirm)[4],
    median_heal= median(cum_heal),P1_heal=fivenum(cum_heal)[2],P3_heal=fivenum(cum_heal)[4],
    median_dead= median(cum_dead),P1_dead=fivenum(cum_dead)[2],P3_dead=fivenum(cum_dead)[4]
    )
basicData 

write.csv(x=basicData,"basicData.csv")

data_last_day_cluster[data_last_day_cluster$`k4$cluster`==1,]$province#2
data_last_day_cluster[data_last_day_cluster$`k4$cluster`==2,]$province#3
data_last_day_cluster[data_last_day_cluster$`k4$cluster`==3,]$province#3

install.packages('httr')
Sys.setenv("TAR" = "internal")
devtools::install_github("chinamap")
remotes::install_github("GuangchuangYu/chinamap")

library(maptools)
#install.packages("maps")
#install.packages("mapdata")
library(maps)
library(mapdata)
map("china")
getwd()
x=readShapePoly('bou2_4p.shp')

plot(x,col=gray(924:0/924));

getColor=function(mapdata,provname,provcol,othercol)
{
  f=function(x,y) ifelse(x %in% y,which(y==x),0);
  colIndex=sapply(mapdata@data$NAME,f,provname);
  col=c(othercol,provcol)[colIndex+1];
  return(col);
}

city = read.csv(text = "城市,jd,wd
               北京,116.4666667,39.9
               上海,121.4833333,31.23333333
               天津,117.1833333,39.15
               重庆,106.5333333,29.53333333
               黑龙,126.6833333,45.75
               吉林,125.3166667,43.86666667
               辽宁,123.4,41.83333333
               内蒙,111.8,40.81666667
               河北,114.4666667,38.03333333
               山西,112.5666667,37.86666667
               山东,117,36.63333333
               河南,113.7,34.8
               陕西,108.9,34.26666667
               甘肃,103.8166667,36.05
               宁夏,106.2666667,38.33333333
               西宁,101.75,36.63333333
               新疆,87.6,43.8
               安徽,117.3,31.85
               江苏,118.8333333,32.03333333
               浙江,120.15,30.23333333
               湖南,113,28.18333333
               江西,115.8666667,28.68333333
               湖北,114.35,30.61666667
               四川,104.0833333,30.65
               贵州,106.7,26.58333333
               福建,119.3,26.08333333
               台湾,121.5166667,25.05
               广东,113.25,23.13333333
               海南,110.3333333,20.03333333
               广西,108.3333333,22.8
               云南,102.6833333,25
               西藏,91.16666667,29.66666667
               香港,114.1666667,22.3
               澳门,113.5,22.2")
city[5,1]="黑龙江"

provname=c("北京市","天津市","河北省","山西省","内蒙古自治区",
           "辽宁省","吉林省","黑龙江省","上海市","江苏省",
           "浙江省","安徽省","福建省","江西省","山东省",
           "河南省","湖北省","湖南省","广东省",
           "广西壮族自治区","海南省","重庆市","四川省","贵州省",
           "云南省","西藏自治区","陕西省","甘肃省","青海省",
           "宁夏回族自治区","新疆维吾尔自治区","台湾省",
           "香港特别行政区","澳门特别行政区");
pop = c(1633, 1115, 6943, 3393, 2405, 4298, 2730, 3824, 1858, 7625,
        5060, 6118, 3581, 4368, 9367, 9360, 5699, 6355, 9449,
        4768, 845, 2816, 8127, 3762, 4514, 284, 3748, 2617,
        552, 610, 2095, 2296, 693,63) * 1.1*10000
sum(pop)*1.1
pop = as.data.frame(cbind(provname,pop))
pop = pop[order(pop$provname),]
provname = provname[order(provname)]
provcol=data_last_day_cluster$`k4$cluster`
cbind(provname,data_last_day_cluster)
provcol [provcol==1] = "yellow"
provcol [provcol==2] = "blue"
provcol [provcol==3] = "red"
provcol [provcol==4] = "green"
tiff(filename = "Rplot%03d.tif",
     width = 4480, height = 4480, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400)
plot(x,col=getColor(x,provname,provcol,"white"),xlab="",ylab="");
points(city$jd,city$wd, pch = 10, col = rgb(0, 0, 0, 0.5)) #注，pch = 20所画出的实心圆点比pch = 19的要小
text(city$jd,city$wd, city$城市, cex =1)
dev.off()

nrow(city)


library(maptools)
x=readShapePoly('bou2_4p.shp')
plot(x,col=gray(924:0/924))
getColor = function(mapdata, provname, provcol, othercol){
  f = function(x, y) ifelse(x %in% y, which(y == x), 0);
  colIndex = sapply(mapdata@data$NAME, f, provname);
  fg = c(othercol, provcol)[colIndex + 1];
  return(fg);
}
provname=c("北京市","天津市","河北省","山西省","内蒙古自治区",
           "辽宁省","吉林省","黑龙江省","上海市","江苏省",
           "浙江省","安徽省","福建省","江西省","山东省",
           "河南省","湖北省","湖南省","广东省",
           "广西壮族自治区","海南省","重庆市","四川省","贵州省",
           "云南省","西藏自治区","陕西省","甘肃省","青海省",
           "宁夏回族自治区","新疆维吾尔自治区","台湾省",
           "香港特别行政区","澳门特别行政区");
plot(x, col = getColor(x, provname, provcol, "white"))


library(ggplot2)


gg <- ggplot(data_cluster, aes(x=time, y=log(cum_confirm), col=cluster)) + 
  geom_point(size=3) +  # Set color to vary based on state categories.
  geom_smooth(method="loess", size=1, se=FALSE) + 
 # coord_cartesian(xlim=c(0, 0.1), ylim=c(0, 1000000)) + 
  labs(title="时间 Vs发病数", y="发病数",
       x="时间")
plot(gg)
head(data_cluster)
ggplot(data_cluster,mapping =aes(x=time, y=log(cum_confirm)))+geom_point()+
  aes(colour=factor(data_cluster$cluster))+stat_smooth() + theme_classic()+
  labs(title="时间 Vs发病数", y="log(发病数)",x="时间")
                                                                                                                                               
ggplot(data_cluster,mapping =aes(x=time, y=log(cum_heal)))+geom_point()+
  aes(colour=factor(data_cluster$cluster))+stat_smooth() + theme_classic()+
  labs(title="时间 Vs治愈数", y="log(治愈数)",x="时间")
data_cluster$cum_dead=as.numeric(data_cluster$cum_dead)
ggplot(data_cluster,mapping =aes(x=time, y=log(cum_dead)))+geom_point()+
  aes(colour=factor(data_cluster$cluster))+stat_smooth() + theme_classic()+
  labs(title="时间 Vs死亡数", y="log(死亡数)",x="时间")


# models <- by_master_index %>% do(mod = lm( RESULTS~order+age, data = .))
# 
# 
# 
data_cluster = merge(data,cluster,by.x="province",by.y="province")
# 
# data_4 = data_cluster[data_cluster$cluster==4,]
# data_3 = data_cluster[data_cluster$cluster==3,]
# 
# plot(data_4$time,data_4$cum_confirm)
# plot(data_3$time,data_3$cum_confirm)
# 
# 
# lm( cum_confirm~time, data =data_4)
# summary(lm( cum_confirm~time, data =data_4))
# summary(lm( cum_confirm~time, data =data_3))

install.packages("deSolve")
library(deSolve) 
library(ggplot2)

seir<-function(time, state, pars){ 
  with(as.list(c(state, pars)),{ 
    dS <-- S * beta * I/N 
    dE <- S * beta * I/N - E * k 
    dI <- E * k - I * (mu + gamma) 
    dR <- I * gamma
    dN <- dS + dE + dI + dR 
    
    list(c(dS,dE,dI,dR,dN)) 
  }) 
} 

N <- 1.9E8 # 总人口
I0 <- 89 # 期初感染数
E0 <- 0 # 期初潜伏数
RM0 <- 0 # 期初移除数
S0 = N - I0 - RM0 # 期初易感人数
init<-c(S = S0, E = E0, I = I0, R = RM0, N = N)	
time <- seq(0, 150, 1) 
pars<-c( 
  beta = 0.55,	#有效接触率
  k = 1,	#潜伏到感染的转化率 
  gamma = 0.2,	#RECOVERY 
  mu=0.02	#感染期死亡率 
) 

res.seir<-as.data.frame(lsoda(y = init, times = time, func = seir, parms = pars)) 

ggplot(res.seir) +
  geom_line(aes(x = time, y = S, col = '2 易感'))+
  geom_line(aes(x = time, y = E, col = '3 潜伏'))+
  geom_line(aes(x = time, y = I, col = '4 感染'))+
  geom_line(aes(x = time, y = R, col = '5 移除'))+
  geom_line(aes(x = time, y = N, col = '1 人口'))+
  theme_light(base_family = 'Kai') +
  scale_colour_manual("",
                      values=c(
                        "2 易感" = "cornflowerblue", "3 潜伏" = "orange",
                        "4 感染" = "darkred", "5 移除" = "forestgreen", 
                        "1 人口" = "black"
                      )
  ) +
  scale_y_continuous('')

ggplot(res.seir) + geom_line(aes(x = time, y = log(I), col = '4 感染'))

#########################################model
data_cluster

by_cluster_4 <- (data_cluster[data_cluster$cluster==1,]) # random problem
by_cluster_4 = by_cluster_4[order(by_cluster_4$time),]
nrow(by_cluster_4)

pop_cluster = cbind(cluster,pop)

N <- as.numeric(as.character(pop_cluster[pop_cluster$cluster==1,]$pop)) # 总人口
I0 <- 41 # 期初感染数,第一个有记录的时间
E0 <- 0 # 期初潜伏数
RM0 <- 0 # 期初移除数
S0 = N - I0 - RM0 # 期初易感人数
init<-c(S = S0, E = E0, I = I0, R = RM0, N = N)	
time <- seq(1,nrow(by_cluster_4), 1) 
#time <- seq(1,200, 1) 
pars<-c( 
  beta = 0.57,	#有效接触率
  k = 1,	#潜伏到感染的转化率 
  gamma = 0.2,	#RECOVERY 
  mu=0.03	#感染期死亡率 
) 

res.seir<-as.data.frame(lsoda(y = init, times = time, func = seir, parms = pars)) 

ggplot() + geom_line(aes(x = by_cluster_4$time, y = I, col = '预测'),data = res.seir)+
  geom_line(aes(x = by_cluster_4$time, y = by_cluster_4$cum_confirm, col = '实际'))+ 
  theme_classic()+
  labs(title="4类城市", y="确诊数量",x="时间")






length(res.seir$I)
length(by_cluster_4$cum_confirm)


by_cluster_3 <- group_by(data_cluster[data_cluster$cluster==2,], time =data_cluster[data_cluster$cluster==2,]$time)

head(by_cluster_3)


by_cluster_3_1 <- by_cluster_3 %>% 
  summarise(
     cum_confirm=sum(cum_confirm) )


N <- sum(as.numeric(as.character(pop_cluster[pop_cluster$cluster==2,]$pop))) # 总人口
I0 <- 41 # 期初感染数,第一个有记录的时间
E0 <- 0 # 期初潜伏数
RM0 <- 0 # 期初移除数
S0 = N - I0 - RM0 # 期初易感人数
init<-c(S = S0, E = E0, I = I0, R = RM0, N = N)	
time <- seq(1,nrow(by_cluster_3_1), 1) 
#time <- seq(1,200, 1) 
pars<-c( 
  beta = 0.49,	#有效接触率
  k = 1,	#潜伏到感染的转化率 
  gamma = 0.2,	#RECOVERY 
  mu=0.01	#感染期死亡率 
) 

res.seir<-as.data.frame(lsoda(y = init, times = time, func = seir, parms = pars)) 

ggplot() + geom_line(aes(x = by_cluster_3_1$time, y = I, col = '预测'),data = res.seir)+
  geom_line(aes(x = by_cluster_3_1$time, y = by_cluster_3_1$cum_confirm, col = '实际'))+ 
  theme_classic()+
  labs(title="3类城市", y="确诊数量",x="时间")


by_cluster_1 <- group_by(data_cluster[data_cluster$cluster==3,], time =data_cluster[data_cluster$cluster==3,]$time)

head(by_cluster_1)


by_cluster_1_1 <- by_cluster_1 %>% 
  summarise(
    cum_confirm=sum(cum_confirm) )



N <- sum(as.numeric(as.character(pop_cluster[pop_cluster$cluster==3,]$pop))) # 总人口
I0 <- 41 # 期初感染数,第一个有记录的时间
E0 <- 0 # 期初潜伏数
RM0 <- 0 # 期初移除数
S0 = N - I0 - RM0 # 期初易感人数
init<-c(S = S0, E = E0, I = I0, R = RM0, N = N)	
time <- seq(1,nrow(by_cluster_1_1), 1) 
#time <- seq(1,200, 1) 
pars<-c( 
  beta = 0.475,	#有效接触率
  k = 1,	#潜伏到感染的转化率 
  gamma = 0.2,	#RECOVERY 
  mu=0.005	#感染期死亡率 
) 

res.seir<-as.data.frame(lsoda(y = init, times = time, func = seir, parms = pars)) 

ggplot() + geom_line(aes(x = by_cluster_1_1$time, y = I, col = '预测'),data = res.seir)+
  geom_line(aes(x = by_cluster_1_1$time, y = by_cluster_1_1$cum_confirm, col = '实际'))+ 
  theme_classic()+
  labs(title="1类城市", y="确诊数量",x="时间")

by_cluster_2 <- group_by(data_cluster[data_cluster$cluster==4,], time =data_cluster[data_cluster$cluster==4,]$time)

head(by_cluster_2)




by_cluster_2_1 <- by_cluster_2 %>% 
  summarise(
    cum_confirm=sum(cum_confirm) )



N <- sum(as.numeric(as.character(pop_cluster[pop_cluster$cluster==4,]$pop))) # 总人口
I0 <- 41 # 期初感染数,第一个有记录的时间
E0 <- 0 # 期初潜伏数
RM0 <- 0 # 期初移除数
S0 = N - I0 - RM0 # 期初易感人数
init<-c(S = S0, E = E0, I = I0, R = RM0, N = N)	
time <- seq(1,nrow(by_cluster_2_1), 1) 
#time <- seq(1,200, 1) 
pars<-c( 
  beta = 0.5,	#有效接触率
  k = 1,	#潜伏到感染的转化率 
  gamma = 0.2,	#RECOVERY 
  mu=0.005	#感染期死亡率 
) 

res.seir<-as.data.frame(lsoda(y = init, times = time, func = seir, parms = pars)) 

ggplot() + geom_line(aes(x = by_cluster_2_1$time, y = I, col = '预测'),data = res.seir)+
  geom_line(aes(x = by_cluster_2_1$time, y = by_cluster_2_1$cum_confirm, col = '实际'))+ 
  theme_classic()+
  labs(title="2类城市", y="确诊数量",x="时间")

2.85

                             