setwd("D:/project/PM2.5和脂代谢")
library(readr)
TCTGLUNGDEMO_IN <- read_csv("D:/project/PM2.5和脂代谢/TGTCLUNG/TCTGLUNGDEMO_IN.csv",col_names = F)
PM2.5=read.csv("D:/project/PM2.5和脂代谢/脂代谢xPM25/pm2.5_.csv")
province=read.csv("D:/project/PM2.5和脂代谢/脂代谢xPM25/province.csv")
clinical_in=TCTGLUNGDEMO_IN[grep(pattern = "肺癌",TCTGLUNGDEMO_IN$X3),]
province=province[,c(1,5)]
province$Province=gsub(" ", "", province$Province)
table(province$Province)
names(PM2.5)[1]="Province"
Pro_PM=merge(province,PM2.5,by="Province")
TCTGLUNGDEMO_Out <- read_csv("D:/project/PM2.5和脂代谢/TGTCLUNG/TCTGLUNGDEMO_Out.csv",col_names = F)
lipid=read.csv("D:/project/PM2.5和脂代谢/TGTCLUNG/Demo_COPDLUNG.csv")
library(tidyr)
tctg_in=merge(clinical_in,lipid,by.x = "X1",by.y ="REG_CODE")
tctg_in=tctg_in[complete.cases(tctg_in$RESULTS),]
library(sqldf)
tctg_in$Age=as.numeric(round((as.Date(substring(tctg_in$REPORT_TIME,1,10))-as.Date(substring(tctg_in$X11,1,10)))/365))
tctg_in=tctg_in[,-c(2:4,6:9,14,16,17,19)]
names(tctg_in)[2]="ID"
library(dplyr)
all_in=left_join(tctg_in,Pro_PM)
names(all_in)[13:30]=substr(colnames(all_in)[13:30],2,5)
all_in=unique(all_in)
all_in$Age=as.numeric(all_in$Age)
all_in$Gender=ifelse(all_in$X12=="女","Female","Male")
all_in$RESULTS=as.numeric(all_in$RESULTS)
tc_in=all_in[all_in$ITEM_ENAME=="TC",]
tg_in=all_in[all_in$ITEM_ENAME=="TG",]
tc_in=arrange(tc_in, X13, REPORT_TIME)
tg_in=arrange(tg_in, X13, REPORT_TIME)
tc_in_du=sqldf("select * from tc_in t where (select count(1) from tc_in where X13 = t.X13 )>1")
tg_in_du=sqldf("select * from tg_in t where (select count(1) from tg_in where X13 = t.X13 )>1")
tc_in_un=tc_in[!duplicated(tc_in$X13),]
tg_in_un=tg_in[!duplicated(tg_in$X13),]
tc_in_un$group=ifelse(tc_in_un$均值>quantile(tc_in_un$均值,.75),"Q4",
                   ifelse(tc_in_un$均值>quantile(tc_in_un$均值,.5),"Q3",
                          ifelse(tc_in_un$均值>quantile(tc_in_un$均值,.25),"Q2","Q1")))
tg_in_un$group=ifelse(tg_in_un$均值>quantile(tg_in_un$均值,.75),"Q4",
                      ifelse(tg_in_un$均值>quantile(tg_in_un$均值,.5),"Q3",
                             ifelse(tg_in_un$均值>quantile(tg_in_un$均值,.25),"Q2","Q1")))
aggregate(tg_in_un[,c(9,38)], by=list(tg_in_un$group), FUN=mean)
library(lattice)
library(MASS)
bwplot(tg_in_un$RESULTS~tg_in_un$group)
model=aov(tg_in_un$RESULTS~tg_in_un$group)
summary(model)
result=TukeyHSD(model)
plot(result)


model=lm(RESULTS~group+X10+X12+Age,tg_in_un)
summary(model)







######################
ggplot(data  = tc_in_un,
       aes(x = 均值,
           y = RESULTS))+
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ 
  geom_smooth(method = lm,
              se     = FALSE, 
              col    = "black",
              size   = .5, 
              alpha  = .8)+ 
  theme_minimal()+
  labs(title    = "PM2.5 vs. TC",
       subtitle = "add regression line")

ggplot(data  = tg_in_un,
       aes(x = group,
           y = RESULTS))+
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ 
  geom_smooth(method = lm,
              se     = FALSE, 
              col    = "black",
              size   = .5, 
              alpha  = .8)+ 
  theme_minimal()+
  labs(title    = "PM2.5 vs. TG",
       subtitle = "add regression line")

plot(tg_in_un$均值,tg_in_un$RESULTS)
library(mgcv) 
summary(lm(tg_in_un$RESULTS~tg_in_un$均值))


write.csv(tc_in_du,"D:/project/PM2.5和脂代谢/result/tc_in_du.csv")
write.csv(tg_in_du,"D:/project/PM2.5和脂代谢/result/tg_in_du.csv")
tc_in_du=read.csv("D:/project/PM2.5和脂代谢/result/tc_in_du.csv")
tc_in_du$Day=as.numeric(tc_in_du$Day)
library(nlme)
m1.nlme = lme(RESULTS ~ 均值,random=~1|X13,data = tc_in_du)
anova(m1.nlme)
library(dplyr)
tc_in_du$X13=as.factor(tc_in_du$X13)
tc_in_du$group=as.factor(tc_in_du$group)
library(ez)
model <- ezANOVA(tc_in_du, dv = RESULTS, wid = X13, within = group, detailed = T)
ezDesign(tc_in_du$RESULTS)
install.packages("remotes")
remotes::install_github("tianshu129/foqat")
library(foqat)



library(tableone)
stable <- CreateTableOne(vars=c("Age","mean","Province","RESULTS","Gender","X10"),
                         strata="group",data=all_in_TG,
                         factorVars=c("Gender","Province","X10"))
stable=print(stable,showAllLevels = TRUE)
stable <- CreateTableOne(vars=c("Age","mean","Province","RESULTS","Gender","X10","group"),
                         data=all_in_TG,factorVars=c("Gender","Province","X10","group"))
stable=print(stable,showAllLevels = TRUE)
stable <- CreateTableOne(vars=c("Age","mean","Province","RESULTS","Gender","X10"),
                         strata="group",data=all_in_TC,
                         factorVars=c("Gender","Province","X10"))
stable=print(stable,showAllLevels = TRUE)
stable <- CreateTableOne(vars=c("Age","mean","Province","RESULTS","Gender","X10","group"),
                         data=all_in_TC,factorVars=c("Gender","Province","X10","group"))
stable=print(stable,showAllLevels = TRUE)

cor.test(x = all_in_TG$RESULTS, y = all_in_TG$mean)
cor.test(x = all_in_TC$RESULTS, y = all_in_TC$mean)
library(mgcv)
model <- gam(RESULTS ~s(mean), data=all_in_TG)
plot(model)
summary(model)
gam.check(model)
plot(all_in_TG$RESULTS,all_in_TG$mean,col = "blue",
     abline(lm(mean~RESULTS,data=all_in_TG)),cex = 1.3,pch = 16,
     xlab = "PM2.5",ylab = "TG")
library(car)
library("psych")
cor_matrix <- corr.test(all_in_TG$RESULTS,all_in_TG$mean)
str(cor_matrix)
