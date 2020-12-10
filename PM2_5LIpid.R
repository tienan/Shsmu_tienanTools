#setwd("D:/project/PM2.5??֬??л")
#install.packages("readr")
library(readr)
TCTGLUNGDEMO_IN <- read_csv("../pm2.5/TCTGLUNGDEMO_IN.csv",col_names = F)
PM2.5=read.csv("../pm2.5/PM2_5")
head(PM2.5)
province=read.csv("../pm2.5/province.csv")
head(province)
clinical_in=TCTGLUNGDEMO_IN[grep(pattern = "",TCTGLUNGDEMO_IN$X3),]
province=province[,c(1,5)]
province$Province=gsub(" ", "", province$Province)
table(province$Province)
names(PM2.5)[1]="Province"
Pro_PM=merge(province,PM2.5,by="Province")
head(Pro_PM)


clinical_out=TCTGLUNGDEMO_Out[grep(pattern = "",TCTGLUNGDEMO_Out$X3),]


TCTGLUNGDEMO_Out <- read_csv("../pm2.5/TCTGLUNGDEMO_Out.csv",col_names = F)


tctg_in=merge(clinical_in,lipid,by.x = "X1",by.y ="REG_CODE")

tctg_out=merge(clinical_out,lipid,by.x = "X1",by.y ="REG_CODE")

tctg_out$X4 = as.Date(tctg_out$X4)



tctg = rbind(tctg_out,tctg_in)
nrow(tctg)


tctg=tctg[complete.cases(tctg$RESULTS),]
#install.packages("sqldf")
library(sqldf)
tctg$Age=as.numeric(round((as.Date(substring(tctg$REPORT_TIME,1,10))-as.Date(substring(tctg$X11,1,10)))/365))
tctg=tctg[,-c(2:4,6:9,14,16,17,19)]
names(tctg)[2]="ID"
#install.packages("dplyr")
library(dplyr)
all=left_join(tctg,Pro_PM)
head(all)
#names(all_in)[13:30]=substr(colnames(all_in)[13:30],2,5)
all=unique(all)
all$Age=as.numeric(all$Age)
all$Gender=ifelse(all$X12=="女","Female","Male")
all$RESULTS=as.numeric(all$RESULTS)
tc=all[all$ITEM_ENAME=="TC",]
tg=all[all$ITEM_ENAME=="TG",]
tc=arrange(tc, X13, REPORT_TIME)
tg=arrange(tg, X13, REPORT_TIME)
tc_du=sqldf("select * from tc t where (select count(1) from tc where X13 = t.X13 )>1")
head(tc_du)
tg_du=sqldf("select * from tg t where (select count(1) from tg where X13 = t.X13 )>1")
head(tg_du)


tc_du$order =c(1:nrow(tc_du))
tc_du = na.omit(tc_du)
tc_du <- group_by(tc_du, X13)


models <-tc_du %>% do(mod = lm( RESULTS~order+Age, data = .))
Sign = sign(as.data.frame(summarise(models, rsq = as.data.frame(summary(mod)$coefficients)[2,1])))

tc_du$sign = Sign


tc_du_1 <- tc_du %>% 
  summarise(median= median(RESULTS),firstValue=first(RESULTS),lastValue= last(RESULTS),
            age=min(Age),ORG_CODE=first(X13),gender=first(Gender),PROVINCE=last(Province),Level = mean(Level))
tc_du_1$sign = Sign

table(tc_du_1$sign)

boxplot(Sign~as.numeric(tc_du_1$Level))


table(tc_du_1$sign$rsq,tc_du_1$PROVINCE)






tg_du$order =c(1:nrow(tg_du))
tg_du = na.omit(tg_du)
tg_du <- group_by(tg_du, X13)


models <-tg_du %>% do(mod = lm( RESULTS~order+Age, data = .))
Sign = sign(as.data.frame(summarise(models, rsq = as.data.frame(summary(mod)$coefficients)[2,1])))

tg_du_1 <- tg_du %>% 
  summarise(median= median(RESULTS),firstValue=first(RESULTS),lastValue= last(RESULTS),
            age=min(Age),ORG_CODE=first(X13),gender=first(Gender),PROVINCE=last(Province),Level = mean(Level))
tg_du_1$sign = Sign

table(tg_du_1$sign)

boxplot(Sign~as.numeric(tc_du_1$Level))


table(tg_du_1$sign$rsq,tg_du_1$PROVINCE)


table(tg_du_1$sign$rsq,tg_du_1$Level)





TCTGLUNGDEMO_IN$X4 = as.Date(TCTGLUNGDEMO_IN$X4)
TCTGLUNGDEMO_Out

TCTGLUNGDEMO = rbind(TCTGLUNGDEMO_Out, TCTGLUNGDEMO_IN)
nrow(TCTGLUNGDEMO )
clinical_in=TCTGLUNGDEMO[grep(pattern = "",TCTGLUNGDEMO_IN$X3),]
tctg=merge(clinical_in,lipid,by.x = "X1",by.y ="REG_CODE")
nrow(tctg)

tctg$Age=as.numeric(round((as.Date(substring(tctg$REPORT_TIME,1,10))-as.Date(substring(tctg$X11,1,10)))/365))
tctg=tctg[,-c(2:4,6:9,14,16,17,19)]


names(tctg)[2]="ID"

all=left_join(tctg,Pro_PM)
head(all)


all=unique(all)
all$Age=as.numeric(all$Age)
all$Gender=ifelse(all$X12=="女","Female","Male")
all$RESULTS=as.numeric(all$RESULTS)
tc=all[all$ITEM_ENAME=="TC",]
tg=all[all$ITEM_ENAME=="TG",]
tc=arrange(tc, X13, REPORT_TIME)
tg=arrange(tg, X13, REPORT_TIME)
tc_du=sqldf("select * from tc t where (select count(1) from tc where X13 = t.X13 )>1")
head(tc_du)
nrow(tc_du)
tg_du=sqldf("select * from tg t where (select count(1) from tg where X13 = t.X13 )>1")
head(tg_du)
nrow(tg_du)

1

basicData <- by_master_index %>% 
  summarise(median= median(RESULTS),firstValue=first(RESULTS),lastValue= last(RESULTS),
            age=min(age),ORG_CODE=first(ORG_CODE),gender=first(sex),PROVINCE=last(PROVINCE))
res = basicData


clinical_in=TCTGLUNGDEMO[grep(pattern = "",TCTGLUNGDEMO_IN$X3),]




lipid=read.csv("../pm2.5/Demo_COPDLUNG.csv")
head(lipid)
#install.packages("tidyr")
library(tidyr)
tctg_in=merge(clinical_in,lipid,by.x = "X1",by.y ="REG_CODE")
nrow(tctg_in)
tctg_in=tctg_in[complete.cases(tctg_in$RESULTS),]
#install.packages("sqldf")
library(sqldf)
tctg_in$Age=as.numeric(round((as.Date(substring(tctg_in$REPORT_TIME,1,10))-as.Date(substring(tctg_in$X11,1,10)))/365))
tctg_in=tctg_in[,-c(2:4,6:9,14,16,17,19)]
names(tctg_in)[2]="ID"
#install.packages("dplyr")
library(dplyr)
all_in=left_join(tctg_in,Pro_PM)
head(all_in)
#names(all_in)[13:30]=substr(colnames(all_in)[13:30],2,5)
all_in=unique(all_in)
all_in$Age=as.numeric(all_in$Age)
all_in$Gender=ifelse(all_in$X12=="女","Female","Male")
all_in$RESULTS=as.numeric(all_in$RESULTS)
tc_in=all_in[all_in$ITEM_ENAME=="TC",]
tg_in=all_in[all_in$ITEM_ENAME=="TG",]
tc_in=arrange(tc_in, X13, REPORT_TIME)
tg_in=arrange(tg_in, X13, REPORT_TIME)
tc_in_du=sqldf("select * from tc_in t where (select count(1) from tc_in where X13 = t.X13 )>1")
head(tc_in_du)
tg_in_du=sqldf("select * from tg_in t where (select count(1) from tg_in where X13 = t.X13 )>1")
head(tg_in_du)



basicData <- by_master_index %>% 
  summarise(median= median(RESULTS),firstValue=first(RESULTS),lastValue= last(RESULTS),
            age=min(age),ORG_CODE=first(ORG_CODE),gender=first(sex),PROVINCE=last(PROVINCE))
res = basicData









tc_in_un=tc_in[!duplicated(tc_in$X13),]
tg_in_un=tg_in[!duplicated(tg_in$X13),]
tc_in_un$group=ifelse(tc_in_un$??ֵ>quantile(tc_in_un$??ֵ,.75),"Q4",
                   ifelse(tc_in_un$??ֵ>quantile(tc_in_un$??ֵ,.5),"Q3",
                          ifelse(tc_in_un$??ֵ>quantile(tc_in_un$??ֵ,.25),"Q2","Q1")))
tg_in_un$group=ifelse(tg_in_un$??ֵ>quantile(tg_in_un$??ֵ,.75),"Q4",
                      ifelse(tg_in_un$??ֵ>quantile(tg_in_un$??ֵ,.5),"Q3",
                             ifelse(tg_in_un$??ֵ>quantile(tg_in_un$??ֵ,.25),"Q2","Q1")))
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
       aes(x = ??ֵ,
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

plot(tg_in_un$??ֵ,tg_in_un$RESULTS)
library(mgcv) 
summary(lm(tg_in_un$RESULTS~tg_in_un$??ֵ))


write.csv(tc_in_du,"D:/project/PM2.5??֬??л/result/tc_in_du.csv")
write.csv(tg_in_du,"D:/project/PM2.5??֬??л/result/tg_in_du.csv")
tc_in_du=read.csv("D:/project/PM2.5??֬??л/result/tc_in_du.csv")
tc_in_du$Day=as.numeric(tc_in_du$Day)
library(nlme)
m1.nlme = lme(RESULTS ~ ??ֵ,random=~1|X13,data = tc_in_du)
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
