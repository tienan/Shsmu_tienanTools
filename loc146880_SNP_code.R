dat = read.table("loc146880_SNP_clinical.txt",header = T,sep = "\t")
head(dat)
dat$lifestatus
dat$survivaltime
dat$stage
?coxph
Surv(dat$survivaltime,dat$lifestatus)
library(survival)
sfit <- coxph(Surv(dat$survivaltime,dat$lifestatus)
               ~ dat$stage+dat$SNP,dat)
summary(sfit)


sfit <- coxph(Surv(dat$survivaltime,dat$lifestatus)
              ~ dat$SNP,dat)
summary(sfit)
#stage 
?table()
table(dat$stage,dat$SNP)

dat$SNP
