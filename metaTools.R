#mate analysis

install.packages("mada")
library("mada")
data(AuditC)
head(AuditC)

# univariate meta-analysis
negLR.DSL <- madauni(AuditC, type = "negLR", method = "DSL")
summary(negLR.DSL)
forest(negLR.DSL)
negLR.MH <- madauni(AuditC, type = "negLR", method = "MH")
summary(negLR.MH)
forest(negLR.MH)
posLR.DSL <- madauni(AuditC, type = "posLR", method = "DSL")
forest(negLR.MH)
summary(posLR.DSL)
forest(posLR.DSL)
posLR.MH <- madauni(AuditC, type = "posLR", method = "MH")
summary(posLR.MH)
forest(posLR.MH)

# bivariate meta-analysis
fit <- reitsma(AuditC)
mcmc_sum <- SummaryPts(fit, n.iter = 10^6)
summary(mcmc_sum)


dat = read.table("../R/METACEUS.txt",header = T,sep = "\t")
dat = dat[order(dat[,1]),]
dat$names=dat[,1]
madad(dat)
forest(madad(dat[,2:6]), type = "spec",snames=dat$ID,xlab = "Specificity",main="")
forest(madad(dat[,2:6]), type = "sens",snames=dat$ID,xlab = "Sensitivity",main="")

ROCellipse(dat[,2:6], pch = "")
points(fpr(dat[,2:6]), sens(dat[,2:6]))

rs <- rowSums(dat[,2:5])
weights <- 4 * rs / max(rs)
crosshair(dat[,2:5], xlim = c(0,1), ylim = c(0,1), col = 1:14, lwd = weights)


fit.DOR.MH <- madauni(dat[,2:5], method = "MH")
forest(fit.DOR.MH)

fit.phm.homo <- phm(dat[,2:5], hetero = FALSE)
summary(fit.phm.homo)
ROCellipse(dat[,2:5], add = F)


fit.reitsma <- reitsma(dat[,2:5])
summary(fit.reitsma)


plot(fit.reitsma, sroclwd = 2,main = "")



rs <- rowSums(dat[,2:5])
weights <- 4 * rs / max(rs)

crosshair(dat[,2:6], xlim = c(0,0.6), ylim = c(0.4,1), col = 1:14, lwd = weights)
ROCellipse(dat[,2:6], pch = "")
points(fpr(dat[,2:5]), sens(dat[,2:5]))

fit.DOR.DSL <- madauni(dat[,2:5])
forest(fit.DOR.DSL)
