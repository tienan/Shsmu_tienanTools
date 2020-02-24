#
install.packages("MatchIt")
library(MatchIt)
?matchit
getwd()
mydata <- read.table("../R/ruijinDepart.txt",sep = "\t",header = T)
attach(mydata)
head(mydata)

mydata = na.omit(mydata)
m.out = matchit ( Group ~ Age+BMI,data = mydata,method ="nearest",ratio = 2)

summary(m.out)

plot (m.out, type = "jitter")
plot (m.out, type = "hist")

m.data1 <- match.data (m.out)

boxplot(m.data1$Age~m.data1$Group)
boxplot(m.data1$BMI~m.data1$Group)
write.csv(m.data1,"BoneOpe.csv")
