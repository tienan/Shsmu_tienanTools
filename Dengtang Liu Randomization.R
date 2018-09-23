# Dengtang Liu Randomization

install.packages("randomizeR")
library("randomizeR")

params <- crPar(100)
rs <- genSeq(params)
rs
getRandList(rs)
plotSeq(rs, plotAllSeq = TRUE)
allSeqs <- getAllSeq(params)



install.packages("randomizr")
library("randomizr")
data(HairEyeColor)
HairEyeColor <- data.frame(HairEyeColor)
hec <- HairEyeColor[rep(1:nrow(HairEyeColor),
                        times = HairEyeColor$Freq), 1:3]

N <- nrow(hec)

# Fix the rownames
rownames(hec) <- NULL


# Set a seed for reproducability
set.seed(343)

# Create untreated and treated outcomes for all subjects
hec <- within(hec,{
  Y0 <- rnorm(n = N,mean = (2*as.numeric(Hair) + -4*as.numeric(Eye) + -6*as.numeric(Sex)), sd = 5)
  Y1 <- Y0 + 6*as.numeric(Hair) + 4*as.numeric(Eye) + 2*as.numeric(Sex)
})

# Calculate true ATE
with(hec, mean(Y1 - Y0))
#> [1] 25


Z <- simple_ra(N = N)
table(Z)

hec$Hair

hec <- within(hec,{
  Z_blocked <- complete_ra(N = N, m_each = c(100, 200, 292),
                           conditions = c("control", "placebo", "treatment"))
  id_var <- 1:nrow(hec)
})



#https://cran.r-project.org/web/packages/randomizr/vignettes/randomizr_vignette.html


#T1: 奥痰平；T2：利培酮；T3：氨磺必利；T4：阿立哌；T5：奋乃静
#7个中心，竞争入组
#1600个患者
?seq
hospital = c(rep(1,time=228),rep(2,time=228),rep(3,time=228),rep(4,time=228),rep(5,time=228),rep(6,time=228),rep(7,time=228))


set.seed(343)
Z <- block_ra(blocks = hospital, num_arms = 5)
table(Z,hospital)
Header = "SMART"
#generate 5 Number

set.seed(343)
a3 = sample(1:10000,length(hospital))
no=paste(Header,"_",a3,sep = "")

write.csv(file = "tableRandom.csv",x = cbind(no,hospital,Z),sep = "/t")


