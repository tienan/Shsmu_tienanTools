#######
r1 <- round(rnorm(20, 10, 4))
r2 <- round(r1 + 10 + rnorm(20, 0, 2))
r3 <- round(r1 + 20 + rnorm(20, 0, 2))
icc(cbind(r1, r2, r3), "twoway")
data(ChickWeight)
ChickWeight

data(anxiety)
icc(anxiety, model="twoway", type="agreement")

install.packages("irr")
library("irr")

data(anxiety)
icc(anxiety, model="twoway", type="agreement")
?icc
r1 <- round(rnorm(20, 10, 4))
r2 <- round(r1 + 10 + rnorm(20, 0, 2))
r3 <- round(r1 + 20 + rnorm(20, 0, 2))
icc(cbind(r1, r1), "oneway","agreement")              # High consistency
icc(cbind(r1, r2, r3), "twoway", "agreement") # Low agreement


r1 <- round(rnorm(20, 10, 4))
r2 <- round(r1 + 10 + rnorm(20, 0, 2))
r3 <- round(r1 + 20 + rnorm(20, 0, 2))
icc(cbind(r1, r2, r3), "twoway")
r = icc(cbind(r1, r1), "oneway")  
r$value
r$p


#########
install.packages("ICC")
library("ICC")
?ICCest
r1 = read.table("xinhua/R1_1person.txt",header = T,sep = "\t")
r2 = read.table("xinhua/R1_2person.txt",header = T,sep = "\t")
r=c()
p=c()
m=c()
s=c()
for (i in 1:ncol(r1)){
  a = icc(cbind(r1[,i],r2[,i]), "oneway")
  r[i] = a$value
  p[i] = a$p.value
  m[i] = mean(r1[,i]-r2[,i])
  s[i] = sd(r1[,i]-r2[,i])
}

res = rbind(r,p,m,s)
source("Rgenetic/basicTools.R")
writeOut(res,"tmp")




r=c()
p=c()
m=c()
s=c()
for (i in 1:ncol(r1)){
  a = cor.test(r1[,i],r2[,i])
  r[i] = a$estimate
  p[i] = a$p.value
  m[i] = mean(r1[,i]-r2[,i])
  s[i] = sd(r1[,i]-r2[,i])
}

res = rbind(r,p,m,s)
source("Rgenetic/basicTools.R")
writeOut(res,"tmp")

#select principle 

sign = read.table("xinhua/tmp.txt",header = T,sep = "\t")

dat = cbind(sign,r2)

# MANOVA test
res.man <- manova(cbind(MinIntensity,                                         
                        MaxIntensity,                                         
                        MedianIntensity,                                      
                        MeanValue,                                            
                        stdDeviation,                                         
                        Variance,                                             
                        VolumeCount,                                          
                        VoxelValueSum,                                        
                        Range,                                                
                        RMS,                                                  
                        MeanDeviation,                                        
                        RelativeDeviation,                                    
                        skewness,                                             
                        kurtosis,                                             
                        uniformity,                                           
                        histogramEnergy,                                      
                        histogramEntropy,                                     
                        FrequencySize,                                        
                        Percentile5,                                          
                        Percentile10,                                         
                        Percentile15,                                         
                        Percentile20,                                         
                        Percentile25,                                         
                        Percentile30,                                         
                        Percentile35,                                         
                        Percentile40,                                         
                        Percentile45,                                         
                        Percentile50,                                         
                        Percentile55,                                         
                        Percentile60,                                         
                        Percentile65,                                         
                        Percentile70,                                         
                        Percentile75,                                         
                        Percentile80,                                         
                        Percentile85,                                         
                        Percentile90,                                         
                        Percentile95,                                         
                        Quantile0.025,                                        
                        Quantile0.25,                                         
                        Quantile0.5,                                          
                        Quantile0.75,                                         
                        Quantile0.975,                                        
                        GLCMEnergy_AllDirection_offset1,                      
                        GLCMEntropy_AllDirection_offset1,                     
                        Inertia_AllDirection_offset1,                         
                        Correlation_AllDirection_offset1,                     
                        InverseDifferenceMoment_AllDirection_offset1,         
                        ClusterShade_AllDirection_offset1,                    
                        ClusterProminence_AllDirection_offset1,               
                        HaralickCorrelation_AllDirection_offset1,             
                        GLCMEnergy_AllDirection_offset1_SD,                   
                        GLCMEntropy_AllDirection_offset1_SD,                  
                        Inertia_AllDirection_offset1_SD,                      
                        Correlation_AllDirection_offset1_SD,                  
                        InverseDifferenceMoment_AllDirection_offset1_SD,      
                        ClusterShade_AllDirection_offset1_SD,                 
                        ClusterProminence_AllDirection_offset1_SD,            
                        HaralickCorrelation_AllDirection_offset1_SD,          
                        GLCMEnergy_angle0_offset1,                            
                        GLCMEntropy_angle0_offset1,                           
                        Inertia_angle0_offset1,                               
                        Correlation_angle0_offset1,                           
                        InverseDifferenceMoment_angle0_offset1,               
                        ClusterShade_angle0_offset1,                          
                        ClusterProminence_angle0_offset1,                     
                        HaralickCorrelation_angle0_offset1,                   
                        GLCMEnergy_angle45_offset1,                           
                        GLCMEntropy_angle45_offset1,                          
                        Inertia_angle45_offset1,                              
                        Correlation_angle45_offset1,                          
                        InverseDifferenceMoment_angle45_offset1,              
                        ClusterShade_angle45_offset1,                         
                        ClusterProminence_angle45_offset1,                    
                        HaralickCorrelation_angle45_offset1,                  
                        GLCMEnergy_angle90_offset1,                           
                        GLCMEntropy_angle90_offset1,                          
                        Inertia_angle90_offset1,                              
                        Correlation_angle90_offset1,                          
                        InverseDifferenceMoment_angle90_offset1,              
                        ClusterShade_angle90_offset1,                         
                        ClusterProminence_angle90_offset1,                    
                        HaralickCorrelation_angle90_offset1,                  
                        GLCMEnergy_angle135_offset1,                          
                        GLCMEntropy_angle135_offset1,                         
                        Inertia_angle135_offset1,                             
                        Correlation_angle135_offset1,                         
                        InverseDifferenceMoment_angle135_offset1,             
                        ClusterShade_angle135_offset1,                        
                        ClusterProminence_angle135_offset1,                   
                        HaralickCorrelation_angle135_offset1,                 
                        GLCMEnergy_AllDirection_offset4,                      
                        GLCMEntropy_AllDirection_offset4,                     
                        Inertia_AllDirection_offset4,                         
                        Correlation_AllDirection_offset4,                     
                        InverseDifferenceMoment_AllDirection_offset4,         
                        ClusterShade_AllDirection_offset4,                    
                        ClusterProminence_AllDirection_offset4,               
                        HaralickCorrelation_AllDirection_offset4,             
                        GLCMEnergy_AllDirection_offset4_SD,                   
                        GLCMEntropy_AllDirection_offset4_SD,                  
                        Inertia_AllDirection_offset4_SD,                      
                        Correlation_AllDirection_offset4_SD,                  
                        InverseDifferenceMoment_AllDirection_offset4_SD,      
                        ClusterShade_AllDirection_offset4_SD,                 
                        ClusterProminence_AllDirection_offset4_SD,            
                        HaralickCorrelation_AllDirection_offset4_SD,          
                        GLCMEnergy_angle0_offset4,                            
                        GLCMEntropy_angle0_offset4,                           
                        Inertia_angle0_offset4,                               
                        Correlation_angle0_offset4,                           
                        InverseDifferenceMoment_angle0_offset4,               
                        ClusterShade_angle0_offset4,                          
                        ClusterProminence_angle0_offset4,                     
                        HaralickCorrelation_angle0_offset4,                   
                        GLCMEnergy_angle45_offset4,                           
                        GLCMEntropy_angle45_offset4,                          
                        Inertia_angle45_offset4,                              
                        Correlation_angle45_offset4,                          
                        InverseDifferenceMoment_angle45_offset4,              
                        ClusterShade_angle45_offset4,                         
                        ClusterProminence_angle45_offset4,                    
                        HaralickCorrelation_angle45_offset4,                  
                        GLCMEnergy_angle90_offset4,                           
                        GLCMEntropy_angle90_offset4,                          
                        Inertia_angle90_offset4,                              
                        Correlation_angle90_offset4,                          
                        InverseDifferenceMoment_angle90_offset4,              
                        ClusterShade_angle90_offset4,                         
                        ClusterProminence_angle90_offset4,                    
                        HaralickCorrelation_angle90_offset4,                  
                        GLCMEnergy_angle135_offset4,                          
                        GLCMEntropy_angle135_offset4,                         
                        Inertia_angle135_offset4,                             
                        Correlation_angle135_offset4,                         
                        InverseDifferenceMoment_angle135_offset4,             
                        ClusterShade_angle135_offset4,                        
                        ClusterProminence_angle135_offset4,                   
                        HaralickCorrelation_angle135_offset4,                 
                        GLCMEnergy_AllDirection_offset7,                      
                        GLCMEntropy_AllDirection_offset7,                     
                        Inertia_AllDirection_offset7,                         
                        Correlation_AllDirection_offset7,                     
                        InverseDifferenceMoment_AllDirection_offset7,         
                        ClusterShade_AllDirection_offset7,                    
                        ClusterProminence_AllDirection_offset7,               
                        HaralickCorrelation_AllDirection_offset7,             
                        GLCMEnergy_AllDirection_offset7_SD,                   
                        GLCMEntropy_AllDirection_offset7_SD,                  
                        Inertia_AllDirection_offset7_SD,                      
                        Correlation_AllDirection_offset7_SD,                  
                        InverseDifferenceMoment_AllDirection_offset7_SD,      
                        ClusterShade_AllDirection_offset7_SD,                 
                        ClusterProminence_AllDirection_offset7_SD,            
                        HaralickCorrelation_AllDirection_offset7_SD,          
                        GLCMEnergy_angle0_offset7,                            
                        GLCMEntropy_angle0_offset7,                           
                        Inertia_angle0_offset7,                               
                        Correlation_angle0_offset7,                           
                        InverseDifferenceMoment_angle0_offset7,               
                        ClusterShade_angle0_offset7,                          
                        ClusterProminence_angle0_offset7,                     
                        HaralickCorrelation_angle0_offset7,                   
                        GLCMEnergy_angle45_offset7,                           
                        GLCMEntropy_angle45_offset7,                          
                        Inertia_angle45_offset7,                              
                        Correlation_angle45_offset7,                          
                        InverseDifferenceMoment_angle45_offset7,              
                        ClusterShade_angle45_offset7,                         
                        ClusterProminence_angle45_offset7,                    
                        HaralickCorrelation_angle45_offset7,                  
                        GLCMEnergy_angle90_offset7,                           
                        GLCMEntropy_angle90_offset7,                          
                        Inertia_angle90_offset7,                              
                        Correlation_angle90_offset7,                          
                        InverseDifferenceMoment_angle90_offset7,              
                        ClusterShade_angle90_offset7,                         
                        ClusterProminence_angle90_offset7,                    
                        HaralickCorrelation_angle90_offset7,                  
                        GLCMEnergy_angle135_offset7,                          
                        GLCMEntropy_angle135_offset7,                         
                        Inertia_angle135_offset7,                             
                        Correlation_angle135_offset7,                         
                        InverseDifferenceMoment_angle135_offset7,             
                        ClusterShade_angle135_offset7,                        
                        ClusterProminence_angle135_offset7,                   
                        HaralickCorrelation_angle135_offset7,                 
                        HaraEntroy,                                           
                        AngularSecondMoment,                                  
                        contrast,                                             
                        HaraVariance,                                         
                        sumAverage,                                           
                        sumVariance,                                          
                        sumEntropy,                                           
                        differenceVariance,                                   
                        differenceEntropy,                                    
                        inverseDifferenceMoment,                              
                        ShortRunEmphasis_AllDirection_offset1,                
                        LongRunEmphasis_AllDirection_offset1,                 
                        GreyLevelNonuniformity_AllDirection_offset1,          
                        RunLengthNonuniformity_AllDirection_offset1,          
                        LowGreyLevelRunEmphasis_AllDirection_offset1,         
                        HighGreyLevelRunEmphasis_AllDirection_offset1,        
                        ShortRunLowGreyLevelEmphasis_AllDirection_offset1,    
                        ShortRunHighGreyLevelEmphasis_AllDirection_offset1,   
                        LongRunLowGreyLevelEmphasis_AllDirection_offset1,     
                        LongRunHighGreyLevelEmphasis_AllDirection_offset1,    
                        ShortRunEmphasis_AllDirection_offset1_SD,             
                        LongRunEmphasis_AllDirection_offset1_SD,              
                        GreyLevelNonuniformity_AllDirection_offset1_SD,       
                        RunLengthNonuniformity_AllDirection_offset1_SD,       
                        LowGreyLevelRunEmphasis_AllDirection_offset1_SD,      
                        HighGreyLevelRunEmphasis_AllDirection_offset1_SD,     
                        ShortRunLowGreyLevelEmphasis_AllDirection_offset1_SD, 
                        ShortRunHighGreyLevelEmphasis_AllDirection_offset1_SD,
                        LongRunLowGreyLevelEmphasis_AllDirection_offset1_SD,  
                        LongRunHighGreyLevelEmphasis_AllDirection_offset1_SD, 
                        ShortRunEmphasis_angle0_offset1,                      
                        LongRunEmphasis_angle0_offset1,                       
                        GreyLevelNonuniformity_angle0_offset1,                
                        RunLengthNonuniformity_angle0_offset1,                
                        LowGreyLevelRunEmphasis_angle0_offset1,               
                        HighGreyLevelRunEmphasis_angle0_offset1,              
                        ShortRunLowGreyLevelEmphasis_angle0_offset1,          
                        ShortRunHighGreyLevelEmphasis_angle0_offset1,         
                        LongRunLowGreyLevelEmphasis_angle0_offset1,           
                        LongRunHighGreyLevelEmphasis_angle0_offset1,          
                        ShortRunEmphasis_angle45_offset1,                     
                        LongRunEmphasis_angle45_offset1,                      
                        GreyLevelNonuniformity_angle45_offset1,               
                        RunLengthNonuniformity_angle45_offset1,               
                        LowGreyLevelRunEmphasis_angle45_offset1,              
                        HighGreyLevelRunEmphasis_angle45_offset1,             
                        ShortRunLowGreyLevelEmphasis_angle45_offset1,         
                        ShortRunHighGreyLevelEmphasis_angle45_offset1,        
                        LongRunLowGreyLevelEmphasis_angle45_offset1,          
                        LongRunHighGreyLevelEmphasis_angle45_offset1,         
                        ShortRunEmphasis_angle90_offset1,                     
                        LongRunEmphasis_angle90_offset1,                      
                        GreyLevelNonuniformity_angle90_offset1,               
                        RunLengthNonuniformity_angle90_offset1,               
                        LowGreyLevelRunEmphasis_angle90_offset1,              
                        HighGreyLevelRunEmphasis_angle90_offset1,             
                        ShortRunLowGreyLevelEmphasis_angle90_offset1,         
                        ShortRunHighGreyLevelEmphasis_angle90_offset1,        
                        LongRunLowGreyLevelEmphasis_angle90_offset1,          
                        LongRunHighGreyLevelEmphasis_angle90_offset1,         
                        ShortRunEmphasis_angle135_offset1,                    
                        LongRunEmphasis_angle135_offset1,                     
                        GreyLevelNonuniformity_angle135_offset1,              
                        RunLengthNonuniformity_angle135_offset1,              
                        LowGreyLevelRunEmphasis_angle135_offset1,             
                        HighGreyLevelRunEmphasis_angle135_offset1,            
                        ShortRunLowGreyLevelEmphasis_angle135_offset1,        
                        ShortRunHighGreyLevelEmphasis_angle135_offset1,       
                        LongRunLowGreyLevelEmphasis_angle135_offset1,         
                        LongRunHighGreyLevelEmphasis_angle135_offset1,        
                        ShortRunEmphasis_AllDirection_offset4,                
                        LongRunEmphasis_AllDirection_offset4,                 
                        GreyLevelNonuniformity_AllDirection_offset4,          
                        RunLengthNonuniformity_AllDirection_offset4,          
                        LowGreyLevelRunEmphasis_AllDirection_offset4,         
                        HighGreyLevelRunEmphasis_AllDirection_offset4,        
                        ShortRunLowGreyLevelEmphasis_AllDirection_offset4,    
                        ShortRunHighGreyLevelEmphasis_AllDirection_offset4,   
                        LongRunLowGreyLevelEmphasis_AllDirection_offset4,     
                        LongRunHighGreyLevelEmphasis_AllDirection_offset4,    
                        ShortRunEmphasis_AllDirection_offset4_SD,             
                        LongRunEmphasis_AllDirection_offset4_SD,              
                        GreyLevelNonuniformity_AllDirection_offset4_SD,       
                        RunLengthNonuniformity_AllDirection_offset4_SD,       
                        LowGreyLevelRunEmphasis_AllDirection_offset4_SD,      
                        HighGreyLevelRunEmphasis_AllDirection_offset4_SD,     
                        ShortRunLowGreyLevelEmphasis_AllDirection_offset4_SD, 
                        ShortRunHighGreyLevelEmphasis_AllDirection_offset4_SD,
                        LongRunLowGreyLevelEmphasis_AllDirection_offset4_SD,  
                        LongRunHighGreyLevelEmphasis_AllDirection_offset4_SD, 
                        ShortRunEmphasis_angle0_offset4,                      
                        LongRunEmphasis_angle0_offset4,                       
                        GreyLevelNonuniformity_angle0_offset4,                
                        RunLengthNonuniformity_angle0_offset4,                
                        LowGreyLevelRunEmphasis_angle0_offset4,               
                        HighGreyLevelRunEmphasis_angle0_offset4,              
                        ShortRunLowGreyLevelEmphasis_angle0_offset4,          
                        ShortRunHighGreyLevelEmphasis_angle0_offset4,         
                        LongRunLowGreyLevelEmphasis_angle0_offset4,           
                        LongRunHighGreyLevelEmphasis_angle0_offset4,          
                        ShortRunEmphasis_angle45_offset4,                     
                        LongRunEmphasis_angle45_offset4,                      
                        GreyLevelNonuniformity_angle45_offset4,               
                        RunLengthNonuniformity_angle45_offset4,               
                        LowGreyLevelRunEmphasis_angle45_offset4,              
                        HighGreyLevelRunEmphasis_angle45_offset4,             
                        ShortRunLowGreyLevelEmphasis_angle45_offset4,         
                        ShortRunHighGreyLevelEmphasis_angle45_offset4,        
                        LongRunLowGreyLevelEmphasis_angle45_offset4,          
                        LongRunHighGreyLevelEmphasis_angle45_offset4,         
                        ShortRunEmphasis_angle90_offset4,                     
                        LongRunEmphasis_angle90_offset4,                      
                        GreyLevelNonuniformity_angle90_offset4,               
                        RunLengthNonuniformity_angle90_offset4,               
                        LowGreyLevelRunEmphasis_angle90_offset4,              
                        HighGreyLevelRunEmphasis_angle90_offset4,             
                        ShortRunLowGreyLevelEmphasis_angle90_offset4,         
                        ShortRunHighGreyLevelEmphasis_angle90_offset4,        
                        LongRunLowGreyLevelEmphasis_angle90_offset4,          
                        LongRunHighGreyLevelEmphasis_angle90_offset4,         
                        ShortRunEmphasis_angle135_offset4,                    
                        LongRunEmphasis_angle135_offset4,                     
                        GreyLevelNonuniformity_angle135_offset4,              
                        RunLengthNonuniformity_angle135_offset4,              
                        LowGreyLevelRunEmphasis_angle135_offset4,             
                        HighGreyLevelRunEmphasis_angle135_offset4,            
                        ShortRunLowGreyLevelEmphasis_angle135_offset4,        
                        ShortRunHighGreyLevelEmphasis_angle135_offset4,       
                        LongRunLowGreyLevelEmphasis_angle135_offset4,         
                        LongRunHighGreyLevelEmphasis_angle135_offset4,        
                        ShortRunEmphasis_AllDirection_offset7,                
                        LongRunEmphasis_AllDirection_offset7,                 
                        GreyLevelNonuniformity_AllDirection_offset7,          
                        RunLengthNonuniformity_AllDirection_offset7,          
                        LowGreyLevelRunEmphasis_AllDirection_offset7,         
                        HighGreyLevelRunEmphasis_AllDirection_offset7,        
                        ShortRunLowGreyLevelEmphasis_AllDirection_offset7,    
                        ShortRunHighGreyLevelEmphasis_AllDirection_offset7,   
                        LongRunLowGreyLevelEmphasis_AllDirection_offset7,     
                        LongRunHighGreyLevelEmphasis_AllDirection_offset7,    
                        ShortRunEmphasis_AllDirection_offset7_SD,             
                        LongRunEmphasis_AllDirection_offset7_SD,              
                        GreyLevelNonuniformity_AllDirection_offset7_SD,       
                        RunLengthNonuniformity_AllDirection_offset7_SD,       
                        LowGreyLevelRunEmphasis_AllDirection_offset7_SD,      
                        HighGreyLevelRunEmphasis_AllDirection_offset7_SD,     
                        ShortRunLowGreyLevelEmphasis_AllDirection_offset7_SD, 
                        ShortRunHighGreyLevelEmphasis_AllDirection_offset7_SD,
                        LongRunLowGreyLevelEmphasis_AllDirection_offset7_SD,  
                        LongRunHighGreyLevelEmphasis_AllDirection_offset7_SD, 
                        ShortRunEmphasis_angle0_offset7,                      
                        LongRunEmphasis_angle0_offset7,                       
                        GreyLevelNonuniformity_angle0_offset7,                
                        RunLengthNonuniformity_angle0_offset7,                
                        LowGreyLevelRunEmphasis_angle0_offset7,               
                        HighGreyLevelRunEmphasis_angle0_offset7,              
                        ShortRunLowGreyLevelEmphasis_angle0_offset7,          
                        ShortRunHighGreyLevelEmphasis_angle0_offset7,         
                        LongRunLowGreyLevelEmphasis_angle0_offset7,           
                        LongRunHighGreyLevelEmphasis_angle0_offset7,          
                        ShortRunEmphasis_angle45_offset7,                     
                        LongRunEmphasis_angle45_offset7,                      
                        GreyLevelNonuniformity_angle45_offset7,               
                        RunLengthNonuniformity_angle45_offset7,               
                        LowGreyLevelRunEmphasis_angle45_offset7,              
                        HighGreyLevelRunEmphasis_angle45_offset7,             
                        ShortRunLowGreyLevelEmphasis_angle45_offset7,         
                        ShortRunHighGreyLevelEmphasis_angle45_offset7,        
                        LongRunLowGreyLevelEmphasis_angle45_offset7,          
                        LongRunHighGreyLevelEmphasis_angle45_offset7,         
                        ShortRunEmphasis_angle90_offset7,                     
                        LongRunEmphasis_angle90_offset7,                      
                        GreyLevelNonuniformity_angle90_offset7,               
                        RunLengthNonuniformity_angle90_offset7,               
                        LowGreyLevelRunEmphasis_angle90_offset7,              
                        HighGreyLevelRunEmphasis_angle90_offset7,             
                        ShortRunLowGreyLevelEmphasis_angle90_offset7,         
                        ShortRunHighGreyLevelEmphasis_angle90_offset7,        
                        LongRunLowGreyLevelEmphasis_angle90_offset7,          
                        LongRunHighGreyLevelEmphasis_angle90_offset7,         
                        ShortRunEmphasis_angle135_offset7,                    
                        LongRunEmphasis_angle135_offset7,                     
                        GreyLevelNonuniformity_angle135_offset7,              
                        RunLengthNonuniformity_angle135_offset7,              
                        LowGreyLevelRunEmphasis_angle135_offset7,             
                        HighGreyLevelRunEmphasis_angle135_offset7,            
                        ShortRunLowGreyLevelEmphasis_angle135_offset7,        
                        ShortRunHighGreyLevelEmphasis_angle135_offset7,       
                        LongRunLowGreyLevelEmphasis_angle135_offset7,         
                        LongRunHighGreyLevelEmphasis_angle135_offset7,        
                        Sphericity,                                           
                        SurfaceArea,                                          
                        VolumeCC,                                             
                        VolumeMM,                                             
                        SurfaceVolumeRatio,                                   
                        Maximum3DDiameter,                                    
                        Compactness1,                                         
                        Compactness2,                                         
                        SphericalDisproportion) ~ T, data = dat)
sig = c()  
a = summary.aov(res.man)
for (i in 1:length(a)){
  
  sig[i] = a[[i]][1,5]
}

writeOut(rbind(c(1:i),sig),"tmp")

res = rbind(r,p,m,s,sig)


choseIndex=r2[,sig<0.05&p<0.05&r>0.75]

ncol(choseIndex)
nrow(choseIndex)
length(dat$T)

save.image("xinhuaIma")


install.packages("e1071")


library(e1071)

testData=cbind(as.factor(dat$T),choseIndex)

head(testData)

tmp = colnames(testData)
tmp[1]="sign"
colnames(testData) = tmp
remove(tmp)

model = svm(sign~. ,data=testData)


x <- subset(testData, select = -sign)
y <- testData$sign
pred <- predict(model,x)
cbind(as.character(pred),as.character(y))

head(dat)

head(testData)


data("iris")
attach(iris)

model = svm(Species~.,data=iris)



dat = read.table("xinhua/R2.txt",header = F, sep = "\t")

dat_1 = t(dat[,-1])

dat_1$

design = cbind(Intercept=1,Group=dat[,1])

dat_1[dat[,1]==1,]=scale(dat_1[dat[,1]==1,])
dat_1[dat[,1]==0,]=scale(dat_1[dat[,1]==0,])


fit <- lmFit(dat_1,design)

tfit <- treat(fit,lfc=log2(1.1))
res = topTreat(tfit,number = 20600,coef=2)
head(res)
writeOut(res,"xinhua/diff.txt")


boxplot(genesExp~group,
        data = data.frame(genesExp = 
                            as.numeric(dat[,385]),
                          group = dat[,1]))
(mean(dat[dat[,1]==1,114])+mean(dat[dat[,1]==0,114]))/2
mean(dat[dat[,1],113])



######norm plot
logi<-lrm(LNM~Size+pT+Location+Ulcer+NerveInvasion+VascularInvasion, x=TRUE, y=TRUE)
nomo<-nomogram(logi, fun=plogis, fun.at=c(.001, .01, .05, seq(.1,.9,by=.1), .95, .99, .999),
               lp=F, funlabel="LNM")
plot(nomo)


