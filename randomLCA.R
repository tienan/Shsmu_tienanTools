install.packages("randomLCA")
library("randomLCA")


randomLCA(myocardial[, 1:4], freq = myocardial$freq,
 nclass = 1)
myocardial.lca2 <- randomLCA(myocardial[, 1:4], freq = myocardial$freq,
                                 nclass = 1)
myocardial.lca2 <- randomLCA(myocardial[, 1:4], freq = myocardial$freq,
                             nclass = 1)

summary(myocardial.lca2)

head(dentistry)

head(symptoms)
nrow(symptoms)
myocardial
?symptoms
?myocardial



myocardial.lca1 <- randomLCA(myocardial[1:2, 1:4], freq = myocardial$freq[c(1:2)],nclass = 1)
myocardial.lca2 <- randomLCA(myocardial[, 1:4], freq = myocardial$freq,nclass = 2)
summary(myocardial.lca2)
outcomeProbs(myocardial.lca2)

plot(myocardial.lca2, type = "b", pch = 1:2, xlab = "Test",
      ylab = "Outcome Probability",
      scales = list(x = list(at = 1:4, labels = names(myocardial)[1:4])),
      key = list(corner = c(0.05, .95), border = TRUE, cex = 1.2,
                   text = list(c("Class 1", "Class 2")),
                   col = trellis.par.get()$superpose.symbol$col[1:2],
                   points = list(pch = 1:2)))

outcomeProbs(myocardial.lca2)

outcomeProbs(myocardial.lca2)

print(postClassProbs(myocardial.lca2), row.names = T)

?postClassProbs


?randomLCA
