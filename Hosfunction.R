scoreHos = function(datReso){
  ncol = ncol(datReso)
  score = data.frame()
  datTmp = data.frame()
  datTmp = datReso
  
  k = 1 # score count
  for (i in 1:ncol){
    tmp = table(datTmp[,i])
    tmpScore=c()
    for (j in 1:length(tmp)){
      tmpScore[j] = (sum(tmp[1:j])/sum(tmp))
    }
    gap = c(length(tmp):1) #
    gap = 1/gap 
    #
    #give the mapping score 
    for (j in 1:length(datTmp[,i])){
      score[j,k] = min(tmpScore[names(tmp)==datTmp[j,i]],gap[names(tmp)==datTmp[j,i]])
    }
    k = k+1
    # tmpScore = tmpScore/max(tmpScore)
  }

    mean = apply(score, 2, mean)
    sd = apply(score,2,sd)

    if(length(mean)>1){
      if ((max(mean) - min(mean)) != 0){
#        weight = (max(mean)-mean)/(max(mean) - min(mean))
        weight = sd
      }else{
        weight=1
      }
      
      #apply(score, 2, sd)
      result=c()
      for (i in 1:nrow(score)){
        result[i] = sum(score[i,] * weight)/sum(weight)
      }
    }else{
      result=c()
      for (i in 1:nrow(score)){
        result[i] = score[i,] 
      }
  }

#  for (i in 1:nrow(score)){
#    cat(result[i])
#   cat("\n")
#  }
  return(result)
}


printRow = function(dat){
  for (i in 1:nrow(dat)){
    cat(dat[i])
    cat("\n")
  }
}


score = function(datReso){
  ncol = ncol(datReso)
  score = data.frame()
  datTmp = data.frame()
  datTmp = datReso
  
  k = 1 # score count
  for (i in 1:ncol){
    tmp = table(datTmp[,i])
    tmpScore=c()
    for (j in 1:length(tmp)){
      tmpScore[j] = (sum(tmp[1:j])/sum(tmp))
    }
    #give the mapping score 
    # tmpScore = tmpScore/max(tmpScore)
    gap = c(length(tmp):1) #
    gap = 1/gap 
    #
    #give the mapping score 
    for (j in 1:length(datTmp[,i])){
      score[j,k] = min(tmpScore[names(tmp)==datTmp[j,i]],gap[names(tmp)==datTmp[j,i]])
    }
    k = k+1
  }
  return(score)
}
