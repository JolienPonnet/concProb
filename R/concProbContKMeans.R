concProbContKMeans <- function(inputDT, nu = 0, nClus, letsTime = TRUE){
  
  checkDT(inputDT, c('observed', 'predicted'))
  checkLogicVec(list(letsTime))
  checkNumOrIntVec(list(nu, nClus))
  checkRanges(list(nu, nClus), list(c('>=', 0), c('>', 0)))
  
  observed <- predicted <- "." <- NULL
  
  if(letsTime)
    ptm <- proc.time()
  
  dataTemp <- inputDT
  
  clusFit <- suppressWarnings(kmeans(dataTemp[, .(predicted, observed)], nClus, iter.max = 20))
  
  clusCents <- clusFit$centers
  clusSize <- clusFit$size
  
  obsCol <- which(colnames(clusCents) == 'observed')
  predCol <- which(colnames(clusCents) == 'predicted')
  
  nConc <- rep(0, nClus)
  nDisc <- rep(0, nClus)
  
  for(iClus in 1:(nClus-1)){
    for(jClus in (iClus+1):nClus){
      if(abs(clusCents[iClus, obsCol] - clusCents[jClus, obsCol]) > nu){
        if( (sign(clusCents[iClus, obsCol] - clusCents[jClus, obsCol]) == sign(clusCents[iClus, predCol] - clusCents[jClus, predCol])) & ((sign(clusCents[iClus, obsCol] - clusCents[jClus, obsCol]) != 0)) ){
          nConc[iClus] <- nConc[iClus] + as.numeric(clusSize[iClus])*as.numeric(clusSize[jClus])
          nConc[jClus] <- nConc[jClus] + as.numeric(clusSize[iClus])*as.numeric(clusSize[jClus])
        } else if ((clusCents[iClus, obsCol] != clusCents[jClus, obsCol]) &( clusCents[iClus, predCol] != clusCents[jClus, predCol])){ #remove ties
          nDisc[iClus] <- nDisc[iClus] + as.numeric(clusSize[iClus])*as.numeric(clusSize[jClus])
          nDisc[jClus] <- nDisc[jClus] + as.numeric(clusSize[iClus])*as.numeric(clusSize[jClus])
        }
      } 
    }
  }
  
  concProb <- sum(nConc)/(sum(nConc) + sum(nDisc))
  
  if (letsTime) {
    time <- proc.time() - ptm 
    return(list(concProbGlobal = concProb, time = time)) 
  } else {
    return(list(concProbGlobal = concProb))
  }
}