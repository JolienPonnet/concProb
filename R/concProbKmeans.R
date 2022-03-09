concProbKMeans <- function(inputDT, nClusMax, letsTime = TRUE){
  checkDT(inputDT, c('observed', 'predicted'))
  checkLogicVec(list(letsTime))
  checkNumOrIntVec(list(nClusMax))
  nClusMax <- round(nClusMax, digits = 0)
  checkRanges(list( nClusMax), list(c('>=', 1)))
  checkLength(list(nClusMax), c(1))
  
  if (letsTime)
    ptm <- proc.time()
  
  observed <- predicted <- NULL
  
  predInputZero <- inputDT[observed == 0, predicted]
  predInputOne <- inputDT[observed == 1, predicted]
  nClus <- min(c(nClusMax, length(predInputZero) - 1, length(predInputOne) - 1))
  zeroClus <- suppressWarnings(kmeans(predInputZero, nClus))
  oneClus <- suppressWarnings(kmeans(predInputOne, nClus))
  predClusZero <- zeroClus$cluster # The cluster to which each prediction belongs
  predClusOne <- oneClus$cluster # The cluster to which each prediction belongs
  nPredZero <- table(predClusZero) # Number of predictions that belong to each cluster
  nPredOne <- table(predClusOne) # Number of predictions that belong to each cluster
  predZero <- zeroClus$centers # The center of each cluster
  predOne <- oneClus$centers # The center of each cluster
  ch <- rep(0, nClus)
  dh <- rep(0, nClus)
  for (iZero in 1:nClus) {
    for (iOne in 1:nClus) {
      if (predZero[iZero] < predOne[iOne]) {
        ch[iZero] <- ch[iZero] + as.numeric(nPredZero[iZero]) * as.numeric(nPredOne[iOne])
      } else if (predZero[iZero] > predOne[iOne]) {
        dh[iZero] <- dh[iZero] + as.numeric(nPredZero[iZero]) * as.numeric(nPredOne[iOne])
      }
    }
  }
  concProb <- sum(ch)/(sum(ch) + sum(dh)) # This equals sum(ch)/(sum(nPredOne)*sum(nPredZero))
  
  if (letsTime) {
    time <- proc.time() - ptm
    return(list(concProbGlobal = concProb, time = time))
  } else {
    return(list(concProbGlobal = concProb))
  }
}
