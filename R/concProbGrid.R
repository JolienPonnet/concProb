concProbGrid <- function (inputDT, quantSplits, letsTime = TRUE){
  checkDT(inputDT, c('observed', 'predicted'))
  checkLogicVec(list(letsTime))
  checkNumOrIntVec(list(quantSplits))
  
  if (letsTime)
    ptm <- proc.time()
  
  observed <- predicted <- NULL
  
  splits <- sort(unique(quantile(inputDT[, predicted], quantSplits))) # Define the split points
  groups <- groupConstruct(inputDT, splits) # Count the number of observations in each interval
  counts <- countConc(groups$zeroGroup, groups$oneGroup) # Determine the elements of the sum in (3.3), without division by the total nb of elements in each group
  
  concProb <- sum(counts$nHigher)/(sum(counts$nHigher) + sum(counts$nLower))
  
  if (letsTime) {
    time <- proc.time() - ptm
    return(list(concProbGlobal = concProb, time = time)) 
  } else {
    return(list(concProbGlobal = concProb)) 
  }
}

groupConstruct <- function(inputDT, splits){
  observed <- predicted <- NULL
  
  X0 <- inputDT[observed == 0, ]
  X0 <- X0[order(predicted), ]
  X1 <- inputDT[observed == 1, ]
  X1 <- X1[order(predicted), ]
  
  zeroGroup = rep(NA,(length(splits) + 1))
  oneGroup = rep(NA,(length(splits) + 1))
  
  vecX0 <- X0$predicted
  vecX1 <- X1$predicted
  
  stepSize <- 100
  
  for(iSplit in 1:(length(splits))){
    seed <- max(round(iSplit/length(splits)*length(vecX0)),1)
    ind0 <- findFirstPos(vecX0, splits[iSplit], seed, stepSize)
    seed <- max(round(iSplit/length(splits)*length(vecX1)),1)
    ind1 <- findFirstPos(vecX1, splits[iSplit], seed, stepSize)
    
    zeroGroup[iSplit] <- ind0 - 1
    oneGroup[iSplit] <- ind1 - 1
  }
  zeroGroup <- c(zeroGroup[1], diff(zeroGroup[-length(zeroGroup)]), length(vecX0) - zeroGroup[length(splits)])
  oneGroup <- c(oneGroup[1], diff(oneGroup[-length(oneGroup)]), length(vecX1) - oneGroup[length(splits)])
  list(zeroGroup = zeroGroup, oneGroup = oneGroup)
}

countConc <- function(zeroGroup, oneGroup){
  nSplits <- length(zeroGroup)
  nLower <- zeroGroup * c(0,cumsum(oneGroup)[-nSplits])
  nHigher <- zeroGroup * c(rev(cumsum(rev(oneGroup)))[-1],0)
  list(nHigher = nHigher, nLower = nLower)
}

# Returns the index of the first number >= targetVal in vec
findFirstPos <- function(vec, targetVal, seed, stepSize){
  vec = as.vector(vec)
  indTemp <- seed
  if(vec[indTemp] <= targetVal){
    while((vec[indTemp] <= targetVal) & (indTemp < length(vec)) ){
      indTemp <- min(indTemp + stepSize, length(vec))
    }
    
    indTemp <- max(indTemp - stepSize + 1, 1)
    while((vec[indTemp] < targetVal) & (indTemp < length(vec)) ){
      indTemp <- indTemp + 1
    }
  } else {
    while((vec[indTemp] > targetVal) & (indTemp > 1)){
      indTemp <- max(indTemp - stepSize, 1)
    }
    
    indTemp <- min(indTemp + stepSize - 1, length(vec))
    while((vec[indTemp] >= targetVal) & indTemp >= 1){
      indTemp <- indTemp - 1
      if(indTemp == 0) break
    }
    indTemp = indTemp + 1
  }
  return(indTemp)
}