

preprocessData <- function(ranges, sampleData){
  startPos <- min(start(ranges))
  endPos <- max(end(ranges))
  nrSamples <- ncol(sampleData)
  ans <- matrix(NA, nrow=nrSamples, ncol=endPos-startPos+1)
  rownames(ans) <- colnames(sampleData)
  colnames(ans) <- startPos:endPos
  for(i in 1:length(ranges)){
    ans[ , as.character(ranges[[i]])] <- sampleData[i, ]
  }
  return(ans)
}


