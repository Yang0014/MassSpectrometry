

plotHDXHeatmap <- function(dataToPlot, interval=300){
  nrRows <- ceiling((max(as.numeric(colnames(dataToPlot))) - 
                     min(as.numeric(colnames(dataToPlot)))) / interval)
  layout(matrix(1:(nrRows+1), ncol=1, nrow=nrRows+1))
  
  ## split the data into nrRows
  newData <- list()
  for(i in 1:nrRows){
    newData[[i]] <- matrix(NA, nrow=nrow(dataToPlot), ncol=interval)
    rownames(newData[[i]]) <- rownames(dataToPlot)
    if(i == nrRows){
      newData[[i]][ , seq(from=1, to=ncol(dataToPlot) - (i-1) * interval)] <-
        dataToPlot[ , seq(from=(i-1) * interval + 1, to=ncol(dataToPlot))]
      colnames(newData[[i]]) <- seq(from=as.integer(colnames(dataToPlot)[(i-1) * interval + 1]), to=as.integer(colnames(dataToPlot)[(i-1) * interval + 1]) + interval-1)
    }else{
      newData[[i]][ , 1:interval] <- dataToPlot[ , seq(from=(i-1)*interval+1,to=i*interval)]
      colnames(newData[[i]]) <- colnames(dataToPlot[ , seq(from=(i-1)*interval+1,to=i*interval)])
    }
  }

  for(i in 1:length(newData)){

  }

}



