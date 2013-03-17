
GLYCANREF <<- list("Glcnac"=203.19, "Man"=162.14, "Fuc"=146.14,
                   "NGNA"=275, "NANA"=291.09)


buildCrossDiffTable = function(inputVector){
  resultTable = matrix(NA, nrow=length(inputVector), ncol=length(inputVector))
  for(i in 1:nrow(resultTable)){
    resultTable[i, ] = inputVector - inputVector[i]
  }
  return(resultTable)
}



