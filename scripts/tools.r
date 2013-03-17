
GLYCANREF <<- list("Glcnac"=203.19, "Man"=162.14, "Fuc"=146.14,
                   "NGNA"=275, "NANA"=291.09)


buildCrossDiffTable = function(inputVector){
  resultTable = matrix(NA, nrow=length(inputVector), ncol=length(inputVector))
  for(i in 1:nrow(resultTable)){
    resultTable[i, ] = inputVector - inputVector[i]
  }
  return(resultTable)
}

solveCombination = function(startMass, currentMass, weightList){
  diffMass = currentMass - startMass
  
}

solveGlycan = function(indexMass, indexFound, mass){
  evolvePath = matrix(0, ncol=length(indexFound), dimnames=list(indexMass, names(indexFound)))
  isFound = TRUE
  i = 1
  while(isFound){
    message("i ", i)
    if(i > nrow(evolvePath)){
      break
    }
    #isFound = FALSE
    indexCurrent = rownames(evolvePath)[i]
    message("indexCurrent: ", indexCurrent)
    for(glycan in names(indexFound)){
      index = indexFound[[glycan]][ ,"row"] == indexCurrent
      if(!any(index)){
        next
      }
     # isFound = TRUE
      indexToAdd = indexFound[[glycan]][index,"col"]
      if(as.character(indexToAdd) %in% rownames(evolvePath)){
        next
      }
      message("indexToAdd: ", indexToAdd)
      rowToAdd = evolvePath[indexCurrent, ]
      rowToAdd[glycan] = rowToAdd[glycan] + 1
      evolvePath = rbind(evolvePath, rowToAdd)
      rownames(evolvePath)[nrow(evolvePath)] = indexToAdd
    }
    i = i + 1
  }
  evolvePath = cbind(evolvePath, "realDiff"=mass[as.integer(rownames(evolvePath))] - mass[indexMass])
  evolvePath = cbind(evolvePath, "theoDiff"=c(evolvePath[ ,names(indexFound)]%*%unlist(GLYCANREF)[names(indexFound)]))
  evolvePath = cbind(evolvePath, "dev"=c(abs(evolvePath[, "realDiff"] - evolvePath[ , "theoDiff"])))
  
}



