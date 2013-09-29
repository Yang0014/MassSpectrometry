GLYCANREF = list("Glcnac"=203.1950, "Man"=162.1424, "Fuc"=146.1430,
                 "Neu5Ac"=291.2579)

glycanRange = function(glycanref, deviationRate=0.015){
  ans = lapply(glycanref, 
               function(glycan){
                 c("min"=(1 - deviationRate)*glycan, 
                   "max"=(1 + deviationRate)*glycan, 
                   "original"=glycan)
               }
               )
  return(ans)
}

buildCrossDiffTable = function(inputVector){
  inputVector = as.numeric(inputVector)
  resultTable = matrix(NA, nrow=length(inputVector), ncol=length(inputVector))
  for(i in 1:nrow(resultTable)){
    resultTable[i, ] = inputVector - inputVector[i]
  }
  return(resultTable)
}

solveCombination = function(startMass, currentMass, weightList){
  diffMass = currentMass - startMass
  
}

evolveGlycan = function(indexMass, indexFound, mass){
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


solveGlycan = function(inputTable, glycanRange, startMass, backBoneMass, startComposition=c("Glcnac"=4, "Man"=3, "Fuc"=1, "NGNA"=0, "NANA"=0)){
  mass = inputTable[ ,2]
  diffTable = buildCrossDiffTable(mass)
  indexFound = list()
  for(glycan in names(glycanRange)){
      indexFound[[glycan]] = which(diffTable >= glycanRange[[glycan]]["min"] & diffTable <= glycanRange[[glycan]]["max"], arr.ind=TRUE)
  }
  massFound = lapply(indexFound, function(indexMatrix){matrix(mass[indexMatrix], ncol=2)})
  massFound = lapply(massFound, function(massMatrix){
                     massMatrix = cbind(massMatrix, abs(massMatrix[ ,1] - massMatrix[ ,2]));
                     colnames(massMatrix) = c("rowMass", "colMass", "diffMass");
                     return(massMatrix)}
                    )
  for(glycan in names(massFound)){
    massFound[[glycan]] = cbind(massFound[[glycan]], 
                                "theoMass"=c(GLYCANREF[[glycan]]))
    massFound[[glycan]] = cbind(massFound[[glycan]], 
                                "devMass"=c(abs(massFound[[glycan]][ ,"theoMass"] - massFound[[glycan]][ ,"diffMass"])))
  }
  indexMass = which(mass == startMass)
  stopifnot(length(indexMass) != 0)
  evolvePath = evolveGlycan(indexMass, indexFound, mass)
  evolvePath[ , names(startComposition)] = sweep(
                                                 evolvePath[ , names(startComposition)],  
                                                 MARGIN=2, startComposition, "+"
                                                 )
  evolvePath[ ,"realDiff"] = evolvePath[ ,"realDiff"] + startMass - backBoneMass
  evolvePath[ ,"theoDiff"] = evolvePath[ ,"theoDiff"] + 
                             sum(startComposition * unlist(GLYCANREF))
  evolvePath[ ,"dev"] = abs(evolvePath[ ,"realDiff"] - evolvePath[ ,"theoDiff"])
  evolvePath = cbind("mass"=mass[as.integer(rownames(evolvePath))], evolvePath)
  #my.write.table(evolvePath, file="43846.40-evole.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
  evolvePath = cbind(evolvePath, 
                     inputTable[as.integer(rownames(evolvePath)), c("Sum.Intensity", "Relative.Abundance", "Measured Average m/z")]
                     )
  evolvePath = evolvePath[order(evolvePath[["mass"]]), ]
  return(evolvePath)
}

