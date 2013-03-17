selfDir = "/Users/gtan/Repositories/MS/scripts"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}

dataFn = "dephosphorylation.txt"
charge = 13
deviationRate = 0.015
glycanRange = lapply(GLYCANREF, function(glycan){c("min"=(1 - deviationRate)*glycan, "max"=(1 + deviationRate)*glycan, "original"=glycan)})

inputTable = read.table(dataFn, header=FALSE)
massOverCharge = inputTable[ ,1]
mass = massOverCharge * charge - charge
diffTable = buildCrossDiffTable(mass)

indexFound = list()
glycan = names(GLYCANREF)[1]
for(glycan in names(GLYCANREF)){
  indexFound[[glycan]] = which(diffTable >= glycanRange[[glycan]]["min"] & diffTable <= glycanRange[[glycan]]["max"], arr.ind=TRUE)
}

massFound = lapply(indexFound, function(indexMatrix){matrix(mass[indexMatrix], ncol=2)})
massFound = lapply(massFound, function(massMatrix){
                   massMatrix = cbind(massMatrix, abs(massMatrix[ ,1] - massMatrix[ ,2]));
                   colnames(massMatrix) = c("rowMass", "colMass", "diffMass");
                   return(massMatrix)}
                  )

glycan = names(massFound)[1]
for(glycan in names(massFound)){
  massFound[[glycan]] = cbind(massFound[[glycan]], "theoMass"=c(GLYCANREF[[glycan]]))
  massFound[[glycan]] = cbind(massFound[[glycan]], "devMass"=c(abs(massFound[[glycan]][ ,"theoMass"] - massFound[[glycan]][ ,"diffMass"])))

}

