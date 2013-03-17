selfDir = "/Users/gtan/Repositories/MS/scripts"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}
selfDir = "/Users/gtan/Repositories/NextGenerationSequencing/scripts"
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

## find the evolve path
startMass = 43846.40
indexMass = which(mass == startMass)
evolvePath = solveGlycan(indexMass, indexFound, mass)
evolvePath[ ,"Glcnac"] = evolvePath[ ,"Glcnac"] + 1
evolvePath[ ,"Man"] = evolvePath[ ,"Man"] + 4
evolvePath[ ,"realDiff"] = evolvePath[ ,"realDiff"] + 43846.40 - 42992.65
evolvePath[ ,"theoDiff"] = evolvePath[ ,"theoDiff"] + 851.75
evolvePath[ ,"dev"] = abs(evolvePath[ ,"realDiff"] - evolvePath[ ,"theoDiff"])

my.write.table(evolvePath, file="43846.40-evole.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)


startMass = 43967.82



