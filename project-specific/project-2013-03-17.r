selfDir = ""

dataFn = "dephosphorylation.txt"
charge = 13
deviationRate = 0.015
glycanRange = lapply(GLYCANREF, function(glycan){c("min"=(1 - deviationRate)*glycan, "max"=(1 + deviationRate)*glycan, "original"=glycan)})

inputTable = read.table(dataFn, header=FALSE)
massOverCharge = inputTable[ ,1]
mass = massOverCharge * charge - charge
diffTable = buildCrossDiffTable(mass)

foundFeature = list()
glycan = names(GLYCANREF)[1]
for(glycan in names(GLYCANREF)){
  foundFeature[[glycan]] = which(diffTable >= glycanRange[[glycan]]["min"] & diffTable <= glycanRange[[glycan]]["max"], arr.ind=TRUE)
}



