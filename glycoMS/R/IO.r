
readXls = function(xlsFn, sheet="Sheet1", chargeState){
  require(XLConnect)
  wb = loadWorkbook(xlsFn)
  sheetData = readWorksheet(wb, sheet=sheet)
  indexMain = which(!is.na(sheetData[ ,"No.."]))
  resultSheet = sheetData[indexMain, ]
  #indexChargeState = which(sheetData[ ,"Average.Mass"] == chargeState)
  indexChargeState = rep(NA, nrow(resultSheet))
  j = 0
  for(i in 1:nrow(sheetData)){
    if(sheetData[i, "Average.Mass"] == "Charge State"){
      j = j + 1
    }
    if(sheetData[i, "Average.Mass"] == as.character(chargeState)){
      indexChargeState[j] = i
    }
  }
  measuredAverage = sheetData[indexChargeState, "Sum.Intensity"]
  stopifnot(length(measuredAverage) == nrow(resultSheet))
  resultSheet = cbind(resultSheet, "Measured Average m/z"=as.numeric(measuredAverage))
  return(resultSheet)
}



