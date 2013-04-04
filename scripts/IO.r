
readXls = function(xlsFn, sheet="Sheet1", chargeState){
  require(XLConnect)
  wb = loadWorkbook(xlsFn)
  sheetData = readWorksheet(wb, sheet=sheet)
  indexMain = which(!is.na(sheetData[ ,"No.."]))
  resultSheet = sheetData[indexMain, ]
  indexChargeState = which(sheetData[ ,"Average.Mass"] == chargeState)
  measuredAverage = sheetData[indexChargeState, "Sum.Intensity"]
  stopifnot(length(measuredAverage) == nrow(resultSheet))
  resultSheet = cbind(resultSheet, "Measured Average m/z"=as.numeric(measuredAverage))
  return(resultSheet)
}



