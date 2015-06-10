remove.packages("glycoMS")

rm glycoMS_0.99.0.tar.gz
R CMD BUILD glycoMS



### install  ###
R CMD INSTALL glycoMS_0.99.0.tar.gz

library(glycoMS)
#### IO.R #######
xlsFn = system.file("extdata/SR_2013-06-26_IgG1_HEK_gly_1_all.xls", package="glycoMS")
inputTable = readProteinDevolutionXls(xlsFn, sheet="Sheet1", chargeState=18)


### match Glycan ###
glycanRange = glycanRange(GLYCANREF, deviationRate=0.015)

startMass = 74585.75
backBoneMass = 73143.890625
glycoTable = solveGlycan(inputTable, glycanRange, startMass=startMass, backBoneMass=backBoneMass, startComposition=c("Glcnac"=4, "Man"=3, "Fuc"=1, "Neu5Ac"=0))

par(mfrow=c(4,1))
xlim = c(74000, 76500)
glycoTable = validateGlyco(glycoTable)
plotPseudoGaussianSpectrum(glycoTable[["mass"]], 100* glycoTable[["Sum.Intensity"]]/max(glycoTable[["Sum.Intensity"]]), labels=NULL, xlim=xlim, sd=3)

glycoTable = dissectGlycoTable(glycoTable, glycos=c("Man"=3))

res = digestGlycoTable(glycoTable, chargeState=18)
plotPseudoGaussianSpectrum(res[["mass"]], 100*res[["Sum.Intensity"]]/max(res[["Sum.Intensity"]]), labels=NULL, sd=3, xlim=xlim)

res = digestGal(glycoTable, chargeState=18)
plotPseudoGaussianSpectrum(res[["mass"]], 100*res[["Sum.Intensity"]]/max(res[["Sum.Intensity"]]), labels=NULL, sd=3, xlim=xlim)

res = digestGalSia(glycoTable, chargeState=18)
plotPseudoGaussianSpectrum(res[["mass"]], 100*res[["Sum.Intensity"]]/max(res[["Sum.Intensity"]]), labels=NULL, sd=3, xlim=xlim)




