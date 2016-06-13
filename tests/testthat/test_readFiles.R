test_that("test_readMZ", {
  mzFns <- system.file(c("threonine/threonine_i2_e35_pH_tree.mzXML",
                         "lockmass/LockMass_test.mzXML"),
                       package = "msdata")
  
  ## Read all the mzXML data
  allIntensities <- readMZ(mzFns)
  expect_equal(length(allIntensities), 2L)
  expect_equal(sapply(allIntensities, nrow), c(353L, 977L),
               check.attributes=FALSE)
  
  starts <- c(50, 70, 80)
  ends <- c(55, 75, 85)
  rangedIntensities <- readMZ(mzFns, starts=starts, ends=ends)
  expect_equal(dim(rangedIntensities), c(2L, 3L))
}
)