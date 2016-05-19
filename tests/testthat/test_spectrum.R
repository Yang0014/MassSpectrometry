test_that("test_normaliseSpectrum", {
  x <- c(50,100,10,200)
  expect_equal(sum(normaliseSpectrum(x, method="sum")), 1)
  expect_equal(max(normaliseSpectrum(x, method="max")), 1)
  expect_equal(sum(normaliseSpectrum(x, method="unit")^2), 1)
}
)