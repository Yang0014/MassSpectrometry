test_that("test_normaliseSpectrum", {
  x <- c(50,100,10,200)
  expect_equal(sum(normaliseSpectrum(x, method="sum")), 1)
  expect_equal(max(normaliseSpectrum(x, method="max")), 1)
  expect_equal(sum(normaliseSpectrum(x, method="unit")^2), 1)
}
)

test_that("test_geometricMF", {
  a <- c(1, 10, 5, 8)
  b <- c(2, 10, 5, 8)
  c <- c(1, 10, 5, 9)
  expect_equal(geometricMF(a, b), 0.9948658, tolerance=1e-7)
  expect_equal(geometricMF(a, c), 0.996804, tolerance=1e-7)
}
)