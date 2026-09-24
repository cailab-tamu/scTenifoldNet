test_that("checkMemory estimates peak memory and the largest gene set that fits", {
  mem <- checkMemory(nGenes = 5000, nNet = 10, nConditions = 1, warn = FALSE)
  expect_named(mem, c("required", "available", "total", "maxGenes"))

  # Calibrated estimate for 5,000 genes and 10 networks (measured: 8.61 GB)
  expect_equal(mem$required / 1e9, 8.61, tolerance = 0.05)

  # Memory grows with the square of the number of genes and with nNet
  expect_gt(checkMemory(2000, warn = FALSE)$required / checkMemory(1000, warn = FALSE)$required, 2.5)
  expect_gt(checkMemory(3000, nNet = 20, warn = FALSE)$required,
            checkMemory(3000, nNet = 5, warn = FALSE)$required)

  # Two conditions need more memory than one
  expect_gt(checkMemory(3000, nConditions = 2, warn = FALSE)$required,
            checkMemory(3000, nConditions = 1, warn = FALSE)$required)

  # maxGenes is the largest gene count whose estimate fits in the available memory
  if (!is.na(mem$available)) {
    expect_lte(checkMemory(mem$maxGenes, warn = FALSE)$required, mem$available)
    expect_gt(checkMemory(mem$maxGenes + 1, warn = FALSE)$required, mem$available)
  }

  expect_error(checkMemory(100, nConditions = 3))
})

test_that("checkMemory warns only when the estimate does not fit", {
  skip_if(is.na(checkMemory(10, warn = FALSE)$available))
  expect_warning(checkMemory(1e6), "at most")
  expect_silent(checkMemory(1e6, warn = FALSE))
  expect_silent(checkMemory(10))
})
