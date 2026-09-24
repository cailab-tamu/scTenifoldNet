test_that("dRegulation does not flag genes whose distance is numerical noise", {
  set.seed(1)
  nGenes <- 50
  coordinates <- matrix(rnorm(nGenes * 3), nGenes, 3)
  genes <- paste0("g", seq_len(nGenes))

  # Identical coordinates up to floating-point noise: nothing changed
  noise <- matrix(rnorm(nGenes * 3, sd = 1e-15), nGenes, 3)
  same <- rbind(coordinates, coordinates + noise)
  rownames(same) <- c(paste0("X_", genes), paste0("Y_", genes))
  expect_warning(dr <- dRegulation(same), "numerical noise")
  expect_true(all(dr$p.value == 1))
  expect_true(all(dr$p.adj == 1))

  # A real change in a few genes is still detected, with no warning
  moved <- coordinates
  moved[1:3, ] <- moved[1:3, ] + 5
  changed <- rbind(coordinates, moved + noise)
  rownames(changed) <- c(paste0("X_", genes), paste0("Y_", genes))
  expect_silent(dr <- suppressMessages(dRegulation(changed)))
  expect_true(all(c("g1", "g2", "g3") %in% dr$gene[dr$p.adj < 0.05]))
  expect_true(all(dr$p.value[!dr$gene %in% c("g1", "g2", "g3")] == 1))
})
