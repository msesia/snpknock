test_that("HMM for haplotypes are loaded correctly", {
  # Load HMM from RData
  hmm.file <- system.file("extdata", "haplotypes_hmm.RData", package = "SNPknock")
  load(hmm.file)
  # Load fastPHASE parameter estimates from text files
  r.file <- system.file("extdata", "haplotypes_rhat.txt", package = "SNPknock")
  alpha.file <- system.file("extdata", "haplotypes_alphahat.txt", package = "SNPknock")
  theta.file <- system.file("extdata", "haplotypes_thetahat.txt", package = "SNPknock")
  char.file <- system.file("extdata", "haplotypes_origchars", package = "SNPknock")
  hmm.compact <- SNPknock.fp.loadFit(r.file, alpha.file, theta.file, char.file)
  hmm.large <-  SNPknock.fp.loadFit_hmm(r.file, alpha.file, theta.file, char.file, phased=T)
  hmm <- c(hmm.compact,hmm.large)
  # Verify that the parameters are identical
  expect_equal(hmm, hmm.hap)
})

test_that("HMM for genotypes are loaded correctly", {
  # Load HMM from RData
  hmm.file <- system.file("extdata", "genotypes_hmm.RData", package = "SNPknock")
  load(hmm.file)
  # Load fastPHASE parameter estimates from text files
  r.file <- system.file("extdata", "genotypes_rhat.txt", package = "SNPknock")
  alpha.file <- system.file("extdata", "genotypes_alphahat.txt", package = "SNPknock")
  theta.file <- system.file("extdata", "genotypes_thetahat.txt", package = "SNPknock")
  char.file <- system.file("extdata", "genotypes_origchars", package = "SNPknock")
  hmm.compact <- SNPknock.fp.loadFit(r.file, alpha.file, theta.file, char.file)
  hmm.large <-  SNPknock.fp.loadFit_hmm(r.file, alpha.file, theta.file, char.file, phased=F)
  hmm <- c(hmm.compact,hmm.large)
  # Verify that the parameters are identical
  expect_equal(hmm, hmm.gen)
})

test_that("Haplotypes are correctly converted into INP format", {
  # Load haplotypes from RData and convert into INP
  hap.file <- system.file("extdata", "haplotypes.RData", package = "SNPknock")
  load(hap.file)
  inp.file <- SNPknock.fp.writeX(H, phased=T)
  # Stored INP file
  inp.file.stored <- system.file("extdata", "haplotypes.inp", package = "SNPknock")
  expect_equal(readLines(inp.file), readLines(inp.file.stored))
})

test_that("Genotypes are correctly converted into INP format", {
  # Load genotypes from RData and convert into INP
  gen.file <- system.file("extdata", "genotypes.RData", package = "SNPknock")
  load(gen.file)
  inp.file <- SNPknock.fp.writeX(X, phased=F)
  # Stored INP file
  inp.file.stored <- system.file("extdata", "genotypes.inp", package = "SNPknock")
  expect_equal(readLines(inp.file), readLines(inp.file.stored))
})
