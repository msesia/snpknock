test_that("Knockoffs for haplotypes are consistent with HMM knockoffs", {
  # Parameters
  p = 100
  n = 1000
  
  # Load HMM
  r_file <- system.file("extdata", "haplotypes_rhat.txt", package = "SNPknock")
  alpha_file <- system.file("extdata", "haplotypes_alphahat.txt", package = "SNPknock")
  theta_file <- system.file("extdata", "haplotypes_thetahat.txt", package = "SNPknock")
  char_file <- system.file("extdata", "haplotypes_origchars", package = "SNPknock")
  hmm.large <- SNPknock.fp.loadFit_hmm(r_file, alpha_file, theta_file, char_file, phased=T)
  hmm.large$Q = hmm.large$Q[1:(p-1),,]
  hmm.large$pEmit = hmm.large$pEmit[1:p,,]
  hmm.small <- SNPknock.fp.loadFit(r_file, alpha_file, theta_file, char_file)
  hmm.small$r = hmm.small$r[1:p]
  hmm.small$alpha = hmm.small$alpha[1:p,]
  hmm.small$theta = hmm.small$theta[1:p,]
  
  # Sample X from this HMM
  X = SNPknock.models.sampleHMM(hmm.large$pInit, hmm.large$Q, hmm.large$pEmit, n=n)
  # Generate simple knockoffs with general algorithm
  Xk <- SNPknock.knockoffHMM(X, hmm.large$pInit, hmm.large$Q, hmm.large$pEmit, display_progress=F)
  # Generate simple knockoffs with specialized algorithm
  Xk.hap <- SNPknock.knockoffHaplotypes(X, hmm.small$r, hmm.small$alpha, hmm.small$theta, display_progress=F)
  
  # Verify that the knockoffs are identical
  expect_equal(Xk, Xk.hap)
})

test_that("Knockoffs for genotypes are consistent with HMM knockoffs", {
  # Parameters
  p = 100
  n = 1000
  
  # Load HMM
  r_file <- system.file("extdata", "haplotypes_rhat.txt", package = "SNPknock")
  alpha_file <- system.file("extdata", "haplotypes_alphahat.txt", package = "SNPknock")
  theta_file <- system.file("extdata", "haplotypes_thetahat.txt", package = "SNPknock")
  char_file <- system.file("extdata", "haplotypes_origchars", package = "SNPknock")
  hmm.large <- SNPknock.fp.loadFit_hmm(r_file, alpha_file, theta_file, char_file, phased=F)
  hmm.large$Q = hmm.large$Q[1:(p-1),,]
  hmm.large$pEmit = hmm.large$pEmit[1:p,,]
  hmm.small <- SNPknock.fp.loadFit(r_file, alpha_file, theta_file, char_file)
  hmm.small$r = hmm.small$r[1:p]
  hmm.small$alpha = hmm.small$alpha[1:p,]
  hmm.small$theta = hmm.small$theta[1:p,]
  
  # Sample X from this HMM
  X = SNPknock.models.sampleHMM(hmm.large$pInit, hmm.large$Q, hmm.large$pEmit, n=n)
  # Generate simple knockoffs with general algorithm
  Xk <- SNPknock.knockoffHMM(X, hmm.large$pInit, hmm.large$Q, hmm.large$pEmit, display_progress=F)
  # Generate simple knockoffs with specialized algorithm
  Xk.gen <- SNPknock.knockoffGenotypes(X, hmm.small$r, hmm.small$alpha, hmm.small$theta, display_progress=F)
  
  # Verify that the knockoffs are identical
  expect_equal(Xk, Xk.gen)
})