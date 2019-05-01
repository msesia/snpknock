test_that("Markov chain group knockoffs are consistent with Markov chain individual knockoffs", {
  #Parameters
  p = 100
  K = 10
  n = 1000
  group.size = 1

  # Generate data
  pInit = rep(1/K,K)
  Q = array(stats::runif((p-1)*K*K),c(p-1,K,K))
  for(j in 1:(p-1)) {
    Q[j,,] = Q[j,,] + diag(rep(1,K))
    Q[j,,] = Q[j,,] / rowSums(Q[j,,])
  }
  X = SNPknock.models.sampleDMC(pInit, Q, n=n)

  # Define groups
  groups = rep(seq(1,p),each=group.size)[1:p]
  # Generate individual knockoffs
  Xk.i = SNPknock.knockoffDMC_single(X, pInit, Q, display_progress=F)
  # Generate group knockoffs
  Xk.g = SNPknock.knockoffDMC(X, pInit, Q, groups=groups, display_progress=F)

  # Verify that the knockoffs are identical
  expect_equal(Xk.i, Xk.g)
})

test_that("HMM group knockoffs are consistent with HMM individual knockoffs", {
  #Parameters
  p = 100
  K = 10
  M = 3
  n = 1000
  group.size = 1

  # Generate data
  pInit = rep(1/K,K)
  Q = array(stats::runif((p-1)*K*K),c(p-1,K,K))
  for(j in 1:(p-1)) {
    Q[j,,] = Q[j,,] + diag(rep(1,K))
    Q[j,,] = Q[j,,] / rowSums(Q[j,,])
  }
  pEmit = array(stats::runif(p*M*K),c(p,M,K))
  for(j in 1:p) { pEmit[j,,] = sweep(pEmit[j,,],2,colSums(pEmit[j,,]),'/') }
  X = SNPknock.models.sampleHMM(pInit, Q, pEmit, n=n)

  # Define groups
  groups = rep(seq(1,p),each=group.size)[1:p]
  # Generate individual knockoffs
  Xk.i = SNPknock.knockoffHMM_single(X, pInit, Q, pEmit, display_progress=F)
  # Generate group knockoffs
  Xk.g = SNPknock.knockoffHMM(X, pInit, Q, pEmit, groups=groups, display_progress=F)

  # Verify that the knockoffs are identical
  expect_equal(Xk.i, Xk.g)
})

test_that("Group knockoffs for haplotypes are consistent with group HMM knockoffs", {
  # Parameters
  p = 200
  n = 1000
  group.size = 11

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

  # Define groups
  groups = sample(p/group.size,p,replace=T)

  # Sample X from this HMM
  X = SNPknock.models.sampleHMM(hmm.large$pInit, hmm.large$Q, hmm.large$pEmit, n=n)
  # Generate group knockoffs with general algorithm
  Xk <- SNPknock.knockoffHMM(X, hmm.large$pInit, hmm.large$Q, hmm.large$pEmit, groups=groups, display_progress=F)
  # Generate group knockoffs with specialized algorithm
  Xk.hap <- SNPknock.knockoffHaplotypes(X, hmm.small$r, hmm.small$alpha, hmm.small$theta, groups=groups, display_progress=F)

  # Verify that the knockoffs are identical
  expect_equal(Xk, Xk.hap)
})

test_that("Group knockoffs for genotypes are consistent with group HMM knockoffs", {
  # Parameters
  p = 200
  n = 1000
  group.size = 11

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

  # Define groups
  groups = sample(p/group.size,p,replace=T)

  # Sample X from this HMM
  X = SNPknock.models.sampleHMM(hmm.large$pInit, hmm.large$Q, hmm.large$pEmit, n=n)
  # Generate group knockoffs with general algorithm
  Xk <- SNPknock.knockoffHMM(X, hmm.large$pInit, hmm.large$Q, hmm.large$pEmit, groups=groups, display_progress=F)
  # Generate group knockoffs with specialized algorithm
  Xk.gen <- SNPknock.knockoffGenotypes(X, hmm.small$r, hmm.small$alpha, hmm.small$theta, groups=groups, display_progress=F)

  # Verify that the knockoffs are identical
  expect_equal(Xk, Xk.gen)
})
