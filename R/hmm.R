#'Forward-backward sampling for hidden Markov models
#'
#' This function samples the latent Markov chain of a hidden Markov model given the observed values.
#'
#' @param X an integer matrix of size n-by-p containing the original variables.
#' @param pInit an array of length K, containing the marginal distribution of the states for the first variable.
#' @param Q an array of size (p-1,K,K), containing a list of p-1 transition matrices between the K states of the Markov chain.
#' @param pEmit an array of size (p,M,K), containing the emission probabilities for each of the M possible emission states,
#' from each of the K hidden states and the p variables.
#' @param seed an integer random seed (default: 123).
#' @param cluster a computing cluster object created by \link[parallel]{makeCluster} (default: NULL).
#' @param display_progress whether to show progress bar (default: FALSE).
#' @return An integer matrix of size n-by-p containing the knockoff variables.
#'
#' @family hmm
#'
#' @details
#' Each element of the matrix X should be an integer value between 0 and M-1.
#' The transition matrices contained in Q are defined with the same convention as in \link{knockoffDMC}.
#' The emission propability matrices contained in pEmit are defined such that \eqn{P[X_{j}=k|H_{j}=l]=\mathrm{pEmit}[j,k,l]},
#' where \eqn{H_j} is the latent variable associated to \eqn{X_j}.
#'
#' @references
#'   \insertRef{sesia2019}{SNPknock}
#'
#' @examples
#' # Generate data
#' p=10; K=5; M=3;
#' pInit = rep(1/K,K)
#' Q = array(stats::runif((p-1)*K*K),c(p-1,K,K))
#' for(j in 1:(p-1)) { Q[j,,] = Q[j,,] / rowSums(Q[j,,]) }
#' pEmit = array(stats::runif(p*M*K),c(p,M,K))
#' for(j in 1:p) { pEmit[j,,] = pEmit[j,,] / rowSums(pEmit[j,,]) }
#' X = sampleHMM(pInit, Q, pEmit, n=20)
#' # Forward-backward sampling
#' H = fwHMM(X, pInit, Q, pEmit)
#'
#' @export
fwHMM <- function(X, pInit, Q, pEmit, seed=123, cluster=NULL, display_progress=FALSE) {
  # Dummy groups
  groups = seq(1, ncol(X))
    
  # Verify dimensions are compatible
  stopifnot(dim(X)[2]==dim(Q)[1]+1)
  stopifnot(length(pInit)==dim(Q)[2])
  stopifnot(dim(pEmit)[3]==dim(Q)[2])
  stopifnot(dim(pEmit)[1]==dim(Q)[1]+1)
  stopifnot(dim(Q)[2]==dim(Q)[3])

  # Verify contents are compatible
  stopifnot(is.integer(X))
  stopifnot(is.numeric(pInit))
  stopifnot(is.numeric(Q))
  stopifnot(is.numeric(pEmit))
  stopifnot(is.numeric(seed))
  stopifnot(seed==floor(seed))
  stopifnot( min(X)>=0 )
  stopifnot( max(X)<dim(pEmit)[2])
  seed = as.integer(seed)
  stopifnot(is.integer(seed))
  stopifnot(is.logical(display_progress))

  # Extract dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  K = length(pInit)
  M = dim(pEmit)[2]
    
  if( (!is.null(cluster)) & (length(cluster)>1) ) {
    if(requireNamespace("doParallel", quietly = TRUE))
    {
      # Count number of workers in the cluster
      ncores = length(cluster)

      # Assign rows to workers
      splits <- cut(1:nrow(X),breaks=ncores,labels=FALSE)

      # Sample knockoffs in parallel
      Xk = do.call(rbind, parallel::parLapply(cluster, 1:ncores, function(i) {
        n.split = sum(splits==i)
        display_progress = (i==1)
        fwHMM_wrapper(X[splits==i,], pInit, Q, pEmit, n.split, ncol(X), K, M, seed+(i-1),
                      groups, display_progress)
      }))
    } else {
      warning("To enable multithreading, please install the `doParallel` package ")
      # Sample knockoffs sequentially
      Xk = fwHMM_wrapper(X, pInit, Q, pEmit, n, p, K, seed, groups, display_progress)
    }
  } else {
    # Sample knockoffs sequentially
    Xk = fwHMM_wrapper(X, pInit, Q, pEmit, n, p, K, M, seed, groups, display_progress)
  }
  storage.mode(Xk) = "integer"
  return(Xk)
}

#'Forward-backward sampling for hidden Markov models for genotypes
#'
#' This function samples the latent Markov chain of a hidden Markov model given the observed values,
#' given an HMM for unphased genotypes.
#'
#' @param X a {0,1,2} matrix of size n-by-p containing the original variables.
#' @param r a vector of length p containing the "r" parameters estimated by fastPHASE.
#' @param alpha a matrix of size p-by-K containing the "alpha" parameters estimated by fastPHASE.
#' @param theta a matrix of size p-by-K containing the "theta" parameters estimated by fastPHASE.
#' @param seed an integer random seed (default: 123).
#' @param cluster a computing cluster object created by \link[parallel]{makeCluster} (default: NULL).
#' @param display_progress whether to show progress bar (default: FALSE).
#' @return A {0,1,2} matrix of size n-by-p containing the knockoff variables.
#'
#' @family hmm
#'
#' @references
#'   \insertRef{sesia2019multi}{SNPknock}
#'
#' @examples
#' # Problem size
#' p = 10
#' n = 100
#' # Load HMM to generate data
#' r_file = system.file("extdata", "haplotypes_rhat.txt", package = "SNPknock")
#' alpha_file = system.file("extdata", "haplotypes_alphahat.txt", package = "SNPknock")
#' theta_file = system.file("extdata", "haplotypes_thetahat.txt", package = "SNPknock")
#' char_file = system.file("extdata", "haplotypes_origchars", package = "SNPknock")
#' hmm.data = loadHMM(r_file, alpha_file, theta_file, char_file, compact=FALSE, phased=FALSE)
#' hmm.data$Q = hmm.data$Q[1:(p-1),,]
#' hmm.data$pEmit = hmm.data$pEmit[1:p,,]
#' # Sample X from this HMM
#' X = sampleHMM(hmm.data$pInit, hmm.data$Q, hmm.data$pEmit, n=n)
#' # Load HMM to generate knockoffs
#' hmm = loadHMM(r_file, alpha_file, theta_file, char_file)
#' hmm$r = hmm$r[1:p]
#' hmm$alpha = hmm$alpha[1:p,]
#' hmm$theta = hmm$theta[1:p,]
#' # Forward-backward sampling
#' H = fwGenotypes(X, hmm$r, hmm$alpha, hmm$theta)

#' @export
fwGenotypes <- function(X, r, alpha, theta, seed=123, cluster=NULL, display_progress=FALSE) {
  # Dummy groups
  groups = seq(1, ncol(X))

  # Verify dimensions are compatible
  stopifnot(dim(X)[2]==length(r))
  stopifnot(dim(X)[2]==dim(alpha)[1])
  stopifnot(dim(X)[2]==dim(theta)[1])
  stopifnot(dim(alpha)[2]==dim(theta)[2])
  stopifnot(dim(X)[2]==length(groups))

  # Verify contents are compatible
  stopifnot(is.integer(X))
  stopifnot(is.numeric(r))
  stopifnot(is.numeric(alpha))
  stopifnot(is.numeric(theta))
  stopifnot(is.numeric(seed))
  stopifnot(seed==floor(seed))
  stopifnot( min(X)>=0 )
  stopifnot( max(X)<=2 )
  seed = as.integer(seed)
  stopifnot(is.integer(seed))
  stopifnot(is.logical(display_progress))
  stopifnot(is.integer(groups))

  # Extract dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  K = dim(alpha)[2]

  # Split non-contigous groups
  groups.cont = rep(1,p)
  group.count = 1
  for(j in 2:p) {
    if(groups[j]!=groups[j-1]) {
      group.count = group.count + 1
    }
    groups.cont[j] = group.count
  }

  if( (!is.null(cluster)) & (length(cluster)>1) ) {
    if(requireNamespace("doParallel", quietly = TRUE))
    {
      # Count number of workers in the cluster
      ncores = length(cluster)

      # Assign rows to workers
      splits <- cut(1:nrow(X),breaks=ncores,labels=FALSE)

      # Sample knockoffs in parallel
      Xk = do.call(rbind, parallel::parLapply(cluster, 1:ncores, function(i) {
        n.split = sum(splits==i)
        display_progress = (i==1)*display_progress
        fwGenotypes_wrapper(X[splits==i,], r, alpha, theta, groups.cont-1, n.split, ncol(X),
                               seed+(i-1), display_progress)
      }))
    } else {
      warning("To enable multithreading, please install the doParallel package ")
      # Sample knockoffs sequentially
      Xk = fwGenotypes_wrapper(X, r, alpha, theta, groups.cont-1, n, p, seed, display_progress)
    }
  } else {
    # Sample knockoffs sequentially
    Xk = fwGenotypes_wrapper(X, r, alpha, theta, groups.cont-1, n, p, seed, display_progress)
  }
  storage.mode(Xk) = "integer"
  return(Xk)
}

#'Forward-backward sampling for hidden Markov models for haplotypes
#'
#' This function samples the latent Markov chain of a hidden Markov model given the observed values,
#' given an HMM for phased haplotypes.
#'
#' @param X a binary matrix of size n-by-p containing the original variables.
#' @param r a vector of length p containing the "r" parameters estimated by fastPHASE.
#' @param alpha a matrix of size p-by-K containing the "alpha" parameters estimated by fastPHASE.
#' @param theta a matrix of size p-by-K containing the "theta" parameters estimated by fastPHASE.
#' @param seed an integer random seed (default: 123).
#' @param cluster a computing cluster object created by \link[parallel]{makeCluster} (default: NULL).
#' @param display_progress whether to show progress bar (default: FALSE).
#' @return A binary matrix of size n-by-p containing the knockoff variables.
#'
#' @family hmm
#'
#' @references
#'   \insertRef{sesia2019multi}{SNPknock}
#'
#' @examples
#' # Problem size
#' p = 10
#' n = 100
#' # Load HMM to generate data
#' r_file = system.file("extdata", "haplotypes_rhat.txt", package = "SNPknock")
#' alpha_file = system.file("extdata", "haplotypes_alphahat.txt", package = "SNPknock")
#' theta_file = system.file("extdata", "haplotypes_thetahat.txt", package = "SNPknock")
#' char_file = system.file("extdata", "haplotypes_origchars", package = "SNPknock")
#' hmm.data = loadHMM(r_file, alpha_file, theta_file, char_file, compact=FALSE, phased=TRUE)
#' hmm.data$Q = hmm.data$Q[1:(p-1),,]
#' hmm.data$pEmit = hmm.data$pEmit[1:p,,]
#' # Sample X from this HMM
#' X = sampleHMM(hmm.data$pInit, hmm.data$Q, hmm.data$pEmit, n=n)
#' # Load HMM to generate knockoffs
#' hmm = loadHMM(r_file, alpha_file, theta_file, char_file)
#' hmm$r = hmm$r[1:p]
#' hmm$alpha = hmm$alpha[1:p,]
#' hmm$theta = hmm$theta[1:p,]
#' # Forward-backward sampling
#' H = fwHaplotypes(X, hmm$r, hmm$alpha, hmm$theta)
#'
#' @export
fwHaplotypes <- function(X, r, alpha, theta, seed=123, cluster=NULL, display_progress=FALSE) {
  # Dummy groups  
  groups = seq(1, ncol(X))

  # Verify dimensions are compatible
  stopifnot(dim(X)[2]==length(r))
  stopifnot(dim(X)[2]==dim(alpha)[1])
  stopifnot(dim(X)[2]==dim(theta)[1])
  stopifnot(dim(alpha)[2]==dim(theta)[2])
  stopifnot(dim(X)[2]==length(groups))

  # Verify contents are compatible
  stopifnot(is.integer(X))
  stopifnot(is.numeric(r))
  stopifnot(is.numeric(alpha))
  stopifnot(is.numeric(theta))
  stopifnot(is.numeric(seed))
  stopifnot(seed==floor(seed))
  stopifnot( min(X)>=0 )
  stopifnot( max(X)<=1 )
  seed = as.integer(seed)
  stopifnot(is.integer(seed))
  stopifnot(is.logical(display_progress))
  stopifnot(is.integer(groups))

  # Extract dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  K = dim(alpha)[2]

  # Split non-contigous groups
  groups.cont = rep(1,p)
  group.count = 1
  for(j in 2:p) {
    if(groups[j]!=groups[j-1]) {
      group.count = group.count + 1
    }
    groups.cont[j] = group.count
  }

  if( (!is.null(cluster)) & (length(cluster)>1) ) {
    if(requireNamespace("doParallel", quietly = TRUE))
    {
      # Count number of workers in the cluster
      ncores = length(cluster)

      # Assign rows to workers
      splits <- cut(1:nrow(X),breaks=ncores,labels=FALSE)

      # Sample knockoffs in parallel
      Xk = do.call(rbind, parallel::parLapply(cluster, 1:ncores, function(i) {
        n.split = sum(splits==i)
        display_progress = (i==1)*display_progress
        fwHaplotypes_wrapper(X[splits==i,], r, alpha, theta, groups.cont-1, n.split, p,
                                seed+(i-1), display_progress)
      }))
    } else {
      warning("To enable multithreading, please install the doParallel package ")
      # Sample knockoffs sequentially
      Xk = fwHaplotypes_wrapper(X, r, alpha, theta, groups.cont-1, n, p, seed, display_progress)
    }
  } else {
    # Sample knockoffs sequentially
    Xk = fwHaplotypes_wrapper(X, r, alpha, theta, groups.cont-1, n, p, seed, display_progress)
  }
  storage.mode(Xk) = "integer"
  return(Xk)
}
