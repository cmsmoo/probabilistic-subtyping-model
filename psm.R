msg <- function(...) message(sprintf(...))

logsumexp <- function(x) {
  ## Compute the sum of values that are represented in log space
  ## without triggering underflow.
  ##
  ## x : Vector of log-space values.
  ##
  ## returns : The sum represented in log-space.
  ##
  m <- max(x)
  m + log(sum(exp(x - m)))
}

soft_assign <- function(d, n) {
  ## Softly assign n objects to each cell of a dimension d grid (sum
  ## of the grid for each object sums to 1).
  ##
  ## d : Dimension of the outcome grid.
  ## n : Number of objects.
  ##
  ## returns : An array with dimension c(d, n)
  ##
  dn <- c(d, n)
  Z  <- array(runif(prod(dn)), dn)
  Zs <- apply(Z, length(dn), sum)
  sweep(Z, length(dn), Zs, "/")
}

lrcovxx <- function(X, K) {
  ## Compute the design matrix variance-covariance statistic for
  ## linear regression.
  ##
  ## X : The design (feature) matrix.
  ## K : A Gram matrix used to weight observations.
  ##
  ## RETURNS
  ##
  ## An ncol(X) x ncol(X) variance-covariance matrix.
  ##
  t(X) %*% solve(K, X)
}

lrcovxy <- function(X, y, K) {
  ## Compute the design-response covariance statistic for linear
  ## regression.
  ##
  ## X : The design (feature) matrix.
  ## y : The response vector.
  ## K : A Gram matrix used to weight observations.
  ##
  ## returns : An ncol(X)-long vector.
  ##
  c(t(X) %*% solve(K, y))
}

fitcoef <- function(bb, ss1, ss2, w, p) {
  ## Compute the coefficients for this basis given a pair of collected
  ## sufficient statistics and their weights.
  ##
  ## bb  : A bspline_basis object.
  ## ss1 : Var.-cov. matrix of design (X' * X).
  ## ss2 : Cov. of design and response (X' * y).
  ## w   : Weights for the sufficient statistics.
  ## p   : The penalty weight.
  ##
  ## returns : A nbases(bb) long vector of coefficients.
  ##
  nb       <- nbases(bb)
  eta1     <- matrix(0, nb, nb)
  eta2     <- numeric(nb)
  eta1[, ] <- p * penalty(bb) + 1e-2 * diag(nb)

  for (i in seq(along = ss1)) {
    eta1 <- eta1 + w[i] * ss1[[i]]
    eta2 <- eta2 + w[i] * ss2[[i]]
  }

  c(solve(eta1, eta2))
}

logp_series <- function(y, B, K, w) {
  ## Compute the log probability of a series given its basis
  ## expansion, covariance matrix, and a set of basis weights.
  ##
  ## y : Series values.
  ## B : Basis expansion.
  ## K : Covariance matrix.
  ## w : Basis weights.
  ##
  ## returns : A singe log probability.
  ##
  yhat <- c(B %*% w)
  mvtnorm::dmvnorm(y, yhat, K, log = TRUE)
}

run_em <- function(t, y, opt) {
  N <- length(t)
  G <- opt$num_subtypes

  bb <- opt$bb
  ## bb <- bspline_basis(opt$knots)

  kfn <- kern_additive(
    kern_diag(opt$noise_v)
  , kern_const(opt$const_v)
  , kern_ou(opt$ou_v, opt$ou_l)
  )

  ## Compute design matrices.
  B <- lapply(t, design, bb = bb)

  ## Compute covariance matrices.
  K <- lapply(t, kfn)

  ## Pre-compute linear regression suff. stats.
  covxx <- mapply(lrcovxx, B, K, SIMPLIFY = FALSE)
  covxy <- mapply(lrcovxy, B, y, K, SIMPLIFY = FALSE)

  ## Initialize model parameters.
  marg <- numeric(G)
  W <- matrix(0, nbases(bb), G)

  ## Initialize expectations.
  Z <- soft_assign(G, N)

  ## Initialize looping.
  maxiter <- opt$maxiter
  tol     <- opt$tol
  logl    <- -Inf

  for (itr in 1:maxiter) {

    ## Estimate subtype marginal.
    marg[] <- rowSums(Z) / sum(Z)

    ## Estimate basis weights.
    for (j in 1:G)
      W[, j] <- fitcoef(bb, covxx, covxy, Z[j, ], 1e-2)

    ## Compute log probability of data.
    lp_marg <- replicate(N, log(marg))
    lp_markers <- matrix(0, G, N)
    for (j in 1:G)
      lp_markers[j, ] <- mapply(
        logp_series
      , y, B, K
      , MoreArgs = list(w = W[, j])
      )

    lp_joint <- lp_marg + lp_markers
    lp <- apply(lp_joint, 2, logsumexp)

    ## Stop if change is small enough.
    old_logl  <- logl
    logl      <- sum(lp)
    delta     <- (logl - old_logl) / abs(logl)
    converged <- delta < tol

    msg("LL=%.02f Convergence=%.06f", logl, delta)
    if (converged) break

    ## Recompute expectations.
    Z <- exp(sweep(lp_joint, 2, lp, "-"))
  }

  ## Rearrange the subtypes in order of increasing severity.
  B0 <- design(bb$lower, bb)
  y0 <- B0 %*% W
  severity <- order(y0, decreasing = TRUE)
  marg <- marg[severity]
  W <- W[, severity]
  Z <- Z[severity, ]

  list(
    nsubtypes = G
  , marg = marg
  , W = W
  , bb = bb
  , kfn = kfn
  , logl = logl
  , Z = Z
  )
}
