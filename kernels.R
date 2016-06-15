kern_diag <- function(v) {
  ## Create a diagonal covariance kernel (uncorrelated).
  ##
  ## PARAM
  ##
  ## v : The noise variance.
  ##
  ## RETURNS
  ##
  ## A function mapping from pairs of vectors to a covariance matrix.
  ##
  force(v)
  function(x1, x2) {

    if (missing(x2))
      diag(v, length(x1))

    else
      matrix(0, length(x1), length(x2))
  }
}

kern_const <- function(v) {
  ## Create a constant covariance kernel (uniform correlation).
  ##
  ## PARAM
  ##
  ## v : The magnitude of the covariance from the mean.
  ##
  ## RETURNS
  ##
  ## A function mapping from pairs of vectors to a covariance matrix.
  ##
  force(v)
  function(x1, x2) {

    if (missing(x2))
      x2 <- x1

    x1 <- c(x1)
    x2 <- c(x2)
    matrix(v, length(x1), length(x2))
  }
}

kern_linear <- function(v) {
  ## Create a linear covariance kernel.
  ##
  ## PARAM
  ##
  ## v : The magnitude of the covariance.
  ##
  ## RETURNS
  ##
  ## A function mapping pairs of vectors to a covariance matrix.
  ##
  force(v)
  function(x1, x2) {

    if (missing(x2))
      x2 <- x1

    x1 <- c(x1)
    x2 <- c(x2)
    v * (x1 %o% x2)
  }
}

kern_expquad <- function(v, l) {
  ## Create an exponentiated quadratic (also known as a squared
  ## exponential) covariance kernel.
  ##
  ## PARAM
  ##
  ## v : The magnitude of the covariance from the mean.
  ## l : The lengthscale of the serial correlation.
  ##
  ## RETURNS
  ##
  ## A function mapping from pairs of vectors to a covariance matrix.
  ##
  force(v)
  force(l)
  function(x1, x2) {

    if (missing(x2))
      x2 <- x1

    x1 <- c(x1)
    x2 <- c(x2)
    d  <- outer(x1, x2, "-") / l
    K  <- v * exp(-0.5 * d ^ 2)

    return(K)
  }
}

kern_ou <- function(v, l) {
  ## Create an Ornstein-Uhlenbeck covariance kernel.
  ##
  ## PARAM
  ##
  ## v : The magnitude of the covariance from the mean.
  ## l : The lengthscale of the serial correlation.
  ##
  ## RETURNS
  ##
  ## A function mapping from pairs of vectors to a covariance matrix.
  ##
  force(v)
  force(l)
  function(x1, x2) {

    if (missing(x2))
      x2 <- x1

    x1 <- c(x1)
    x2 <- c(x2)
    d  <- outer(x1, x2, "-") / l
    K  <- v * exp(- abs(d))

    return(K)
  }
}

kern_spline <- function(bb, lambda = 1) {
  force(bb)
  force(lambda)
  function(x1, x2) {

    if (missing(x2))
      x2 <- x1

    x1 <- c(x1)
    x2 <- c(x2)
    X1 <- design(x1, bb)
    X2 <- design(x2, bb)
    P <- penalty(bb) + lambda * diag(nbases(bb))
    X1 %*% solve(P, t(X2))
  }
}

kern_additive <- function(diag_cov, ...) {
  ## Create a kernel formed by the sumof component kernels.
  ##
  ## PARAM
  ##
  ## diag_cov : A diagonal covariance kernel.
  ## ...      : Component covariance kernels.
  ##
  ## RETURNS
  ##
  ## A function mapping from pairs of vectors to a covariance matrix.
  ##
  force(diag_cov)
  components <- list(...)
  function(x1, x2) {

    if (missing(x2)) {
      K <- diag_cov(x1)
      x2 <- x1

    } else {
      K <- matrix(0, length(x1), length(x2))
    }

    for (kfn in components)
      K <- K + kfn(x1, x2)

    return(K)
  }
}
