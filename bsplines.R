bspline_basis <- function(xl, xr, ndx, bdeg) {
  ## Create a new B-spline basis object.
  ##
  ## PARAM
  ##
  ## xl   : Domain lower bound.
  ## xr   : Domain upper bound.
  ## ndx  : Number of interior pieces.
  ## bdeg : Degree of the B-spline components.
  ##
  ## RETURNS
  ##
  ## A bspline_basis object; used to compute features or penalties.
  ##
  dx    <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, dx)

  bb <- list(
    lower = xl
  , upper = xr
  , ndx   = ndx
  , bdeg  = bdeg
  , dx    = dx
  , knots = knots
  )

  class(bb) <- c("bspline_basis", class(bb))

  return(bb)
}

nbases <- function(bb)
  ## Compute the number of basis functions in a design matrix produced
  ## by this bspline_basis.
  ##
  ## PARAM
  ##
  ## bb : A bspline_basis object.
  ##
  ## RETURNS
  ##
  ## The number of bases (a scalar).
  ##
  bb$ndx + bb$bdeg + 1

design <- function(x, bb) {
  ## Compute a design (feature) matrix using the basis bb at the
  ## points in the vector x.
  ##
  ## PARAM
  ##
  ## bb : A bspline_basis object.
  ## x  : A vector of points in the domain of the basis.
  ##
  ## RETURNS
  ##
  ## A matrix with length(x) rows and nbases(bb) columns.
  ##
  B <- splines::spline.des(bb$knots, x, bb$bdeg + 1, 0 * x)$design
  X <- cbind(1, B)
  return(X)
}

penalty <- function(bb) {
  ## Compute a first-order differences penalty matrix for the basis
  ## coefficients to be included in a Tikhonov-regularized linear fit.
  ##
  ## PARAM
  ##
  ## bb: A bspline_basis object.
  ##
  ## RETURNS
  ##
  ## A nbases(bb) x nbases(bb) matrix.
  ##
  P <- diag(1e-4, nbases(bb))
  D <- diag(nbases(bb) - 1)
  D <- diff(D)
  D <- t(D) %*% D
  P[-1, -1] <- D
  P
}
