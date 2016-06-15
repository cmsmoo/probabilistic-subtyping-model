## The MIT License (MIT)

## Copyright (c) 2016 Peter Schulam

## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:

## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.

## -----------------------------------------------------------------------------

## CITATION: If you use this software, please cite us:

## @inproceedings{schulam2015clustering,
##   title={Clustering Longitudinal Clinical Marker Trajectories from
##          Electronic Health Data: Applications to Phenotyping and
##          Endotype Discovery.},
##   author={Schulam, Peter and Wigley, Fredrick and Saria, Suchi},
##   booktitle={AAAI},
##   pages={2956--2964},
##   year={2015}
## }


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
