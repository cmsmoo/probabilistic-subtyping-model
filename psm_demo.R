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


source("bsplines.R")
source("kernels.R")
source("psm.R")

library(dplyr)
library(readr)

set.seed(20150610)

bb <- bspline_basis(0, 23, 3, 2)

opt <- list(
  num_subtypes = 9     ## How many subpopulations to use.
, bb           = bb    ## The B-spline basis object.
, noise_v      = 1.0   ## The random noise variance.
, const_v      = 16.0  ## The magnitude of the long-term kernel component.
, ou_v         = 64.0  ## The magnitude of the short-term kernel component.
, ou_l         = 4.0   ## The length-scale of the short-term kernel component.
, maxiter      = 100   ## The max. number of EM iterations to run.
, tol          = 1e-4  ## The relative improvement threshold to use as a stopping point.
  )

benchmark <- read_csv("your_long_format_data.csv")
dataset <- benchmark %>% group_by(patient_id) %>% do(t = .$timestamp, y = .$observed_value)

tt <- dataset$t
yt <- dataset$y

nobs <- vapply(tt, length, numeric(1))
minobs <- vapply(tt, min, numeric(1))
valid <- minobs <= 2 & nobs > 2

fit <- run_em(tt[valid], yt[valid], opt)

saveRDS(fit, "model.rds")
