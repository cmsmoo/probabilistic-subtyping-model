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
