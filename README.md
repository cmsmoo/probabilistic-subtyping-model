# Probabilistic Subtyping Model (PSM)

This is an R implementation of the probabilistic subtyping model (PSM). The code implements a procedure that clusters "unbalanced" longitudinal trajectories (i.e. the number of observations and the timing of observations are different across trajectories). The clusters will be invariant to certain additive transformations (see the associated [paper](http://pschulam.com/papers/schulam+wigley+saria_aaai_2015.pdf)). By default, these are vertical shifts and transient spikes, which are both substantively motivated in the paper. If you use this code, please cite us:

```
@inproceedings{schulam2015clustering,
  title={Clustering Longitudinal Clinical Marker Trajectories from Electronic Health Data: Applications to Phenotyping and Endotype Discovery.},
  author={Schulam, Peter and Wigley, Fredrick and Saria, Suchi},
  booktitle={AAAI},
  pages={2956--2964},
  year={2015}
}
```

# Usage

The following code shows how to use the model (this is copied from the file `psm_demo.R`). The data is expected to be in a "long format" csv. In a long format dataset, there are three columns and each row records an observation made on one patient at one point in time. The three columns are: the patient ID (some unique key identifying each individual), the timestamp (some real number representing time since an alignment event), and the observed value itself. See [this blog post](http://www.r-bloggers.com/managing-longitudinal-data-conversion-between-the-wide-and-the-long/) for more information about wide vs. long format longitudinal datasets.

```{r}
source("bsplines.R")
source("kernels.R")
source("psm.R")

library(dplyr)
library(readr)

set.seed(1)

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
```

# Bugs

Please report any bugs to Peter Schulam (author and maintaner) or submit them as issues in this repository.
