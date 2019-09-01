## Welcome to PoSI site

This page is to demonstrate simulations comparing the assumption-lean PoSI with various post-selection inference methods. Details will be updated soon on [Valid Post-selection Inference in Assumption-lean Linear Regression](https://arxiv.org/abs/1806.04119). Our package will be also up soon.

### Sample generation scheme

The following code generates samples using setup specified in `opt`. 

| Parameter | Description 												| 
| --------- | ------------------------------------------------------ 	| 
| xmat		| Sample setup                                             	|
|           | `a`: orthogonal design; `b`: exchangeable design; `c`: worst-case design                           |
| nrow		| Sample size												|
| ncol		| Number of covariates 										|
| maxk 		| Maximum model size 										|
| seed_beta | Random seed for X 										|
| seed_eps  | Random seed for error										|
| conf_level| Confidence level 											| 
| nboot 	| Bootstrap sample size										| 
| method 	| Model selection methods 									|
|           | `fs`: forward selection; `lar`: LARS; `bic`: model with the smallest BIC | 


```r
library(pracma)
library(matrixStats)

source("utilities.R")

# Sample setup
opt <- NULL
opt$xmat <- "a"			# sample setup
opt$nrow <- 200		# sample size
opt$ncol <- 15			# number of covariates
opt$maxk <- 5			# max model size
opt$seed_beta <- 123	# random seed for X, beta
opt$seed_eps <- 100		# random seed for error
opt$conf_level <- .95	# confidence level
opt$nboot <- 200		# bootstrap sample size
opt$method <- "fs"		# model selection method

# Generate sample
data <- Generate(opt)
xx <- data$x
yy <- data$y

```

### PoSI vs [Berk et al.](https://projecteuclid.org/euclid.aos/1369836961)

The following chunk computes the proposed `PoSI`, projected PoSI `PoSIBox` and `Berk` et al. PoSI. 

```r
if (!require("tmax")) install.packages("tmax_1.0.tar.gz", repos=NULL, dependencies=T)
require("tmax")

# selected model
M <- c(1,2)

# PoSI
posi_fit <- fixedx_posi(xx, yy, alpha = 1-opt$conf_level, Nboot = opt$nboot)
posi_ret <- posi(posi_fit, M)

# Berk
## it might take a while
berk_fit <- maxt_posi(xx, yy, maxk = opt$maxk, sandwich = FALSE, 
											alpha = 1-opt$conf_level, Nboot = opt$nboot)
berk_ret <- posi(berk_fit, M, sigma = 1)		# assume sigma to be known here

```


### PoSI vs [selectiveInference](https://projecteuclid.org/euclid.aos/1460381681)

The following chunk computes confidence regions using PoSI and selective inference for the first `opt$maxk` steps of forward stepwise `opt$method="fs"` or LARS `opt$method="lar"`. 

```r
library(selectiveInference)
# selectiveInference
if(opt$method == "fs") {
    fit <- fs(xx, yy, maxsteps = opt$maxk, intercept = F)
    fit_si_inf <- fsInf(fit, type = "active")
} else if(opt$method == "lar") {
    fit <- lar(xx, yy, maxsteps = opt$maxk+1, intercept = F)
    fit_si_inf <- larInf(fit, k = opt$maxk, type = "active")
}
si_box <- fit_si_inf$ci

# PoSI
posi_fit <- fixedx_posi(xx, yy, alpha = 1-opt$conf_level, Nboot = opt$nboot)
posi_ret <- posi(posi_fit, fit_si_inf$vars)

```


### PoSI vs [splitSample](https://arxiv.org/abs/1611.05401)

The following chunk computes confidence regions using PoSI and split sample method. 
For forward stepwise `opt$method="fs"` and LARS `opt$method="lar"`, we compute the confidence regions for the model at `opt$maxk` step. For `opt$method="bic"`, we compute the confidence regions for the model with smallest BIC after `opt$maxk` steps of forward stepwise selection. We use Bonferroni correction to achieve simultaneous coverage for the split sample method.

```r
library(selectiveInference)
library(leaps)

# split sample
sample_idx <- sample(opt$nrow, opt$nrow/2)
sample_ci_idx <- setdiff(1:opt$nrow, sample_idx)

# use half of the sample to select model 
if(opt$method == "fs") {
  fit <- fs(xx[sample_idx,], yy[sample_idx], maxsteps = opt$maxk, intercept = F)
  fit_si_inf <- fsInf(fit, k = opt$maxk)
  selected_vars <- fit_si_inf$vars
} else if(opt$method == "lar") {
  fit <- lar(xx[sample_idx,], yy[sample_idx], maxsteps = opt$maxk+1, intercept = F)
  fit_si_inf <- larInf(fit, k = opt$maxk)
  selected_vars <- fit_si_inf$vars
} else if(opt$method == "bic") {
  fit <- regsubsets(xx[sample_idx,], yy[sample_idx], nvmax = opt$maxk, method = "forward",  intercept = F)
  fit.s <- summary(fit)
  selected_vars <- which(fit.s$which[which.min(fit.s$bic),])
}
# use the other half to produce inference
fit <- lm(yy[sample_ci_idx] ~ xx[sample_ci_idx, selected_vars] - 1) 
split_box <- confint(fit, level = 1-0.05/length(selected_vars)) 	# Bonferroni correction


# PoSI
if(opt$method == "fs") {
  fit <- fs(xx, yy, maxsteps = opt$maxk, intercept = F)
  fit_si_inf <- fsInf(fit, k = opt$maxk)
  selected_vars <- fit_si_inf$vars
} else if(opt$method == "lar") {
  fit <- lar(xx, yy, maxsteps = opt$maxk+1, intercept = F)
  fit_si_inf <- larInf(fit, k = opt$maxk)
  selected_vars <- fit_si_inf$vars
} else if(opt$method == "bic") {
  fit <- regsubsets(xx, yy, nvmax = opt$maxk, method = "forward",  intercept = F)
  fit.s <- summary(fit)
  selected_vars <- which(fit.s$which[which.min(fit.s$bic),])
}

posi_fit <- fixedx_posi(xx, yy, alpha = 1-opt$conf_level, Nboot = opt$nboot)
posi_ret <- posi(posi_fit, selected_vars)

```
