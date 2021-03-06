## Welcome to PoSI site

This page is to demonstrate simulations comparing the assumption-lean PoSI with various post-selection inference methods. Please [download](https://github.com/post-selection-inference/R/archive/master.zip) or clone this repo and install the packages if necessary. Details of the simulation setup will be updated soon on [Valid Post-selection Inference in Assumption-lean Linear Regression](https://arxiv.org/abs/1806.04119). Our package will be also up soon.

The above mentioned paper provides valid confidence regions post-variable selection in the context of linear regression. Suppose <img src="https://latex.codecogs.com/gif.latex?(X_i,&space;Y_i),&space;1\le&space;i\le&space;n" title="(X_i, Y_i), 1\le i\le n" /> represent the regression data. The OLS estimator for constructed based on <img src="https://latex.codecogs.com/gif.latex?(X_{i,M},&space;Y_i)" title="(X_{i,M}, Y_i)" /> for a subset <img src="https://latex.codecogs.com/gif.latex?M\subseteq\{1,2,\ldots,p\}" title="M\subseteq\{1,2,\ldots,p\}" /> is given by

<img src="https://latex.codecogs.com/gif.latex?\hat{\beta}_{M}&space;:=&space;\left(\frac{1}{n}\sum_{i=1}^n&space;X_{i,M}X_{i,M}^{\top}\right)^{-1}\left(\frac{1}{n}\sum_{i=1}^n&space;X_{i,M}Y_i\right)\in\mathbb{R}^{|M|}." title="\hat{\beta}_{M} := \left(\frac{1}{n}\sum_{i=1}^n X_{i,M}X_{i,M}^{\top}\right)^{-1}\left(\frac{1}{n}\sum_{i=1}^n X_{i,M}Y_i\right)\in\mathbb{R}^{|M|}." />

(This is constructed based only on the covariates with indices in <img src="https://latex.codecogs.com/gif.latex?M" title="M" />.) The target of this OLS estimator is given by

<img src="https://latex.codecogs.com/gif.latex?{\beta}_{M}&space;:=&space;\left(\frac{1}{n}\sum_{i=1}^n&space;\mathbb{E}\left[X_{i,M}X_{i,M}^{\top}\right&space;]\right)^{-1}\left(\frac{1}{n}\sum_{i=1}^n&space;\left[X_{i,M}Y_i&space;\right&space;]\right)\in\mathbb{R}^{|M|}." title="{\beta}_{M} := \left(\frac{1}{n}\sum_{i=1}^n \mathbb{E}\left[X_{i,M}X_{i,M}^{\top}\right ]\right)^{-1}\left(\frac{1}{n}\sum_{i=1}^n \left[X_{i,M}Y_i \right ]\right)\in\mathbb{R}^{|M|}." />

The reason for calling this the target of <img src="https://latex.codecogs.com/gif.latex?\hat{\beta}_M" title="\hat{\beta}_M" /> is shown in the paper. For the case of fixed covariates, the expectation is only with respect to <img src="https://latex.codecogs.com/gif.latex?Y_i" title="Y_i" />'s. The proposed confidence regions for <img src="https://latex.codecogs.com/gif.latex?\beta_{\hat{M}}" title="\beta_{\hat{M}}" /> for a randomly selected model <img src="https://latex.codecogs.com/gif.latex?\hat{M}" title="\hat{M}" /> (in case of fixed covariates) is given by

<img src="https://latex.codecogs.com/gif.latex?\hat{\mathcal{R}}_{n,\hat{M}}&space;:=&space;\left\{\theta\in\mathbb{R}^{|\hat{M}|}:\,\|\hat{\Sigma}_{n,\hat{M}}(\hat{\beta}_{n,\hat{M}}&space;-&space;\theta)\|_{\infty}&space;\le&space;C_n^{\Gamma}(\alpha)\right\}," title="\hat{\mathcal{R}}_{n,\hat{M}} := \left\{\theta\in\mathbb{R}^{|\hat{M}|}:\,\|\hat{\Sigma}_{n,\hat{M}}(\hat{\beta}_{n,\hat{M}} - \theta)\|_{\infty} \le C_n^{\Gamma}(\alpha)\right\}," />

where <img src="https://latex.codecogs.com/gif.latex?\textstyle\hat{\Sigma}_{n,\hat{M}}&space;:=&space;n^{-1}\sum_{i=1}^n&space;X_{i,\hat{M}}X_{i,\hat{M}}^{\top}," title="\textstyle\hat{\Sigma}_{n,\hat{M}} := n^{-1}\sum_{i=1}^n X_{i,\hat{M}}X_{i,\hat{M}}^{\top}," />
and <img src="https://latex.codecogs.com/gif.latex?C_n^{\Gamma}(\alpha)" title="C_n^{\Gamma}(\alpha)" /> represents the <img src="https://latex.codecogs.com/gif.latex?(1-\alpha)" title="(1-\alpha)" />-th quantile of <img src="https://latex.codecogs.com/gif.latex?\textstyle\|n^{-1}\sum_{i=1}^n&space;\{X_{i}Y_i&space;-&space;\mathbb{E}[X_iY_i]\}\|_{\infty}" title="\textstyle\|n^{-1}\sum_{i=1}^n \{X_{i}Y_i - \mathbb{E}[X_iY_i]\}\|_{\infty}" />.

### Sample generation scheme

The following code generates samples using setup specified in `opt`. 

| Parameter | Description        												| 
| --------- | ------------------------------------------------------ 	| 
| xmat		  | Sample setup                              |
|           | `a`: orthogonal design; `b`: exchangeable design; `c`: worst-case design |
| nrow		  | Sample size												        |
| ncol		  | Number of covariates 										  |
| maxk      | Maximum model size                        |
| seed_beta | Random seed for X 										    |
| seed_eps  | Random seed for error								   		|
| conf_level| Confidence level 											    | 
| nboot     | Bootstrap sample size									  	| 
| method 	  | Model selection methods 									|
|           | `fs`: forward selection; `lar`: LARS; `bic`: model with the smallest BIC | 


```r
library(pracma)
library(matrixStats)

source("utilities.R")
# This file contains the functions Generate()
# and fixedx_posi().

# Sample setup
opt <- NULL
opt$xmat <- "a"				# sample setup
opt$nrow <- 200				# sample size
opt$ncol <- 15				# number of covariates
opt$maxk <- 5					# max model size
opt$seed_beta <- 123	# random seed for X, beta
opt$seed_eps <- 100		# random seed for error
opt$conf_level <- .95	# confidence level
opt$nboot <- 200			# bootstrap sample size
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
