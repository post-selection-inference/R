# generate possible models given all variables
gen_submodel_indices_k <- function(num_vars, total_vars, force_in=NULL) {
  vars_to_choose <- setdiff(1:total_vars, force_in)
  num_vars_to_choose <- num_vars - length(force_in)
  ind_mat <- t(combn(vars_to_choose, num_vars_to_choose))
  if(is.null(force_in)) {return(ind_mat)}
  return(t(apply(ind_mat, 1, function(i) sort(c(i, force_in)))))
}

# Generate sample based on setting specified in opt
# beta0 defaults to 0
Generate <- function(opt, beta0 = NULL){
  n <- opt$nrow; p <- opt$ncol;
  setting <- opt$xmat
  if(is.null(beta0)) {
    set.seed(opt$seed_beta)
    beta0 <- runif(opt$ncol, -1, 1)*0 
  }
  set.seed(opt$seed_eps)
  eps <- rnorm(opt$nrow, 0, 1)
  
  set.seed(opt$seed_beta)
  Xpre <- matrix(rnorm(n*p), nrow=n, ncol = p)
  X <- pracma::gramSchmidt(Xpre)$Q
  ## X satisfies t(X)%*%X = I_p.
  if(setting == "a"){
    ## Setting a means alpha = 0; orthogonal design 	
    x <- X
  }
  if(setting == "b"){
    ## Setting b means alpha = -1/(p+2); t(x)%*%x = I_p + alpha*1_p1_p^t.
    lambda_2 <- 1
    lambda_1 <- (1 - p/(p+2)) ## 32 should have been p+2.
    SigmaHalf <- sqrt(lambda_1)*matrix(1/p,ncol = p, nrow = p) + 
      sqrt(lambda_2)*(diag(rep(1,p)) - matrix(1/p,nrow = p,ncol = p))
    x <- X%*%SigmaHalf
  }
  if(setting == "c"){
    ## Setting c means x^tx = (e_1,...,e_{p-1},X_p(c)).	
    SigmaHalf <- cbind(rbind(diag(rep(1,p-1)),rep(0,p-1)), 
                       c(rep(1/sqrt(2*(p-1)), p-1),1/sqrt(2)))
    x <- X%*%SigmaHalf
  }
  x <- x*sqrt(n)
  y <- x%*%beta0 + eps
  return(list(y = y, x = x, xtx = t(x)%*%x, opt = opt))
}

# PoSI constant for fixed X
fixedx_posi <- function(x, y, alpha = 0.05, Nboot = 1000){
  # x = full regressor matrix
  # y = response vector
  # alpha = alpha level
  # Nboot = bootstrap sample
  x <- as.matrix(x)
  y <- as.vector(y)
  namex <- colnames(x)
  n <- nrow(x)
  p <- ncol(x)
  if(n != length(y))
    stop("Number of rows of x should match length of y!")
  xx <- t(x)%*%x/n   # Gram matrix of regressors
  Y <- matrix(y,nrow = n)
  xy <- t(x)%*%Y/n   # 'Gamma vector'
  e.mat <- matrix(rnorm(n*Nboot), nrow=n) # Multiplier
  r1.vec <- colMaxs(abs(t(x)%*%sweep(e.mat, 1, y, "*")-outer(c(xy), colSums(e.mat), "*")))/n
  c1.alpha <- unname(quantile(r1.vec, prob = 1 - alpha))
  
  ret <- list(x = x, y = y, C1Alpha = c1.alpha, Nboot = Nboot)
  class(ret) <- "PoSI"
  
  return(ret)
}

# Berk et al. PoSI
maxt_posi <- function(x, y, maxk = ncol(x), sandwich = TRUE, 
                      alpha = 0.05, Nboot = 200) {
  # x = full regressor matrix 
  # y = response vector
  # maxk = max model size
  # sandwich = indicator of using sandwich estimator
  # alpha = alpha level
  # Nboot = bootstrap sample
  n <- nrow(x); p <- ncol(x)
  maxt_0 <- rep(-Inf, Nboot)
  
  for(k in 1:maxk) {
    varlist <- gen_submodel_indices_k(k, p)
    set.seed(123)
    EE <- matrix(rnorm(n*Nboot), nrow = n, ncol = Nboot)
    tmax_ret <- tmax::max_t_mul_boot_by_k(x, y, EE, sandwich = sandwich, opt$nboot, varlist)
    maxt_0 <- pmax(maxt_0, colMax(tmax_ret$max_t))
  }
  
  k0 <- quantile(maxt_0, 1 - alpha)
  ret <- list(x = x, y = y, k = k0, Nboot = Nboot, sandwich = sandwich)
  class(ret) <- "PoSI_Berk"
  
  return(ret)
}

# generic function for PoSI and PoSI_Berk
posi <- function(posi_obj, M, sigma = NULL) {
  # posi_obj = PoSI object
  # M = model as integer vector into columns of x
  # sigma 
  UseMethod("posi", posi_obj)
}

# PoSI box and PoSI parallelpiped volume 
posi.PoSI_Berk <- function(posi_obj, M, sigma = NULL){
  # posi_obj = PoSI object
  # M = model as integer vector into columns of x
  # sigma 
  if(class(posi_obj) != "PoSI_Berk") stop("Please input PoSI object")
  x <- posi_obj$x; y <- posi_obj$y; k <- posi_obj$k
  n <- nrow(x); p <- length(M)
  mod.lm <- lm(y ~ x[,M] - 1)
  beta.M <- unname(mod.lm$coeff)
  Sigma <- t(x)%*%x/n
  Sigma.M <- Sigma[M,M,drop=FALSE] # Covariance matrix of model M
  if(!is.null(sigma)) {
    sehat <- sigma*sqrt(diag(solve(Sigma.M))/n)
  } else if(posi_obj$sandwich) {
    tmp_coef <- lm.fit(x[,M,drop=F], diag(mod.lm$residuals))$coeff
    sehat <- sqrt(diag(tmp_coef %*% t(tmp_coef)))
  } else {
    sehat <- sqrt(sum((mod.lm$residuals)^2)/(n-p-1)*diag(solve(Sigma.M))/n)
  }
  box <- cbind(beta.M - k*sehat, beta.M + k*sehat)
  vol.M <- prod(2*k*sehat)
  colnames(box) <- c("Lower", "Upper")
  rownames(box) <- colnames(x)[M]
  return(list(beta = beta.M, intervals.M = box, Volume.M = vol.M))
}

# PoSI box and PoSI parallelpiped volume 
posi.PoSI <- function(posi_obj, M){
  # posi_obj = PoSI object
  # M = model as integer vector into columns of x
  if(class(posi_obj) != "PoSI") stop("Please input PoSI object")
  x <- posi_obj$x; y <- posi_obj$y; C1Alpha <- posi_obj$C1Alpha
  Sigma <- t(x)%*%x/length(y)
  Sigma.M <- Sigma[M,M,drop=FALSE] # Covariance matrix of model M
  mod.lm <- lm(y ~ x[,M] - 1)
  beta.M <- unname(mod.lm$coeff)
  l1.Sigma.inv <- rowSums(abs(solve(Sigma.M))) * C1Alpha
  box <- cbind(beta.M - l1.Sigma.inv, beta.M + l1.Sigma.inv)
  vol.M <- (2*C1Alpha)^nrow(Sigma.M)/det(Sigma.M)
  colnames(box) <- c("Lower", "Upper")
  rownames(box) <- colnames(x)[M]
  return(list(beta = beta.M, intervals.M = box, 
              Volume.M = vol.M, Volume.box = prod(2*l1.Sigma.inv)))
}
