# generate possible models given all variables
gen_submodel_indices_k <- function(num_vars, total_vars, force_in=NULL, force_out=NULL) {
  vars_to_choose <- setdiff(1:total_vars, c(force_in, force_out))
  num_vars_to_choose <- num_vars - length(force_in)
  ind_mat <- t(combn(vars_to_choose, num_vars_to_choose))
  if(is.null(force_in)) {return(ind_mat)}
  return(t(apply(ind_mat, 1, function(i) sort(c(i, force_in)))))
}

gen_submodel_indices <- function(max_vars, total_vars, force_in=NULL) {
  vars_to_choose <- setdiff(1:total_vars, force_in)
  num_vars_to_choose <- max_vars - length(force_in)
  # if (max_vars>num_vars_to_choose) {max_vars<-num_vars_to_choose}
  # if (max_vars<num_vars_to_choose) {num_vars_to_choose <- max_vars}
  submodels <- vector("list")
  i <- 1
  for (num_of_vars in 1:num_vars_to_choose) {
    combs <- combn(vars_to_choose, num_of_vars)
    for (comb in 1:ncol(combs)) {
      submodels[[i]] <- combs[, comb]
      i <- i+1
    }
  }
  if(is.null(force_in)) {return(submodels)}
  lapply(1:length(submodels), function(i) sort(c(submodels[[i]], force_in)))
}

colMax <- function(data) apply(data, 2, function(x) max(x, na.rm = TRUE))

which.colMax <- function(data) apply(data, 2, function(x) which.max(x))

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
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
    SigmaHalf <- cbind(rbind(diag(rep(1,p-1)),rep(0,p-1)), c(rep(1/sqrt(2*(p-1)), p-1),1/sqrt(2)))
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
  # alpha = level of 
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

# Legacy
# return maxt, beta, se, by model size opt$k
max_t_mul_boot_sample_C <- function(xx, yy, opt, sandwich=TRUE, return_sample=TRUE, force_in=NULL) {
  nn <- nrow(xx)
  p <- ncol(xx)
  
  # get varlist
  if (!opt$k) {
    stop("Please input number of variables to select")
  } else {varlist <- gen_submodel_indices_k(opt$k, p, force_in = force_in)}
  num_submodels <- nrow(varlist)
  
  if(return_sample) {
    zz <- rep(0, p)
    all_submodels <- apply(varlist, 1, function(k) {zz[k]=1;paste0(zz, collapse = "")})
  }
  
  set.seed(opt$seed_eps)
  EE <- matrix(rnorm(nn*opt$nboot), nrow = nn, ncol = opt$nboot)
  boot_sample_C <- tmax::max_t_mul_boot_by_k(xx, yy, EE, sandwich, opt$nboot, varlist)
  
  ret <- list(X = xx, Y = yy, X_setting = opt$xmat, NumVars = opt$k, 
              beta = boot_sample_C$beta, se = boot_sample_C$se,
              BootSample = boot_sample_C$max_t, 
              AllModelNames = all_submodels, BootstrapSize = opt$nboot)
  
  ret$call <- match.call()
  class(ret) <- "MultiplierBootstrapSample"
  
  return(ret)
}

#' Return Max-T multiplier bootstrap sample of model of size $k$
#' 
#' @param xx The whole model matrix.
#' @param yy Response vector.
#' @param k Numeric. Model size.
#' @param sandwhich Logical. If TRUE, the variance is obtained by the sandwich estimator; if FALSE, the variance is obtained by the model-based variance. TODO: bootstrap based se
#' @param return_sample Logical. If TRUE, return all the bootstrap sample; if false, otherwise.
#' @param force_in Index vector. Indices of variables to force in.
#' @param Nboot Numeric. The number of bootstrap sample.
#' @param intercept Logical. If TRUE, add an intercept to the model matrix.
#' @param individual Index vector. Indices of variables that we want the individual Max-T bootstrap sample. (posi-1).
#' @return model matrix xx, response vector yy, NumVars model size, coefficient estimate beta, standard eror of \hat\beta se, 
## the vector of Max-T bootstrap sample BootSample, the vector of Max-T bootstrap sample ranking BootRank, 
## the matrix of Max-T bootstrap sample BootSample of selected individual variables BootSample1, 
## the matrix of Max-T bootstrap sample BootSample ranking of selected individual variables BootRank1, 
## matrix of all models considered AllModelNames, bootstrap sample BootstrapSize
max_t_mul_boot_k <- function(xx, yy, k, 
                             sandwich=TRUE, return_sample=TRUE, force_in=NULL, 
                             Nboot = 200, intercept = T, individual = NULL) {
  if(is.data.frame(xx)) xx <- as.matrix(xx)
  if(intercept) {
    xx <- cbind(1, xx); force_in = c(1, force_in+1); k <- k + 1
    if(!is.null(individual)) individual <- individual + 1
  }
  
  nn <- nrow(xx)
  p <- ncol(xx)
  
  # get varlist
  if (!k) {
    stop("Please input number of variables to select")
  } else {
    varlist <- gen_submodel_indices_k(k, p, force_in = force_in)
  }
  num_submodels <- nrow(varlist)
  
  if(return_sample) {
    zz <- rep(0, p)
    all_submodels <- apply(varlist, 1, function(k) {zz[k]=1;paste0(zz, collapse = "")})
  }
  
  #set.seed(123)
  EE <- matrix(rnorm(nn*Nboot), nrow = nn, ncol = Nboot)
  boot_sample_C <- max_t_mul_boot_by_k(xx, yy, EE, sandwich, Nboot, varlist, individual)
  
  boot_rank <- t(apply(boot_sample_C$max_t, 1, rank))
  
  if(!is.null(individual)) {
    boot_rank1 <- t(apply(boot_sample_C$max_t1, 1, rank))
    ret <- list(X = xx, Y = yy, NumVars = k, 
                beta = boot_sample_C$beta, se = boot_sample_C$se,
                BootSample = boot_sample_C$max_t, 
                BootRank = boot_rank,
                BootSample1 = boot_sample_C$max_t1, 
                BootRank1 = boot_rank1,
                AllModelNames = all_submodels, BootstrapSize = Nboot)
  } else {
    ret <- list(X = xx, Y = yy, NumVars = k, 
                beta = boot_sample_C$beta, se = boot_sample_C$se,
                BootSample = boot_sample_C$max_t, 
                BootRank = boot_rank,
                AllModelNames = all_submodels, BootstrapSize = Nboot)
  }
  
  ret$call <- match.call()
  class(ret) <- "MultiplierBootstrapSample"
  
  return(ret)
}


#' Return Max-T multiplier bootstrap sample of model M
#' 
#' @param xx The whole model matrix.
#' @param yy Response vector.
#' @param sandwhich Logical. If TRUE, the variance is obtained by the sandwich estimator; if FALSE, the variance is obtained by the model-based variance. TODO: bootstrap based se
#' @param return_sample Logical. If TRUE, return all the bootstrap sample; if false, otherwise.
#' @param M Index vector. The indices of variables correspond to the model of interest.
#' @param Nboot Numeric. The number of bootstrap sample.
#' @param intercept Logical. If TRUE, add an intercept to the model matrix.
#' @param individual Index vector. Indices of variables that we want the individual Max-T bootstrap sample. (posi-1).
#' @return model matrix xx, response vector yy, model index M, coefficient estimate beta, standard eror of \hat\beta se, 
## the matrix of T bootstrap sample BootSample, the vector of Max-T bootstrap sample BootSample BootSample1, bootstrap sample BootstrapSize
max_t_mul_boot_M <- function(xx, yy, sandwich=TRUE, return_sample=TRUE,
                             M, Nboot = 200, intercept, individual=TRUE) {
  if(is.data.frame(xx)) xx <- as.matrix(xx)
  nn <- nrow(xx)
  
  #set.seed(123)
  EE <- matrix(rnorm(nn*Nboot), nrow = nn, ncol = Nboot)
  xx <- xx[, M]
  if(intercept) xx <- cbind(1, xx)
  boot_sample_C <- tmax::max_t_mul_boot(xx, yy, EE, sandwich, Nboot)
  
  # TODO: update max_t_mul_boot
  boot_sample_C$se <- sqrt(diag(boot_sample_C$se))
  boot_sample_C$max_t1 <- abs(sweep(boot_sample_C$bootcoef, 1, boot_sample_C$se, "/"))
  boot_sample_C$max_t <- colMax(boot_sample_C$max_t1)
  
  if(intercept) boot_sample_C$max_t1 <- boot_sample_C$max_t1[2:nrow(boot_sample_C$max_t1),]
  
  ret <- list(X = xx, Y = yy, M = M,
              beta = boot_sample_C$coef, se = boot_sample_C$se,
              BootSample = boot_sample_C$max_t,
              BootSample1 = boot_sample_C$max_t1, 
              BootstrapSize = Nboot)
  
  ret$call <- match.call()
  class(ret) <- "MultiplierBootstrapSample"
  
  return(ret)
}

# Berk et al. PoSI
maxt_posi <- function(xx, yy, maxk = ncol(xx), sandwich = TRUE, alpha = 0.05, Nboot = 200) {
  n <- nrow(xx); p <- ncol(xx)
  maxt_0 <- rep(-Inf, Nboot)
  
  for(k in 1:maxk) {
    varlist <- gen_submodel_indices_k(k, p)
    #set.seed(123)
    EE <- matrix(rnorm(n*Nboot), nrow = n, ncol = Nboot)
    tmax_ret <- tmax::max_t_mul_boot_by_k(xx, yy, EE, sandwich = sandwich, Nboot, varlist)
    maxt_0 <- pmax(maxt_0, colMax(tmax_ret$max_t))
  }
  
  k0 <- quantile(maxt_0, 1 - alpha)
  ret <- list(x = xx, y = yy, k = k0, Nboot = Nboot, sandwich = sandwich)
  class(ret) <- "PoSI_Berk"
  
  return(ret)
}

# generic function for PoSI and PoSI_Berk
posi <- function(posi_obj, M, sigma = NULL) {
  UseMethod("posi", posi_obj)
}

# PoSI box and PoSI parallelpiped volume 
posi.PoSI_Berk <- function(posi_obj, M, sigma = NULL){
  # posi_obj = PoSI object
  # M = model as integer vector into columns of x
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

# PoSI box and PoSI parallelpiped volume 
verextreme.posix <- function(x, y, C1Alpha, M){
  Sigma <- t(x)%*%x/length(y)
  Sigma.M <- Sigma[M,M,drop=FALSE] # Covariance matrix of model M
  mod.lm <- lm(y ~ x[,M] - 1)
  beta.M <- unname(mod.lm$coeff)
  l1.Sigma.inv <- rowSums(abs(solve(Sigma.M))) * C1Alpha
  box <- cbind(beta.M - l1.Sigma.inv, beta.M + l1.Sigma.inv)
  vol.M <- (2*C1Alpha)^nrow(Sigma.M)/det(Sigma.M)
  colnames(box) <- c("Lower", "Upper")
  rownames(box) <- colnames(x)[M]
  return(list(intervals.M = box, Volume.M = vol.M))
}


get_T <- function(xx = NULL, 
                  yy = NULL, 
                  M = NULL,
                  transformation = NULL,
                  Hm, 
                  maxk = nrow(Hm), 
                  Nboot = 200, 
                  alpha = 0.05,
                  intercept = T,
                  adjust = F,
                  force_in = NULL) {
  # TODO: force.in
  Hstar <- apply(Hm, 2, cummax)
  dim(Hstar) <- dim(Hm)
  
  if(adjust) {
    if(class(M) == "numeric" | class(M) == "integer") {
      tmp <- M
      M <- rep(F, ncol(xx))
      M[tmp] <- T
    }
    
    Hsmax <- rep(0, Nboot)
    for(i in 1:maxk){
      Hscdf <- ecdf(Hstar[i,])
      Hsmax <- pmax(Hsmax, Hscdf(Hm[i,]))
    }
    K <- quantile(Hsmax, 1-alpha)
    
    if(is.null(transformation)) {
      ss <- sum(M)
    } else {
      ss <- transformation
    }
    
    Hmhat <- max_t_mul_boot_M(xx, yy, M=M, Nboot=Nboot, intercept = intercept)$BootSample  
    Hmhat <- Hmhat[order(Hmhat)]
    return(Hmhat[quantile(Hstar[ss,], K)])
    
    # Hmhat <- ecdf(max_t_mul_boot_M(xx, yy, M=M, Nboot=Nboot, intercept = intercept)$BootSample)
    # return(quantile(Hmhat, Hs(K)) )
  } else {
    K <- quantile(Hstar[maxk, ], 1-alpha, type=1)
    return(K)
  }
}

get_p <- function(tval,
                  xx = NULL, 
                  yy = NULL, 
                  M = NULL,
                  transformation = NULL,
                  Hm,
                  maxk = nrow(Hm), 
                  Nboot = 200, 
                  alpha = 0.05,
                  intercept = T,
                  adjust = F,
                  force_in = NULL) {
  # TODO: force.in
  # fit <- lm(yy ~ xx[, M])
  # se <- diag(sqrt(vcovHC(fit, "HC0")))
  # tval <- (coef(fit)/se)[-1]
  
  Hstar <- apply(Hm, 2, cummax)
  dim(Hstar) <- dim(Hm)
  
  if(adjust) {
    if(class(M) == "numeric" | class(M) == "integer") {
      tmp <- M
      M <- rep(F, ncol(xx))
      M[tmp] <- T
    }
    
    Hsmax <- rep(0, Nboot)
    for(i in 1:maxk){
      Hscdf <- ecdf(Hstar[i,])
      Hsmax <- pmax(Hsmax, Hscdf(Hm[i,]))
    }
    
    pval <- seq(0,1,0.001)
    K <- quantile(Hsmax, pval)

    if(is.null(transformation)) {
      ss <- sum(M)
    } else {
      ss <- transformation
    }

    Hmhat <- max_t_mul_boot_M(xx, yy, M=M, Nboot=Nboot, intercept = intercept)$BootSample
    # Hmhat_cdf <- ecdf(Hmhat)
    # Hm_tval <- (Hmhat_cdf(abs(tval)))
    # # Hmm_tval <- quantile(Hm[5,], Hm_tval)
    # Hstar_tval <- quantile(Hstar[ss, ], Hm_tval)
    # Hsmax <- Hsmax[order(Hsmax)]
    # 1-Hsmax[Hstar_tval]
    
    Hmhat <- Hmhat[order(Hmhat)]
    Hmhat_cdf <- ecdf(Hmhat[quantile(Hstar[ss,], K)])
    return(1-Hmhat_cdf(abs(tval)))
    
  } else {
    Hstar_cdf <- ecdf(Hstar[maxk, ])
    return(1-Hstar_cdf(tval))
  }
}

### TODO: buggy
# get_T1 <- function(xx = NULL, 
#                   yy = NULL, 
#                   M = NULL, 
#                   transformation = NULL,
#                   Hmp, 
#                   maxk = nrow(Hm), 
#                   Nboot = 200, 
#                   alpha = 0.05,
#                   intercept = T,
#                   adjust = F,
#                   force_in = NULL) {
#   # TODO: force.in
#   Hmp_ <- Hmp[,,M,drop=F]
#   Hstar <- apply(Hmp_, c(2,3), cummax)
#   dim(Hstar) <- dim(Hmp_)
#   
#   if(adjust) {
#     if(class(M) == "numeric" | class(M) == "integer") {
#       tmp <- M
#       M <- rep(F, ncol(xx))
#       M[tmp] <- T
#     }
#     
#     if(is.null(transformation)) {
#       ss <- sum(M)
#     } else {
#       ss <- transformation
#     }
#     Hsmax <- matrix(0, Nboot*ss, nrow = sum(M), ncol = Nboot)
#     for(j in 1:sum(M)) {
#       for(i in 1:maxk){
#         Hscdf <- ecdf(Hstar[i,,j])
#         Hsmax[j,] <- pmax(Hsmax[j,], Hscdf(Hmp_[i,,j]))
#       }
#     }
#       
#     K <- apply(Hsmax, 1, quantile, 1-alpha)
#     
#     Hmhat <- max_t_mul_boot_M(xx, yy, M=M, Nboot=Nboot, 
#                               intercept = intercept, individual = T)$BootSample1
#     Hmhat <- t(apply(matrix(Hmhat, ncol=Nboot), 1, sort))
#     
#     if(is.null(transformation)) {
#       # Kalpha <- mapply(function(Mj, k) quantile(Hstar[ss,,Mj], k), 1:ss, K)
#       # mapply(function(Mj, k) Hmhat[Mj, k], 1:ss, Kalpha)
#       return(mapply(function(Mj, k) Hmhat[Mj, quantile(Hstar[ss,,Mj], k)], 1:ss, K)) 
#     } else {
#       return(mapply(function(Mj, k) Hmhat[Mj, quantile(Hstar[ss,,Mj], k)], which(M), K)) 
#     }
#     
#     # Hmhat <- ecdf(max_t_mul_boot_M(xx, yy, M=M, Nboot=Nboot, intercept = intercept)$BootSample)
#     # return(quantile(Hmhat, Hs(K)) )
#   } else {
#     K <- apply(t(Hmp_[maxk,,]), 1, quantile, 1-alpha)
#     return(K)
#   }
# }
