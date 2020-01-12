#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <iostream>

using namespace Rcpp;

int geqrf(int m, int n, double* A, int lda, double *tau);
int ormqr(char side, char trans, int m, int n, int k, 
          double *A, int lda, double *tau, double* C, int ldc);
int trtrs(char uplo, char trans, char diag, 
          int n, int nrhs, 
          double* A, int lda, double* B, int ldb);
void dgemv(char* TRANS, const int* M, const int* N,
           double* alpha, double* A, const int* LDA, double* X,
           const int* INCX, double* beta, double* C, const int* INCY);
void dgemm(char* TRANSA, char* TRANSB, const int* M,
           const int* N, const int* K, double* alpha, double* A,
           const int* LDA, double* B, const int* LDB, double* beta,
           double* C, const int* LDC);
void print_matrix(double* Xr, int n, int d);
void print_vector(double* Xr, int n);

void inverse(double* A, int N)
{
  int *IPIV = new int[N+1];
  int LWORK = N*N;
  double *WORK = new double[LWORK];
  int INFO;
  
  dgetrf_(&N,&N,A,&N,IPIV,&INFO);
  dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
  
  delete[] IPIV;
  delete[] WORK;
}

int geqrf(int m, int n, 
          double* A, int lda, double *tau) {
  int info = 0, lwork = -1;
  double iwork;
  
  dgeqrf_(&m, &n, A, &lda, tau, 
          &iwork, &lwork, &info);
  lwork = (int)iwork;
  double* work = new double[lwork];
  
  dgeqrf_(&m, &n, A, &lda, tau, 
          work, &lwork, &info);
  delete[] work;
  return info;
}

int ormqr(char side, char trans, int m, int n, int k, 
          double *A, int lda, double *tau, double* C, int ldc) {
  int info = 0, lwork = -1;
  double iwork;
  
  dormqr_(&side, &trans, &m, &n, &k, 
          A, &lda, tau, C, &ldc, &iwork, &lwork, &info);
  lwork = (int)iwork;
  double* work = new double[lwork];
  
  dormqr_(&side, &trans, &m, &n, &k, 
          A, &lda, tau, C, &ldc, work, &lwork, &info);
  
  delete[] work;
  return info;
}

int trtrs(char uplo, char trans, char diag, 
          int n, int nrhs, 
          double* A, int lda, double* B, int ldb) {
  int info = 0;
  dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, 
          A, &lda, B, &ldb, &info);
  return info;
}

// alpha A*b + beta y
void matvecprod(char no, double* A, double* b, double* y, 
                double alpha, double beta, 
                int M, int N) { 
  int m= M, n= N, lda= M, incx= 1, incy= 1;
  
  dgemv_(&no, &m, &n, &alpha, 
         A, &lda, b, &incx, &beta, y, &incy);
}

// alpha A*B + beta C
int matmatprod(char transA, char transB, 
               int m, int n, int k,
               double alpha, double* A, int lda,
               double* B, int ldb, double beta, 
               double* C, int ldc) {
  int info = 0;
  dgemm_(&transA, &transB, &m, &n, &k, 
         &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
  return info;
}


// [[Rcpp::export]]
List lm_qr (NumericMatrix& Xr, NumericMatrix &y) {
  int n = Xr.nrow(), d = Xr.ncol(), dy = y.ncol();
  NumericMatrix A = clone(Xr), b = clone(y);
  NumericVector beta(d), tau(d);
  
  geqrf(n, d, A.begin(), n, tau.begin());
  
  ormqr('L', 'T', n, dy, d, A.begin(), n, tau.begin(), b.begin(), n);
  
  trtrs('U', 'N', 'N', d, dy, A.begin(), n, b.begin(), n);
  
  // for (int i=0; i < d; i++) {beta[i] = b[i];}
  return List::create(_("coefficient") = b,
                      _("qr")          = A,
                      _("qraux")       = tau);
  
}

// [[Rcpp::export]]
double max_t (NumericMatrix& Xr, NumericVector &y) {
  int n = Xr.nrow(), d = Xr.ncol(), dy = 1, dresid = n;
  double max_t = -1;
  NumericMatrix A = clone(Xr), resid_mat(n, n), se(d, d);
  NumericVector tau(d), tstats(d), b = clone(y), resid = clone(y), rn = rnorm(n);
  
  // b = argmin|| y - X beta ||
  geqrf(n, d, A.begin(), n, tau.begin());
  
  ormqr('L', 'T', n, dy, d, A.begin(), n, tau.begin(), b.begin(), n);
  
  trtrs('U', 'N', 'N', d, dy, A.begin(), n, b.begin(), n);
  
  // resid = y - X b
  matvecprod('N', Xr.begin(), b.begin(), resid.begin(), -1.0, 1.0, n, d);
  resid_mat = diag(resid);
  
  // SE
  A = clone(Xr);
  geqrf(n, d, A.begin(), n, tau.begin());
  
  ormqr('L', 'T', n, dresid, d, A.begin(), n, tau.begin(), resid_mat.begin(), n);
  
  trtrs('U', 'N', 'N', d, dresid, A.begin(), n, resid_mat.begin(), n);
  
  matmatprod('N', 'T', d, d, n, 1.0, 
             resid_mat.begin(), n, resid_mat.begin(), n, 0.0,
             se.begin(), d);
  
  // bb = argmin|| rn*resid - X bb || 
  A = clone(Xr);
  for (int i = 0; i < n; i++) {resid[i] *= rn[i];}
  geqrf(n, d, A.begin(), n, tau.begin());
  
  ormqr('L', 'T', n, dy, d, A.begin(), n, tau.begin(), resid.begin(), n);
  
  trtrs('U', 'N', 'N', d, dy, A.begin(), n, resid.begin(), n);
  
  // max t
  for (int j = 0; j < d; j++) {tstats[j] = resid[j]/sqrt(se(j,j));}
  max_t = max(abs(tstats));
  
  return max_t;
}

// // [[Rcpp::export]]
// NumericVector max_t_mul_boot (NumericMatrix& Xr, NumericVector &y, int nboot) {
//     int n = Xr.nrow(), d = Xr.ncol(), dy = 1, dresid = n;
//     NumericMatrix A = clone(Xr), resid_mat(n, n), se(d, d);
//     NumericVector tau(d), tstats(d), max_t(nboot), rn(n),
//         b = clone(y), resid = clone(y);

//     // b = argmin|| y - X beta ||
//     geqrf(n, d, A.begin(), n, tau.begin());

//     ormqr('L', 'T', n, dy, d, A.begin(), n, tau.begin(), b.begin(), n);

//     trtrs('U', 'N', 'N', d, dy, A.begin(), n, b.begin(), n);

//     // resid = y - X b
//     matvecprod('N', Xr.begin(), b.begin(), resid.begin(), -1.0, 1.0, n, d);
//     resid_mat = diag(resid);

//     // SE
//     A = clone(Xr);
//     geqrf(n, d, A.begin(), n, tau.begin());

//     ormqr('L', 'T', n, dresid, d, A.begin(), n, tau.begin(), resid_mat.begin(), n);

//     trtrs('U', 'N', 'N', d, dresid, A.begin(), n, resid_mat.begin(), n);

//     matmatprod('N', 'T', d, d, n, 1.0, 
//         resid_mat.begin(), n, resid_mat.begin(), n, 0.0,
//         se.begin(), d);

//     // bootstrap bb = argmin|| rn*resid - X bb || 
//     for (int iboot = 0; iboot < nboot; iboot++) {
//         rn = rnorm(n);
//         A = clone(Xr);
//         for (int i = 0; i < n; i++) {resid[i] *= rn[i];}
//         geqrf(n, d, A.begin(), n, tau.begin());

//         ormqr('L', 'T', n, dy, d, A.begin(), n, tau.begin(), resid.begin(), n);

//         trtrs('U', 'N', 'N', d, dy, A.begin(), n, resid.begin(), n);

//         // max t
//         for (int j = 0; j < d; j++) {tstats[j] = resid[j]/sqrt(se(j,j));}
//         max_t[iboot] = max(abs(tstats));
//     }

//     return max_t;
// }


// [[Rcpp::export]]
List max_t_mul_boot (NumericMatrix& Xr, NumericVector &y, 
                     NumericMatrix& normmat, int sandwich, int nboot) {
  int n = Xr.nrow(), d = Xr.ncol(), dy = 1, dresid = n;
  NumericMatrix A = clone(Xr), resid_mat(n, n), se(d, d), ee = clone(normmat);
  NumericVector tau(d), tstats(d), max_t(nboot),
  b = clone(y), resid = clone(y);
  
  // b = argmin|| y - X beta ||
  geqrf(n, d, A.begin(), n, tau.begin());
  
  ormqr('L', 'T', n, dy, d, A.begin(), n, tau.begin(), b.begin(), n);
  
  trtrs('U', 'N', 'N', d, dy, A.begin(), n, b.begin(), n);
  
  // resid = y - X b
  matvecprod('N', Xr.begin(), b.begin(), resid.begin(), -1.0, 1.0, n, d);
  
  if(sandwich) {
    resid_mat = diag(resid);
    
    // SE
    A = clone(Xr);
    geqrf(n, d, A.begin(), n, tau.begin());
    ormqr('L', 'T', n, dresid, d, A.begin(), n, tau.begin(), resid_mat.begin(), n);
    trtrs('U', 'N', 'N', d, dresid, A.begin(), n, resid_mat.begin(), n);
    matmatprod('N', 'T', d, d, n, 1.0, 
               resid_mat.begin(), n, resid_mat.begin(), n, 0.0,
               se.begin(), d);
  } else {
    matmatprod('T', 'N', d, d, n, 1.0, 
               Xr.begin(), n, Xr.begin(), n, 0.0,
               se.begin(), d);
    inverse(se.begin(), d);
  }
  
  // boot
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < nboot; j++) {
      ee(i,j) *= resid[i];
    }
  }
  A = clone(Xr);
  geqrf(n, d, A.begin(), n, tau.begin());
  
  ormqr('L', 'T', n, nboot, d, A.begin(), n, tau.begin(), ee.begin(), n);
  
  trtrs('U', 'N', 'N', d, nboot, A.begin(), n, ee.begin(), n);
  
  return List::create(_("coef")   = b[Range(0,d-1)],
                      _("resid")  = resid,
                      _("se")     = se,
                      _("bootcoef")  = ee(Range(0,d-1), _));
}

// [[Rcpp::export]]
List max_t_mul_boot_by_k (NumericMatrix& Xr, NumericVector &y, 
                          NumericMatrix& normmat, 
                          int sandwich, int nboot, 
                          NumericMatrix& ind_mat,
                          Nullable<NumericVector> individual_ = R_NilValue) {
  int n = Xr.nrow(), dy = 1, dresid = n, n_mod = ind_mat.nrow(), k = ind_mat.ncol();
  NumericMatrix maxt(n_mod, nboot), XXr(n, k);
  NumericMatrix beta_mat(n_mod, k), se_mat(n_mod, k);
  NumericVector individual; NumericMatrix maxt1(Xr.ncol(), nboot);
  if(individual_.isNotNull()) {
    individual = individual_;
    std::fill(maxt1.begin(), maxt1.end(), 0.0);
  }

  // std::cout << "dev" << std::endl;
  
  for(int i_mod = 0; i_mod < n_mod; i_mod++) {
    // subset matrix X
    for (int i=0; i<k; i++) {
      XXr(_,i) = Xr(_, ind_mat(i_mod, i)-1);
    }
    
    NumericMatrix A = clone(XXr), resid_mat(n, n), se(k, k), ee = clone(normmat);
    NumericVector tau(k), tstats(k), b = clone(y), resid = clone(y);

    // QR
    geqrf(n, k, A.begin(), n, tau.begin());

    // b = argmin|| y - X beta ||
    ormqr('L', 'T', n, dy, k, A.begin(), n, tau.begin(), b.begin(), n);
    trtrs('U', 'N', 'N', k, dy, A.begin(), n, b.begin(), n);
    beta_mat(i_mod, _) = b[Range(0,k-1)];
    
    // resid = y - X b
    matvecprod('N', XXr.begin(), b.begin(), resid.begin(), -1.0, 1.0, n, k);
    
    if(sandwich) {
      resid_mat = diag(resid);
      
      // SE
      // A = clone(XXr);
      // geqrf(n, k, A.begin(), n, tau.begin());
      ormqr('L', 'T', n, dresid, k, A.begin(), n, tau.begin(), resid_mat.begin(), n);
      trtrs('U', 'N', 'N', k, dresid, A.begin(), n, resid_mat.begin(), n);
      matmatprod('N', 'T', k, k, n, 1.0, 
                 resid_mat.begin(), n, resid_mat.begin(), n, 0.0,
                 se.begin(), k);
    } else {
      matmatprod('T', 'N', k, k, n, 1.0, 
                 XXr.begin(), n, XXr.begin(), n, 0.0,
                 se.begin(), k);
      inverse(se.begin(), k);
    }
    
    for(int i = 0; i < k; i++) {
      se_mat(i_mod, i) = sqrt(se(i,i));
    }
    
    // boot
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nboot; j++) {
        ee(i,j) *= resid[i];
      }
    }
    // A = clone(XXr);
    // geqrf(n, k, A.begin(), n, tau.begin());
    ormqr('L', 'T', n, nboot, k, A.begin(), n, tau.begin(), ee.begin(), n);
    trtrs('U', 'N', 'N', k, nboot, A.begin(), n, ee.begin(), n);  
    
    for(int iboot = 0; iboot < nboot; iboot++) {
      for(int i = 0; i < k; i++) {
        // maxt(i,iboot) = bootcoef(i,iboot)/se(i,i);
        double tmp = std::abs(ee(i, iboot)/sqrt(se(i, i)));
        if(maxt(i_mod, iboot) < tmp) maxt(i_mod, iboot) = tmp;
      }
    }

    if(individual_.isNotNull()) {
      for(int i = 0; i < individual.length(); ++i) {
        int ind1 = -1;
        for(int j = 0; j < k; ++j) {
          if(ind_mat(i_mod, j) == individual(i)) {ind1 = j; break;}
        }
        if(ind1 > -1) {
          for(int iboot = 0; iboot < nboot; ++iboot) {
            double tmp = std::abs(ee(ind1, iboot)/sqrt(se(ind1, ind1)));
            if(maxt1(i, iboot) < tmp) maxt1(i, iboot) = tmp;
          }
        }
      }
    }
    
    // coef_se[[i]] <- rbind(ret$coef[1:opt$k], sqrt(diag(ret$se)))
    // scaledt <- ret$bootcoef[1:opt$k,]/coef_se[[i]][2,]
    // max_scaled_t <- apply(as.matrix(abs(scaledt)), 2, max)
  }
  
  if(individual_.isNotNull()) {
    return List::create(_("max_t")   = maxt,
                        _("max_t1")   = maxt1(Range(0, individual.length()-1), _),
                      _("beta")    = beta_mat,
                      _("se")      = se_mat);
  } else {
    return List::create(_("max_t")   = maxt,
                      _("beta")    = beta_mat,
                      _("se")      = se_mat);
  }
}

// [[Rcpp::export]]
NumericMatrix inverse_C (NumericMatrix& xx, int N){
  inverse(xx.begin(), N);
  return xx;
}

// [[Rcpp::export]]
NumericMatrix cov_mat_C (NumericMatrix& Xr){
  int n = Xr.nrow(), d = Xr.ncol();
  NumericMatrix se(d, d);
  
  matmatprod('T', 'N', d, d, n, 1.0, 
             Xr.begin(), n, Xr.begin(), n, 0.0,
             se.begin(), d);
  
  return se;
}


// [[Rcpp::export]]
List true_max_t (NumericMatrix& xtx, NumericMatrix &Xr, NumericVector &y, 
                 NumericVector &beta_true, NumericMatrix& ind_mat) {
  int n = Xr.nrow(), d = Xr.ncol(), n_mod = ind_mat.nrow(), k = ind_mat.ncol();
  NumericMatrix XXr(n, k), xtxr(k, k);
  NumericVector resid = clone(y), tstat(n_mod);
  
  matvecprod('N', Xr.begin(), beta_true.begin(), resid.begin(), -1.0, 1.0, n, d);
  
  for(int i_mod = 0; i_mod < n_mod; i_mod++) {
    Rcout << "Model: " << i_mod << "\n";
    // subset covariance mat
    for (int i=0; i<k; i++) {
      for (int j=0; j<k; j++) {
        xtxr(i,j) = xtx(ind_mat(i_mod,i)-1, ind_mat(i_mod,j)-1); 
      }
    }
    Rcout << xtxr << "\n";
    
    // subset matrix X
    for (int i=0; i<k; i++) {
      XXr(_,i) = Xr(_, ind_mat(i_mod, i)-1);
    }
    
    NumericMatrix A = clone(xtxr);
    NumericVector tau(k), xproj(k);  
    
    matvecprod('T', XXr.begin(), resid.begin(), xproj.begin(), -1.0, 0, k, n);
    
    Rcout << "xproj: " <<  xproj << "\n";
    
    // b = argmin|| y - X beta ||
    NumericVector b = clone(xproj);
    geqrf(k, k, A.begin(), k, tau.begin());
    ormqr('L', 'T', k, 1, k, A.begin(), k, tau.begin(), b.begin(), k);
    trtrs('U', 'N', 'N', k, 1, A.begin(), k, b.begin(), k);

    Rcout << "b: " << b << "\n";
    
    inverse(xtxr.begin(), k);
    double tmp_max = 0, tmp = 0;
    for(int i=0; i<k; i++) {
       tmp = sqrt(xtxr(i,i)) * b[i];
      if(tmp_max < tmp) tmp_max = tmp;
    }
    tstat[i_mod] = tmp;
  }
  
  return List::create(_("max_true")   = tstat);
}

// [[Rcpp::export]]
NumericVector residual(NumericMatrix& Xr, NumericVector &beta, NumericVector &y) {
  int n = Xr.nrow(), d = Xr.ncol();
  NumericVector ret(n), yy = clone(y);
  
  matvecprod('N', Xr.begin(), beta.begin(), yy.begin(), -1.0, 1.0, n, d);
  
  return yy;
}

// [[Rcpp::export]]
List QR(NumericMatrix& Xr) {
  int n = Xr.nrow(), d = Xr.ncol();
  NumericMatrix A = clone(Xr);
  NumericVector tau(d);
  // double* A = new double[n*d];
  // NumericMatrix ret(n, d);
  
  //Lapack has column-major order
  // for(int col=0, D1_idx=0; col<d; ++col)
  // {
  //     for(int row = 0; row<n; ++row)
  //     {
  //         // Lapack uses column major format
  //         A[D1_idx++] = Xr(row, col);
  //     }
  // }
  
  geqrf(n, d, A.begin(), n, tau.begin());
  
  // for (int i = 0; i < n*d; i++) {ret[i] = A[i];}
  
  return List::create(_("R")  = A,
                      _("tau")  = tau);
}

// [[Rcpp::export]]
NumericVector compute_qb(NumericMatrix& Xr, NumericVector& tau, NumericVector& b) {
  int n = Xr.nrow(), d = Xr.ncol();
  
  ormqr('L', 'T', n, 1, d, Xr.begin(), n, tau.begin(), b.begin(), n);
  
  return b;
}

// [[Rcpp::export]]
NumericVector compute_b(NumericMatrix& Xr, NumericVector& b) {
  int n = Xr.nrow(), d = Xr.ncol();
  
  trtrs('U', 'N', 'N', d, 1, Xr.begin(), n, b.begin(), n);
  
  return b;
}


void print_matrix(double* Xr, int n, int d) {
  int i, j;
  
  for (i = 0; i < n; i++) {
    for (j = 0; j < d; j++) {
      Rcout << Xr[i+j*n] << "\t";
    }
    Rcout << std::endl;
  }
}


void print_vector(double* Xr, int n) {
  int i;
  
  for (i = 0; i < n; i++) {
    if (i % 10 == 0 ) {Rcout << std::endl;}
    Rcout << Xr[i] << "\t";
  }
}

