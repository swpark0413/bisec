// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;




// [[Rcpp::export]]
arma::mat innerprod(const arma::mat A, const arma::mat B) {
  arma::mat m = A.t() * B;
  return m;
}

// [[Rcpp::export]]
arma::mat arma_matmul(const arma::mat A, const arma::mat B){
  return A * B;
}


// [[Rcpp::export]]
List splitList1(const List L){
  int niter = L.size();
  List A, B, C, subL;
  for(int i = 0; i < niter; i++){
    subL = L[i];
    A.push_back(as<NumericVector>(subL[0]));
    B.push_back(as<NumericVector>(subL[1]));
    C.push_back(as<NumericVector>(subL[2]));
  }
  return List::create(Named("cov") = A,
                      Named("evals") = B,
                      Named("evecs") = C);
}



// [[Rcpp::export]]
List splitList2(const List L){
  int niter = L.size();
  List A, B, C, D, E, F, subL;
  for(int i = 0; i < niter; i++){
    subL = L[i];
    A.push_back(as<NumericVector>(subL[0]));
    B.push_back(as<NumericVector>(subL[1]));
    C.push_back(as<NumericVector>(subL[2]));
    D.push_back(as<NumericMatrix>(subL[3]));
    E.push_back(as<NumericMatrix>(subL[4]));
    F.push_back(as<NumericVector>(subL[5]));
  }
  return List::create(Named("beta") = A,
                      Named("eta") = B,
                      Named("evals") = C,
                      Named("evecs") = D,
                      Named("sxy") = E,
                      Named("syy") = F);
}


// [[Rcpp::export]]
List splitList3(const List L){
  List A;
  int niter = L.size();
  NumericMatrix temp = as<NumericMatrix>(L[0]);
  int nrow = temp.nrow();
  int ncol = temp.ncol();
  for(int j = 0; j < ncol; j++){
    NumericMatrix res = no_init_matrix(nrow,niter);
    for(int i = 0; i < niter; i++){
      temp = as<NumericMatrix>(L[i]);
      res(_, i) = temp(_,j);
    }
    A.push_back(res);
  }
  return A;
}

