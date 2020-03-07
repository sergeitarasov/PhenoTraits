#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//' Matrix Exponential (scalar implementation)
//'
//' @param t time
//' @param A rate matrix
//'
//' @return probability matrix
//' @export
// [[Rcpp::export]]
arma::mat ExpMat_C(double t, arma::mat A) {
  arma::mat B = expmat(A*t);
  return(B);
}


//' Matrix Exponential (vector implementation)
//'
// [[Rcpp::export]]
arma::mat vExpMat_C(NumericVector t, arma::mat Q) {

  int tsize = t.size();
  int out_rows = pow(Q.n_cols, 2);

  //intitilize matrix to keep results
  arma::mat OUT(out_rows, tsize, arma::fill::zeros);
  for(int i = 0; i < tsize; ++i) {
    arma::mat Qout = expmat(Q*t[i]);
    arma::vec V = vectorise(Qout);
    OUT.col(i) = V;
  }
  return(OUT);
}


// [[Rcpp::export]]
arma::vec Exp_2Q_C(double t1, double t2, double Tau, arma::mat Q1, arma::mat Q2) {
  arma::mat B = expmat(Q1*(t2-t1))*expmat(Q2*(Tau-t2));
  arma::vec V = vectorise(B);
  return(V);
}


// [[Rcpp::export]]
arma::mat vExp_2Q_t1_C(NumericVector t1, double t2, double Tau, arma::mat Q1, arma::mat Q2) {
  int t1size = t1.size();
  int out_rows = pow(Q1.n_cols, 2);
  //intitilize matrix to keep results
  arma::mat OUT(out_rows, t1size, arma::fill::zeros);
  for(int i = 0; i < t1size; ++i) {
    OUT.col(i) = Exp_2Q_C(t1[i], t2, Tau, Q1, Q2);
  }
  return(OUT);
}



// [[Rcpp::export]]
arma::vec Exp_2Q_scen3_C(double t1, double Tau, arma::mat Q1, arma::mat Q2) {
  arma::mat B = expmat(Q1*t1) * expmat(Q2*(Tau-t1));
  arma::vec V = vectorise(B);
  return(V);
}


// [[Rcpp::export]]
arma::mat vExp_2Q_scen3_C(NumericVector t1, double Tau, arma::mat Q1, arma::mat Q2) {

  int t1size = t1.size();
  int out_rows = pow(Q1.n_cols, 2);

  //intitilize matrix to keep results
  arma::mat OUT(out_rows, t1size, arma::fill::zeros);
  for(int i = 0; i < t1size; ++i) {
    OUT.col(i) = Exp_2Q_scen3_C(t1[i], Tau, Q1, Q2);
  }
  return(OUT);
}



