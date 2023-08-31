#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//' @title Fast Iterative Reweighed Least Square algorithm for Binomials
//' @description IRLS to estimate network score of Binomial nodes.
//' @keywords internal
//' @export
// [[Rcpp::export]]

Rcpp::List irls_binomial_cpp_fast(arma::mat A, arma::vec b, double maxit, double tol)
{
//Def
arma::vec x;
x.zeros(A.n_cols,1);
arma::vec xold;
arma::mat varmatrix;

double nobs;
nobs = A.n_rows;
double df;
df = A.n_cols;
//double n;

double ll;
double aic;
double bic;
double mdl;

arma::vec W(nobs);
arma::vec unit(nobs);
unit.ones(nobs);
arma::vec eta(nobs);
arma::vec g(nobs);
arma::vec gprime(nobs);
arma::vec z(nobs);
//mod
//arma::vec bprime(nobs);

//int k;

for (int i = 0; i < maxit; ++i) {
  eta = A * x;

  for (int j=0; j < nobs; ++j) {
    double e = exp(-1.0 * eta[j]);
    g[j] = 1.0 / (1.0 + e);
    gprime[j] = e / ((1.0 + e) * (1.0 + e));

    //mod
    //bprime[j] = (b[j]+(df/nobs)*(0.5))/(1+(df/nobs));
    //bprime[j] = b[j]+(sum(b)/nobs)*(0.5);
  }

  //z = eta+(b-g)/gprime;
  z = eta+(b-g)/gprime;

  W = (gprime % gprime);
  //mod
  //W %= (unit+(unit*(sum(b)/nobs)));
  W /= (g % (unit-g));
  //W += unit;
  xold = x;

  //coefficients
  //x = arma::solve(A.t()*(W % A.each_col()), A.t()*(W % z), arma::solve_opts::no_approx);
  varmatrix = A.t()*(W % A.each_col());
  //varmatrix %=(g-0.5*unit);
  x = arma::solve(varmatrix, A.t()*(W % z), arma::solve_opts::no_approx + arma::solve_opts::fast + arma::solve_opts::likely_sympd);
  //k = i;

if(sqrt(arma::dot(x-xold,x-xold)) < tol){
 break;
}}

//n = A.n_rows;

//arma::vec e;
//double ssr;
//e = (b - A*x);
//e = (bprime - A*x);
//ssr = accu(e.t()*e);

//scores

//ll = arma::accu(-arma::dot(b,log(unit + exp(-(A*x)))) - arma::dot((unit-b),log(unit + exp(A*x))));
ll = arma::accu(-arma::dot(b,log(unit + exp(-(eta)))) - arma::dot((unit-b),log(unit + exp(eta))));

aic = - 2 * ll + 2 * df;

bic = - 2 * ll + log(nobs) * df;

//mdl

// arma::mat xz;
// xz.zeros(size(x));
//
// arma::vec ez;
// double ssrz;
// double ssrtot;
// double RR;
// double F;
// double mdl;
// arma::vec yaverage(n);
//
// ez = (b - A*xz);
// ssrz = accu(ez.t()*ez);
// F = (((ssrz - ssr)/df)/(ssr/((n-(df + 1)))));
//
// for (int j=0; j < n; ++j) {
// yaverage[j] = b[j] - arma::mean(b);
// }
//
// ssrtot = accu(yaverage.t()*yaverage);
//
// RR = 1-(ssr/ssrtot);
//
// if (RR > (df/n)) {
// mdl = (n/2) * log(ssr/(n-df)) + (df/2) * log(F) + log(n);
// } else {
// mdl = (n/2) * log((accu(b.t()*b))/n) + 0.5 * log(n);
// }

mdl = 1;

//return
return Rcpp::List::create(
  Rcpp::Named("loglik") = ll,
  Rcpp::Named("aic") = aic,
  Rcpp::Named("bic") = bic,
  Rcpp::Named("mdl") = mdl
);

}
