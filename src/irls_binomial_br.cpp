#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//' @title BR Iterative Reweighed Least Square algorithm for Binomials
//' @description IRLS to estimate network score of Binomial nodes.
//' @keywords internal
//' @export
// [[Rcpp::export]]

Rcpp::List irls_binomial_cpp_br(arma::mat A, arma::vec b, double maxit, double tol)
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
arma::vec e;
double ssr;

arma::vec W(nobs);
arma::vec unit(nobs);
unit.ones(nobs);
arma::vec eta(nobs);
arma::vec g(nobs);
arma::vec gprime(nobs);
arma::vec z(nobs);
//mod
arma::vec bprime(nobs);

//int k;

for (int i = 0; i < maxit; ++i) {
  eta = A * x;

  for (int j=0; j < nobs; ++j) {
    g[j] = 1.0 / (1.0 + exp(-1.0 * eta[j]));
    gprime[j] = exp (-1.0 * eta[j]) / ((1.0 + exp (-1.0 * eta[j])) * (1.0 + exp (-1.0 * eta[j])));
        //mod
    bprime[j] = (b[j]+(df/nobs)*(0.5))/(1+(df/nobs));

  }

  //z = eta+(b-g)/gprime;
  z = eta+(bprime-g)/gprime;

  W = (gprime % gprime);
  W /= (g % (unit-g));
  //W += unit*(df/nobs);
  xold = x;

  //coefficients
  //x = arma::solve(A.t()*(W % A.each_col()), A.t()*(W % z), arma::solve_opts::no_approx);
  varmatrix = A.t()*(W % A.each_col());
  x = arma::solve(varmatrix, A.t()*(W % z), arma::solve_opts::no_approx);
  //k = i;

if(sqrt(arma::dot(x-xold,x-xold)) < tol){
 break;
}}

//loglik

ll = arma::accu(-arma::dot(b,log(unit + exp(-(A*x)))) - arma::dot((unit-b),log(unit + exp(A*x))));

//aic

aic = - 2 * ll + 2 * df;

//bic

bic = - 2 * ll + log(nobs) * df;

//sse
e = (b - A*x);
ssr = accu(e.t()*e);

mdl = 1;

//return
return Rcpp::List::create(
  Rcpp::Named("coefficients") = x,
  Rcpp::Named("loglik") = ll,
  Rcpp::Named("aic") = aic,
  Rcpp::Named("bic") = bic,
  Rcpp::Named("mdl") = mdl,
  Rcpp::Named("sse") = ssr,
  Rcpp::Named("varcov") = varmatrix
);

}
