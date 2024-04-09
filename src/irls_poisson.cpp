#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace R;
//' @title Factorial
//' @description Calculate the factorial in C##
//' @keywords internal
//' @returns a double
//' @export
// [[Rcpp::export]]
double factorial(double n)
{
  return (n == 1.0 || n == 0.0) ? 1 : factorial(n - 1.0) * n;
}


//' @title Iterative Reweighed Least Square algorithm for Poissons
//' @description IRLS to estimate network score of Poisson nodes.
//' @keywords internal
//' @returns a list
//' @export
// [[Rcpp::export]]

Rcpp::List irls_poisson_cpp(arma::mat A, arma::vec b, double maxit, double tol)
{
//Def
arma::vec x;
x.zeros(A.n_cols,1);
arma::vec xold;
arma::mat varmatrix;

double nobs;
nobs = A.n_rows;
double ll;
arma::vec e;
double ssr;
double aic;
double bic;
double df;
double mdl;


arma::vec W(nobs);
arma::vec unit(nobs);
unit.ones(nobs);
arma::vec eta(nobs);
arma::vec g(nobs);
arma::vec f(nobs);
arma::vec gprime(nobs);
arma::vec z(nobs);

for (int i = 0; i < maxit; ++i) {
  eta = A * x;

  W = exp(eta);

  z = eta+(b-W)/W;
  xold = x;

  //coefficients
  varmatrix = A.t()*(W % A.each_col());
  x = arma::solve(varmatrix, A.t()*(W % z), arma::solve_opts::no_approx);

if(sqrt(pow(arma::norm(x-xold), 2)) < tol){
 break;
}}

df = A.n_cols;

//loglik

    for (int j = 0; j < nobs; ++j) {
      f[j] = log(factorial(1.0 * b[j]));
      }

    ll = arma::accu(b % (eta) - exp(eta) - f);

aic = - 2 * ll + 2 * df;

bic = - 2 * ll + log(nobs) * df;

//sse
e = (b - eta);
ssr = arma::dot(e, e);

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



