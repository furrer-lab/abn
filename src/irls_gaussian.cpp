#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//' @title Iterative Reweighed Least Square algorithm for Gaussians
//' @description IRLS to estimate network score of Gaussian nodes.
//' @keywords internal
//' @export
// [[Rcpp::export]]

Rcpp::List irls_gaussian_cpp(arma::mat A, arma::mat b, double maxit, double tol)
{

//Def
arma::mat x;
x.zeros(A.n_cols,1);
//arma::mat varmatrix;

arma::mat gprime;
gprime.ones(A.n_rows,1);

arma::vec e;
double ssr;
double ll;
double n;
double aic;
double bic;
double mdl;
double df;
df = A.n_cols + 1;

//coefficients
x = arma::solve(A.t() * (gprime % A.each_col()), A.t() * (gprime % b));

//varmatrix = arma::solve(A,A,arma::solve_opts::no_approx);

//loglik
e = (b - A*x);
ssr = accu(e.t()*e);
n=A.n_rows;

ll = 0.5 * ( - n * (log(2 * arma::datum::pi) + 1 - log(n) + log(ssr)));

aic = - 2 * ll + 2 * df;

bic = - 2 * ll + log(n) * df;

mdl = 1;

return Rcpp::List::create(
  Rcpp::Named("coefficients") = x,
  Rcpp::Named("loglik") = ll,
    Rcpp::Named("aic") = aic,
  Rcpp::Named("bic") = bic,
  Rcpp::Named("mdl") = mdl,
  Rcpp::Named("sse") = ssr
  //Rcpp::Named("varcov") = varmatrix
);
}
