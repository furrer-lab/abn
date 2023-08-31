#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//' @title Fast Iterative Reweighed Least Square algorithm for Gaussians
//' @description IRLS to estimate network score of Gaussian nodes.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List irls_gaussian_cpp_fast(arma::mat A, arma::vec b, double maxit, double tol)
{

//Def
arma::mat x;
x.zeros(A.n_cols,1);

arma::mat gprime;
gprime.ones(A.n_rows,1);

arma::vec e;
double ssr;

double ll;
double n;
double df;
double aic;
double bic;
double mdl;

//coefficients
x = arma::solve(A.t() * (gprime % A.each_col()), A.t() * (gprime % b));

//loglik
e = (b - A*x);
ssr = accu(e.t()*e);
n = A.n_rows;
df = A.n_cols + 1;

//scores

ll = 0.5 * ( - n * (log(2 * arma::datum::pi) + 1 - log(n) + log(ssr)));

aic = - 2 * ll + 2 * df;

bic = - 2 * ll + log(n) * df;


//mdl



// arma::mat xz;
// xz.zeros(size(x));
// arma::vec ez;
//
//
// double ssrz;
// double ssrtot;
// double RR;
// double F;
//
//
//
// arma::vec yaverage(n);
//
// ez = (b - A*xz);
// ssrz = accu(ez.t()*ez);
// F = (((ssrz - ssr)/df)/(ssr/((n-(df + 1)))));
//
// for (int j=0; j < (n); ++j) {
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

return Rcpp::List::create(
  Rcpp::Named("loglik") = ll,
  Rcpp::Named("aic") = aic,
  Rcpp::Named("bic") = bic,
  Rcpp::Named("mdl") = mdl
);
}
