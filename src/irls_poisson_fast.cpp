#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace R;

//' @title Fast Iterative Reweighed Least Square algorithm for Poissons
//' @description IRLS to estimate network score of Poisson nodes.
//' @keywords internal
//' @export
// [[Rcpp::export]]
double factorial_fast(double n)
{
  return (n == 1.0 || n == 0.0) ? 1 : factorial_fast(n - 1.0) * n;
}

// [[Rcpp::export]]

Rcpp::List irls_poisson_cpp_fast(arma::mat A, arma::vec b, double maxit, double tol)
{
//Def
  arma::vec x;
  x.zeros(A.n_cols,1);
  arma::vec xold;
  arma::mat varmatrix;

  double nobs;
  nobs = A.n_rows;
  double ll;
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
    x = arma::solve(varmatrix, A.t()*(W % z), arma::solve_opts::likely_sympd);

    if(sqrt(pow(arma::norm(x-xold), 2)) < tol){
      break;
    }
  }

  arma::vec e;
  e = (b - eta);

  df = A.n_cols;

//scores
  for (int j = 0; j < nobs; ++j) {
    f[j] = log(factorial_fast(1.0 * b[j]));
  }
  ll = arma::accu(b % (eta) - exp(eta) - f);

  aic = - 2 * ll + 2 * df;
  bic = - 2 * ll + log(nobs) * df;
  mdl = 1;

//mdl

// arma::mat xz;
// xz.zeros(size(x));
//
// arma::vec ez;
// double ssrz;
// double ssrtot;
// double RR;
// double F;
//
// arma::vec yaverage(nobs);
//
// ez = (b - A*xz);
// ssrz = accu(ez.t()*ez);
// F = (((ssrz - ssr)/df)/(ssr/((nobs-(df + 1)))));
//
// for (int j=0; j < nobs; ++j) {
//   yaverage[j] = b[j] - arma::mean(b);
// }
//
// ssrtot = accu(yaverage.t()*yaverage);
//
// RR = 1-(ssr/ssrtot);
//
// if (RR > (df/nobs)) {
//   mdl = (nobs/2) * log(ssr/(nobs-df)) + (df/2) * log(F) + log(nobs);
// } else {
//   mdl = (nobs/2) * log((accu(b.t()*b))/nobs) + 0.5 * log(nobs);
// }


//return
return Rcpp::List::create(
  Rcpp::Named("loglik") = ll,
  Rcpp::Named("aic") = aic,
  Rcpp::Named("bic") = bic,
  Rcpp::Named("mdl") = mdl
  );
}


// note that valgrind might issue some possible memory loss.
// IMO false positive, see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=36298
