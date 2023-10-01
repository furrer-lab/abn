#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace R;

//' @title Fast Factorial
//' @description Calculate the factorial in C##
//' @keywords internal
//' @export
// [[Rcpp::export]]
double factorial_fast(double n)
{
  return (n == 1.0 || n == 0.0) ? 1 : factorial_fast(n - 1.0) * n;
}

//' @title Fast Iterative Reweighed Least Square algorithm for Poissons
//' @description IRLS to estimate network score of Poisson nodes.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List irls_poisson_cpp_fast(arma::mat A, arma::vec b, double maxit, double tol)
{
  // // print out A
  // Rcpp::Rcout << "A: " << A << std::endl;

  // // print out b
  // Rcpp::Rcout << "b: " << b << std::endl;

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
    // // print out i
    // Rcpp::Rcout << "i: " << i << std::endl;

    eta = A * x;
    // // print out eta
    // Rcpp::Rcout << "eta: " << eta << std::endl;

    W = exp(eta);
    // // print out W
    // Rcpp::Rcout << "W: " << W << std::endl;

    z = eta+(b-W)/W;
    // // print out z
    // Rcpp::Rcout << "z: " << z << std::endl;

    xold = x;
    // // print out xold
    // Rcpp::Rcout << "xold: " << xold << std::endl;

    //coefficients
    varmatrix = A.t()*(W % A.each_col());
    // // print out varmatrix
    // Rcpp::Rcout << "varmatrix: " << varmatrix << std::endl;

    // if varmat has no nan solve it else raise warning and break
    if (varmatrix.has_nan() || varmatrix.has_inf()) {
      Rcpp::warning("irls_poisson_fast.cpp: varmatrix has nan. Stopping early at iteration %d", i);
      break;
    } else {
      // Rcpp::Rcout << "varmatrix has no nan" << std::endl;
      x = arma::solve(varmatrix, A.t()*(W % z), arma::solve_opts::likely_sympd);
      // // print out x
      // Rcpp::Rcout << "x: " << x << std::endl;
    }

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
