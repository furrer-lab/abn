#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include <Rmath.h> // lgammafn

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace R;

//' @title Fast Iterative Reweighed Least Square algorithm for Poissons
//' @description IRLS to estimate network score of Poisson nodes.
//' @keywords internal
//' @returns a list
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
  arma::vec f(nobs);

  // Pre-calculate the log-factorial term using a stable function
  for (int j = 0; j < nobs; ++j) {
    f[j] = R::lgammafn(b[j] + 1.0); // Equivalent to lfactorial(b[j]) in R
  }

  // Initial values calculation to avoid division by zero
  arma::vec mu = b + 0.1;
  arma::vec eta_init = arma::log(mu);
  arma::vec W_init = arma::exp(eta_init);
  arma::vec z_init = eta_init + (b - W_init) / W_init;

  // Check if initial values are finite and stable
  if (!W_init.is_finite() || !z_init.is_finite()) {
    Rcpp::warning("Initial values for W or z are not finite. Starting with x = 0.");
    x.zeros(A.n_cols, 1);
  } else {
    arma::mat varmatrix_init = A.t() * (W_init % A.each_col());
    if (varmatrix_init.has_nan() || varmatrix_init.has_inf()) {
      Rcpp::warning("Initial varmatrix is unstable. Starting with x = 0.");
      x.zeros(A.n_cols, 1);
    } else {
      x = arma::solve(varmatrix_init, A.t() * (W_init % z_init), arma::solve_opts::likely_sympd);
      if (!x.is_finite()) {
        Rcpp::warning("Solving for initial x resulted in non-finite values. Starting with x = 0.");
        x.zeros(A.n_cols, 1);
      }
    }
  }

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
