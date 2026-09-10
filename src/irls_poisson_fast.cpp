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
Rcpp::List irls_poisson_cpp_fast(arma::mat A, arma::vec b, double maxit, double tol) {
  // Definitions
  arma::vec x;
  arma::vec xold;
  arma::mat varmatrix;

  double nobs = A.n_rows;
  double ll = -std::numeric_limits<double>::infinity();
  double aic, bic, df, mdl;

  arma::vec W(nobs);
  arma::vec eta(nobs);
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

  // Iterative Reweighted Least Squares
  for (int i = 0; i < maxit; ++i) {
    eta = A * x;
    W = exp(eta);

    // Check for non-finite weights which can cause division by zero
    if (!W.is_finite()) {
      Rcpp::warning("Non-finite weights encountered at iteration %d. Stopping.", i);
      break;
    }
    z = eta + (b - W) / W;
    xold = x;

    // Coefficients update
    double current_ll = arma::accu(b % eta - exp(eta) - f);
    if (!std::isfinite(current_ll)) {
      current_ll = -std::numeric_limits<double>::infinity();
    }

    // Step-Halving to avoid overly aggressive steps leading to instability
    varmatrix = A.t() * (W % A.each_col());
    double lambda = 1e-8;
    varmatrix.diag() += lambda;

    // Check if varmatrix is stable before solving
    if (varmatrix.has_nan() || varmatrix.has_inf()) {
      Rcpp::warning("varmatrix became unstable at iteration %d. Stopping.", i);
      break;
    }

    // Calculate the update direction
    arma::vec x_update_direction = arma::solve(varmatrix, A.t() * (W % z), arma::solve_opts::likely_sympd) - xold;

    double step = 1.0;
    int max_half_steps = 10;
    bool update_accepted = false;

    // Step-halving loop to find a valid update
    for (int j = 0; j < max_half_steps; ++j) {
      arma::vec x_new = xold + step * x_update_direction;

      if (!x_new.is_finite()) {
        step /= 2.0;
        continue; // Try a smaller step
      }

      // Calculate new eta and log-likelihood
      arma::vec eta_new = A * x_new;
      double new_ll = arma::accu(b % eta_new - exp(eta_new) - f);

      // Check if the new log-likelihood is finite and improves upon the current one
      if (std::isfinite(new_ll) && new_ll >= current_ll) {
        x = x_new; // Accept the update
        ll = new_ll;
        update_accepted = true;
        break;
      }

      // If not accepted, halve the step size and try again
      step /= 2.0; // Reduce step size
    }

    if (!update_accepted) {
      Rcpp::warning("Step halving failed to find a valid step at iteration %d.", i);
      x = xold; // Revert to old coefficients if no improvement found
    }
    // --- End of Step-Halving ---

    // Check convergence
    if (sqrt(pow(arma::norm(x - xold), 2)) < tol) {
      break;
    }
  }

  // Final calculations
  eta = A * x;
  df = A.n_cols;
  ll = arma::accu(b % eta - exp(eta) - f);
  aic = -2 * ll + 2 * df;
  bic = -2 * ll + log(nobs) * df;
  mdl = 1;
  // Calculate the variance-covariance matrix
  varmatrix = A.t() * (exp(eta) % A.each_col());
  if (varmatrix.has_nan() || varmatrix.has_inf()) {
    Rcpp::warning("Variance-covariance matrix is unstable. Returning NaN.");
    varmatrix.fill(NAN);
  }
  // Calculate the sum of squared errors
  arma::vec e = (b - eta);
  double ssr = arma::dot(e, e);

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
