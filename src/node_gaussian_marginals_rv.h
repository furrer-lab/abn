void calc_gaussian_marginal_rv_R(network *dag, datamatrix *obsdata, int nodeid,  int errverbose, int trace, 
                                datamatrix *designmatrix, const double priormean, const double priorsd, const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs,double epsabs_inner, int maxiters_inner, double finitestepsize, int verbose,
				 double h_guess, double h_epsabs, int maxiters_hessian,
			       double *denom_modes, int paramid, double betafixed, double mlik, double *posterior,
				 double max_hessian_error,double myfactor_brent, int maxiters_hessian_brent, double num_intervals_brent);
				
double g_outer_gaus_marg_R (int Rn, double *betaincTauDBL, void *params);

void rv_dg_outer_gaus_marg_R (int n, double *betaDBL, double *dgvaluesDBL,void *params);

double compute_mlik_gaus_marg_nm(const gsl_vector *finitestepsize_vec, void *params);

int rv_hessg_outer_gaus_marg( gsl_vector* beta, void* params, gsl_matrix* hessgvalues,double h, gsl_matrix* hessgvalues3pt);

double compute_mlik_gaus_marg_brent(double finitestepsize, void *params);

double get_best_stepsize_gaus_marg(double delta,double lower,double upper,int maxiters_hessian, struct fnparams *gparams,
			 double (* compute_mlik_nm_brent) (double finitestepsize, void *params), gsl_min_fminimizer *s1, double *finitestepsize,double *saverror);
			 
