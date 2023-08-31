#include <gsl/gsl_deriv.h>

void calc_node_Score_pois_rv_R( network *dag,  datamatrix *obsdata, int nodeid, int errverbose, int trace,
                                datamatrix *designmatrix, const double priormean, const double priorsd,const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs, int storeModes, double epsabs_inner, int maxiters_inner, double finitestepsize, int verbose,
				   double h_guess, double h_epsabs, int maxiters_hessian, int ModesONLY,
				   double max_hessian_error,double myfactor_brent, int maxiters_hessian_brent, double num_intervals_brent);
				   
double g_pois_outer_R (int Rn, double *betaincTauDBL, void *params);

void rv_dg_pois_outer_R (int n, double *betaDBL, double *dgvaluesDBL,void *params);

double compute_pois_mlik_nm(const gsl_vector *finitestepsize_vec, void *params);

int rv_hessg_pois_outer( gsl_vector* beta, void* params, gsl_matrix* hessgvalues,double h, gsl_matrix* hessgvalues3pt);

double compute_mlik_pois_brent(double finitestepsize, void *params);

double get_best_stepsize_pois(double delta,double lower,double upper,int maxiters_hessian, struct fnparams *gparams,
			 double (* compute_mlik_nm_brent) (double finitestepsize, void *params), gsl_min_fminimizer *s1, double *finitestepsize,double *saverror, int errverbose);
