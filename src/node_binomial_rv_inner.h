
void build_designmatrix_rv(network *dag,datamatrix *obsdata, double priormean, double priorsd,const double priorgamshape, const double priorgamscale,
			   datamatrix *designmatrix, int nodeid, int storeModes);

double g_inner( gsl_vector *beta, const datamatrix *designdata, int groupid,double epsabs,int maxiters, int verbose);

int rv_g_inner (const gsl_vector *epsilonvec, void *params, double *gvalue);

int rv_dg_inner (const gsl_vector *epsilon, void *params, gsl_vector *dgvalues);

int rv_hessg_inner (const gsl_vector *epsilon, void *params, gsl_matrix *hessgvalues);

int wrapper_rv_fdf_inner (const gsl_vector *beta, void *gparams,
                     gsl_vector *dgvalues, gsl_matrix *hessgvalues);

int generate_rv_inits(gsl_vector *myBeta,struct fnparams *gparams);

double g_outer_single (double x, void *params);

double rv_partial_deriv (const gsl_vector *beta, void *params, int index);

double rv_partial_deriv_wrapper(double x, void *params);
