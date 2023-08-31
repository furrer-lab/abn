				
void build_designmatrix_gaus_rv(network *dag,datamatrix *obsdata, double priormean, double priorsd,const double priorgamshape, const double priorgamscale, 
			   datamatrix *designmatrix, int nodeid, int storeModes);
			   
int generate_gaus_rv_inits(gsl_vector *myBeta,struct fnparams *gparams);

double g_inner_gaus( gsl_vector *beta, const datamatrix *designdata, int groupid, double epsabs, int maxiters, int verbose);

int rv_dg_inner_gaus (const gsl_vector *epsilonvec, void *params, gsl_vector *dgvalues);

int rv_hessg_inner_gaus (const gsl_vector *epsilonvec, void *params,gsl_matrix *hessgvalues);

int wrapper_rv_fdf_inner_gaus (const gsl_vector *beta, void *gparams,gsl_vector *dgvalues, gsl_matrix *hessgvalues);

int rv_g_inner_gaus (const gsl_vector *epsilonvec, void *params, double *gvalue);

double g_outer_gaus_single (double x, void *params);




