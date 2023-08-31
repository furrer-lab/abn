				
void build_designmatrix_pois_rv(network *dag,datamatrix *obsdata, double priormean, double priorsd,const double priorgamshape, const double priorgamscale, 
			   datamatrix *designmatrix, int nodeid, int storeModes);
			   
int generate_pois_rv_inits(gsl_vector *myBeta,struct fnparams *gparams);

double g_pois_inner( gsl_vector *beta, const datamatrix *designdata, int groupid, double epsabs, int maxiters, int verbose);

int rv_dg_pois_inner (const gsl_vector *epsilonvec, void *params, gsl_vector *dgvalues);

int rv_hessg_pois_inner (const gsl_vector *epsilonvec, void *params,gsl_matrix *hessgvalues);

int wrapper_rv_fdf_pois_inner (const gsl_vector *beta, void *gparams, gsl_vector *dgvalues, gsl_matrix *hessgvalues);

int rv_g_pois_inner (const gsl_vector *epsilonvec, void *params, double *gvalue);

double g_outer_pois_single (double x, void *params);
			   
