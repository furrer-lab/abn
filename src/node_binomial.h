void calc_node_Score_binary( network *dag,  datamatrix *obsdata, int nodeid, int verbose,
                                datamatrix *designmatrix, const double priormean, const double priorsd,
                                const int maxiters, const double epsabs, int storeModes);
				
void build_designmatrix(network *dag,datamatrix *obsdata, double priormean, double priorsd,datamatrix *designmatrix, int nodeid, int storeModes);

int laplace_g (const gsl_vector *beta, void *params,double *gvalue);

int laplace_dg (const gsl_vector *beta, void *params, gsl_vector *dgvalues);

int laplace_hessg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues);

int wrapper_fdf (const gsl_vector *beta, void *gparams,gsl_vector *dgvalues, gsl_matrix *hessgvalues);

int generate_inits_n(gsl_vector *myBeta,struct fnparams *gparams);

void calc_binary_marginal(network *dag, datamatrix *obsdata, int nodeid,  int verbose,
                                datamatrix *designmatrix, const double priormean, const double priorsd,
                                const int maxiters, const double epsabs, double *denom_modes,int paramid, double betafixed,
			        double mlik, double *posterior);
				
int laplace_g_marg (const gsl_vector *beta, void *params,double *gvalue);

int laplace_dg_marg (const gsl_vector *beta, void *params, gsl_vector *dgvalues);

int laplace_hessg_marg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues);

int wrapper_fdf_marg (const gsl_vector *beta, void *gparams, gsl_vector *dgvalues, gsl_matrix *hessgvalues);

		     
