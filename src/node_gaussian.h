void calc_node_Score_gaus( network *dag,  datamatrix *obsdata, int nodeid, int verbose,
                                datamatrix *designmatrix, const double priormean, const double priorsd,const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs,int storeModes);
				
void build_designmatrix_gaus(network *dag,datamatrix *obsdata, double priormean, double priorsd,double priorgamshape, double priorgamscale,datamatrix *designmatrix, int nodeid,
			     int storeModes);
                                
int laplace_gaus_g (const gsl_vector *beta, void *params,double *gvalue);

int laplace_gaus_dg (const gsl_vector *beta, void *params, gsl_vector *dgvalues);

int laplace_gaus_hessg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues);

int wrapper_gaus_fdf (const gsl_vector *beta, void *gparams,gsl_vector *dgvalues, gsl_matrix *hessgvalues);

int generate_gaus_inits(gsl_vector *myBeta,struct fnparams *gparams, int errverbose);

void calc_gaussian_marginal(network *dag, datamatrix *obsdata, int nodeid,  int verbose,datamatrix *designmatrix, 
		                const double priormean, const double priorsd,const double priorgamshape,const double priorgamscale,
                                const int maxiters, const double epsabs, double *denom_modes,int paramid,
			        double betafixed, double mlik, double *posterior);
				
 int laplace_gaus_g_marg (const gsl_vector *beta, void *params,double *gvalue);

int laplace_gaus_dg_marg (const gsl_vector *beta, void *params, gsl_vector *dgvalues);

int laplace_gaus_hessg_marg (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues);

int wrapper_gaus_fdf_marg (const gsl_vector *beta, void *gparams, gsl_vector *dgvalues, gsl_matrix *hessgvalues);

int generate_gaus_inits_marg(gsl_vector *myBeta,struct fnparams *gparams); 
		     
