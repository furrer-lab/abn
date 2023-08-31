void calc_parameter_marginal(network *dag,datamatrix *obsdata, datamatrix *designmatrix,
				const double priormean, const double priorsd, const double priorgamshape, const double priorgamscale,
			     const int maxiters, const double epsabs, int verbose, const int errverbose, int trace,
			      double *denom_modes, int childid, int paramid,  
			     double epsabs_inner, int maxiters_inner, double finitestepsize,
			     double h_guess, double h_epsabs,int maxiters_hessian,
			     double betafixed, double mlik, double *posterior,
			     double max_hessian_error,double myfactor_brent, int maxiters_hessian_brent, double num_intervals_brent);
			     
			       
			       
