/** **********************************************************************************************************************/
/** **********************************************************************************************************************/
/** README. */
/** The top level function in this file is calc_node_Score_binary_rv(), which computes the log marginal likelihood for a */
/** binary logit model with one level of random effects to adjust for within group correlation. Due to the numerics are  */
/** far more involved than for glm's - in particular the Laplace approx to the mlik involves summing over Laplace approxs*/
/** - one for each group of observations. Hence there is a single outer Laplace calc. and then many "inner" Laplace calcs*/
/** key differences with non-mixed models - the design matrix contains the local rv term (epsilon, say). There are a     */ 
/** number of different accuracy parameters to deal with e.g. step sizes for numerical derivs etc.                       */
/** init values - as for binary use glm least square ests but do NOT use Gaussian variance value but rather the use      */
/** p(1-p) where p is the intercept term in the mean e.g. p=)(exp(b0)/(1+exp(b0))                                        */

/** top level - calc_node_Score_binary_rv(): uses a minimiser (bfgs2) rather than root finder and 
                                  g_outer(): the objective function -1/n*log(f(D|theta)f(thera))
                              rv_dg_outer(): the derivative of the objective function - note this is a numeric derivative
                                             using central finite differences
                           rv_hessg_outer(): numerical derivative using central finite differences
                           
    second level -                g_outer(): this is basically just a wrapper which computes the prior for the means and
                                             calls a function which computes a Laplace approx for each separate data group 
                                  g_inner(): this does the actual Laplace approx calc for each group and uses root 
                                             finding with analytical hessian
    
    third level                   g_inner(): does root finding using next two functions and then computes the Laplace value 
                              rv_dg_inner(): first derviative 
                           rv_hessg_inner(): second derivatives
                               rv_g_inner(): value of the -1/n*log(f(D|theta)f(theta) for each group of data D            
    
    others
                           g_outer_single(): as g_outer() but where one variable is allowed to vary - used in num derivative
                           
                                                                                                                         */
/** **********************************************************************************************************************/
/** **********************************************************************************************************************/
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h" 
#include "node_binomial_rv_inner.h"
#include "node_binomial_rv.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_deriv.h>
#include "node_gaussian_rv.h"
#include "node_gaussian_rv_inner.h"
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#define NEWHESS
#define NOPRIOR1 /**turn on for removing effect of priors - this is purely to enable a check of loglike values with lme4 */
/** ****************************************************************************************************
 ***** calc an individual logistic regression model 
 *******************************************************************************************************/
void calc_node_Score_gaus_rv_R(network *dag, datamatrix *obsdata, int nodeid,  int errverbose, int trace,
                                datamatrix *designmatrix, const double priormean, const double priorsd, const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs,int storeModes, double epsabs_inner, int maxiters_inner, double finitestepsize, int verbose,
				 double h_guess, double h_epsabs, int maxiters_hessian, int ModesONLY,
			         double max_hessian_error,double myfactor_brent, int maxiters_hessian_brent, double num_intervals_brent)
{
#ifdef NOPRIOR
Rprintf("############ Warning - Priors turned off - use only for checking mlik value! ################\n");
#endif
  
  int i,status,sss,index=0,iter;
  /*int j;*/
  gsl_vector *myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,*localbeta,*localbeta2,/* *dgvalues,*/ *finitefactors,/* *factorindexes,*/ *finitestepsize_vec=0,*nmstepsize=0;/** this will hold the parameter point estimates + 1 is for precision of rv's */
  struct fnparams gparams;/** for passing to the gsl zero finding functions */
  double gvalue;double nm_size=0.0;
  gsl_matrix *mattmp2,*mattmp3,*mattmp4,*hessgvalues,*hessgvalues3pt;
  double mydet=0.0,logscore=0.0;/*,logscore3pt=0.0;*/
  gsl_permutation *initsperm;
  gsl_permutation *perm=0; 
  const gsl_multimin_fminimizer_type *T;
       gsl_multimin_fminimizer *s;
     gsl_multimin_function F; 
     
  double lower,upper,lower_f,upper_f;int found=0; double delta=0.0,new_f_min=0.0, increLogscale=0.0, best_Error=0.0,best_h=0.0;
  const gsl_min_fminimizer_type *T1;
  gsl_min_fminimizer *s1;  
  int n,m;
  int nDim;/** dimension of optim problem */
  int *nbd;/** nbd is an integer array of dimension nDim.
	                                      On entry nbd represents the type of bounds imposed on the variables, and must be specified as follows:
	                                      nbd(i)=0 if x(i) is unbounded,
		                              1 if x(i) has only a lower bound,
		                              2 if x(i) has both lower and upper bounds, and
		                              3 if x(i) has only an upper bound.
	                                      On exit nbd is unchanged.*/
  
  double *lowerbounds,*upperbounds; /* h_gvalue;*//*,lowestHesserror,beststepsize;*/
  int failcode;/** check code see R ?optim - if non-zero then a problem **/
  double factr=1e-07;/** error size scaler - this is the default value*/
  double pgtol=1e-07;/** default value is zero - this is the gradient tolerance - mmm what does that actually mean? */
  int fncount,grcount;/** hold number of evaluations */
  char msg[60];/** error message */
  /* int trace=errverbose;  */ /** like verbose */
  int nREPORT=1000;/** report freq*/
  int lmm=5;/** see R ?optim - number of function evals to store - default is 5 */
  /** want to find the modes of the function g(betas) where betas=b_0,b_1,,...,tau, the latter being precision */
  /** g(betas) is not differentiable as it contains integrals (which themselves need Laplace estimates **/
  
    
  /** SETUP things which are the same across all data chunks - groups  */
  /** build design matrix which is designmatrix->datamatrix=X, designmatrix->Y=Y plus priors designmatrix->priorsd, designmatrix->priormean **/
  /** NOTE: design matrix here does include the random effect term **/
  /** note - numparams does NOT include either precision term - numpars +1 */
  
  build_designmatrix_gaus_rv(dag,obsdata,priormean, priorsd,priorgamshape,priorgamscale,designmatrix,nodeid,storeModes);
  
  
  nDim=designmatrix->numparams+1+1;/** +1 for group precision +1 for residual precision */
  lowerbounds=(double *)R_alloc(nDim,sizeof(double*));
  upperbounds=(double *)R_alloc(nDim,sizeof(double*));
  nbd=(int *)R_alloc(nDim,sizeof(int*));
  for(i=0;i<nDim-1;i++){lowerbounds[i]=-DBL_MAX;
                        upperbounds[i]=DBL_MAX;
			nbd[i]=0;}
			nbd[nDim-2]=1;lowerbounds[nDim-2]=0.001;/** lower bound for group precision */
			nbd[nDim-1]=1;lowerbounds[nDim-1]=0.001;/** lower bound for residual precision */
  finitefactors = gsl_vector_alloc(7);/** used to change stepsize in hessian estimate **/			
  gsl_vector_set(finitefactors,0,1.0E-03);gsl_vector_set(finitefactors,1,1.0E-02);gsl_vector_set(finitefactors,2,1.0E-01);
  gsl_vector_set(finitefactors,3,1.0);gsl_vector_set(finitefactors,4,1.0E+01);gsl_vector_set(finitefactors,5,1.0E+02);
  gsl_vector_set(finitefactors,6,1.0E+03);
  
  /*factorindexes = gsl_vector_alloc(7);*//** used to change stepsize in hessian estimate **/			
  /*for(i=0;i<7;i++){gsl_vector_set(factorindexes,i,i);}*/
  
  /** change finite.step.size by 0.1,1, and 10 factors respectively **/
  
  /*dgvalues = gsl_vector_alloc (designmatrix->numparams+1+1);*//** inc rv precision & inc residual precision */
  
  myBeta = gsl_vector_alloc (designmatrix->numparams+1+1);/** inc rv precision & inc residual precision */
  vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
  vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
  mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
  mattmp3 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
  mattmp4 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
  initsperm = gsl_permutation_alloc (designmatrix->numparams);/** for use with initial guesses */
  vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
  vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
  localbeta = gsl_vector_alloc (designmatrix->numparams);/** scratch space in later functions - excl. all precision **/
  localbeta2 = gsl_vector_alloc (designmatrix->numparams+1+1);/** scratch space in later functions - inc. both group and residual precision **/
  
  hessgvalues = gsl_matrix_alloc (designmatrix->numparams+1+1,designmatrix->numparams+1+1);
  hessgvalues3pt = gsl_matrix_alloc (designmatrix->numparams+1+1,designmatrix->numparams+1+1);
  
  gparams.designdata=designmatrix;
  
   gparams.vectmp1=vectmp1;
   gparams.vectmp2=vectmp2;
   gparams.mattmp2=mattmp2;
   gparams.mattmp3=mattmp3;
   gparams.mattmp4=mattmp4;
   gparams.perm=initsperm;
   gparams.vectmp1long=vectmp1long;
   gparams.vectmp2long=vectmp2long;
   gparams.beta=localbeta;
   gparams.betaincTau=localbeta2;
   gparams.epsabs_inner=epsabs_inner;
   gparams.maxiters_inner=maxiters_inner;
   gparams.verbose=verbose;
   gparams.finitestepsize=finitestepsize;
   
   dag->nodeScoresErrCode[nodeid]=0;/** reset error code to no error **/
   
   /*status=GSL_SUCCESS;*/
   generate_gaus_rv_inits(myBeta,&gparams);/** get initial values, n.b. classical variance est split equally between group and residual variance **/
 
   
  for(i=0;i<finitefactors->size;i++){/** go forwards through the factors so start with SMALLEST STEPSIZE */
   /*Rprintf("step size iteration %d\n",i);*/
     failcode=0;/** reset*/
    gparams.finitestepsize=gsl_vector_get(finitefactors,i)*finitestepsize;
  
     lbfgsb(nDim, lmm, myBeta->data, lowerbounds, upperbounds, nbd, &gvalue, &g_outer_gaus_R,
                      &rv_dg_outer_gaus_R, &failcode, 
	              &gparams,
	              factr,
                      pgtol, &fncount, &grcount,
                      maxiters, msg, trace, nREPORT);
		      
    if(!failcode){dag->nodeScoresErrCode[nodeid]=0;/*bestsize=gparams.finitestepsize;*/break;}/** break out of for loop if no error as we are done **/	     
   
   } /** end of for loop so now have mode estimates */
 
   
   /*g_outer_gaus_R(nDim,myBeta->data,&gparams);*/
   
   /*Rprintf("loglike value=%f\n",gvalue);*/
   /*for(i=0;i<4;i++){Rprintf("best est.=%e ",gsl_vector_get(myBeta,i));}Rprintf("\n");*/
   
   if(failcode){Rprintf("%s at node %d\n",msg,nodeid+1);/** notify if there is an error and set final error code **/
		     dag->nodeScoresErrCode[nodeid]=1;
   } 
     
    gparams.finitestepsize=finitestepsize;/** reset */
    if(storeModes){/** keep a copy of the parameter modes found for use later in other function calls etc**/
	 index=0;    /*Rprintf("size of beta=%d %f %f\n",myBeta->size, gsl_vector_get(myBeta,0),gsl_vector_get(myBeta,1));*/
		     for(i=0;i<dag->numNodes+3;i++){/** roll myBeta into dag->modes into the appropriate columns**/
		       if(gsl_matrix_get(dag->modes,nodeid,i)!=DBL_MAX){
			 gsl_matrix_set(dag->modes,nodeid,i,gsl_vector_get(myBeta,index++));}} 
                   /*for(i=0;i<dag->numNodes+3;i++){Rprintf("%e ",gsl_matrix_get(dag->modes,nodeid,i));}Rprintf("\n");*/
		   
		   }  
		   
  if(!ModesONLY){/** only want modes so can skip the rest **/
     
   /** now compute the hessian at the step size with lowest error **/
   /*Rprintf("starting hessian estimation\n");*/
   n=obsdata->numDataPts;
   m=designmatrix->numparams+1+1;/** inc group precision and residual precision*/
   perm = gsl_permutation_alloc (m);
 
   /** just re-use as much of existing gparams as possible - so names are meaningless e.g. betafixed is actually gvalue */
   gparams.betaincTau=myBeta;
   gparams.nDim=n;
   gparams.mDim=m;
   gparams.perm=perm;
   gparams.mattmp2=hessgvalues;
   gparams.mattmp3=hessgvalues3pt;
   gparams.gvalue=gvalue;
   
   /** optimisation is to find the best step_size h, e.g. stepsize with the smallest local error between 5pt and 3pt rules */
    F.f = &compute_mlik_gaus_nm;
    F.params = &gparams;
    F.n = 1;
   
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 1);
   
    finitestepsize_vec = gsl_vector_alloc (1);
    gsl_vector_set (finitestepsize_vec, 0, h_guess);
    nmstepsize = gsl_vector_alloc (1);
    gsl_vector_set_all (nmstepsize, h_guess); 
    gsl_multimin_fminimizer_set (s, &F, finitestepsize_vec, nmstepsize);
    status = GSL_SUCCESS;
    
     iter=0;
   
    do
         {
           iter++;
           status = gsl_multimin_fminimizer_iterate (s);
     
           if (status) 
             break;
	   
	   nm_size = gsl_multimin_fminimizer_size (s);
           status = gsl_multimin_test_size (nm_size, h_epsabs);
     /*
           if (status == GSL_SUCCESS)
             {
               Rprintf ("converged to minimum at\n");
             }
     */
           /*Rprintf ("iter=%5d error in mlik=%3.5e using fin.diff step= %3.2e nmsize=%3.2e\n", iter,s->fval,gsl_vector_get (s->x, 0),nm_size);*/
    
         }
       while (status == GSL_CONTINUE && iter < maxiters_hessian);
       if( (status != GSL_SUCCESS)){/*actual_status=status;*//** copy for use later **/
                                    status=GSL_FAILURE;} /** solution failed to achieve a value below h_epsabs **/                                                               
	 
    finitestepsize=gsl_vector_get(s->x,0);/** get best fin.diff stepsize **/
    
    dag->hessianError[nodeid]= s->fval;/** get fin.diff error **/
    
    gsl_multimin_fminimizer_free (s);
    
    /*if(dag->hessianError[nodeid]==DBL_MAX){Rprintf("############################## got DBL_MAX\n");}*/
    
    if(dag->hessianError[nodeid]!=DBL_MAX && dag->hessianError[nodeid]>max_hessian_error){Rprintf("Error in mlik = %e > tolerance of %e so continuing optimisation using Brent initial guess h=%e\n",
                                                   dag->hessianError[nodeid],max_hessian_error,finitestepsize); 
     /* Rprintf("stepsize after NM= %e\n",finitestepsize);*/
  
     T1 = gsl_min_fminimizer_brent;
     s1 = gsl_min_fminimizer_alloc (T1);
	 
      /** must find lower and upper such that f(lower)<f(finitestepsize)<f(upper) **/ 
      /** use an interval of lower=finitestepsize/FACTOR, upper=finitestepsize*FACTOR and then start at the lower end and travel up until
      find suitable endpoints - seems to work but not exactly fast!**/
      best_Error=dag->hessianError[nodeid];/** original error from nelder */
      best_h=finitestepsize;               /** original stepsize from nelder */
      found=0;/** flag for found good result */
      lower=finitestepsize/myfactor_brent;
      upper=myfactor_brent*finitestepsize;
      lower_f=compute_mlik_gaus_brent(lower, &gparams);/** value at lower point **/
      upper_f=compute_mlik_gaus_brent(upper, &gparams);/** value at higher point **/
      increLogscale=(gsl_sf_log(upper)-gsl_sf_log(lower))/num_intervals_brent;/** on a log scale */
      for(delta=gsl_sf_log(lower)+increLogscale;delta<gsl_sf_log(upper);delta+=increLogscale){/** linear increments on a log scale **/
	R_CheckUserInterrupt();/** allow an interupt from R console */ 
	/** find a point which has f(x) lower than f(lower) and f(upper) **/
	 new_f_min=compute_mlik_gaus_brent(gsl_sf_exp(delta), &gparams); 
	 if (errverbose) 	 Rprintf("GN: lower=%e, delta=%e, upper=%e\n",lower,gsl_sf_exp(delta),upper);
        if(lower_f>new_f_min && new_f_min<upper_f  && get_best_stepsize_gaus(gsl_sf_exp(delta),lower,upper,maxiters_hessian_brent,&gparams, &compute_mlik_gaus_brent,
									     s1,&finitestepsize,&(dag->hessianError[nodeid]), errverbose )<=max_hessian_error){/** have an interval suitable for bracketing **/
	                                                           /** above is address so can store error withouth rerunning function */
	  /*finitestepsize=delta;*/
	  found=1;
	  status=GSL_SUCCESS;
	  break;/** break out of delta - so have found new x_min **/
	} else {/** have not got a good enough error but save the best error and stepsize so far found **/
	        if(dag->hessianError[nodeid]<best_Error){best_Error=dag->hessianError[nodeid];
	                                                best_h=finitestepsize;}
	        }
      } /** end of search for interval and good error **/
         
      if(!found){/** have not found a suitably small error but may have found a better error than nelder mead **/
        
       /** best_Error will either be the original nelder mean value or better, and best_h is the corresponding stepsize**/
	                                 dag->hessianError[nodeid]=best_Error;
					 finitestepsize=best_h;
        /** reset back to nelder-mead estimate **/
	status=GSL_FAILURE;/** set to failure since we did not achieve the lower error asked for */
      Rprintf("failed to meet tolerance of %e and using best error estimate found of %e\n",max_hessian_error,dag->hessianError[nodeid]);}

    gsl_min_fminimizer_free (s1);
   
   } /** end of error being too large **/
 
   if(dag->hessianError[nodeid]==DBL_MAX){/** in this case nelder mead could not estimate the hessian error so abort as something is probably
                                               very wrong here */
                                          error("");}/** use the R tryCatch rather than the switch for status below **/
      /*Rprintf("GVALUE=%e nodeid=%d\n",gvalue,nodeid+1);*/
       /*finitestepsize=1E-04;
       h_guess=finitestepsize;
*/
       switch(status){  /** choose which type of node we have */
                     case GSL_SUCCESS:{    
		                     /** successful finite diff so now do final computation with the optimal step size **/
                                     /*Rprintf("search for optimal step size : status = %s at nodeid %d\n", gsl_strerror (status),nodeid+1);*/
                                     rv_hessg_outer_gaus(myBeta,&gparams, hessgvalues,finitestepsize,hessgvalues3pt);/** EDIT BACK to "finitestepsize" start with LARGEST STEPSIZE **/
				    /* Rprintf("HESSIAN using stepsize=%e\n",finitestepsize);
				     for(i1=0;i1<hessgvalues->size1;i1++){
				        for(i2=0;i2<hessgvalues->size2;i2++){Rprintf("%e ",gsl_matrix_get(hessgvalues,i1,i2));}Rprintf("\n");}   */
                                     status=gsl_linalg_LU_decomp(hessgvalues,perm,&sss);
                                     mydet=gsl_linalg_LU_lndet(hessgvalues);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
                                     logscore= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
				     if(gsl_isnan(logscore)){logscore= R_NaN;dag->nodeScoresErrCode[nodeid]=2;}
				     dag->nodeScores[nodeid]=logscore;
				       
		                      break;  
		     }
       
		     case GSL_FAILURE: {/** the minimiser did not find a minimum meeting the accuracy requirements and so may be unreliable **/
		                       Rprintf ("-- ERROR! -- search for optimal step size error: status = %s at nodeid %d\n", gsl_strerror (status),nodeid+1);
                                       rv_hessg_outer_gaus(myBeta,&gparams, hessgvalues,finitestepsize,hessgvalues3pt);/** start with LARGEST STEPSIZE **/
                                       /*Rprintf("HESSIAN using stepsize=%e\n",finitestepsize);
				       for(i1=0;i1<hessgvalues->size1;i1++){
				        for(i2=0;i2<hessgvalues->size2;i2++){Rprintf("%e ",gsl_matrix_get(hessgvalues,i1,i2));}Rprintf("\n");} */
				        
				       status=gsl_linalg_LU_decomp(hessgvalues,perm,&sss);
                                       mydet=gsl_linalg_LU_lndet(hessgvalues);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
                                       logscore= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
				       dag->nodeScoresErrCode[nodeid]=4;
				       if(gsl_isnan(logscore)){logscore= R_NaN;dag->nodeScoresErrCode[nodeid]=2;}
				       dag->nodeScores[nodeid]=logscore;
				       
		                       
				       break; 
		     }
		     
		     default:{Rprintf("got case %s\n",gsl_strerror (status)); error("in default switch in calc_node_Score_binary_rv_R() - should never get here!");}  
		     
          }
          
        
   /** try the bounded search for h stepsize rather than one-dim min which needs bound specified **/     
   } /** end of ModesONLY **/     
  
   /** now free up allocated memory **/
   for(i=0;i<designmatrix->numUnqGrps;i++){gsl_matrix_free(designmatrix->array_of_designs[i]);
                                           gsl_vector_free(designmatrix->array_of_Y[i]);}
   gsl_vector_free(designmatrix->priormean);
   gsl_vector_free(designmatrix->priorsd);
   gsl_vector_free(designmatrix->priorgamshape);
   gsl_vector_free(designmatrix->priorgamscale);
   gsl_vector_free(designmatrix->Y);
   gsl_matrix_free(designmatrix->datamatrix_noRV);
   /*gsl_vector_free(dgvalues);*/
   gsl_vector_free(myBeta); 
   gsl_vector_free(vectmp1);
   gsl_vector_free(vectmp2);
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(mattmp3);
   gsl_matrix_free(mattmp4);
   gsl_permutation_free(initsperm);
   gsl_vector_free(vectmp1long);
   gsl_vector_free(vectmp2long);
   gsl_vector_free(localbeta);
   gsl_vector_free(localbeta2);
   gsl_matrix_free(hessgvalues);
   gsl_matrix_free(hessgvalues3pt);
   gsl_vector_free(finitefactors);
   /*gsl_vector_free(factorindexes);*/
   
   if(!ModesONLY){/** didn't allocate these so don't unallocate! */
    gsl_permutation_free(perm);
    gsl_vector_free(finitestepsize_vec);
    gsl_vector_free(nmstepsize);}
   
   /*if(!failcode){*//*}*/
   
   /*dag->nodeScores[nodeid]=logscore;*/

  /*error("end of calc_node_Score_gaus_rv_R\n");*/

}

/** *************************************************************************************************************/
/** **************************************************************************************************************/
/** g_outer(y) = -(1/n)* log( prod[g_inner()].prior )                                                            */ 
/** this function is the product (or sum of logs) of the likelihoods for each individual groups of data          */
/** where these individual likelihood have integrals in them - to integrat out the latent rv effect. g_outer     */
/** then calls g_inner() which then deals with each of the individual likelihoods in turn                        */
/** **************************************************************************************************************/
/** **************************************************************************************************************/
double g_outer_gaus_R (int Rn, double *betaincTauDBL, void *params) /*typedef double optimfn(int n, double *par, void *ex);*/
{
  int i,j;
  double term1=0.0,singlegrp=0.0;
  const datamatrix *designdata = ((struct fnparams *) params)->designdata;/** all design data inc Y and priors **/

const gsl_vector *priormean = designdata->priormean;
  const gsl_vector *priorsd   = designdata->priorsd;
  const gsl_vector *priorgamshape   = designdata->priorgamshape;
  const gsl_vector *priorgamscale   = designdata->priorgamscale;
   gsl_vector *beta   = ((struct fnparams *) params)->beta;/** does not include precision */
  gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
  gsl_vector *vectmp2 =((struct fnparams *) params)->vectmp2;/** numparams long*/
  gsl_vector *betaincTau=((struct fnparams *) params)->betaincTau;/** to copy betaincTauDBL into **/
  double epsabs_inner=((struct fnparams *) params)->epsabs_inner;/** absolute error in internal laplace est */
  int maxiters_inner=((struct fnparams *) params)->maxiters_inner;/** number of steps for inner root finder */
  int verbose=((struct fnparams *) params)->verbose;/**  */
  
  int n_betas= (designdata->datamatrix_noRV)->size2;/** number of mean terms excl rv and precision **/
  int n=(designdata->datamatrix_noRV)->size1;/** total number of obs **/
  
  double term2=0.0,term3=0.0,term4=0.0,gval=0.0, term5=0.0;
  /*Rprintf("%d %d\n",n_betas,betaincTau->size);*/
  double tau_rv,tau_resid;
  
  for(i=0;i<betaincTau->size;i++){gsl_vector_set(betaincTau,i,betaincTauDBL[i]);} /** copy R double array into gsl vect inc. group and residual precision **/
  
  /*Rprintf("got = %f %f %f\n",gsl_vector_get(betaincTau,0),gsl_vector_get(betaincTau,1),gsl_vector_get(betaincTau,2));*/
    
  tau_rv=gsl_vector_get(betaincTau,betaincTau->size-2);/** extract the tau-precision from *beta - last entry */
  /*Rprintf("g_outer_rv tau=%f\n",tau_rv);*/
  tau_resid=gsl_vector_get(betaincTau,betaincTau->size-1);/** extract the tau-precision from *beta - last entry */
  /*Rprintf("g_outer_resid tau=%f\n",tau_resid);*/
  
  if(tau_rv<0.0){Rprintf("tau_rv negative =%e in g_outer_gaus!\n",tau_rv);error("");}
  if(tau_resid<0.0){Rprintf("tau_resid negative=%e in g_outer_gaus!\n",tau_resid);error("");}
  
  /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last two entries are group precision then residual precision */
       for(i=0;i<n_betas;i++){gsl_vector_set(beta,i,gsl_vector_get(betaincTau,i));/*Rprintf("passed beta=%f\n",gsl_vector_get(beta,i));*/
       }
     
  /** part 1 - the integrals over each group of observations - use laplace for this and that is taken care of in g_inner */ 
  /** first we want to evaluate each of the integrals for each data group **/ 
       for(j=0;j<designdata->numUnqGrps;j++){/** for each data group **/
	/*j=0;*/
	/*Rprintf("processing group %d\n",j+1);*/
	 /*Rprintf("tau in loop=%f\n",gsl_vector_get(betaincTau,n_betas));*/
	  singlegrp=g_inner_gaus(betaincTau,designdata,j, epsabs_inner,maxiters_inner,verbose);
	  /*Rprintf("loglike for group data = %e\n",singlegrp);*/
	 if(gsl_isnan(singlegrp)){error("nan in g_outer_gaus_R\n");}
	  term1+= singlegrp;
      }
    
    /** NOTE: uncomment next line as useful for debugging as this should be the same as logLik value from lme4 */
   /** Rprintf("total loglike=%e\n",term1);**/

  /** part 2 the priors for the means **/
  term2=0; for(i=0;i<n_betas;i++){term2+=-log(sqrt(2.0*M_PI)*gsl_vector_get(priorsd,i));}
  /** Calc this in parts: R code "term3<- sum( (-1/(2*sd.loc*sd.loc))*(mybeta-mean.loc)*(mybeta-mean.loc) );" **/
  gsl_vector_memcpy(vectmp1,beta);/** copy beta to temp vec */
  gsl_vector_memcpy(vectmp2,priormean);
  gsl_vector_scale(vectmp2,-1.0);
  gsl_vector_add(vectmp1,vectmp2);/** vectmp1= beta-mean**/
  gsl_vector_memcpy(vectmp2,vectmp1);/** copy vectmp1 to vectmp2 **/
  gsl_vector_mul(vectmp2,vectmp1);/** square all elements in vectmp1 and store in vectmp2 */
  gsl_vector_memcpy(vectmp1,priorsd);
  gsl_vector_mul(vectmp1,priorsd);/** square all elements in priorsd and store in vectmp1 */
  gsl_vector_div(vectmp2,vectmp1);/** vectmp2/vectmp1 and store in vectmp2 **/
  gsl_vector_scale(vectmp2,-0.5); /** scale by -1/2 */
  gsl_vector_set_all(vectmp1,1.0); /** ones vector */
  gsl_blas_ddot (vectmp2, vectmp1, &term3);/** DOT product simply to calcu sum value */
  
  
  /** part 3 the prior for the group precision tau_rv **/
  term4=  -gsl_vector_get(priorgamshape,0)*log(gsl_vector_get(priorgamscale,0))
             -gsl_sf_lngamma(gsl_vector_get(priorgamshape,0)) 
	     +(gsl_vector_get(priorgamshape,0)-1)*log(tau_rv)
	     -(tau_rv/gsl_vector_get(priorgamscale,0));
   
  /** part 4 the prior for the residual precision tau_resid **/
  term5=  -gsl_vector_get(priorgamshape,0)*log(gsl_vector_get(priorgamscale,0))
             -gsl_sf_lngamma(gsl_vector_get(priorgamshape,0)) 
	     +(gsl_vector_get(priorgamshape,0)-1)*log(tau_resid)
	     -(tau_resid/gsl_vector_get(priorgamscale,0));
   
	     
  gval=(-1.0/n)*(term1+term2+term3+term4+term5);
  
   /** NO PRIOR */
  /* Rprintf("WARNING - NO PRIOR\n");*/
  #ifdef NOPRIOR
  gval=(-1.0/n)*(term1);
  #endif
   if(gsl_isnan(gval)){error("g_outer_gaus_R\n");}

/*Rprintf("g_outer_final=%f term1=%f term2=%f term3=%f term4=%f term5=%f total=%f %d\n",gval,term1,term2,term3,term4,term5,term1+term2+term3+term4,n);*/	
	return(gval);/** negative since its a minimiser */
}
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/ 
void rv_dg_outer_gaus_R (int n, double *betaDBL, double *dgvaluesDBL,void *params)
{
  struct fnparams *gparams = ((struct fnparams *) params);
  
  /** fixed covariate and precision terms **/
  gsl_function F;int i;int haveTau=0;
  double result, abserr;
  double h=((struct fnparams *) gparams)->finitestepsize;
  /*double h_adj=0.0;*//** a new h is existing h is too small **/
  gsl_vector *betaincTau=((struct fnparams *) gparams)->betaincTau;/** scratch space to copy betaincTauDBL into **/
 
  for(i=0;i<betaincTau->size;i++){gsl_vector_set(betaincTau,i,betaDBL[i]);} /** copy into gsl_vector **/
   
  if(betaDBL[n-1]<0.0){error("negative tau_rv=%e in rv_dg_outer_gaus_R\n",betaDBL[n-1]);}
  if(betaDBL[n-2]<0.0){error("negative tau_resid=%e in rv_dg_outer_gaus_R\n",betaDBL[n-2]);}
  /** if get tau which is negative */
  
       F.function = &g_outer_gaus_single;
       F.params = gparams;
   
  for(i=0;i<n;i++){ /** each of the partial derivatives for beta terms, precision_rv and precision_resid */    
  gparams->fixed_index=i;
  if(i== n-1 || i== n-2){haveTau=1;} else {haveTau=0;}
    
  /** readme - evaluating f() at a negative value e.g. tau-h **/
  if(!haveTau){gsl_deriv_central (&F, betaDBL[i], h, &result, &abserr);/*Rprintf("fixed=%d val=%f\n",i,result);*/
  } else { /** first try central and if this goes into negative tau area then revert to forwards **/
           gsl_deriv_central(&F, betaDBL[i], h, &result, &abserr); 
	   if(gsl_isnan(abserr)){gsl_deriv_forward(&F, betaDBL[i], h, &result, &abserr);}
  }
  
  dgvaluesDBL[i]=result;
  }
  
/* Rprintf("dg_outer\n"); 
 for(i=0;i<n;i++){Rprintf("%f ",dgvaluesDBL[i]);}Rprintf("\n");
  */
}
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
int rv_hessg_outer_gaus( gsl_vector* beta, void* params, gsl_matrix* hessgvalues,double h, gsl_matrix* hessgvalues3pt)
{
  
  struct fnparams *gparams = ((struct fnparams *) params);
  int i,j;
  int haveTau=0;
  gsl_function F;
  gparams->betaincTau=beta;  
  
   
  /** if get tau which is negative */
  F.function = &g_outer_gaus_single;
  F.params = gparams;
  
 if(gsl_vector_get(beta,beta->size-1)<0.0){Rprintf("negative tau_rv in hess %e\n",gsl_vector_get(beta,beta->size-1));error("");}
 if(gsl_vector_get(beta,beta->size-2)<0.0){Rprintf("negative tau_resid in hess %e\n",gsl_vector_get(beta,beta->size-2));error("");}
 
  /** diagnonal terms d^2f/dx^2 - finite differences - difference between two */
  for(i=0;i<hessgvalues->size1;i++){
   for(j=0;j<hessgvalues->size2;j++){
        if(i>=j){/** hessian is symmetric **/
	        /** d^2f/ dbeta[i] dbeta[j] */
	 gparams->fixed_index=i;/** d^2f/dbeta1 dbeta0 = central diff of beta1 when beta0=beta0+h/2 and similar when beta0=beta0-h/2  -rest are fixed */
         if(i== hessgvalues->size1-1 || i== hessgvalues->size1-2){haveTau=1;} else {haveTau=0;}
   
         
	 gsl_matrix_set(hessgvalues,i,j,get_second_deriv_5pt(gparams,i,j,h,haveTau,&F));/** last args: havetau=0 **/
	 gsl_matrix_set(hessgvalues3pt,i,j,get_second_deriv_3pt(gparams,i,j,h,haveTau,&F));/** last args: havetau=0 **/
	} /** end i>j **/
	
   }
  }
  
  /** finally fill out matrix on upper diagnonal */
  for(i=0;i<hessgvalues->size1;i++){
   for(j=0;j<hessgvalues->size2;j++){
     if(i>=j){gsl_matrix_set(hessgvalues,j,i,gsl_matrix_get(hessgvalues,i,j));}}}
 
 /** finally fill out matrix on upper diagnonal */
  for(i=0;i<hessgvalues3pt->size1;i++){
   for(j=0;j<hessgvalues3pt->size2;j++){
     if(i>=j){gsl_matrix_set(hessgvalues3pt,j,i,gsl_matrix_get(hessgvalues3pt,i,j));}}}
 
  return GSL_SUCCESS;
  
}
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
double compute_mlik_gaus_nm(const gsl_vector *finitestepsize_vec, void *params)
{
   struct fnparams *gparams = ((struct fnparams *) params);
   gsl_vector *myBeta=gparams->betaincTau;
   int n=gparams->nDim;
   int m=gparams->mDim;
   gsl_permutation *perm=gparams->perm;
   gsl_matrix *hessgvalues=gparams->mattmp2;
   gsl_matrix *hessgvalues3pt=gparams->mattmp3;
   double gvalue=gparams->gvalue;
   int sss;
   double mydet;
   double logscore,logscore3pt;
   double finitestepsize=gsl_vector_get(finitestepsize_vec,0);
   double error_val=0.0;
   /*int i,j;
   
   Rprintf("got h=%e n=%d m=%d gvalue=%e\n",finitestepsize,n,m,gvalue);
  for(i=0;i<myBeta->size;i++){Rprintf("beta= %f ",gsl_vector_get(myBeta,i));}Rprintf("\n");
   */
   
   rv_hessg_outer_gaus(myBeta,gparams, hessgvalues,finitestepsize,hessgvalues3pt);/** start with LARGEST STEPSIZE **/
   
   /*for(i=0;i<hessgvalues3pt->size1;i++){for(j=0;j<hessgvalues3pt->size2;j++){Rprintf("%e ",gsl_matrix_get(hessgvalues3pt,i,j));}Rprintf("\n");}*/
   
   /*status=*/gsl_linalg_LU_decomp(hessgvalues,perm,&sss);
   mydet=gsl_linalg_LU_lndet(hessgvalues);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
   logscore= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
  
   /*status=*/gsl_linalg_LU_decomp(hessgvalues3pt,perm,&sss);
   mydet=gsl_linalg_LU_lndet(hessgvalues3pt);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
   logscore3pt= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
   
   error_val=fabs(logscore-logscore3pt);
   
   /*gparams->logscore=logscore;
   gparams->logscore3pt=logscore3pt;*/
   /*Rprintf("logscore=%e logscore3pt=%e\n",logscore,logscore3pt);*/
   
   if(gsl_isnan(error_val) || gsl_isinf(error_val)){return(DBL_MAX);/*error("Non-finite value in mlik error estimation");*/}
   return(error_val);
 
}	
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
/** *************************************************************************************************************************/
double compute_mlik_gaus_brent(double finitestepsize, void *params)
{
   struct fnparams *gparams = ((struct fnparams *) params);
   gsl_vector *myBeta=gparams->betaincTau;
   int n=gparams->nDim;
   int m=gparams->mDim;
   gsl_permutation *perm=gparams->perm;
   gsl_matrix *hessgvalues=gparams->mattmp2;
   gsl_matrix *hessgvalues3pt=gparams->mattmp3;
   double gvalue=gparams->gvalue;
   int sss;
   double mydet;
   double logscore,logscore3pt;
   double error_val=0.0;
   /*double finitestepsize=gsl_vector_get(finitestepsize_vec,0);*/
   
  /* int i,j;
   
   Rprintf("got h=%e n=%d m=%d gvalue=%e\n",finitestepsize,n,m,gvalue);
  for(i=0;i<myBeta->size;i++){Rprintf("beta= %f ",gsl_vector_get(myBeta,i));}Rprintf("\n");
  */
   
   rv_hessg_outer_gaus(myBeta,gparams, hessgvalues,finitestepsize,hessgvalues3pt);/** start with LARGEST STEPSIZE **/
   
   /*for(i=0;i<hessgvalues3pt->size1;i++){for(j=0;j<hessgvalues3pt->size2;j++){Rprintf("%e ",gsl_matrix_get(hessgvalues3pt,i,j));}Rprintf("\n");}*/
   
   /*status=*/gsl_linalg_LU_decomp(hessgvalues,perm,&sss);
   mydet=gsl_linalg_LU_lndet(hessgvalues);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
   logscore= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
  
   /*status=*/gsl_linalg_LU_decomp(hessgvalues3pt,perm,&sss);
   mydet=gsl_linalg_LU_lndet(hessgvalues3pt);/** compute determinant but this might be a nan - overflow? gsl_linalg_LU_lndet*/
   logscore3pt= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/** this is the final value */
 
   error_val=fabs(logscore-logscore3pt);
   /* Rprintf("error_val=%e\n",error_val);*/
   /*gparams->logscore=logscore;
   gparams->logscore3pt=logscore3pt;*/
   /*Rprintf("logscore=%e logscore3pt=%e\n",logscore,logscore3pt);*/
   
   if(gsl_isnan(error_val) || gsl_isinf(error_val)){return(DBL_MAX);/*error("Non-finite value in mlik error estimation");*/}
   return(error_val); 
 
}
/** ***********************************************************************************************/
/** ***********************************************************************************************/
double get_best_stepsize_gaus(double delta,double lower,double upper,int maxiters_hessian, struct fnparams *gparams,
			      double (* compute_mlik_nm_brent) (double finitestepsize, void *params), gsl_min_fminimizer *s1, double *finitestepsize,double *saverror, int errverbose)
{
  gsl_function F1; 
  /*const gsl_min_fminimizer_type *T1;
  gsl_min_fminimizer *s1;*/ 
  int status=GSL_SUCCESS;
  int iter;double myerror=0.0;
  *finitestepsize=delta;/** --copy-- current delta to finitestepsize as finitstepsize will be changed in here **/
  /** could probably do the memory withouth alloc free in here but outside of function to avoid extra mallocs */
   F1.function = compute_mlik_gaus_brent;
   F1.params = gparams;
   gsl_min_fminimizer_set (s1, &F1, *finitestepsize,lower,upper);
   iter=0;
        do
         {
           iter++; 
	 /*Rprintf("iter=%d\n",iter);*/
	   status = gsl_min_fminimizer_iterate (s1);
           *finitestepsize = gsl_min_fminimizer_x_minimum (s1);
           lower=gsl_min_fminimizer_x_lower (s1);
	   upper=gsl_min_fminimizer_x_upper (s1);
          status = gsl_min_test_interval (lower,upper, 0.00001, 0.0);
	  /*Rprintf("[%e,%e] min=%e\n",lower,upper,finitestepsize);*/
           /*if (status==GSL_SUCCESS)*/ 
             /*break;*/
         }
        while (status == GSL_CONTINUE && iter < maxiters_hessian); 
     
     /** now get the error in hessian **/
     
     /*compute_mlik_nm_brent(finitestepsize, &gparams);*/
    
     /*gsl_min_fminimizer_free (s1);*/
 
  myerror=compute_mlik_gaus_brent(*finitestepsize, gparams);
  *saverror=myerror;/** this is a pointer to the dag->HessianError[nodeid] */
  if (errverbose>0) { 
    Rprintf("Gaussian node: Brent minimiser: error in mlik=%e in [%e,%e] with best h=%e\n",myerror,lower,upper,*finitestepsize);
  }
return(myerror);  
}	
