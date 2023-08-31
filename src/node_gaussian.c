#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h" 
#include "node_gaussian.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

/** ****************************************************************************************************
 ***** calc an individual logistic regression model 
 *******************************************************************************************************/
void calc_node_Score_gaus(network *dag, datamatrix *obsdata, int nodeid,  int verbose,
                                datamatrix *designmatrix, 
			    const double priormean, const double priorsd, const double priorgamshape, const double priorgamscale,
                                const int maxiters, const double epsabs,int storeModes)
{
 int i,ss,status,index;
 /*unsigned int k,j;*/
 int iter=0;
 /*unsigned int numparents=0;*/
 double logscore=0.0;
 double mydet=0.0;
 const gsl_multiroot_fdfsolver_type *T;
 gsl_multiroot_fdfsolver *s;
 gsl_multiroot_function_fdf FDF;
 gsl_vector *myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,*term1,*term2,*term3,*vectmp3long,*localbeta;
 gsl_matrix *hessgvalue,*mattmp2,*mattmp3,*mattmp4;
 struct fnparams gparams;/** for passing to the gsl zero finding functions */
 double gvalue,n,m;
 gsl_permutation *perm=0;
 gsl_permutation *initsperm;
 
 /** build design matrix which is designmatrix->defn, designmatrix->Y plus priors designmatrix->priorsd, designmatrix->priormean **/
 build_designmatrix_gaus(dag,obsdata,priormean, priorsd,priorgamshape,priorgamscale,designmatrix,nodeid, storeModes); 
 
   vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
   vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
   vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
   vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
   vectmp3long = gsl_vector_alloc (obsdata->numDataPts);
   term1 = gsl_vector_alloc (designmatrix->numparams);
   term2 = gsl_vector_alloc (designmatrix->numparams);
   term3 = gsl_vector_alloc (designmatrix->numparams);
   hessgvalue = gsl_matrix_alloc (designmatrix->numparams+1,designmatrix->numparams+1);/** will hold hessian matrix NOTE: NOT used in solver. the +1 is due to precision term
                                                                                            as it allocates its own storage based on the problem dimension**/
   mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
   mattmp3 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
   mattmp4 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
   localbeta = gsl_vector_alloc (designmatrix->numparams);/** scratch space in later functions **/
   initsperm = gsl_permutation_alloc (designmatrix->numparams);/** for use with initial guesses */
   myBeta = gsl_vector_alloc (designmatrix->numparams+1);/** this will hold the parameter point estimates INC tau-precision*/
   
   /** now store in gparams for sending to gsl solver routines */
   gparams.Y=designmatrix->Y;
   gparams.X=designmatrix->datamatrix;
   gparams.priormean=designmatrix->priormean;
   gparams.priorsd  =designmatrix->priorsd;
   gparams.priorgamshape=designmatrix->priorgamshape;
   gparams.priorgamscale  =designmatrix->priorgamscale;
   
   gparams.vectmp1=vectmp1;
   gparams.vectmp2=vectmp2;
   gparams.vectmp1long=vectmp1long;
   gparams.vectmp2long=vectmp2long;
   gparams.vectmp3long=vectmp3long;
   gparams.term1=term1;
   gparams.term2=term2;
   gparams.term3=term3;
   
   gparams.mattmp2=mattmp2;
   gparams.mattmp3=mattmp3;
   gparams.mattmp4=mattmp4;
   
   gparams.beta=localbeta;
   
   gparams.perm=initsperm;
 
   dag->nodeScoresErrCode[nodeid]=0;/** reset error code to no error **/
 
   /** now we need to solve system defined in laplace_dg()=0 */
           
    FDF.f = &laplace_gaus_dg;
    FDF.df = &laplace_gaus_hessg;
    FDF.fdf = &wrapper_gaus_fdf;
    FDF.n = designmatrix->numparams+1;/** +1 for the precision term */
    FDF.params = &gparams;  
   /** ******************** FIRST TRY for a root using hybridsj  *******************************************************/
    iter=0;
    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc (T, designmatrix->numparams+1); /** +1 for the precision term */
    /*status_inits=*/generate_gaus_inits(myBeta,&gparams, verbose);
    status=GSL_FAILURE;/** just set it to something not equal to GSL_SUCCESS */
    
    gsl_multiroot_fdfsolver_set (s, &FDF, myBeta);
      
    iter=0; 
       do
         {
           iter++;
       
           status = gsl_multiroot_fdfsolver_iterate (s);
           #ifdef PRINTGSL
           print_state (iter, s);
           #endif
          if (status)
             break;
     
           status = gsl_multiroot_test_residual (s->f, epsabs);
         }
       while (status == GSL_CONTINUE && iter < maxiters);
     if(status != GSL_SUCCESS){dag->nodeScoresErrCode[nodeid]=1;} else {dag->nodeScoresErrCode[nodeid]=0;}/** reset error code to no error **/
       if( (status != GSL_SUCCESS) && verbose){Rprintf ("Zero finding warning: status = %s at nodeid %d\n", gsl_strerror (status),nodeid);}
       gsl_vector_memcpy(myBeta,s->x);
     gsl_multiroot_fdfsolver_free(s);
       
       /** ******************** END of FIRST TRY for a root using hybridsj  *******************************************************/
   if(status!=GSL_SUCCESS){/** try other solver **/
    /*Rprintf("gaussian: using hybridsj failed so re-trying with unscaled hybridj: nodeid=\n",nodeid+1);*/ 
    iter=0; 
    T = gsl_multiroot_fdfsolver_hybridj;
    s = gsl_multiroot_fdfsolver_alloc (T, designmatrix->numparams+1); /** +1 for the precision term */
    /*status_inits=*/generate_gaus_inits(myBeta, &gparams, verbose);
    status=GSL_FAILURE;/** just set it to something not equal to GSL_SUCCESS */
    
    gsl_multiroot_fdfsolver_set (s, &FDF, myBeta);
      
    iter=0; 
       do
         {
           iter++;
       
           status = gsl_multiroot_fdfsolver_iterate (s);
           #ifdef PRINTGSL
           print_state (iter, s);
           #endif
          if (status)
             break;
     
           status = gsl_multiroot_test_residual (s->f, epsabs);
         }
       while (status == GSL_CONTINUE && iter < maxiters);
     if(status != GSL_SUCCESS){dag->nodeScoresErrCode[nodeid]=1;} else {dag->nodeScoresErrCode[nodeid]=0;}/** reset error code to no error **/
       if( (status != GSL_SUCCESS) && verbose){Rprintf ("Zero finding warning: status = %s at nodeid %d\n", gsl_strerror (status),nodeid);}
       gsl_vector_memcpy(myBeta,s->x);
     gsl_multiroot_fdfsolver_free(s);
    
   }
     
       if(storeModes){/** keep a copy of the parameter modes found for use later in other function calls etc**/
	 index=0;    /*Rprintf("size of beta=%d %f %f\n",myBeta->size, gsl_vector_get(myBeta,0),gsl_vector_get(myBeta,1));*/
		     for(i=0;i<dag->numNodes+2;i++){/** roll myBeta into dag->modes into the appropriate columns**/
		       if(gsl_matrix_get(dag->modes,nodeid,i)!=DBL_MAX){
			 gsl_matrix_set(dag->modes,nodeid,i,gsl_vector_get(myBeta,index++));}} 
                   /*for(i=0;i<dag->numNodes;i++){Rprintf("%e ",gsl_matrix_get(dag->modes,nodeid,i));}Rprintf("\n");*/
		   
		   }        
             

  /** ************************************/         
  /** we now have all the individual parts so put it together to the laplace approx */
 /* if(status != GSL_SUCCESS){logscore= -DBL_MAX;Rprintf("gaussian no root at node %d\n",nodeid+1);*//** root finding failed so discard model e.g. set fit to worst possible */
  /*} else { */
     
  laplace_gaus_g(myBeta,&gparams, &gvalue);
  laplace_gaus_hessg(myBeta,&gparams, hessgvalue);
  
   n=obsdata->numDataPts;
   m=designmatrix->numparams+1;
   perm = gsl_permutation_alloc (m);
   /*status2=*/gsl_linalg_LU_decomp(hessgvalue,perm,&ss);
   /*if(status2!= GSL_SUCCESS){logscore= -DBL_MAX;Rprintf("gaussian no inversion at node %d\n",nodeid+1);
   } else {*/
     mydet=gsl_linalg_LU_lndet(hessgvalue);/** compute determinant but this might be a nan - overflow?*/
     /*if(gsl_isnan(mydet)){Rprintf("no det - gaussian node at node %d\n",nodeid+1);logscore= -DBL_MAX;*/ /** is nan so return default mlik error value */
     /*} else {*/ /** all is ok so now compute the actual laplace value */
     logscore= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n); /** this is the final value */
         /* }*/
     if(gsl_isnan(logscore)){logscore= R_NaN;
                           dag->nodeScoresErrCode[nodeid]=2;} 
   /*}*/
    /*Rprintf("node=%d logscore=%f\n",nodeid,logscore);*/
    /*** Last Step before return - free all the gsl vectors, matrices, other etc **/
   gsl_vector_free(designmatrix->Y);
   gsl_matrix_free(designmatrix->datamatrix);
   gsl_vector_free(designmatrix->priormean);
   gsl_vector_free(designmatrix->priorsd);
   gsl_vector_free(designmatrix->priorgamshape);
   gsl_vector_free(designmatrix->priorgamscale);
   
   gsl_vector_free(myBeta);
   gsl_vector_free(localbeta);
   gsl_vector_free(vectmp1);
   gsl_vector_free(vectmp2);
   gsl_vector_free(vectmp1long);
   gsl_vector_free(vectmp2long);
   gsl_vector_free(term1);
   gsl_vector_free(term2);
   gsl_vector_free(term3);
   gsl_vector_free(vectmp3long);
   gsl_matrix_free(hessgvalue);
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(mattmp3);
   gsl_matrix_free(mattmp4);
   
   gsl_permutation_free(initsperm);
   /*if(status == GSL_SUCCESS){*/gsl_permutation_free(perm);/*}*/ /** only allocate this is status==GSL_SUCCESS */
   /*gsl_multiroot_fdfsolver_free (s);*/
    
   /*return(logscore);*/ 
   dag->nodeScores[nodeid]=logscore;

}/** ****************************************************************************************************
 *******************************************************************************************************/
/** ****************************************************************************************************
 *******************************************************************************************************/
void calc_gaussian_marginal(network *dag, datamatrix *obsdata, int nodeid,  int verbose,
                                datamatrix *designmatrix, 
		                const double priormean, const double priorsd,const double priorgamshape,const double priorgamscale,
                                const int maxiters, const double epsabs, double *denom_modes,int paramid, double betafixed, double mlik, double *posterior)
{			        
int i,j,ss,status/* ,status2,status_inits,index*/;
 /* unsigned int k,j;*/
 int iter=0;
 /* unsigned int numparents=0;*/ 
 double logscore=0.0;
 const gsl_multiroot_fdfsolver_type *T;
 gsl_multiroot_fdfsolver *s;
 gsl_multiroot_function_fdf FDF;
 gsl_vector *myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,*term1,*term2,*term3,*vectmp3long,*localbeta,*betafull,*dgvaluesfull;
 gsl_matrix *hessgvalue,*mattmp2,*mattmp3,*mattmp4,*hessgvaluefull;
 struct fnparams gparams;/** for passing to the gsl zero finding functions */
 double gvalue,n,m;
 gsl_permutation *perm=0;
 gsl_permutation *initsperm;
 double val=0.0;
 
 /** build design matrix which is designmatrix->defn, designmatrix->Y plus priors designmatrix->priorsd, designmatrix->priormean **/
 build_designmatrix_gaus(dag,obsdata,priormean, priorsd,priorgamshape,priorgamscale,designmatrix,nodeid,0); 

   /** down to here is as for the network score calc which is an integral across all parameters - we now adjust this to that its across all parameters
      minus one, where this one is fixed at values across a grid **/
  /** ALWAYS HAVE AT LEAST TWO PARAMS: mean+precision if only a model with a single parameter then no integration required just evaluation of (-1/n)*g() **/
  /** ********************************************************************************************************************/
  switch(designmatrix->numparams+1){
    case 1:{ error("must always have at least two parameters - a mean term and a precision/variance!\n");break;}

   /** ********************************************************************************************************************/
    /** THIS IS THE MAIN CASE AND REST OF FUNCTION IS IN HERE***************************************************************/
    default:{/** **/
    /** allocate only once here since same dimension within one node for the marginals**/
   /** GENERAL IDEA - keep the same dimensions as in the full margLik calc but drop off terms at the end if needed */  
  
    vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
    vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
    vectmp3long = gsl_vector_alloc (obsdata->numDataPts);
    term1 = gsl_vector_alloc (designmatrix->numparams);
    term2 = gsl_vector_alloc (designmatrix->numparams);
    term3 = gsl_vector_alloc (designmatrix->numparams);
    hessgvalue = gsl_matrix_alloc (designmatrix->numparams+1-1,designmatrix->numparams+1-1);/** will hold hessian matrix NOTE: NOT used in solver
                                                                                            as it allocates its own storage based on the problem dimension**/
    mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    mattmp3 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    mattmp4 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    localbeta = gsl_vector_alloc (designmatrix->numparams);/** scratch space in later functions **/
    initsperm = gsl_permutation_alloc (designmatrix->numparams);/** for use with initial guesses */
    
    betafull = gsl_vector_alloc (designmatrix->numparams+1);/** this will hold the re-build full vector of all parameters */
    dgvaluesfull = gsl_vector_alloc (designmatrix->numparams+1);/** this will hold the re-build full vector of all parameters derivs */
    hessgvaluefull = gsl_matrix_alloc (designmatrix->numparams+1,designmatrix->numparams+1);/**  will hold hessian matrix **/
    
   /* iter=0; */      
    FDF.f = &laplace_gaus_dg_marg;
    FDF.df = &laplace_gaus_hessg_marg;
    FDF.fdf = &wrapper_gaus_fdf_marg;
    FDF.n = designmatrix->numparams+1-1;/** +1 for the precision term -1 as we drop off a single term as marginal post. estimation*/
    FDF.params = &gparams;  
    myBeta = gsl_vector_alloc (designmatrix->numparams+1-1);/** this will hold the parameter point estimates INC tau-precision -1=marginal*/
    T = gsl_multiroot_fdfsolver_hybridj;
    s = gsl_multiroot_fdfsolver_alloc (T, designmatrix->numparams+1-1); /** +1 for the precision term -1=marginal */
   
   /** now send */ 
   gparams.Y=designmatrix->Y;
   gparams.X=designmatrix->datamatrix;
   gparams.priormean=designmatrix->priormean;
   gparams.priorsd  =designmatrix->priorsd;
   gparams.priorgamshape=designmatrix->priorgamshape;
   gparams.priorgamscale  =designmatrix->priorgamscale;
   
   gparams.vectmp1=vectmp1;
   gparams.vectmp2=vectmp2;
   gparams.vectmp1long=vectmp1long;
   gparams.vectmp2long=vectmp2long;
   gparams.vectmp3long=vectmp3long;
   gparams.term1=term1;
   gparams.term2=term2;
   gparams.term3=term3;
 
   gparams.mattmp2=mattmp2;
   gparams.mattmp3=mattmp3;
   gparams.mattmp4=mattmp4;
   
   gparams.beta=localbeta;
   gparams.perm=initsperm;
   gparams.betafull=betafull;
   gparams.dgvaluesfull=dgvaluesfull;
   gparams.hessgvalues=hessgvaluefull;
   gparams.betafixed=0.0;/** these will be changed in loop below*/
   gparams.betaindex=paramid;/** this is fixed - the variable for which the posterior is calculated **/
   
   generate_gaus_inits_marg(myBeta,&gparams);/** generate initial estimates for the remaining variable - not the posterior variable **/
   
    /** README. To avoid the user needing to specify a range of evaluation we run all the integrations twice, first to find the two end points close to zeros
   and then one final time and actual save the results. We start from the mode and run iterations to the left (increasing) and to the right (decreasing) and
   then finally run the integration between these end points and save the result. The code looks a lot but its largely just repetition */
  
   /** set guess for other parameters to modes */
   j=0;
   for(i=0;i<designmatrix->numparams+1;i++){if(i!= paramid){gsl_vector_set(myBeta,j++,denom_modes[i]);}} /** use modes as initial values ignoring current mode which is to be marginalised **/ 

     gparams.betafixed=betafixed;
    gsl_multiroot_fdfsolver_set (s, &FDF, myBeta);
    
    iter=0; 
    do
      {
	iter++;
	status = gsl_multiroot_fdfsolver_iterate (s);
	if (status)
	  break;
	
	status = gsl_multiroot_test_residual (s->f, epsabs);
      }
    while (status == GSL_CONTINUE && iter < maxiters);
    
    if(status != GSL_SUCCESS){
      if( verbose>0 ) {
	Rprintf ("Zero finding error: status = %s at x=%f\n", gsl_strerror (status),gparams.betafixed);/*exit(1);*/
      }
    }
    gsl_vector_memcpy(myBeta,s->x);
    /** we now have all the individual parts so put it together to the laplace approx */
    
    
    /*for(iter=0;iter< myBeta->size;iter++){Rprintf("%f ",gsl_vector_get(myBeta,iter));}Rprintf("\n");*/
    
    laplace_gaus_g_marg(myBeta,&gparams, &gvalue);
    laplace_gaus_hessg_marg(myBeta,&gparams, hessgvalue);
    n=obsdata->numDataPts;
    m=designmatrix->numparams+1-1;/** -1 since a marginal calc. one less dimension in the integral **/
    perm = gsl_permutation_alloc (m);
    gsl_linalg_LU_decomp(hessgvalue,perm,&ss);
    logscore= -n*gvalue-0.5*gsl_linalg_LU_lndet(hessgvalue)+(m/2.0)*log((2.0*M_PI)/n); /** this is the final value */ 
    val=exp(logscore-mlik);
    /*Rprintf("got betafixed=%f mlik=%f and gval=%f value=%f %f\n",betafixed,mlik,logscore,val,gsl_linalg_LU_lndet(hessgvalue));*/
    /*Rprintf("got betafixed=%f mlik=%f logscore=%f exp(logscore-mlik)=%f gvalue=%f\n",betafixed,mlik,logscore,val,gvalue);*/
    if(gsl_isnan(val)){*posterior=R_NaN;
    } else {*posterior=val;}
    
    
    
   gsl_vector_free(myBeta);
   gsl_vector_free(vectmp1);
   gsl_vector_free(vectmp2);
   gsl_vector_free(vectmp1long);
   gsl_vector_free(vectmp2long);
   gsl_vector_free(dgvaluesfull); 
   gsl_vector_free(term1);
   gsl_vector_free(term2);
   gsl_vector_free(term3);
   gsl_vector_free(vectmp3long);
   gsl_vector_free(betafull);
   gsl_vector_free(localbeta);
   gsl_matrix_free(hessgvalue);
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(mattmp3);
   gsl_matrix_free(mattmp4);
   gsl_matrix_free(hessgvaluefull);
   gsl_permutation_free(initsperm);
   gsl_permutation_free(perm);
   gsl_multiroot_fdfsolver_free (s);
   
   gsl_vector_free(designmatrix->Y);
   gsl_matrix_free(designmatrix->datamatrix);
   gsl_vector_free(designmatrix->priormean);
   gsl_vector_free(designmatrix->priorsd);
   gsl_vector_free(designmatrix->priorgamshape);
   gsl_vector_free(designmatrix->priorgamscale);
   }
   
    } /** end of switch **/

   
}   

/** **************************************************************************************************************/
/** **************************************************************************************************************/
/** build the design matrix - plus other associated things *******************************************************/
void build_designmatrix_gaus(network *dag,datamatrix *obsdata, double priormean, double priorsd,double priorgamshape, double priorgamscale,datamatrix *designmatrix, 
			     int nodeid, int storeModes)
{
  
 int i,j,k;
 int numparents=0;
 gsl_vector *Y,*vecpriormean,*vecpriorsd,*vecpriorgamshape,*vecpriorgamscale;
 gsl_vector_int *parentindexes=0;
 gsl_matrix *datamat;

 if(dag->maxparents>0){
 parentindexes=gsl_vector_int_alloc(dag->maxparents);
 /** collect parents of this node **/
 for(j=0;j<dag->numNodes;j++){
              if(   dag->defn[nodeid][j]==1    /** got a parent so get its index **/
                 && numparents<dag->maxparents /** if numparents==dag->maxparents then we are done **/
                ){
		        gsl_vector_int_set(parentindexes,numparents++,j);/** store index of parent **/
                  }
		}
 } /** end maxparents=0 check */
/** this part is new and just for posterior param est - it does not affect Laplace approx in any way****/
  /** setup matrix where each non DBL_MAX entry in a row is for a parameter to be estimated and the col is which param
      first col is for the intercept */
  if(storeModes){
    for(k=0;k<dag->numNodes+3;k++){gsl_matrix_set(dag->modes,nodeid,k,DBL_MAX);} /** initialise row to DBL_MAX **/
    gsl_matrix_set(dag->modes,nodeid,0,1);/** the intercept term - always have an intercept - but not in dag.m definition */  
    for(k=0;k<numparents;k++){gsl_matrix_set(dag->modes,nodeid,gsl_vector_int_get(parentindexes,k)+1,1);} /** offset is 1 due to intercept */
    gsl_matrix_set(dag->modes,nodeid,dag->numNodes+1,1);/** the precision term put at end of other params */  
  }
  
  datamat=gsl_matrix_alloc(obsdata->numDataPts,numparents+1); /** only +1 since this is design matrix for the mean */
  designmatrix->datamatrix=datamat;
  Y=gsl_vector_alloc(obsdata->numDataPts);
  designmatrix->Y=Y;
  vecpriormean=gsl_vector_alloc(numparents+1);
  designmatrix->priormean=vecpriormean;	
  vecpriorsd=gsl_vector_alloc(numparents+1);
  designmatrix->priorsd=vecpriorsd;
  vecpriorgamshape=gsl_vector_alloc(1); /** only 1 of these per node */
  designmatrix->priorgamshape=vecpriorgamshape;
  vecpriorgamscale=gsl_vector_alloc(1); /** only 1 of these per node */
  designmatrix->priorgamscale=vecpriorgamscale;
  
  /** create design matrix - copy relevant cols from the observed data **/
 /** int** designmatrix is just used as storage space, fill up from left cols across until as far as needed */
 for(i=0;i<obsdata->numDataPts;i++){/** for each observed data point **/
   /*designmatrix->data[i][0]=1;*//** intercept term **/
   gsl_matrix_set(designmatrix->datamatrix,i,0,1.0); 
  /** copy values at node - response values - into vector: */
  gsl_vector_set(designmatrix->Y,i,obsdata->defn[i][nodeid]);
   
   for(k=0;k<numparents;k++){/* now build design matrice of explanatories */
	           gsl_matrix_set(designmatrix->datamatrix,i,k+1,obsdata->defn[i][gsl_vector_int_get(parentindexes,k)]); 
                            } /** end of explanatories **/
                         
   } /** end of data point loop */   
                        
   designmatrix->numparams=numparents+1;/** +1 for intercept - excludes precision **/
   /** now set the priormean and priorsd vector - choose the correct prior values */
   
   for(k=0;k<designmatrix->numparams;k++){
                                          gsl_vector_set(designmatrix->priormean,k,priormean);
                                          gsl_vector_set(designmatrix->priorsd,k,priorsd);
					  }
   gsl_vector_set(designmatrix->priorgamshape,0,priorgamshape);
   gsl_vector_set(designmatrix->priorgamscale,0,priorgamscale);
   
   gsl_vector_int_free(parentindexes);/** finished with this **/
  
 /*  Rprintf("##########################################\n");
  Rprintf("got %d parents\n",numparents);
for(i=0;i<obsdata->numDataPts;i++){
   Rprintf("Y=%f\t",gsl_vector_get(designmatrix->Y,i));  
   for(k=0;k<numparents+1;k++){Rprintf("%f\t",gsl_matrix_get(designmatrix->datamatrix,i,k));}Rprintf("\n");}
   */
   /*Rprintf("mean=%f sd=%f gam1=%f gam2=%f",gsl_vector_get(designmatrix->priormean,0),
	                                   gsl_vector_get(designmatrix->priorsd,0),
	                                   gsl_vector_get(designmatrix->priorgamshape,0),
	                                   gsl_vector_get(designmatrix->priorgamscale,0));*/
     
     


}  


/** ***************************************************************************************************************
******************************************************************************************************************* 
** laplace method = int^b_a exp(-lambda g(y)) h(y) dy = exp(-lambda g(y*)) h(y*) (2PI/lambda)^(d/2) det(hess)^(1/2)
** lambda = sample size n, g(y) = -(1/n)* log( f(D|betas)f(betas) ) e.g. -(1/n)* log (like*prior)
*******************************************************************************************************************
******************************************************************************************************************/
/** GAUSSIAN NODE ************************************************************************************************/  
/** **************************************************************************************************************/
/** g(y) = -(1/n)* log( f(D|betas)f(betas) */ 
/** **************************************************************************************************************/
int laplace_gaus_g (const gsl_vector *betaincTau, void *params,double *gvalue)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;
	     /*gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;*/
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
       const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale;
       gsl_vector *beta = ((struct fnparams *) params)->beta;
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** number of coefficients excluding tau-precision */
       double term2=0;
       double term3=0;
       double term4=0;
       double term5=0;
       double term6=0;
       double storedbl1,storedbl2,storedbl3;
      
       double tau=gsl_vector_get(betaincTau,m);/** extract the tau-precision from *beta */
       /*Rprintf("passed tau=%f\n",tau);
       Rprintf("Y=%f\n",gsl_vector_get(Y,0));
       Rprintf("X[0]=%f %f %f\n",gsl_matrix_get(X,0,0),gsl_matrix_get(X,0,1),gsl_matrix_get(X,0,2));
       Rprintf("X[10]=%f %f %f\n",gsl_matrix_get(X,10,0),gsl_matrix_get(X,10,1),gsl_matrix_get(X,10,2));*/
       int i=0;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betaincTau,i));/*Rprintf("passed beta=%f\n",gsl_vector_get(beta,i));*/
       }
     
     /** same as in logistic model */
     term2=0; for(i=0;i<m;i++){term2+=-log(sqrt(2.0*M_PI)*gsl_vector_get(priorsd,i));}
     /*Rprintf("term2 (Rterm3)=%f\n",term2);*/
     
     /** same as in logistic model */
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
     /*Rprintf("term3 (Rterm4)=%f\n",term3); */    
     
     /** Need -0.5tau*(Y%*%Y+ (X%*%myBeta)%*%(X%*%myBeta)-2Y%*%X%*%myBeta) */
     
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_ddot (Y, vectmp1long, &storedbl1);/** storedbl1 holds Y%*%(X%*%mybeta)**/
     storedbl1= -2.0*storedbl1;/** now gives 2Y%*%X%*%myBeta */
     gsl_blas_ddot (vectmp1long, vectmp1long, &storedbl2);/** storebdl2 is (X%*%myBeta)%*%(X%*%myBeta) **/
     gsl_blas_ddot (Y, Y, &storedbl3); /** storebdl3 is Y%*%Y**/
     term4= -(tau/2.0)*(storedbl1+storedbl2+storedbl3);
     /*Rprintf("term4 (Rterm2)=%f\n",term4);*/
     
     term5= (n/2.0)*log(tau/(2.0*M_PI));
     /*Rprintf("term5 (Rterm1)=%f\n",term5);*/
     term6=  -gsl_vector_get(priorgamshape,0)*log(gsl_vector_get(priorgamscale,0))
             -gsl_sf_lngamma(gsl_vector_get(priorgamshape,0)) 
	     +(gsl_vector_get(priorgamshape,0)-1)*log(tau)
	     -(tau/gsl_vector_get(priorgamscale,0));
     /*Rprintf("term6 (Rterm5)=%5.10f\n",term6);*/
     *gvalue=(-1.0/n)*(term2+term3+term4+term5+term6);

     /*Rprintf("loglike=%f\n",term5+term4);*/

       return GSL_SUCCESS;
     }
        
     
/** **************************************************************************************************************/
/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_gaus_dg (const gsl_vector *betaincTau, void *params, gsl_vector *dgvalues)
{      
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
        gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;/** numobs long **/
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
       const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale;
       gsl_vector *beta = ((struct fnparams *) params)->beta;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       gsl_vector *term2 = ((struct fnparams *) params)->term2;
       /*gsl_vector *term3 = ((struct fnparams *) params)->term3;*/
       /*gsl_matrix *mattmp1 = ((struct fnparams *) params)->mattmp1;*/
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** number of coefficients excluding tau-precision */
       double tau=gsl_vector_get(betaincTau,m);/** extract the tau-precision from *beta */
       int i=0; 
       double storedbl1/*,storedbl2,storedbl3,tmptau*/;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betaincTau,i));}
 
      /** do dg_beta terms first - store these in term1 vec*/ 
      /** -(beta_j - mu_j)/sd_j^2" **/
     gsl_vector_memcpy(vectmp1,beta);/** copy beta to temp vec */
     gsl_vector_memcpy(vectmp2,priormean);
     gsl_vector_scale(vectmp2,-1.0);
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1= beta-mean**/
     gsl_vector_memcpy(vectmp2,priorsd);
     gsl_vector_mul(vectmp2,priorsd);/** square all elements in priorsd and store in vectmp2 */
     gsl_vector_div(vectmp1,vectmp2);
     gsl_vector_scale(vectmp1,-1.0); 
     gsl_vector_memcpy(term1,vectmp1);/** the prior term in dg_dbeta **/
  
     /**  tau *sum(y_i X_ij - X_i beta X_ij) */
     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);/** each entry in vectmp1 is sum(y_i x_ij) for a fixed j**/
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp1long, 0.0, vectmp2);/** each entry in vectmp2 is sum(x_ij*sum(x_i*beta)) **/
     gsl_vector_scale(vectmp2,-1.0);/** need negative e.g. sum(-X_i beta X_ij) */
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1 contains sum(y_i X_ij - X_i beta X_ij) **/
     gsl_vector_scale(vectmp1,tau);/** mult by tau **/
     gsl_vector_memcpy(term2,vectmp1);/** the remaining term in dg_beta **/
     /** put the two parts of the dg_dbeta terms together */
     gsl_vector_add(term1,term2);
     gsl_vector_scale(term1,(-1.0/n));
     /** now store the dg_dbeta terms */
     for(i=0;i<m;i++){gsl_vector_set(dgvalues,i,gsl_vector_get(term1,i));}
     
     /** dg_dtau **/
     gsl_vector_scale(vectmp1long,-1.0);/** from above vectmp1long is X*beta so this is -X*beta */
     gsl_vector_add(vectmp1long,Y);/** Y-X*beta **/
     gsl_vector_memcpy(vectmp2long,vectmp1long);
     gsl_blas_ddot (vectmp2long, vectmp1long, &storedbl1);/** get sum((Y_i-X_i*beta)^2) using dot product */
     storedbl1= storedbl1*(-0.5);
     gsl_vector_set(dgvalues,m,(-1.0/n)*(
                                          (n/(2.0*tau))+
				          storedbl1 + 
				          (gsl_vector_get(priorgamshape,0)-1.0)/tau
				          - (1.0/gsl_vector_get(priorgamscale,0)))
					  
                                          );   
					 
 return GSL_SUCCESS;
     }


/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_gaus_hessg (const gsl_vector *betaincTau, void *params, gsl_matrix *hessgvalues)
{      
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
       gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
       gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
       /*gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;*//** numobs long **/
       /*const gsl_vector *priormean = ((struct fnparams *) params)->priormean; */
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
       /*const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale;*/
       gsl_vector *beta = ((struct fnparams *) params)->beta;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       /*gsl_vector *term2 = ((struct fnparams *) params)->term2;*/
       /*gsl_vector *term3 = ((struct fnparams *) params)->term3;*/
       gsl_matrix *mattmp3 = ((struct fnparams *) params)->mattmp3;
       gsl_matrix *mattmp2 = ((struct fnparams *) params)->mattmp2;
       int n=Y->size;/** no. observations **/
       int m=X->size2;/** number of coefficients excluding tau-precision */
       /*Rprintf("M=%d %f\n",m, gsl_matrix_get(X,2,2));*/
       double tau=gsl_vector_get(betaincTau,m);/** extract the tau-precision from *beta */
       int i=0;
       int j,k;
       double tmp1;
      /* double storedbl1,storedbl2,storedbl3,tmptau; */
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betaincTau,i));}
 
  /** matrix X%*%t(X) - this is symmetrical so dg_bj_b_k=dg_b_k_b_j **/
  gsl_matrix_memcpy(mattmp2,X);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                       1.0, X, mattmp2,
                       0.0, mattmp3);
  
  /*for(i=0;i<100;i++){Rprintf("%f \n",gsl_matrix_get(X,i,0));}*/
  
   /** matrix of hess excluding tau derivatives for the moment - note this is symmetrical*/
     for(j=0;j<m;j++){
       for(k=0;k<m;k++){
                    if(j!=k){/** second derivatives */
                      tmp1= (-1.0/n)*(-tau*gsl_matrix_get(mattmp3,j,k));
		      *gsl_matrix_ptr(hessgvalues,j,k)=tmp1;
                    } else {
		      tmp1= (-1.0/n)*(-tau*gsl_matrix_get(mattmp3,j,k)-1.0/(gsl_vector_get(priorsd,j)*gsl_vector_get(priorsd,j)));
		      *gsl_matrix_ptr(hessgvalues,j,k)=tmp1;}
                     }
                     }
   
  
  /** now for dg_dtau second deriv **/
  tmp1=(-1.0/n)*( (-n/(2.0*tau*tau)) - ( (gsl_vector_get(priorgamshape,0)-1)/(tau*tau)) );
  *gsl_matrix_ptr(hessgvalues,m,m)=tmp1;
  /*Rprintf("final cell is %d %d %f\n",m,m,*gsl_matrix_ptr(hessgvalues,m,m));*/
  
  /** now for dg_dtau_dbeta sum(y_i X_ij - X_i beta X_ij) - last ROW of hessian */
     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);/** each entry in vectmp1 is sum(y_i x_ij) for a fixed j**/
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp1long, 0.0, vectmp2);/** each entry in vectmp2 is sum(x_ij*sum(x_i*beta)) **/
     /*Rprintf("beta=%10.15f %10.15f\n",gsl_vector_get(beta,0),gsl_vector_get(beta,1));
     Rprintf("==%f %f\n",gsl_vector_get(vectmp2,0),gsl_vector_get(vectmp2,1));
     Rprintf("==%f %f\n",gsl_vector_get(vectmp1,0),gsl_vector_get(vectmp1,1));*/
     gsl_vector_scale(vectmp2,-1.0);/** need negative e.g. sum(-X_i beta X_ij) */
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1 contains sum(y_i X_ij - X_i beta X_ij) **/
     gsl_vector_memcpy(term1,vectmp1);/** the remaining term in dg_beta **/
     gsl_vector_scale(term1,(-1.0/n));
     
     /** last row in hessian **/
     for(j=0;j<m;j++){*gsl_matrix_ptr(hessgvalues,m,j)=gsl_vector_get(term1,j);}
 
 /** dg_dbeta_dtau is same as dg_dtau_dbeta - note symmetry here **/
     
     /** last col in hessian **/
     for(j=0;j<m;j++){*gsl_matrix_ptr(hessgvalues,j,m)=gsl_vector_get(term1,j);}
 
  return GSL_SUCCESS;
} 
 
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int generate_gaus_inits(gsl_vector *myBeta, struct fnparams *gparams, int errverbose){

    /** just hard coded for moment */
           
    /** use least squares estimates to get started **/
    /** beta_hat= (X^T X)^{-1} X^T y **/
       int haveError=0;
       const gsl_vector *Y = gparams->Y;/** design matrix **/
       const gsl_matrix *X = gparams->X;/** response variable **/
       gsl_vector *vectmp1= gparams->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = gparams->vectmp2;
       gsl_vector *vectmp1long = gparams->vectmp1long;/** numobs long **/
       gsl_vector *vectmp2long = gparams->vectmp2long;/** numobs long **/
      /* const gsl_vector *priormean = gparams->priormean;
       const gsl_vector *priorsd   = gparams->priorsd;
       const gsl_vector *priorgamshape = gparams->priorgamshape;
       const gsl_vector *priorgamscale = gparams->priorgamscale; */
       /*gsl_vector *beta = gparams->beta;*/
       /*gsl_vector *term1 = gparams->term1;
       gsl_vector *term2 = gparams->term2;
       gsl_vector *term3 = gparams->term3;*/
       gsl_matrix *mattmp2 = gparams->mattmp2;/** same dim as X*/
       gsl_matrix *mattmp3 = gparams->mattmp3;/** p x p **/
       gsl_matrix *mattmp4 = gparams->mattmp4;/** p x p **/
       gsl_permutation *perm = gparams->perm;
      double n=Y->size;/** no. observations **/
     double m=X->size2;/** number of coefficients excluding tau-precision */
     unsigned int i;
     int ss;
     double variance=0.0;
     
    /*Rprintf("X: %d %d %d %d %d %d\n",X->size1,X->size2,mattmp2->size1,mattmp2->size2,mattmp3->size1,mattmp3->size2); */
    gsl_matrix_memcpy(mattmp2,X);
    gsl_blas_dgemm (CblasTrans, CblasNoTrans,    /** mattmp3 is p x p matrix X^T X **/
                       1.0, X, mattmp2,
                       0.0, mattmp3);
    gsl_permutation_init(perm);/** reset - might not be needed */                   
    gsl_linalg_LU_decomp(mattmp3,perm,&ss);      
    gsl_set_error_handler_off();/**Turning off GSL Error handler as this may fail as mattmp3 may be singular */     
    haveError=gsl_linalg_LU_invert (mattmp3, perm, mattmp4);/** mattmp4 is now inv (X^T X) */                   
    if(!haveError){ /** don't have error */
     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1); /** X^T Y */
     gsl_blas_dgemv (CblasNoTrans, 1.0, mattmp4, vectmp1, 0.0, vectmp2); 
     for(i=0;i<myBeta->size-1;i++){gsl_vector_set(myBeta,i,gsl_vector_get(vectmp2,i));} /** set to Least squares estimate */
    } else {
      if (errverbose > 0) {
	Rprintf ("caught gsl error - singular matrix in initial guess estimates\n");
      }
      /** must have a singular matrix so try all zeros for initial values*//*Rprintf("failed on invert\n");*/
      for(i=0;i<myBeta->size-1;i++){gsl_vector_set(myBeta,i,0.01);}
    }
    gsl_set_error_handler (NULL);/** restore the error handler*/
    /*for(i=0;i<myBeta->size-1;i++){Rprintf("%f ",gsl_vector_get(myBeta,i));}Rprintf("\n"); */
    /** now for variance estimate */
    /** first get y_hat estimate */
    gsl_blas_dgemv (CblasNoTrans, 1.0, X, vectmp2, 0.0, vectmp1long); /** vectmp1 is y_hat */ 
    gsl_vector_scale(vectmp1long,-1.0);/** - y_hat */
    gsl_vector_add(vectmp1long,Y);/** now have Y-y_hat (or -y_hat + Y) */
    gsl_vector_memcpy(vectmp2long,vectmp1long);
    gsl_blas_ddot (vectmp1long, vectmp2long, &variance);/** got sum((Y-Y_hat)^2) */
    variance=variance/(n-m);/** unbiased estimator using denominator n-#term in regression equation **/
    
    gsl_vector_set(myBeta,myBeta->size-1,1.0/variance);/** as are using precision parameterization not variance **/ 
   
    return GSL_SUCCESS;
}   
 

/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int wrapper_gaus_fdf (const gsl_vector *beta, void *gparams,
                     gsl_vector *dgvalues, gsl_matrix *hessgvalues)
     {
       laplace_gaus_dg(beta, gparams, dgvalues);
       laplace_gaus_hessg(beta, gparams, hessgvalues);
     
       return GSL_SUCCESS;
     }  
     
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/     
  /** **************************************************************************************************************/ 
/**** Gaussian marginal functions down here *********************************************************************/
/****************************************************************************************************************/
/** **************************************************************************************************************/
int laplace_gaus_g_marg (const gsl_vector *betashort, void *params,double *gvalue)
{      
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;
	    /* gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long; */
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
       const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale;
       gsl_vector *beta = ((struct fnparams *) params)->beta;/** only long enough for beta terms - no precision */
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** number of coefficients excluding tau-precision */
       double term2=0;
       double term3=0;
       double term4=0;
       double term5=0;
       double term6=0;
       double storedbl1,storedbl2,storedbl3;
       double tau;      
       int i=0;
       /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *betafull = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
            
        /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(betafull,0,betafixed);
                     for(i=1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(betafull->size-1)){gsl_vector_set(betafull,betafull->size-1,betafixed);
                     for(i=0;i<betafull->size-1;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}}
       
     if(betaindex>0 && betaindex<(betafull->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(betafull,betaindex,betafixed);
	 for(i=betaindex+1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}
     }	 
       
       tau=gsl_vector_get(betafull,m);/** extract the tau-precision from *beta */
      
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betafull,i));/*Rprintf("passed beta=%f\n",gsl_vector_get(beta,i));*/}
     
     /** same as in logistic model */
     term2=0; for(i=0;i<m;i++){term2+=-log(sqrt(2.0*M_PI)*gsl_vector_get(priorsd,i));}
    /* Rprintf("term2 (Rterm3)=%f\n",term2);*/
     
     /** same as in logistic model */
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
    /* Rprintf("term3 (Rterm4)=%f\n",term3);*/     
     
     /** Need -0.5tau*(Y%*%Y+ (X%*%myBeta)%*%(X%*%myBeta)-2Y%*%X%*%myBeta) */
     
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_ddot (Y, vectmp1long, &storedbl1);/** storedbl1 holds Y%*%(X%*%mybeta)**/
     storedbl1= -2.0*storedbl1;/** now gives 2Y%*%X%*%myBeta */
     gsl_blas_ddot (vectmp1long, vectmp1long, &storedbl2);/** storebdl2 is (X%*%myBeta)%*%(X%*%myBeta) **/
     gsl_blas_ddot (Y, Y, &storedbl3); /** storebdl3 is Y%*%Y**/
     term4= -(tau/2.0)*(storedbl1+storedbl2+storedbl3);
     /*Rprintf("term4 (Rterm2)=%f\n",term4);*/
     
     term5= (n/2.0)*log(tau/(2.0*M_PI));
     /*Rprintf("term5 (Rterm1)=%f\n",term5);*/
     term6=  -gsl_vector_get(priorgamshape,0)*log(gsl_vector_get(priorgamscale,0))
             -gsl_sf_lngamma(gsl_vector_get(priorgamshape,0)) 
	     +(gsl_vector_get(priorgamshape,0)-1)*log(tau)
	     -(tau/gsl_vector_get(priorgamscale,0));
    /* Rprintf("term6 (Rterm5)=%5.10f\n",term6);*/
     *gvalue=(-1.0/n)*(term2+term3+term4+term5+term6);

     /*Rprintf("loglike=%f\n",term5+term4);*/

     
     return GSL_SUCCESS;
     }
        
     
/** **************************************************************************************************************/
/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_gaus_dg_marg (const gsl_vector *betashort, void *params, gsl_vector *dgvaluesshort)
{      
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
        gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;/** numobs long **/
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
       const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale;
       gsl_vector *beta = ((struct fnparams *) params)->beta;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       gsl_vector *term2 = ((struct fnparams *) params)->term2;
       /*gsl_vector *term3 = ((struct fnparams *) params)->term3;*/
      /* gsl_matrix *mattmp1 = ((struct fnparams *) params)->mattmp1; */
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** number of coefficients excluding tau-precision */
       double tau;/*=gsl_vector_get(betaincTau,m);*//** extract the tau-precision from *beta */
       int i=0;int col;double tmp;
       double storedbl1/*,storedbl2,storedbl3,tmptau*/;
      /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *betafull = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
       gsl_vector *dgvaluesfull = ((struct fnparams *) params)->dgvaluesfull;
       
       /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(betafull,0,betafixed);
                     for(i=1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(betafull->size-1)){gsl_vector_set(betafull,betafull->size-1,betafixed);
                     for(i=0;i<betafull->size-1;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}}
       
     if(betaindex>0 && betaindex<(betafull->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(betafull,betaindex,betafixed);
	 for(i=betaindex+1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}
     }	 
       
       tau=gsl_vector_get(betafull,m);/** extract the tau-precision from *beta */
      
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betafull,i));/*Rprintf("passed beta=%f\n",gsl_vector_get(beta,i));*/} 
       
      /** do dg_beta terms first - store these in term1 vec*/ 
      /** -(beta_j - mu_j)/sd_j^2" **/
     gsl_vector_memcpy(vectmp1,beta);/** copy beta to temp vec */
     gsl_vector_memcpy(vectmp2,priormean);
     gsl_vector_scale(vectmp2,-1.0);
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1= beta-mean**/
     gsl_vector_memcpy(vectmp2,priorsd);
     gsl_vector_mul(vectmp2,priorsd);/** square all elements in priorsd and store in vectmp2 */
     gsl_vector_div(vectmp1,vectmp2);
     gsl_vector_scale(vectmp1,-1.0); 
     gsl_vector_memcpy(term1,vectmp1);/** the prior term in dg_dbeta **/
  
     /**  tau *sum(y_i X_ij - X_i beta X_ij) */
     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);/** each entry in vectmp1 is sum(y_i x_ij) for a fixed j**/
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp1long, 0.0, vectmp2);/** each entry in vectmp2 is sum(x_ij*sum(x_i*beta)) **/
     gsl_vector_scale(vectmp2,-1.0);/** need negative e.g. sum(-X_i beta X_ij) */
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1 contains sum(y_i X_ij - X_i beta X_ij) **/
     gsl_vector_scale(vectmp1,tau);/** mult by tau **/
     gsl_vector_memcpy(term2,vectmp1);/** the remaining term in dg_beta **/
     /** put the two parts of the dg_dbeta terms together */
     gsl_vector_add(term1,term2);
     gsl_vector_scale(term1,(-1.0/n));
     /** now store the dg_dbeta terms */
     for(i=0;i<m;i++){gsl_vector_set(dgvaluesfull,i,gsl_vector_get(term1,i));}
     
     /** dg_dtau **/
     gsl_vector_scale(vectmp1long,-1.0);/** from above vectmp1long is X*beta so this is -X*beta */
     gsl_vector_add(vectmp1long,Y);/** Y-X*beta **/
     gsl_vector_memcpy(vectmp2long,vectmp1long);
     gsl_blas_ddot (vectmp2long, vectmp1long, &storedbl1);/** get sum((Y_i-X_i*beta)^2) using dot product */
     storedbl1= storedbl1*(-0.5);
     gsl_vector_set(dgvaluesfull,m,(-1.0/n)*(
                                          (n/(2.0*tau))+
				          storedbl1 + 
				          (gsl_vector_get(priorgamshape,0)-1.0)/tau
				          - (1.0/gsl_vector_get(priorgamscale,0)))
					  
                                          ); 
     
     /** need to drop one cell in term1 before copying back */
     /** create an adjusted term1 which contains the term1 without the  re-inserted at the correct place **/
    col=0;
     for(i=0;i<betafull->size;i++){
       if(i!=betaindex){/** unless fixed variable then **/
	 tmp=gsl_vector_get(dgvaluesfull,i);
	 col=i;
	 if(i>betaindex){col=i-1;} 
                               gsl_vector_set(dgvaluesshort,col,tmp);}
	}
 
 /*Rprintf("end gaus_dg\n");
 for(i=0;i<betafull->size-1;i++){Rprintf("=%f\n",gsl_vector_get(dgvaluesshort,i));}	
 */
   return GSL_SUCCESS;
     }


/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */ 
/** **************************************************************************************************************/
int laplace_gaus_hessg_marg (const gsl_vector *betashort, void *params, gsl_matrix *hessgvaluesshort)
{      
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
       gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
       gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
      /* gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;*//** numobs long **/
      /* const gsl_vector *priormean = ((struct fnparams *) params)->priormean; */
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       const gsl_vector *priorgamshape = ((struct fnparams *) params)->priorgamshape;
      /* const gsl_vector *priorgamscale = ((struct fnparams *) params)->priorgamscale; */
       gsl_vector *beta = ((struct fnparams *) params)->beta;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       /*gsl_vector *term2 = ((struct fnparams *) params)->term2;
       gsl_vector *term3 = ((struct fnparams *) params)->term3;*/
       gsl_matrix *mattmp3 = ((struct fnparams *) params)->mattmp3;
       gsl_matrix *mattmp2 = ((struct fnparams *) params)->mattmp2;
       int n=Y->size;/** no. observations **/
       int m=X->size2;/** number of coefficients excluding tau-precision */
       double tau;
       /*double tau=gsl_vector_get(betaincTau,m);*//** extract the tau-precision from *beta */
       int i=0;int row,col;double tmp;
       int j,k;
       double tmp1;
      /* double storedbl1,storedbl2,storedbl3,tmptau; */
      /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *betafull = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
       gsl_matrix *hessgvaluefull = ((struct fnparams *) params)->hessgvalues;/** will hold "full hessian matrix" **/
       
       /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(betafull,0,betafixed);
                     for(i=1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(betafull->size-1)){gsl_vector_set(betafull,betafull->size-1,betafixed);
                     for(i=0;i<betafull->size-1;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}}
       
     if(betaindex>0 && betaindex<(betafull->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(betafull,betaindex,betafixed);
	 for(i=betaindex+1;i<betafull->size;i++){gsl_vector_set(betafull,i,gsl_vector_get(betashort,i-1));}
     }	 
       
       tau=gsl_vector_get(betafull,m);/** extract the tau-precision from *beta */       
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */
       /** STOP - NEED TO copy betaincTau into shorter beta since last entry is tau = precision */
       for(i=0;i<m;i++){gsl_vector_set(beta,i,gsl_vector_get(betafull,i));/*Rprintf("passed beta=%f\n",gsl_vector_get(beta,i));*/} 
       
  /** matrix X%*%t(X) - this is symmetrical so dg_bj_b_k=dg_b_k_b_j **/
  /*Rprintf("%d %d %d %d hesssize=%d %d\n",X->size1,X->size2,mattmp2->size1,mattmp2->size2,hessgvalues->size1,hessgvalues->size2);*/
  gsl_matrix_memcpy(mattmp2,X);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                       1.0, X, mattmp2,
                       0.0, mattmp3);
  
  /*for(i=0;i<100;i++){Rprintf("%f \n",gsl_matrix_get(X,i,0));}*/
  
   /** matrix of hess excluding tau derivatives for the moment - note this is symmetrical*/
     for(j=0;j<m;j++){
       for(k=0;k<m;k++){
                    if(j!=k){/** second derivatives */
                      tmp1= (-1.0/n)*(-tau*gsl_matrix_get(mattmp3,j,k));
		      *gsl_matrix_ptr(hessgvaluefull,j,k)=tmp1;
                    } else {
		      tmp1= (-1.0/n)*(-tau*gsl_matrix_get(mattmp3,j,k)-1.0/(gsl_vector_get(priorsd,j)*gsl_vector_get(priorsd,j)));
		      *gsl_matrix_ptr(hessgvaluefull,j,k)=tmp1;}
                     }
                     }
                    
  /** now for dg_dtau second deriv **/
  tmp1=(-1.0/n)*( (-n/(2.0*tau*tau)) - ( (gsl_vector_get(priorgamshape,0)-1)/(tau*tau)) );
  *gsl_matrix_ptr(hessgvaluefull,m,m)=tmp1;
  /*Rprintf("final cell is %d %d %f\n",m,m,*gsl_matrix_ptr(hessgvaluefull,m,m));*/
  
  /** now for dg_dtau_dbeta sum(y_i X_ij - X_i beta X_ij) - last ROW of hessian */
     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);/** each entry in vectmp1 is sum(y_i x_ij) for a fixed j**/
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/  
     gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp1long, 0.0, vectmp2);/** each entry in vectmp2 is sum(x_ij*sum(x_i*beta)) **/
     gsl_vector_scale(vectmp2,-1.0);/** need negative e.g. sum(-X_i beta X_ij) */
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1 contains sum(y_i X_ij - X_i beta X_ij) **/
     gsl_vector_memcpy(term1,vectmp1);/** the remaining term in dg_beta **/
     gsl_vector_scale(term1,(-1.0/n));
     
     /** last row in hessian **/
     for(j=0;j<m;j++){*gsl_matrix_ptr(hessgvaluefull,m,j)=gsl_vector_get(term1,j);}
 
 /** dg_dbeta_dtau is same as dg_dtau_dbeta - note symmetry here **/
     
     /** last col in hessian **/
     for(j=0;j<m;j++){*gsl_matrix_ptr(hessgvaluefull,j,m)=gsl_vector_get(term1,j);}
     
 /** now remove the row and col for the fixed variable */
 /** need to drop a row and drop a col **/
     row=0;
     col=0;
     for(i=0;i<betafull->size;i++){
        for(j=0;j<betafull->size;j++){
       if(i!=betaindex && j!=betaindex){/** unless fixed variable then **/
	 tmp=gsl_matrix_get(hessgvaluefull,i,j);
	 row=i;col=j;
	 if(i>betaindex){row=i-1;} 
	 if(j>betaindex){col=j-1;}
                               gsl_matrix_set(hessgvaluesshort,row,col,tmp);}
	}
       }
 /*Rprintf("end gaus_dgHESS\n");
 for(i=0;i<betafull->size-1;i++){Rprintf("=%f\n",gsl_matrix_get(hessgvaluesshort,i,i));}*/	
 
  return GSL_SUCCESS;
} 
 
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int generate_gaus_inits_marg(gsl_vector *betashort,struct fnparams *gparams){
    
    /** betashort is all the beta's minus the one to be marginalised */
    /** general idea: find all the estimates and then simply drop the one that is fixed */    
    
    /** use least squares estimates to get started **/
    /** beta_hat= (X^T X)^{-1} X^T y **/
    
    const gsl_vector *Y = gparams->Y;/** design matrix **/
       const gsl_matrix *X = gparams->X;/** response variable **/
       gsl_vector *vectmp1= gparams->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = gparams->vectmp2;
       gsl_vector *vectmp1long = gparams->vectmp1long;/** numobs long **/
       gsl_vector *vectmp2long = gparams->vectmp2long;/** numobs long **/

       gsl_matrix *mattmp2 = gparams->mattmp2;/** same dim as X*/
       gsl_matrix *mattmp3 = gparams->mattmp3;/** p x p **/
       gsl_matrix *mattmp4 = gparams->mattmp4;/** p x p **/
       gsl_permutation *perm = gparams->perm;
      double n=Y->size;/** no. observations **/
     double m=X->size2;/** number of coefficients excluding tau-precision */
     int ss,i,j;
     double variance=0.0;
     gsl_vector *beta = gparams->betafull;/** will hold "full beta vector inc precision" **/   
     int betaindex = gparams->betaindex;
     
    /*Rprintf("X: %d %d %d %d %d %d\n",X->size1,X->size2,mattmp2->size1,mattmp2->size2,mattmp3->size1,mattmp3->size2);*/ 
    gsl_matrix_memcpy(mattmp2,X);
    gsl_blas_dgemm (CblasTrans, CblasNoTrans,    /** mattmp3 is p x p matrix X^T X **/
                       1.0, X, mattmp2,
                       0.0, mattmp3);
    gsl_permutation_init(perm);/** reset - might not be needed */                   
    gsl_linalg_LU_decomp(mattmp3,perm,&ss);
    gsl_linalg_LU_invert (mattmp3, perm, mattmp4);/** mattmp4 is now inv (X^T X) */                   
    
    gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1); /** X^T Y */
    gsl_blas_dgemv (CblasNoTrans, 1.0, mattmp4, vectmp1, 0.0, vectmp2); 
    /*Rprintf("rrr\n");*/
   for(i=0;i<beta->size-1;i++){gsl_vector_set(beta,i,gsl_vector_get(vectmp2,i));} /** set to Least squares estimate */
   /*Rprintf("rrr\n");*/
    /** now for variance estimate */
    /** first get y_hat estimate */
    gsl_blas_dgemv (CblasNoTrans, 1.0, X, vectmp2, 0.0, vectmp1long); /** vectmp1 is y_hat */ 
    gsl_vector_scale(vectmp1long,-1.0);/** - y_hat */
    gsl_vector_add(vectmp1long,Y);/** now have Y-y_hat (or -y_hat + Y) */
    gsl_vector_memcpy(vectmp2long,vectmp1long);
    gsl_blas_ddot (vectmp1long, vectmp2long, &variance);/** got sum((Y-Y_hat)^2) */
    variance=variance/(n-m);/** unbiased estimator using denominator n-#term in regression equation **/
    
    gsl_vector_set(beta,beta->size-1,1.0/variance);/** as are using precision parameterization not variance **/ 
    
    /** now copy into smaller vector and drop the variable which is fixed **/
    i=0;
    for(j=0;j<beta->size;j++){
      if(j!=betaindex){gsl_vector_set(betashort,i,gsl_vector_get(beta,j)); i++;}
    } /** set to Least squares estimate */
    
    /*Rprintf("passed varnum %d\n",betaindex);
    for(i=0;i<beta->size;i++){Rprintf("beta=%f\n",gsl_vector_get(beta,i));}
    for(i=0;i<betashort->size;i++){Rprintf("betashort=%f\n",gsl_vector_get(betashort,i));}
    */
    return GSL_SUCCESS;
}   
 

/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/          
int wrapper_gaus_fdf_marg (const gsl_vector *beta, void *gparams,
                     gsl_vector *dgvalues, gsl_matrix *hessgvalues)
     { 
       laplace_gaus_dg_marg(beta, gparams, dgvalues);
       laplace_gaus_hessg_marg(beta, gparams, hessgvalues);
     
       return GSL_SUCCESS;
     }     
