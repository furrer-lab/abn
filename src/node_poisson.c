#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "utility.h"
#include "node_poisson.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

#define PRINTGSL1
/** ****************************************************************************************************
 ***** calc an individual poisson regression model
 *******************************************************************************************************/
void calc_node_Score_pois(network *dag, datamatrix *obsdata, int nodeid,  int verbose,
                                datamatrix *designmatrix, const double priormean, const double priorsd,
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
 gsl_vector *myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,*dgvalue,*term1,*term2,*term3,*vectmp3long;
 gsl_matrix *hessgvalue,*mattmp1,*mattmp2,*mattmp3,*mattmp4;
 struct fnparams gparams;/** for passing to the gsl zero finding functions */
 double gvalue,n,m;
 gsl_permutation *perm=0;
 gsl_permutation *initsperm;
 double inits[2];

 /** build design matrix which is designmatrix->defn, designmatrix->Y plus priors designmatrix->priorsd, designmatrix->priormean **/
 build_designmatrix_pois(dag,obsdata,priormean, priorsd,designmatrix,nodeid,storeModes);

 vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
 vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
 vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
 vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
 myBeta = gsl_vector_alloc (designmatrix->numparams);/** this will hold the parameter point estimates */
 mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
 mattmp3 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
 mattmp4 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
 initsperm = gsl_permutation_alloc (designmatrix->numparams);/** for use with initial guesses */
 term1 = gsl_vector_alloc (designmatrix->numparams);
 term2 = gsl_vector_alloc (designmatrix->numparams);
 term3 = gsl_vector_alloc (designmatrix->numparams);
 dgvalue = gsl_vector_alloc (designmatrix->numparams);/** will hold partial derivates **/
 vectmp3long = gsl_vector_alloc (obsdata->numDataPts);
 mattmp1 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
 hessgvalue = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);/** will hold hessian matrix **/

 /** now store in gparams for sending to gsl solver routines */
   gparams.Y=designmatrix->Y;
   gparams.X=designmatrix->datamatrix;
   gparams.priormean=designmatrix->priormean;
   gparams.priorsd  =designmatrix->priorsd;

   gparams.vectmp1=vectmp1;
   gparams.vectmp2=vectmp2;
   gparams.vectmp1long=vectmp1long;
   gparams.vectmp2long=vectmp2long;

   gparams.mattmp2=mattmp2;
   gparams.mattmp3=mattmp3;
   gparams.mattmp4=mattmp4;
   gparams.perm=initsperm;

   gparams.term1=term1;
   gparams.term2=term2;
   gparams.term3=term3;

   gparams.vectmp3long=vectmp3long;
   gparams.mattmp1=mattmp1;
   inits[0]=0.1;
   inits[1]=1.0;/** hard coded initial guess offsets */

   dag->nodeScoresErrCode[nodeid]=0;/** reset error code to no error **/

   /** now we need to solve system defined in laplace_dg()=0 */
    FDF.f = &laplace_dg_pois;
    FDF.df = &laplace_hessg_pois;
    FDF.fdf = &wrapper_fdf_pois;
    FDF.n = designmatrix->numparams;
    FDF.params = &gparams;

    /** ******************** FIRST TRY for a root using hybridsj  *******************************************************/
    for(i=0;i<2;i++){/** try two different init conditions **/
    gparams.inits_adj=inits[i];
    iter=0;
    T = gsl_multiroot_fdfsolver_hybridj;
    s = gsl_multiroot_fdfsolver_alloc (T, designmatrix->numparams);
    status=GSL_FAILURE;/** just set it to something not equal to GSL_SUCCESS */
    /*status_inits=*/generate_inits_n_pois(myBeta,&gparams);

    gsl_multiroot_fdfsolver_set (s, &FDF, myBeta);

   #ifdef PRINTGSL
   if(nodeid==11){print_state (iter, s);}
   #endif
    iter=0;
       do
         {
           iter++;

           status = gsl_multiroot_fdfsolver_iterate (s);
           #ifdef PRINTGSL
           if(nodeid==11){ print_state (iter, s);}
           #endif
          if (status)
             break;

           status = gsl_multiroot_test_residual (s->f, epsabs);
         }
       while (status == GSL_CONTINUE && iter < maxiters);
      if(status != GSL_SUCCESS){dag->nodeScoresErrCode[nodeid]=1;} else {dag->nodeScoresErrCode[nodeid]=0;}/** reset error code to no error **/
      gsl_vector_memcpy(myBeta,s->x);
      gsl_multiroot_fdfsolver_free(s);
      if(status==GSL_SUCCESS){break;}
    }
   /** ******************** END of FIRST TRY for a root using hybridsj  *******************************************************/
  if(status!=GSL_SUCCESS){/** try other solver **/
    /*Rprintf("poisson: using hybridj failed so re-trying with scaled hybridsj\n");*/
    for(i=0;i<2;i++){/** try two different init conditions **/
    gparams.inits_adj=inits[i];
    iter=0;
    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc (T, designmatrix->numparams);
    status=GSL_FAILURE;/** just set it to something not equal to GSL_SUCCESS */
    /*status_inits=*/generate_inits_n_pois(myBeta,&gparams);
    /*if(nodeid==11){gsl_vector_set(myBeta,0,-1.007);gsl_vector_set(myBeta,0,1.649);gsl_vector_set(myBeta,0,-0.05411);gsl_vector_set(myBeta,0,2.4);gsl_vector_set(myBeta,0,0.008);gsl_vector_set(myBeta,0,-0.05);}
   */
    gsl_multiroot_fdfsolver_set (s, &FDF, myBeta);

   #ifdef PRINTGSL
   if(nodeid==11){print_state (iter, s);}
   #endif
    iter=0;
       do
         {
           iter++;

           status = gsl_multiroot_fdfsolver_iterate (s);
           #ifdef PRINTGSL
           if(nodeid==11){print_state (iter, s);}
           #endif
          if (status)
             break;

           status = gsl_multiroot_test_residual (s->f, epsabs);
         }
       while (status == GSL_CONTINUE && iter < maxiters);
       if(status != GSL_SUCCESS){dag->nodeScoresErrCode[nodeid]=1;} else {dag->nodeScoresErrCode[nodeid]=0;}/** reset error code to no error **/
       gsl_vector_memcpy(myBeta,s->x);
      gsl_multiroot_fdfsolver_free(s);
      if(status==GSL_SUCCESS){break;}
    } /** inits cond loop */
    }
   /** ******************** END of SECOND TRY for a root  *******************************************************/
   if(storeModes){/** keep a copy of the parameter modes found for use later in other function calls etc**/
	 index=0;
		     for(i=0;i<dag->numNodes+1;i++){/** roll myBeta into dag->modes into the appropriate columns**/
		       if(gsl_matrix_get(dag->modes,nodeid,i)!=DBL_MAX){
			 gsl_matrix_set(dag->modes,nodeid,i,gsl_vector_get(myBeta,index++));}}
                   /*for(i=0;i<dag->numNodes;i++){Rprintf("%10.10e ",gsl_matrix_get(dag->modes,nodeid,i));}Rprintf("\n");*/

		   }
  /** we now have all the individual parts so put it together to the laplace approx */
  /*if(status != GSL_SUCCESS){logscore= -DBL_MAX;Rprintf("poisson no root at node %d\n",nodeid+1); *//** root finding failed so discard model by setting fit to worst possible */
 /* } else {*/
  laplace_g_pois(myBeta,&gparams, &gvalue);
  laplace_hessg_pois(myBeta,&gparams, hessgvalue);
   n=obsdata->numDataPts;
   m=designmatrix->numparams;
   perm = gsl_permutation_alloc (m);
   /*status2=*/gsl_linalg_LU_decomp(hessgvalue,perm,&ss);
   /*if(status2 != GSL_SUCCESS){logscore= -DBL_MAX; Rprintf("poisson no inversion at node %d\n",nodeid+1);
   } else {*/
     mydet=gsl_linalg_LU_lndet(hessgvalue);/** compute determinant but this might be a nan - overflow?*/
     /*if(gsl_isnan(mydet)){Rprintf("no det - poisson node at node %d\n",nodeid+1);logscore= -DBL_MAX;*/ /** is nan so return default mlik error value */
     /*} else {*/ /** all is ok so now compute the actual laplace value */

     logscore= -n*gvalue-0.5*mydet+(m/2.0)*log((2.0*M_PI)/n);/*}*/ /** this is the final value */
         /* }*/
     if(gsl_isnan(logscore)){logscore= R_NaN;
                           dag->nodeScoresErrCode[nodeid]=2;}
   /*}*/


   gsl_vector_free(myBeta);
   gsl_vector_free(vectmp1);
   gsl_vector_free(vectmp2);
   gsl_vector_free(vectmp1long);
   gsl_vector_free(vectmp2long);
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(mattmp3);
   gsl_matrix_free(mattmp4);
   gsl_permutation_free(initsperm);
   gsl_vector_free(dgvalue);
   gsl_vector_free(term1);
   gsl_vector_free(term2);
   gsl_vector_free(term3);
   gsl_matrix_free(mattmp1);
   gsl_matrix_free(hessgvalue);
   gsl_vector_free(vectmp3long);

   gsl_vector_free(designmatrix->Y);
   gsl_matrix_free(designmatrix->datamatrix);
   gsl_vector_free(designmatrix->priormean);
   gsl_vector_free(designmatrix->priorsd);

   /*if(status == GSL_SUCCESS){*/gsl_permutation_free(perm);/*}*/ /** only allocate this is status==GSL_SUCCESS */
   /*gsl_multiroot_fdfsolver_free (s);*/

   /*return(logscore);*/
   dag->nodeScores[nodeid]=logscore;

}
/** ****************************************************************************************************
 ***** marginal distribution for logistic regression model
 *******************************************************************************************************/
void calc_poisson_marginal(network *dag, datamatrix *obsdata, int nodeid,  int verbose,
                                datamatrix *designmatrix, const double priormean, const double priorsd,
                                const int maxiters, const double epsabs, double *denom_modes,int paramid, double betafixed, double mlik, double *posterior)
{

 int i,j,ss,status/*,status2,status_inits,index*/;
 /*unsigned int k,j;*/
 int iter=0;
 /*unsigned int numparents=0;*/
 double logscore=0.0;
 const gsl_multiroot_fdfsolver_type *T;
 gsl_multiroot_fdfsolver *s;
 gsl_multiroot_function_fdf FDF;
 gsl_vector *myBeta,*vectmp1,*vectmp2,*vectmp1long,*vectmp2long,*dgvalue,*term1,*term2,*term3,*vectmp3long, *betafull;
 gsl_matrix *hessgvalue,*mattmp1,*mattmp2,*mattmp3,*mattmp4,*hessgvaluefull;
 struct fnparams gparams;/** for passing to the gsl zero finding functions */
 double gvalue,n,m;
 gsl_permutation *perm=0;
 gsl_permutation *initsperm;
 double val=0.0;

 /** build design matrix which is designmatrix->defn, designmatrix->Y plus priors designmatrix->priorsd, designmatrix->priormean **/
 build_designmatrix_pois(dag,obsdata,priormean, priorsd,designmatrix,nodeid,0);
 /** down to here is as for the network score calc which is an integral across all parameters - we now adjust this to that its across all parameters
      minus one, where this one is fixed at values across a grid **/
  /** SPECIAL CASE if only a model with a single parameter then no integration required just evaluation of (-1/n)*g() **/

  /** ********************************************************************************************************************/
  switch(designmatrix->numparams){
    case 1:{/** only a constant term **/
    vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
    vectmp2long = gsl_vector_alloc (obsdata->numDataPts);

   /** now store in gparams for sending to gsl solver routines */
   gparams.Y=designmatrix->Y;
   gparams.X=designmatrix->datamatrix;
   gparams.priormean=designmatrix->priormean;
   gparams.priorsd  =designmatrix->priorsd;

   gparams.vectmp1=vectmp1;
   gparams.vectmp2=vectmp2;
   gparams.vectmp1long=vectmp1long;
   gparams.vectmp2long=vectmp2long;

   myBeta = gsl_vector_alloc (designmatrix->numparams);
   n=obsdata->numDataPts;

      gsl_vector_set(myBeta,0,betafixed);
      laplace_g_pois(myBeta,&gparams, &gvalue);
      logscore= -n*gvalue;
      val=exp(logscore-mlik);
    /* Rprintf("got betafixed=%f mlik=%f and value=%f\n",betafixed,mlik,val); */
    *posterior=val;

    /** Last Step before return - free all the gsl vectors, matrices, other etc **/
   gsl_vector_free(designmatrix->Y);
   gsl_matrix_free(designmatrix->datamatrix);
   gsl_vector_free(designmatrix->priormean);
   gsl_vector_free(designmatrix->priorsd);

   gsl_vector_free(myBeta);
   gsl_vector_free(vectmp1);
   gsl_vector_free(vectmp2);
   gsl_vector_free(vectmp1long);
   gsl_vector_free(vectmp2long);

   break;}

   default:{
     /** allocate only once here since same dimension within one node for the marginals**/
  /** GENERAL IDEA - keep the same dimensions as in the full margLik calc but drop off terms at the end if needed */
    vectmp1 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp2 = gsl_vector_alloc (designmatrix->numparams);/** scratch space **/
    vectmp1long = gsl_vector_alloc (obsdata->numDataPts);/** scratch space **/
    vectmp2long = gsl_vector_alloc (obsdata->numDataPts);
    mattmp2 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    mattmp3 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    mattmp4 = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);
    initsperm = gsl_permutation_alloc (designmatrix->numparams);/** for use with initial guesses */
    term1 = gsl_vector_alloc (designmatrix->numparams);
    term2 = gsl_vector_alloc (designmatrix->numparams);
    term3 = gsl_vector_alloc (designmatrix->numparams);
    dgvalue = gsl_vector_alloc (designmatrix->numparams);/** will hold partial derivates **/
    vectmp3long = gsl_vector_alloc (obsdata->numDataPts);
    mattmp1 = gsl_matrix_alloc (obsdata->numDataPts,designmatrix->numparams);
    betafull = gsl_vector_alloc (designmatrix->numparams);/** this will hold the re-build full vector of all parameters */
    hessgvaluefull = gsl_matrix_alloc (designmatrix->numparams,designmatrix->numparams);/**  will hold hessian matrix **/

    n=obsdata->numDataPts;
    m=designmatrix->numparams-1;/** IMPORTANT: -1 since now a marginal calculation **/
    FDF.f = &laplace_dg_pois_marg;
    FDF.df = &laplace_hessg_pois_marg;
    FDF.fdf = &wrapper_fdf_pois_marg;
    FDF.n = designmatrix->numparams-1;
    FDF.params = &gparams;
    myBeta = gsl_vector_alloc (designmatrix->numparams-1);/** this will hold the parameter point estimates */
    hessgvalue = gsl_matrix_alloc (designmatrix->numparams-1,designmatrix->numparams-1);/**  IMPORTANT: -1 since now a marginal calculation will hold hessian matrix **/
    T = gsl_multiroot_fdfsolver_hybridj;
    s = gsl_multiroot_fdfsolver_alloc (T, designmatrix->numparams-1);
    perm = gsl_permutation_alloc (m);

 /** now store in gparams for sending to gsl solver routines */
    /** now send */
   gparams.Y=designmatrix->Y;
   gparams.X=designmatrix->datamatrix;
   gparams.vectmp1=vectmp1;
   gparams.vectmp2=vectmp2;
   gparams.vectmp1long=vectmp1long;
   gparams.vectmp2long=vectmp2long;
   gparams.vectmp3long=vectmp3long;
   gparams.term1=term1;
   gparams.term2=term2;
   gparams.term3=term3;
   gparams.priormean=designmatrix->priormean;
   gparams.priorsd  =designmatrix->priorsd;
   gparams.mattmp1=mattmp1;
   gparams.mattmp2=mattmp2;
   gparams.mattmp3=mattmp3;
   gparams.mattmp4=mattmp4;
   gparams.perm=initsperm;
   gparams.betafull=betafull;
   gparams.hessgvalues=hessgvaluefull;
   gparams.betafixed=0.0;/** these will be changed in loop below*/
   gparams.betaindex=paramid;/** this is fixed - the variable for which the posterior is calculated **/

   gparams.inits_adj=0.1;  /* added RF 21042022, addressing valgrind issue */


   /** generate initial estimates for the remaining variable - not the posterior variable **/
   generate_inits_n_pois(myBeta,&gparams);

      j=0;
      for(i=0;i<designmatrix->numparams;i++){if(i!= paramid){gsl_vector_set(myBeta,j++,denom_modes[i]);}} /** use modes as initial values ignoring current mode which is to be marginalised **/
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

       if(status != GSL_SUCCESS){       if (verbose>0) {
			Rprintf ("Zero finding error: status = %s at x=%f\n", gsl_strerror (status),gparams.betafixed);/*exit(1);*/}
       }
      gsl_vector_memcpy(myBeta,s->x);
      /** we now have all the individual parts so put it together to the laplace approx */
      laplace_g_pois_marg(myBeta,&gparams, &gvalue);
      laplace_hessg_pois_marg(myBeta,&gparams, hessgvalue);
      gsl_linalg_LU_decomp(hessgvalue,perm,&ss);
      logscore= -n*gvalue-0.5*gsl_linalg_LU_lndet(hessgvalue)+(m/2.0)*log((2.0*M_PI)/n); /** this is the final value */
      val=exp(logscore-mlik);
      /*Rprintf("got betafixed=%f mlik=%f and value=%f\n",betafixed,mlik,val);*/
   *posterior=val;
   /*** Last Step before return - free all the gsl vectors, matrices, other etc **/

   gsl_vector_free(myBeta);
   gsl_vector_free(vectmp1);
   gsl_vector_free(vectmp2);
   gsl_vector_free(vectmp1long);
   gsl_vector_free(vectmp2long);
   gsl_vector_free(term1);
   gsl_vector_free(term2);
   gsl_vector_free(term3);
   gsl_vector_free(vectmp3long);

   gsl_vector_free(betafull);
   gsl_matrix_free(hessgvalue);
   gsl_matrix_free(mattmp1);
   gsl_matrix_free(mattmp2);
   gsl_matrix_free(mattmp3);
   gsl_matrix_free(mattmp4);
   gsl_vector_free(dgvalue);

   gsl_matrix_free(hessgvaluefull);
   gsl_permutation_free(perm);
   gsl_permutation_free(initsperm);
   gsl_multiroot_fdfsolver_free (s);

   gsl_vector_free(designmatrix->Y);
   gsl_matrix_free(designmatrix->datamatrix);
   gsl_vector_free(designmatrix->priormean);
   gsl_vector_free(designmatrix->priorsd);

   } /** end of default case */

  } /** end of switch control */

}

/** **************************************************************************************************************/
/** **************************************************************************************************************/
/** build the design matrix - plus other associated things *******************************************************/
void build_designmatrix_pois(network *dag,datamatrix *obsdata, double priormean, double priorsd,datamatrix *designmatrix, int nodeid, int storeModes)
{

 int i,j,k;
 int numparents=0;
 gsl_vector *Y,*vecpriormean,*vecpriorsd;
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
 } /** maxparent=0 check */
  /** this part is new and just for posterior param est - it does not affect Laplace approx in any way****/
  /** setup matrix where each non DBL_MAX entry in a row is for a parameter to be estimated and the col is which param
      first col is for the intercept */
  if(storeModes){
    for(k=0;k<dag->numNodes+3;k++){gsl_matrix_set(dag->modes,nodeid,k,DBL_MAX);} /** initialise row to DBL_MAX n.b. +2 here is need in fitabn.R part**/
    gsl_matrix_set(dag->modes,nodeid,0,1);/** the intercept term - always have an intercept - but not in dag.m definition */
    for(k=0;k<numparents;k++){gsl_matrix_set(dag->modes,nodeid,gsl_vector_int_get(parentindexes,k)+1,1);} /** offset is 1 due to intercept */
  }
  /** ****************************************************************************************************/

  datamat=gsl_matrix_alloc(obsdata->numDataPts,numparents+1);
  designmatrix->datamatrix=datamat;
  Y=gsl_vector_alloc(obsdata->numDataPts);
  designmatrix->Y=Y;
  vecpriormean=gsl_vector_alloc(numparents+1);
  designmatrix->priormean=vecpriormean;
  vecpriorsd=gsl_vector_alloc(numparents+1);
  designmatrix->priorsd=vecpriorsd;

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

   designmatrix->numparams=numparents+1;/** +1 for intercept**/
   /** now set the priormean and priorsd vector - choose the correct prior values */

   for(k=0;k<designmatrix->numparams;k++){
                                          gsl_vector_set(designmatrix->priormean,k,priormean);
                                          gsl_vector_set(designmatrix->priorsd,k,priorsd);
   }

   gsl_vector_int_free(parentindexes);/** finished with this **/

 /*    Rprintf("##########################################\n");
  Rprintf("got %d parents\n",numparents);
for(i=0;i<obsdata->numDataPts;i++){
   Rprintf("Y=%f\t",gsl_vector_get(designmatrix->Y,i));
   for(k=0;k<numparents+1;k++){Rprintf("%f\t",gsl_matrix_get(designmatrix->datamatrix,i,k));}Rprintf("\n");}
   */


}

/** ***************************************************************************************************************
*******************************************************************************************************************
** laplace method = int^b_a exp(-lambda g(y)) h(y) dy = exp(-lambda g(y*)) h(y*) (2PI/lambda)^(d/2) det(hess)^(1/2)
** lambda = sample size n, g(y) = -(1/n)* log( f(D|betas)f(betas) ) e.g. -(1/n)* log (like*prior)
*******************************************************************************************************************
******************************************************************************************************************/

/** **************************************************************************************************************/
/** g(y) = -(1/n)* log( f(D|betas)f(betas) */
/** **************************************************************************************************************/
int laplace_g_pois (const gsl_vector *beta, void *params,double *gvalue)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;
	gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       double n=Y->size;/** no. observations **/
       double m=X->size2;
       double term1=0;
       double term2=0;
       double term3=0;
       double storedbl1;
       unsigned int i=0;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */

     /** DO IN THREE PARTS - term1, term2, term3 */
     /** R code "term2<-sum( log(1/(sqrt(2*pi)*sd.loc)) );" **/
     term2=0; for(i=0;i<m;i++){term2+=-log(sqrt(2.0*M_PI)*gsl_vector_get(priorsd,i));}

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

     /** Rcode  Y%*%(X%*%mybeta)-sum( log(Y!) + exp(X%*%mybeta) );  */
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/
     gsl_blas_ddot (Y, vectmp1long, &storedbl1);/** storedbl1 holds Y%*%(X%*%mybeta)**/
     term1+=storedbl1;

     for(i=0;i<vectmp1long->size;i++){
       gsl_vector_set(vectmp2long,i,-gsl_sf_lnfact(gsl_vector_get(Y,i))-exp(gsl_vector_get(vectmp1long,i)));
                               } /** vectmp2 holds -log(Y!)-exp(X%*%mybeta)) */
     gsl_vector_set_all(vectmp1long,1.0); /** ones vector */
     gsl_blas_ddot (vectmp2long, vectmp1long, &storedbl1);/** DOT product simply to calc -sum( log(Y!)+exp(X%*%mybeta))) */
     term1+=storedbl1;

     *gvalue=(-1.0/n)*(term1+term2+term3);

     /** Rprintf("Poisson raw loglike=%f\n", term1);**/
       return GSL_SUCCESS;
     }

/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */
/** **************************************************************************************************************/
int laplace_dg_pois (const gsl_vector *beta, void *params, gsl_vector *dgvalues)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
        gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;/** numobs long **/
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       gsl_vector *term2 = ((struct fnparams *) params)->term2;
       gsl_vector *term3 = ((struct fnparams *) params)->term3;
       double n=Y->size;/** no. observations **/

       unsigned int i=0;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */

     /** DO IN THREE PARTS - term1, term2, term3 */
      /** term1 (beta_j - mu_j)/sd_j^2" **/
     gsl_vector_memcpy(vectmp1,beta);/** copy beta to temp vec */
     gsl_vector_memcpy(vectmp2,priormean);
     gsl_vector_scale(vectmp2,-1.0);
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1= beta-mean**/
     gsl_vector_memcpy(vectmp2,priorsd);
     gsl_vector_mul(vectmp2,priorsd);/** square all elements in priorsd and store in vectmp2 */
     gsl_vector_div(vectmp1,vectmp2);
     gsl_vector_scale(vectmp1,-1.0);
     gsl_vector_memcpy(term1,vectmp1);

     /**     */
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/
     for(i=0;i<vectmp1long->size;i++){
       gsl_vector_set(vectmp2long,i,-exp(gsl_vector_get(vectmp1long,i)) );
                                      } /** vectmp2long holds -exp(X%*%mybeta) */

     /*Rprintf("=%d %d %d %d\n",X->size1, X->size2, beta->size,vectmp2long->size);*/
     gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp2long, 0.0, vectmp1);/** vectmp1 hold X%*% -exp(X%*%mybeta) **/
     gsl_vector_memcpy(term2,vectmp1);

     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);/** sum(Y_i*X_ij) over all j - one entry per j */
     gsl_vector_memcpy(term3,vectmp1);

     /*Rprintf("==%f %f %f\n",gsl_vector_get(term1,0),gsl_vector_get(term2,0),gsl_vector_get(term3,0));*/

     gsl_vector_add(term1,term2);/** add term 2 to term 1 */
     gsl_vector_add(term1,term3);/** add term 3 to term 1 */
     gsl_vector_scale(term1,-1.0/n);

     gsl_vector_memcpy(dgvalues,term1);


 return GSL_SUCCESS;
     }

/** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */
/** **************************************************************************************************************/
int laplace_hessg_pois (const gsl_vector *beta, void *params, gsl_matrix *hessgvalues)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
        gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;/** numobs long **/
        gsl_vector *vectmp3long = ((struct fnparams *) params)->vectmp3long;/** numobs long **/
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       gsl_vector *term2 = ((struct fnparams *) params)->term2;
       gsl_matrix *mattmp1 = ((struct fnparams *) params)->mattmp1;
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** no. params to estimate*/

       unsigned int i=0;unsigned int j=0;unsigned int k=0;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */

     /** do in multiple parts - need to do element operations first */
     /** first exp(Xb) */
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/
     for(i=0;i<vectmp1long->size;i++){
	   gsl_vector_set(vectmp2long,i,exp(gsl_vector_get(vectmp1long,i))/n); /** vettmp2long hold exp(X%*%mybeta)/n */
                                      } /** vectmp2long holds the main complicated term */

     gsl_matrix_memcpy(mattmp1,X);/** make a copy of X*/
     gsl_matrix_mul_elements (mattmp1, X);/** calc X^2 is in mattmp1*/

     /*Rprintf("=%d %d %d %d\n",X->size1, X->size2, beta->size,vectmp2long->size);*/
     gsl_blas_dgemv (CblasTrans, 1.0, mattmp1, vectmp2long, 0.0, vectmp1);/** vecttmp2long hold Xij^2*complicated **/
     gsl_vector_memcpy(term1,vectmp1);


     gsl_vector_set_all(term2,0.0); /** zeros vector */

     gsl_vector_memcpy(vectmp1,priorsd);/** copy priorsd in vectmp1 **/
     gsl_vector_mul(vectmp1,priorsd);/** square priorsd */
     gsl_vector_scale(vectmp1,n);/** now have n*sigma^2 **/
     gsl_vector_set_all(vectmp2,1.0); /** ones vector */
     gsl_vector_div(vectmp2,vectmp1);/** get 1/(n*sigma^2) into vectmp2 **/
     gsl_vector_add(term2,vectmp2);/** add to term2*/

     gsl_vector_add(term1,term2);

     /*Rprintf("hess[1,1] at beta=%5.10f is %5.10f\n",gsl_vector_get(beta,0),gsl_vector_get(term1,0)); */


     /** STILL TO DO OFF DIAGONAL ELEMENTS - check for triangular?*/
     for(j=0;j<m;j++){
       for(k=0;k<m;k++){
                    if(j!=k){/** dealt with j==k case above */
                          /** NOTE - vectmp2long is the SAME as in the j==k case so can use this directly **/
                          /** need X[,j]*X[,k] - element wise mult **/
                            gsl_matrix_get_col(vectmp1long,X,j); /** get col j in X **/
                            gsl_matrix_get_col(vectmp3long,X,k); /** get col k in X **/
                            gsl_vector_mul(vectmp1long,vectmp3long); /** element by element multiplication - result in vecttmp1long **/

                            gsl_blas_ddot (vectmp1long, vectmp2long, gsl_matrix_ptr(hessgvalues,j,k));/** DOT product simply to calcu sum value */
                    } else {*gsl_matrix_ptr(hessgvalues,j,k)=gsl_vector_get(term1,j);}
                     }
                     }

 return GSL_SUCCESS;
     }

/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/
int generate_inits_n_pois(gsl_vector *myBeta,struct fnparams *gparams){

    /** this is the SAME CODE as in the Gaussian case  */

    /** beta_hat= (X^T X)^{-1} X^T y **/

       const gsl_vector *Y = gparams->Y;/** design matrix **/
       const gsl_matrix *X = gparams->X;/** response variable **/
       gsl_vector *vectmp1= gparams->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = gparams->vectmp2;
       gsl_matrix *mattmp2 = gparams->mattmp2;/** same dim as X*/
       gsl_matrix *mattmp3 = gparams->mattmp3;/** p x p **/
       gsl_matrix *mattmp4 = gparams->mattmp4;/** p x p **/
       gsl_vector *vectmp1long = gparams->vectmp1long;/** scratch space **/
       gsl_permutation *perm = gparams->perm;
       // double adj=gparams->inits_adj;/** help with initial guess - an offset */

     unsigned int i;
     int ss;
     int haveError;

    /*Rprintf("X: %d %d %d %d %d %d\n",X->size1,X->size2,mattmp2->size1,mattmp2->size2,mattmp3->size1,mattmp3->size2); */
    gsl_matrix_memcpy(mattmp2,X);
    gsl_blas_dgemm (CblasTrans, CblasNoTrans,    /** mattmp3 is p x p matrix X^T X **/
                       1.0, X, mattmp2,
                       0.0, mattmp3);
    gsl_permutation_init(perm);/** reset - might not be needed */
    gsl_linalg_LU_decomp(mattmp3,perm,&ss);
    gsl_set_error_handler_off();/**Turning off GSL Error handler as this may fail as mattmp3 may be singular */
    haveError=gsl_linalg_LU_invert (mattmp3, perm, mattmp4);/** mattmp4 is now inv (X^T X) */
    if(!haveError){/** no error */
      /** copy Y into vectmp1long and +1 and take logs since poisson has log link - this is a fudge */

      /* 042022 valgrind temp try
      Rprintf("-> got in inits=%f %f %d\n",gsl_vector_get(Y, 0), adj, vectmp1long->size);
      */

      /*orig    for(i=0;i<Y->size;i++){gsl_vector_set(vectmp1long,i,log(gsl_vector_get(Y,i)+adj));} */ /** NOTE -  +1.0 or 0.1 give different reliablity*/
      for(i=0;i<Y->size;i++){gsl_vector_set(vectmp1long,i,log(gsl_vector_get(Y,i)+0.1));} /** NOTE -  +1.0 or 0.1 give different reliablity*/
    /*    for(i=0;i<vectmp1long->size;i++){gsl_vector_set(vectmp1long,i,log(gsl_vector_get(Y,i)+adj));} * NOTE -  +1.0 or 0.1 give different reliablity*/


    gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp1long, 0.0, vectmp1); /** X^T Y */
    gsl_blas_dgemv (CblasNoTrans, 1.0, mattmp4, vectmp1, 0.0, vectmp2);
             for(i=0;i<myBeta->size;i++){
	                   /*if(gsl_vector_get(vectmp2,i)>0.0){*/
			   gsl_vector_set(myBeta,i,gsl_vector_get(vectmp2,i));
			   /*} else {Rprintf("got in inits=%f %d\n",gsl_vector_get(vectmp2,i),i);
			           error("mmmm....this should not happen...negative initial guess in generate_init_n_pois()\n");}
			  */ }
    } else {
      Rprintf ("caught gsl error - singular matrix in initial guess estimates\n");/** singular to set initial values all to zero **/
      for(i=0;i<myBeta->size;i++){gsl_vector_set(myBeta,i,0.0);}}
    gsl_set_error_handler (NULL);/** restore the error handler*/
    /*Rprintf("inits\n");for(i=0;i<myBeta->size;i++){Rprintf("%10.15e ",gsl_vector_get(myBeta,i));} Rprintf("\n");*//** set to Least squares estimate */


    return GSL_SUCCESS;
}
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/
int wrapper_fdf_pois (const gsl_vector *beta, void *gparams,
                     gsl_vector *dgvalues, gsl_matrix *hessgvalues)
     {
       laplace_dg_pois(beta, gparams, dgvalues);
       laplace_hessg_pois(beta, gparams, hessgvalues);

       return GSL_SUCCESS;
     }
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/
/** **************************************************************************************************************/
/** g(y) = -(1/n)* log( f(D|betas)f(betas) */
/** **************************************************************************************************************/
int laplace_g_pois_marg (const gsl_vector *betashort, void *params,double *gvalue)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
       gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;
       gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
       gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;
       gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       double n=Y->size;/** no. observations **/
       double m=X->size2;
        /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *beta = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
       double term1=0;
       double term2=0;
       double term3=0;
       double storedbl1;
       int i=0;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */

      /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(beta,0,betafixed);
                     for(i=1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(beta->size-1)){gsl_vector_set(beta,beta->size-1,betafixed);
                     for(i=0;i<beta->size-1;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}}

     if(betaindex>0 && betaindex<(beta->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(beta,betaindex,betafixed);
	 for(i=betaindex+1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}
     }

     /** DO IN THREE PARTS - term1, term2, term3 */
     /** R code "term2<-sum( log(1/(sqrt(2*pi)*sd.loc)) );" **/
     term2=0; for(i=0;i<m;i++){term2+=-log(sqrt(2.0*M_PI)*gsl_vector_get(priorsd,i));}
    /** DO IN THREE PARTS - term1, term2, term3 */
     /** R code "term2<-sum( log(1/(sqrt(2*pi)*sd.loc)) );" **/
     term2=0; for(i=0;i<m;i++){term2+=-log(sqrt(2.0*M_PI)*gsl_vector_get(priorsd,i));}

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

     /** Rcode  Y%*%(X%*%mybeta)-sum( log(Y!) + exp(X%*%mybeta) );  */
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/
     gsl_blas_ddot (Y, vectmp1long, &storedbl1);/** storedbl1 holds Y%*%(X%*%mybeta)**/
     term1+=storedbl1;

     for(i=0;i<vectmp1long->size;i++){
       gsl_vector_set(vectmp2long,i,-gsl_sf_lnfact(gsl_vector_get(Y,i))-exp(gsl_vector_get(vectmp1long,i)));
                               } /** vectmp2 holds -log(Y!)-exp(X%*%mybeta)) */
     gsl_vector_set_all(vectmp1long,1.0); /** ones vector */
     gsl_blas_ddot (vectmp2long, vectmp1long, &storedbl1);/** DOT product simply to calc -sum( log(Y!)+exp(X%*%mybeta))) */
     term1+=storedbl1;

     *gvalue=(-1.0/n)*(term1+term2+term3);

       return GSL_SUCCESS;

     }

 /** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */
/** **************************************************************************************************************/
int laplace_dg_pois_marg (const gsl_vector *betashort, void *params, gsl_vector *dgvaluesshort)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
        gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;/** numobs long **/
       const gsl_vector *priormean = ((struct fnparams *) params)->priormean;
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       gsl_vector *term2 = ((struct fnparams *) params)->term2;
       gsl_vector *term3 = ((struct fnparams *) params)->term3;
       /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *beta = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
       double n=Y->size;/** no. observations **/

       int i=0; double tmp=0;int col;
       /** beta are the parameters values at which the function is to be evaluated **/
     /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(beta,0,betafixed);
                     for(i=1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(beta->size-1)){gsl_vector_set(beta,beta->size-1,betafixed);
                     for(i=0;i<beta->size-1;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}}


     if(betaindex>0 && betaindex<(beta->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(beta,betaindex,betafixed);
	 for(i=betaindex+1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}
     }

     /** DO IN THREE PARTS - term1, term2, term3 */
      /** term1 (beta_j - mu_j)/sd_j^2" **/
     gsl_vector_memcpy(vectmp1,beta);/** copy beta to temp vec */
     gsl_vector_memcpy(vectmp2,priormean);
     gsl_vector_scale(vectmp2,-1.0);
     gsl_vector_add(vectmp1,vectmp2);/** vectmp1= beta-mean**/
     gsl_vector_memcpy(vectmp2,priorsd);
     gsl_vector_mul(vectmp2,priorsd);/** square all elements in priorsd and store in vectmp2 */
     gsl_vector_div(vectmp1,vectmp2);
     gsl_vector_scale(vectmp1,-1.0);
     gsl_vector_memcpy(term1,vectmp1);

     /**     */
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/
     for(i=0;i<vectmp1long->size;i++){
       gsl_vector_set(vectmp2long,i,-exp(gsl_vector_get(vectmp1long,i)) );
                                      } /** vectmp2long holds -exp(X%*%mybeta) */

     /*Rprintf("=%d %d %d %d\n",X->size1, X->size2, beta->size,vectmp2long->size);*/
     gsl_blas_dgemv (CblasTrans, 1.0, X, vectmp2long, 0.0, vectmp1);/** vectmp1 hold X%*% -exp(X%*%mybeta) **/
     gsl_vector_memcpy(term2,vectmp1);

     gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);/** sum(Y_i*X_ij) over all j - one entry per j */
     gsl_vector_memcpy(term3,vectmp1);

     /*Rprintf("==%f %f %f\n",gsl_vector_get(term1,0),gsl_vector_get(term2,0),gsl_vector_get(term3,0));*/

     gsl_vector_add(term1,term2);/** add term 2 to term 1 */
     gsl_vector_add(term1,term3);/** add term 3 to term 1 */
     gsl_vector_scale(term1,-1.0/n);

     /** need to drop one cell in term1 before copying back */
     /** create an adjusted term1 which contains the term1 without the  re-inserted at the correct place **/
    col=0;
     for(i=0;i<beta->size;i++){
       if(i!=betaindex){/** unless fixed variable then **/
	 tmp=gsl_vector_get(term1,i);
	 col=i;
	 if(i>betaindex){col=i-1;}
                               gsl_vector_set(dgvaluesshort,col,tmp);}
	}


 return GSL_SUCCESS;
     }

  /** **************************************************************************************************************/
/** partial_g(y)/partial_beta vector of first derivatives                                                        */
/** **************************************************************************************************************/
int laplace_hessg_pois_marg (const gsl_vector *betashort, void *params, gsl_matrix *hessgvaluesshort)
{
       const gsl_vector *Y = ((struct fnparams *) params)->Y;/** design matrix **/
       const gsl_matrix *X = ((struct fnparams *) params)->X;/** response variable **/
        gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
        gsl_vector *vectmp2 = ((struct fnparams *) params)->vectmp2;
        gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;/** numobs long **/
        gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;/** numobs long **/
        gsl_vector *vectmp3long = ((struct fnparams *) params)->vectmp3long;/** numobs long **/
       /*const gsl_vector *priormean = ((struct fnparams *) params)->priormean;*/
       const gsl_vector *priorsd   = ((struct fnparams *) params)->priorsd;
       gsl_vector *term1 = ((struct fnparams *) params)->term1;
       gsl_vector *term2 = ((struct fnparams *) params)->term2;
       gsl_matrix *mattmp1 = ((struct fnparams *) params)->mattmp1;
       /** this is extra stuff to deal with the fixed beta **/
       gsl_vector *beta = ((struct fnparams *) params)->betafull;/** will hold "full beta vector" **/
       double betafixed = ((struct fnparams *) params)->betafixed;/** the fixed beta value passed through**/
       int betaindex = ((struct fnparams *) params)->betaindex;
       gsl_matrix *hessgvalues = ((struct fnparams *) params)->hessgvalues;
       double n=Y->size;/** no. observations **/
       double m=X->size2;/** no. params to estimate*/

       int i=0;int j=0;int k=0;int row,col;double tmp;
       /** beta are the parameters values at which the function is to be evaluated **/
       /** gvalue is the return value - a single double */

       /** create an adjusted beta which contains the FIXED beta re-inserted at the correct place **/
     if(betaindex==0){gsl_vector_set(beta,0,betafixed);
                     for(i=1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}}
     if(betaindex==(beta->size-1)){gsl_vector_set(beta,beta->size-1,betafixed);
                     for(i=0;i<beta->size-1;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}}


     if(betaindex>0 && betaindex<(beta->size-1)){
         for(i=0;i<betaindex;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i));}
         gsl_vector_set(beta,betaindex,betafixed);
	 for(i=betaindex+1;i<beta->size;i++){gsl_vector_set(beta,i,gsl_vector_get(betashort,i-1));}
     }

    /** do in multiple parts - need to do element operations first */
     /** first exp(Xb) */
     gsl_blas_dgemv (CblasNoTrans, 1.0, X, beta, 0.0, vectmp1long);/** vectmp1long hold X%*%mybeta **/
     for(i=0;i<vectmp1long->size;i++){
	   gsl_vector_set(vectmp2long,i,exp(gsl_vector_get(vectmp1long,i))/n); /** vettmp2long hold exp(X%*%mybeta)/n */
                                      } /** vectmp2long holds the main complicated term */

     gsl_matrix_memcpy(mattmp1,X);/** make a copy of X*/
     gsl_matrix_mul_elements (mattmp1, X);/** calc X^2 is in mattmp1*/

     /*Rprintf("=%d %d %d %d\n",X->size1, X->size2, beta->size,vectmp2long->size);*/
     gsl_blas_dgemv (CblasTrans, 1.0, mattmp1, vectmp2long, 0.0, vectmp1);/** vecttmp2long hold Xij^2*complicated **/
     gsl_vector_memcpy(term1,vectmp1);


     gsl_vector_set_all(term2,0.0); /** zeros vector */

     gsl_vector_memcpy(vectmp1,priorsd);/** copy priorsd in vectmp1 **/
     gsl_vector_mul(vectmp1,priorsd);/** square priorsd */
     gsl_vector_scale(vectmp1,n);/** now have n*sigma^2 **/
     gsl_vector_set_all(vectmp2,1.0); /** ones vector */
     gsl_vector_div(vectmp2,vectmp1);/** get 1/(n*sigma^2) into vectmp2 **/
     gsl_vector_add(term2,vectmp2);/** add to term2*/

     gsl_vector_add(term1,term2);

     /*Rprintf("hess[1,1] at beta=%5.10f is %5.10f\n",gsl_vector_get(beta,0),gsl_vector_get(term1,0)); */


     /** STILL TO DO OFF DIAGONAL ELEMENTS - check for triangular?*/
     for(j=0;j<m;j++){
       for(k=0;k<m;k++){
                    if(j!=k){/** dealt with j==k case above */
                          /** NOTE - vectmp2long is the SAME as in the j==k case so can use this directly **/
                          /** need X[,j]*X[,k] - element wise mult **/
                            gsl_matrix_get_col(vectmp1long,X,j); /** get col j in X **/
                            gsl_matrix_get_col(vectmp3long,X,k); /** get col k in X **/
                            gsl_vector_mul(vectmp1long,vectmp3long); /** element by element multiplication - result in vecttmp1long **/

                            gsl_blas_ddot (vectmp1long, vectmp2long, gsl_matrix_ptr(hessgvalues,j,k));/** DOT product simply to calcu sum value */
                    } else {*gsl_matrix_ptr(hessgvalues,j,k)=gsl_vector_get(term1,j);}
                     }
                     }

     /** need to drop a row and drop a col **/
     row=0;
     col=0;
     for(i=0;i<beta->size;i++){
        for(j=0;j<beta->size;j++){
       if(i!=betaindex && j!=betaindex){/** unless fixed variable then **/
	 tmp=gsl_matrix_get(hessgvalues,i,j);
	 row=i;col=j;
	 if(i>betaindex){row=i-1;}
	 if(j>betaindex){col=j-1;}
                               gsl_matrix_set(hessgvaluesshort,row,col,tmp);}
	}
       }


 return GSL_SUCCESS;
     }

 /** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/
int wrapper_fdf_pois_marg (const gsl_vector *beta, void *gparams,
                     gsl_vector *dgvalues, gsl_matrix *hessgvalues)
     {
       laplace_dg_pois_marg(beta, gparams, dgvalues);
       laplace_hessg_pois_marg(beta, gparams, hessgvalues);

       return GSL_SUCCESS;
     }
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/
