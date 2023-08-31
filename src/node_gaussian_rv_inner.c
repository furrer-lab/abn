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
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "node_gaussian_rv_inner.h"
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
#define PRINTGSL1
#define NOPRIOR1 /**turn on for removing effect of priors - this is purely to enable a check of loglike values with lme4 */
/** ****************************************************************************************************
 ***** calc an individual logistic regression model
 *******************************************************************************************************/

/** **************************************************************************************************************/
/** **************************************************************************************************************/
/** build the design matrix - plus other associated things *******************************************************/
/** NOTE:key feature of this design matrix is that the last column is full of 1's and this is for the random effect   */
/**    epsilon term so it can easily be included in matrix operations. this is different from the non-rv case    */
/** **************************************************************************************************************/
void build_designmatrix_gaus_rv(network *dag,datamatrix *obsdata, double priormean, double priorsd,const double priorgamshape, const double priorgamscale,datamatrix *designmatrix, int nodeid, int storeModes)
{

 int i,j,k;
 int numparents=0;
 gsl_vector_int *parentindexes=0;
 int num_unq_grps=0;
 int *groupcnts;
 int *curindex;
 gsl_matrix **array_of_designs;
 gsl_vector **array_of_Y;

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
 } /** check for maxparent=0 */
  /** this part is new and just for posterior param est - it does not affect Laplace approx in any way****/
  /** setup matrix where each non DBL_MAX entry in a row is for a parameter to be estimated and the col is which param
      first col is for the intercept */
 if(storeModes){
    for(k=0;k<dag->numNodes+3;k++){gsl_matrix_set(dag->modes,nodeid,k,DBL_MAX);} /** initialise row to DBL_MAX n.b. +2 here is need in fitabn.R part**/
    gsl_matrix_set(dag->modes,nodeid,0,1);/** the intercept term - always have an intercept - but not in dag.m definition */
    for(k=0;k<numparents;k++){gsl_matrix_set(dag->modes,nodeid,gsl_vector_int_get(parentindexes,k)+1,1);} /** offset is 1 due to intercept */
    gsl_matrix_set(dag->modes,nodeid,dag->numNodes+1,1);/** the residual precision **/
    gsl_matrix_set(dag->modes,nodeid,dag->numNodes+2,1);/** the group level precision term put at end of other params */
 }
  /** ****************************************************************************************************/

  designmatrix->datamatrix=gsl_matrix_alloc(obsdata->numDataPts,numparents+1+1);/** +1=intercept +1=rv_precision - note this is just for the mean so no extra term for gaussian node **/
  designmatrix->Y=gsl_vector_alloc(obsdata->numDataPts);
  designmatrix->priormean=gsl_vector_alloc(numparents+1);
  designmatrix->priorsd=gsl_vector_alloc(numparents+1);
  designmatrix->priorgamshape=gsl_vector_alloc(1); /** only 1 of these per node - NOTE: use same prior for group precision and overall precision */
  designmatrix->priorgamscale=gsl_vector_alloc(1); /** only 1 of these per node - NOTE: use same prior for group precision and overall precision */

  designmatrix->datamatrix_noRV=gsl_matrix_alloc(obsdata->numDataPts,numparents+1);/** drop the last col - used for initial value estimation only**/

  /** create design matrix - ALL DATA POINTS - copy relevant cols from the observed data **/
 /** int** designmatrix is just used as storage space, fill up from left cols across until as far as needed */
 for(i=0;i<obsdata->numDataPts;i++){/** for each observed data point **/
   gsl_matrix_set(designmatrix->datamatrix,i,0,1.0); /** set first column - intercept -  to 1's **/
   gsl_matrix_set(designmatrix->datamatrix_noRV,i,0,1.0);/** build matrix same as datamatrix just without the last col (which contains 1's for epsilon rv term */

   gsl_matrix_set(designmatrix->datamatrix,i,(designmatrix->datamatrix)->size2-1,1.0);/** set last column - rv precision to 1.0 **/

   gsl_vector_set(designmatrix->Y,i,obsdata->defn[i][nodeid]);/** copy values at node - response values - into vector Y */

   for(k=0;k<numparents;k++){/** now build design matrix of explanatories other than intercept*/

     gsl_matrix_set(designmatrix->datamatrix,i,k+1,obsdata->defn[i][gsl_vector_int_get(parentindexes,k)]);
     gsl_matrix_set(designmatrix->datamatrix_noRV,i,k+1,obsdata->defn[i][gsl_vector_int_get(parentindexes,k)]);
                            } /** end of explanatories **/
   } /** end of data point loop */

   designmatrix->numparams=numparents+1;/** +1 for intercept - excludes precisions **/
   /** now set the priormean and priorsd vector */
   for(k=0;k<designmatrix->numparams;k++){/** num params does NOT include precision term **/
                                          gsl_vector_set(designmatrix->priormean,k,priormean);
                                          gsl_vector_set(designmatrix->priorsd,k,priorsd);
   }
   /** set prior for precision **/
   gsl_vector_set(designmatrix->priorgamshape,0,priorgamshape);/** prior for precision term */
   gsl_vector_set(designmatrix->priorgamscale,0,priorgamscale);/** prior for precision term */

   gsl_vector_int_free(parentindexes);/** finished with this **/

   /** ***********************************************************************************************************************/
   /** ***********************************************************************************************************************/
   /** DOWN HERE is splitting the single single design matrix and Y into separate chunks *************************************/
   /** we now want to split designmatrix->datamatrix and designmatrix->Y into grouped blocks                                **/
   /** ***********************************************************************************************************************/

   /** get number of unique groups - equal to max int since using R factors **/
   num_unq_grps=0;for(i=0;i<obsdata->numDataPts;i++){if(obsdata->groupIDs[i]>num_unq_grps){num_unq_grps=obsdata->groupIDs[i];}}
   groupcnts=(int *)R_alloc(num_unq_grps,sizeof(int));/** will hold n_j's e.g. counts of how many obs in each group **/
   curindex=(int *)R_alloc(num_unq_grps,sizeof(int));
   for(i=0;i<num_unq_grps;i++){groupcnts[i]=0;curindex[i]=0;}

   for(i=0;i<num_unq_grps;i++){/** for each unique group of data **/
     for(j=0;j<obsdata->numDataPts;j++){/** for each observation **/
         if( (obsdata->groupIDs[j]-1)==i){groupcnts[i]++;/** increment count **/}
     }
   }

  /** create an array of gsl_matrix where each one is the design matrix for a single group of data **/
  array_of_designs=(gsl_matrix **)R_alloc(num_unq_grps,sizeof(gsl_matrix*));/** a list of design matrix, one for each group */
  array_of_Y=(gsl_vector **)R_alloc(num_unq_grps,sizeof(gsl_vector*)); /** a list of Y vectors,, one for each group */
  for(i=0;i<num_unq_grps;i++){array_of_designs[i]=gsl_matrix_alloc(groupcnts[i],(designmatrix->datamatrix)->size2);
                              array_of_Y[i]=gsl_vector_alloc(groupcnts[i]);}


  /** now loop through group j; for fixed j loop through each record in total design matrix and copying group members into new group design matrix **/
  for(j=0;j<num_unq_grps;j++){/** for each group **/
    for(i=0;i<obsdata->numDataPts;i++){/** for each data point **/
	   if( (obsdata->groupIDs[i]-1)==j){/** if current data point is for group j then store **/
	     for(k=0;k<(designmatrix->datamatrix)->size2;k++){/** for each member of the row in design matrix **/
             gsl_matrix_set(array_of_designs[j],curindex[j],k,gsl_matrix_get(designmatrix->datamatrix,i,k));}
             gsl_vector_set(array_of_Y[j],curindex[j],gsl_vector_get(designmatrix->Y,i));/** copy relevant Ys for fixed j */
	     curindex[j]++;
	   }
    }
  }

  /** uncomment to print out array of design matrices - one for each data group **/
 /* Rprintf("no cols=%d\n",array_of_designs[0]->size2);
   for(j=0;j<num_unq_grps;j++){Rprintf("-------group %d------\n",j);
     for(i=0;i<array_of_designs[j]->size1;i++){
       Rprintf("Y=%f\t",gsl_vector_get(array_of_Y[j],i));
       for(k=0;k<array_of_designs[j]->size2;k++){
       Rprintf("%f ",gsl_matrix_get(array_of_designs[j],i,k));
       }
       Rprintf("\n");
     }
   }
 */
 /** down to here we now have a split the design matrix and Y up into separate matrices and vectors, one for each observational group */
 /** so we can free the previous datamatrix as this is not needed, also the previous Y **/
 gsl_matrix_free(designmatrix->datamatrix);
 /*gsl_vector_free(designmatrix->Y);*/ /** need to keep y*/

 /** copy addresses */
 designmatrix->numUnqGrps=num_unq_grps;
 designmatrix->array_of_designs=array_of_designs;
 designmatrix->array_of_Y=array_of_Y;



}

/** *******************************************************************************************************************************************/
/** *******************************************************************************************************************************************/
int generate_gaus_rv_inits(gsl_vector *myBeta,struct fnparams *gparams){

    /** this is the SAME CODE as in the Gaussian case  */

    /** beta_hat= (X^T X)^{-1} X^T y **/
    const datamatrix *designdata = ((struct fnparams *) gparams)->designdata;/** all design data inc Y and priors **/

       const gsl_vector *Y = designdata->Y;/** response vector **/
       const gsl_matrix *X = designdata->datamatrix_noRV ;/** design matrix - with one too few cols! **/
       gsl_vector *vectmp1= gparams->vectmp1;/** numparams long*/
       gsl_vector *vectmp2 = gparams->vectmp2;/** numparams long*/
       gsl_matrix *mattmp2 = gparams->mattmp2;/** same dim as X*/
       gsl_matrix *mattmp3 = gparams->mattmp3;/** p x p **/
       gsl_matrix *mattmp4 = gparams->mattmp4;/** p x p **/
       gsl_vector *vectmp1long = gparams->vectmp1long;/** scratch space **/
       gsl_vector *vectmp2long = gparams->vectmp2long;/** scratch space **/
       gsl_permutation *perm = gparams->perm;
     unsigned int i;
     int ss;
     int haveError;
     double variance=0.0;
     double n=Y->size;/** no. observations **/
     double m=X->size2;/** number of coefficients excluding tau-precision */

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
      /*for(i=0;i<vectmp1long->size;i++){gsl_vector_set(vectmp1long,i,log(gsl_vector_get(Y,i)+DBL_MIN)/(log(1-gsl_vector_get(Y,i)+DBL_MIN)));}  */
    /*for(i=0;i<vectmp1long->size;i++){gsl_vector_set(vectmp1long,i,log(gsl_vector_get(Y,i)+1)/(log(1-gsl_vector_get(Y,i)+1)));} */

    gsl_blas_dgemv (CblasTrans, 1.0, X, Y, 0.0, vectmp1);  /** X^T Y */
    gsl_blas_dgemv (CblasNoTrans, 1.0, mattmp4, vectmp1, 0.0, vectmp2);
    for(i=0;i<myBeta->size-2;i++){gsl_vector_set(myBeta,i,gsl_vector_get(vectmp2,i));} /** size myBeta->size-2 as last two entries are precisions **/
    } else {/** singular to set initial values all to zero **/
           Rprintf("caught gsl error - singular matrix in initial guess estimates\n");
           for(i=0;i<myBeta->size;i++){gsl_vector_set(myBeta,i,0.01);}}
    gsl_set_error_handler (NULL);/** restore the error handler*/
   /*Rprintf("inits\n");for(i=0;i<myBeta->size;i++){Rprintf("%10.15e ",gsl_vector_get(myBeta,i));} Rprintf("\n");*//** set to Least squares estimate */
     /** now for variance estimate */
    /** first get y_hat estimate */

    gsl_blas_dgemv (CblasNoTrans, 1.0, X, vectmp2, 0.0, vectmp1long); /** vectmp1 is y_hat */
    /*for(i=0;i<vectmp1long->size;i++){Rprintf("y_hat=%f\n",gsl_vector_get(vectmp1long,i));}*/
    /*error("");*/
    gsl_vector_scale(vectmp1long,-1.0);/** - y_hat */
    gsl_vector_add(vectmp1long,Y);/** now have Y-y_hat (or -y_hat + Y) */
    /*for(i=0;i<vectmp1long->size;i++){gsl_vector_set(vectmp1long,i,fabs(gsl_vector_get(vectmp1long,i)));}
    for(i=0;i<vectmp1long->size;i++){Rprintf("y_hat=%f\n",gsl_vector_get(vectmp1long,i));}
    for(i=0;i<vectmp1long->size;i++){gsl_vector_set(vectmp1long,i,log(gsl_vector_get(vectmp1long,i))/log(1-gsl_vector_get(vectmp1long,i)));}*/ /** errors on logit scale **/
    /*gsl_vector_set_all(vectmp2long,1);*/
    gsl_vector_memcpy(vectmp2long,vectmp1long);
    gsl_blas_ddot (vectmp1long, vectmp2long, &variance);/** got sum((Y-Y_hat)^2) */
    variance=variance/(n-m);/** unbiased estimator using denominator n-#term in regression equation **/
   /* Rprintf("variance estimator=%f precision=%f\n",variance,1/variance);*/
  /* variance=0.086;*/
    /*variance=exp(gsl_vector_get(myBeta,0))/(1+exp(gsl_vector_get(myBeta,0)));*/
    /** variance here is for plain glm - just split this 50/50 between residual error and group level variance **/

    gsl_vector_set(myBeta,myBeta->size-2,1.0/(0.5*variance));/** - estimate for rv precision **/
    gsl_vector_set(myBeta,myBeta->size-1,1.0/(0.5*variance));/** - estimate for residual precision **/
    /*gsl_vector_set(myBeta,myBeta->size-1,1.0/variance); */
    /*gsl_vector_set(myBeta,0,0.9 );gsl_vector_set(myBeta,1,0.9);gsl_vector_set(myBeta,2,1.5);*/

    #ifdef junk
    Rprintf("------------ TEMP: Using Fixed initial values from LME4----------------\n");
    gsl_vector_set(myBeta,0,0.062044233);/** intercept */
    gsl_vector_set(myBeta,1,-0.1229382094322);/** slope g2 */
    gsl_vector_set(myBeta,2,1.0/0.1570366587829);/** group level precision */
    gsl_vector_set(myBeta,3,1.0/0.8628565204966);/** residual precision */
    #endif
    /*Rprintf("inits\n");for(i=0;i<myBeta->size;i++){Rprintf("%10.15e ",gsl_vector_get(myBeta,i));} Rprintf("\n");*//** set to Least squares estimate */




    return GSL_SUCCESS;
}
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/
double g_inner_gaus( gsl_vector *beta, const datamatrix *designdata, int groupid, double epsabs, int maxiters, int verbose){

   /** this function perform a Laplace approx on a single data group given fixed beta, so only integrate over single term epsilon **/

/* const gsl_multiroot_fdfsolver_type *T;
 gsl_multiroot_fdfsolver *s;
 gsl_multiroot_function_fdf FDF;*/
 struct fnparams gparams;/** for passing to the gsl zero finding functions */
 /*double epsilon=0;*//** the variable we want to find the root of **/
 gsl_vector *epsilon = gsl_vector_alloc (1);
 gsl_vector *dgvalues = gsl_vector_alloc (1);
 gsl_matrix *hessgvalue = gsl_matrix_alloc (1,1);
 /*int iter=0;*/
 /*int status;*/
 /*double epsabs=1e-5;
 int maxiters=100;*/
 /*int verbose=1;*/

 gsl_vector *vectmp1 = gsl_vector_alloc (designdata->numparams+1);/** scratch space same length as number of params inc precision **/
 gsl_vector *vectmp1long = gsl_vector_alloc ( ((designdata->array_of_Y)[groupid])->size);/** scratch space same length as number of obs in group j**/
 gsl_vector *vectmp2long = gsl_vector_alloc ( ((designdata->array_of_Y)[groupid])->size);
 double logscore;
 double gvalue;int n,m;

 /*for(i=0;i<beta->size;i++){Rprintf("g_inner_gaus=%f\n",gsl_vector_get(beta,i));}*/

 /*Rprintf("I HAVE epsabs_inner=%f maxiters_inner=%d verbose=%d\n",epsabs,maxiters,verbose);*/

 /*
 FDF.f = &rv_dg_inner_gaus;
 FDF.df = &rv_hessg_inner_gaus;
 FDF.fdf = &wrapper_rv_fdf_inner_gaus;
 FDF.n = 1;
 FDF.params = &gparams;
 */
 gparams.Y=designdata->array_of_Y[groupid];
 gparams.X=designdata->array_of_designs[groupid];
 gparams.beta=beta;/** inc group and residual precision **/
 /*Rprintf("tau in g_inner=%f\n",gsl_vector_get(beta,beta->size-1));
 if(gsl_vector_get(beta,beta->size-1)<0.0){Rprintf("got negative tau!!=%f\n",gsl_vector_get(beta,beta->size-1));error("");}*/
 gparams.vectmp1=vectmp1;/** same length as beta but used as scratch space */
 gparams.vectmp1long=vectmp1long;
 gparams.vectmp2long=vectmp2long;

  /** ******************** FIRST TRY for a root using hybridsj  *******************************************************/
 #ifdef NO
  iter=0;
    /*T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc (T);*/
    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc (T, 1);
    status=GSL_FAILURE;/** just set it to something not equal to GSL_SUCCESS */
    /*status_inits=generate_inits_rv_n(x,&gparams);*/
    gsl_vector_set(epsilon,0,0.0);/** initial guess */
    /*gsl_root_fdfsolver_set (s, &FDF, epsilon);*/
   gsl_multiroot_fdfsolver_set (s, &FDF, epsilon);

    /*Rprintf ("using %s method\n",
               gsl_root_fdfsolver_name (s));

       Rprintf ("%-5s %10s %10s %10s\n",
               "iter", "root", "err", "err(est)");
   */

   /*print_state (iter, s);*/

    iter=0;
       do
         {
           iter++;

           status = gsl_multiroot_fdfsolver_iterate (s);

           /*print_state (iter, s);*/

          if (status)
             break;

           status = gsl_multiroot_test_residual (s->f, epsabs);
         }
       while (status == GSL_CONTINUE && iter < maxiters);
       if( status != GSL_SUCCESS){Rprintf ("Zero finding warning: internal--- epsilon status = %s\n", gsl_strerror (status));
                                  /*for(i=0;i<s->x->size;i++){Rprintf("0epsilon=%f ",gsl_vector_get(s->x,i));}Rprintf("\n");*/}
       gsl_vector_memcpy(epsilon,s->x);
       Rprintf("modes: %f\n",gsl_vector_get(epsilon,0));
     gsl_multiroot_fdfsolver_free(s);
    /*Rprintf("x=%5.10f f=%5.10f\n",gsl_root_fdfsolver_root(s),rv_dg_inner(gsl_root_fdfsolver_root(s),&gparams));*/

  /* if(status != GSL_SUCCESS){*//*error("no root\n");*//*Rprintf("binary no root at node %d\n",groupid+1);*//*logscore= DBL_MAX;*/ /** root finding failed so discard model by setting fit to worst possible */
  /*} else {*/
    /*gsl_vector_set(epsilon,0,0.3);*/
 #endif

  rv_dg_inner_gaus(epsilon,&gparams, dgvalues);/** value is returned in dgvalues - first entry **/
  gsl_vector_memcpy(epsilon,dgvalues);/** copy value dgvalues into epsilon */
  /*Rprintf("mode for epsilon=%f\n",gsl_vector_get(epsilon,0));*/

  rv_g_inner_gaus(epsilon,&gparams, &gvalue);/*Rprintf("==>g()=%e %f tau=%f\n",gvalue,gsl_vector_get(epsilon,0),gsl_vector_get(beta,2));*/
  /*if(status != GSL_SUCCESS){Rprintf("1epsilon=%f %f\n",gsl_vector_get(epsilon,0), gvalue);}*/
  rv_hessg_inner_gaus(epsilon,&gparams, hessgvalue);

                  /* Rprintf("node=%d hessian at g\n",nodeid+1);
		   for(j=0;j<myBeta->size;j++){Rprintf("%f ",gsl_vector_get(myBeta,j));}Rprintf("\n");
                   for(j=0;j<hessgvalue->size1;j++){
                   for(k=0;k<hessgvalue->size2;k++){Rprintf("%f ",gsl_matrix_get(hessgvalue,j,k));} Rprintf("\n");}*/
   /*Rprintf("epsilon=%f\n",epsilon);*/
   n=((designdata->array_of_designs)[groupid])->size1;/** number of obs in group */
   m=1;/** number of params */
   /*Rprintf("gvalue in g_inner=|%f| n=|%d| |%f|\n",gvalue,n,-n*gvalue);*/
   /*if(status != GSL_SUCCESS){Rprintf("2epsilon=%f %f\n",gsl_vector_get(epsilon,0), gvalue);}*/
   logscore= -n*gvalue-0.5*log(gsl_matrix_get(hessgvalue,0,0))+(m/2.0)*log((2.0*M_PI)/n); /** this is the final value */
   if(gsl_isnan(logscore)){error("nan in g_inner hessmat=%f epsilon=%f gvalue=%f\n",gsl_matrix_get(hessgvalue,0,0),gsl_vector_get(epsilon,0),gvalue);}
     /*}*/

   /*Rprintf("group=%d logscore=%f\n",groupid+1,logscore);*/


      gsl_vector_free(dgvalues);
      gsl_vector_free(epsilon);
      gsl_matrix_free(hessgvalue);
      gsl_vector_free(vectmp1);
      gsl_vector_free(vectmp1long);
      gsl_vector_free(vectmp2long);

      return(logscore);



}
/** *************************************************************************************
*****************************************************************************************/
/** Note this function does not evaluate the derivate but rather give the turning point
solution directly **/
/** ***************************************************************************************/
int rv_dg_inner_gaus (const gsl_vector *epsilonvec, void *params, gsl_vector *dgvalues)
{

   /*double epsilon=0.3; */
   double epsilon=gsl_vector_get(epsilonvec,0);
   const gsl_vector *Y = ((struct fnparams *) params)->Y;/** response variable **/
   const gsl_matrix *X = ((struct fnparams *) params)->X;/** design matrix INC epsilon col **/
   const gsl_vector *beta = ((struct fnparams *) params)->beta;/** fixed covariate and precision terms **/
   gsl_vector *vectmp1 = ((struct fnparams *) params)->vectmp1;
   gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;
   gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;

   double tau_rv = gsl_vector_get(beta,beta->size-2);/** inc the precision terms - second last entries */
   double tau_resid = gsl_vector_get(beta,beta->size-1);/** last entry - residual precision */
   double n = (double)(Y->size);/** number of observations */
   int i;

   double /*term3,*/term2;

   /*term3 = (tau_rv*epsilon)/n;*/ /** correct sign */

   /** now for the more complex term */
   /** the design matrix does not include precisions but does include epsilon, beta includes precisions but not epsilon. To use matrix operations
       we make a copy of beta and replace one precision value with value for epsilon - copy into vectmp1 */
   for(i=0;i<beta->size-2;i++){gsl_vector_set(vectmp1,i,gsl_vector_get(beta,i));} /** copy **/
   gsl_vector_set(vectmp1,beta->size-2,epsilon); /** last entry in vectmp1 is not precision but epsilon **/

   /*for(i=0;i<vectmp1->size;i++){Rprintf("=>%f\n",gsl_vector_get(vectmp1,i));} */

   /** get X%*%beta where beta = (b0,b1,...,epsilon) and so we get a vector of b0*1+b1*x1i+b2*x2i+epsilon*1 for each obs i */
    gsl_blas_dgemv (CblasNoTrans, 1.0, X, vectmp1, 0.0, vectmp1long);/** vectmp1long hold X%*%vectmp1 = X%*%mybeta **/
    /*for(i=0;i<vectmp1long->size;i++){Rprintf("=%f\n",gsl_vector_get(vectmp1long,i));}*/

    /*Rprintf("---\n");for(i=0;i<X->size1;i++){for(j=0;j<X->size2;j++){Rprintf("%f ",gsl_matrix_get(X,i,j));}Rprintf("\n");}Rprintf("---\n");*/

   /*for(i=0;i<vectmp2long->size;i++){Rprintf(">%f\n",gsl_vector_get(vectmp2long,i));}*/

   gsl_vector_scale(vectmp1long,-1.0);/** multiple each entry by -1 **/
   gsl_vector_memcpy(vectmp2long,Y);/** vectmp2long becomes Y **/
   gsl_vector_add(vectmp2long,vectmp1long);/** vectmp2long becomes Y-XB **/
   gsl_vector_set_all(vectmp1long,1.0);/** reset each value to unity **/
   gsl_blas_ddot (vectmp2long, vectmp1long, &term2);/** just to get the sum of vectmp2long */

  /* Rprintf("analytical solution=%f\n",(tau_resid*term2)/(tau_resid*n + tau_rv));*/

   /** This derivative can be solved analytically so no need to return value of d_g/d_epsilon
   term2 *= -tau_resid/n;

     gsl_vector_set(dgvalues,0,term2+term3);
     if(gsl_isnan(gsl_vector_get(dgvalues,0))){error("rv_dg_inner is nan %f %f %f\n",term2,term3);}
   **/

   gsl_vector_set(dgvalues,0,(tau_resid*term2)/(tau_resid*n + tau_rv)); /** solves dg/d_epsilon=0 */

 return GSL_SUCCESS;
}


/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/
int rv_hessg_inner_gaus (const gsl_vector *epsilonvec, void *params,gsl_matrix *hessgvalues)
{

   /*double epsilon=0.3;*/
   /*double epsilon=gsl_vector_get(epsilonvec,0);*/
   const gsl_vector *Y = ((struct fnparams *) params)->Y;/** response variable **/

   const gsl_vector *beta = ((struct fnparams *) params)->beta;/** fixed covariate and precision terms **/

   double tau_rv = gsl_vector_get(beta,beta->size-2);/** inc the precision terms - second last entries */
   double tau_resid = gsl_vector_get(beta,beta->size-1);/** last entry - residual precision */
   double n = (double)(Y->size);/** number of observations */

   gsl_matrix_set(hessgvalues,0,0,(tau_rv/n)+tau_resid);

   if(gsl_isnan(gsl_matrix_get(hessgvalues,0,0))){error("rv_hess_inner_gaus is nan\n");}

   return GSL_SUCCESS;

}

/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/
int wrapper_rv_fdf_inner_gaus (const gsl_vector *beta, void *gparams,gsl_vector *dgvalues, gsl_matrix *hessgvalues)
     {
       rv_dg_inner_gaus(beta, gparams, dgvalues);
       rv_hessg_inner_gaus(beta, gparams, hessgvalues);

       return GSL_SUCCESS;
     }

/** **************************************************************************************************************/
/** **************************************************************************************************************/
int rv_g_inner_gaus (const gsl_vector *epsilonvec, void *params, double *gvalue)
{

  double epsilon=gsl_vector_get(epsilonvec,0);
   const gsl_vector *Y = ((struct fnparams *) params)->Y;/** response variable **/
   const gsl_matrix *X = ((struct fnparams *) params)->X;/** design matrix INC epsilon col **/
   const gsl_vector *beta = ((struct fnparams *) params)->beta;/** fixed covariate and precision terms **/
   gsl_vector *vectmp1 = ((struct fnparams *) params)->vectmp1;
   gsl_vector *vectmp1long = ((struct fnparams *) params)->vectmp1long;
   gsl_vector *vectmp2long = ((struct fnparams *) params)->vectmp2long;

   double tau_rv = gsl_vector_get(beta,beta->size-2);/** inc the precision terms - second last entries */
   double tau_resid = gsl_vector_get(beta,beta->size-1);/** last entry - residual precision */
   double n = (double)(Y->size);/** number of observations */
   int i;

   double term1,term2;


   /** easy terms collected together - no Y,X, or betas **/
   term1 = (n/2.0)*log(tau_resid/(2.0*M_PI)) - (tau_rv/2.0)*epsilon*epsilon + 0.5*log(tau_rv/(2.0*M_PI));


   /** now for the more complex term */
   /** the design matrix does not include precisions but does include epsilon, beta includes precisions but not epsilon. To use matrix operations
       we make a copy of beta and replace one precision value with value for epsilon - copy into vectmp1 */
   for(i=0;i<beta->size-2;i++){gsl_vector_set(vectmp1,i,gsl_vector_get(beta,i));} /** copy **/
   gsl_vector_set(vectmp1,beta->size-2,epsilon); /** last entry in vectmp1 is not precision but epsilon **/

   /*for(i=0;i<vectmp1->size;i++){Rprintf("=>%f\n",gsl_vector_get(vectmp1,i));} */

   /** get X%*%beta where beta = (b0,b1,...,epsilon) and so we get a vector of b0*1+b1*x1i+b2*x2i+epsilon*1 for each obs i */
    gsl_blas_dgemv (CblasNoTrans, 1.0, X, vectmp1, 0.0, vectmp1long);/** vectmp1long hold X%*%vectmp1 = X%*%mybeta **/
    /*for(i=0;i<vectmp1long->size;i++){Rprintf("=%f\n",gsl_vector_get(vectmp1long,i));}*/

    /*Rprintf("---\n");for(i=0;i<X->size1;i++){for(j=0;j<X->size2;j++){Rprintf("%f ",gsl_matrix_get(X,i,j));}Rprintf("\n");}Rprintf("---\n");*/

   /*for(i=0;i<vectmp2long->size;i++){Rprintf(">%f\n",gsl_vector_get(vectmp2long,i));}*/

   gsl_vector_scale(vectmp1long,-1.0);/** multiple each entry by -1 **/
   gsl_vector_memcpy(vectmp2long,Y);/** vectmp2long becomes Y **/
   gsl_vector_add(vectmp2long,vectmp1long);/** vectmp2long becomes Y-XB **/

   /*for(i=0;i<vectmp2long->size;i++){Rprintf("> %f\n",gsl_vector_get(vectmp2long,i));}*/

   /** need sum of (Y-XB)^2 so just do a dot product **/
   gsl_vector_memcpy(vectmp1long,vectmp2long);/** copy vectmp2long into vectmp1long */
   gsl_blas_ddot (vectmp2long, vectmp1long, &term2);/** just to get the sum of (Y-XB)^2 */
   term2 *= -(tau_resid/2.0);

   /*Rprintf("term2=%f epsilon=%f tau_resid=%f\n",term2,epsilon,tau_resid);*/

  *gvalue = (-1.0/n)*(term1 + term2);
   /*Rprintf("\n----value of term1 %f %f %f----\n",((storedbl1+storedbl2)*(-1/n)),term2,term3); */
  if(gsl_isnan(*gvalue)){error("\n oops - got an NAN! in g_rv_g_inner_gaus-----\n");}

  return GSL_SUCCESS;
}
/** **************************************************************************************************************/
/** **************************************************************************************************************/
double g_outer_gaus_single (double x, void *params)
{
  int i,j;
  double term1=0.0;
  const datamatrix *designdata = ((struct fnparams *) params)->designdata;/** all design data inc Y and priors **/
  gsl_vector *betaincTau   = ((struct fnparams *) params)->betaincTau;/** include precision */
  int fixed_beta =((struct fnparams *) params)->fixed_index;/** which parameter is to be treated as fixed */
  const gsl_vector *priormean = designdata->priormean;
  const gsl_vector *priorsd   = designdata->priorsd;
  const gsl_vector *priorgamshape   = designdata->priorgamshape;
  const gsl_vector *priorgamscale   = designdata->priorgamscale;
  gsl_vector *beta   = ((struct fnparams *) params)->beta;/** does not include precision */
  gsl_vector *vectmp1= ((struct fnparams *) params)->vectmp1;/** numparams long*/
  gsl_vector *vectmp2 =((struct fnparams *) params)->vectmp2;/** numparams long*/
  double epsabs_inner=((struct fnparams *) params)->epsabs_inner;/** absolute error in internal laplace est */
  int maxiters_inner=((struct fnparams *) params)->maxiters_inner;/** number of steps for inner root finder */
  int verbose=((struct fnparams *) params)->verbose;/**  */

  int n_betas= (designdata->datamatrix_noRV)->size2;/** number of mean terms excl rv and precision **/
  int n=(designdata->datamatrix_noRV)->size1;/** total number of obs **/

  double term2=0.0,term3=0.0,term4=0.0,gval=0.0, term5=0.0;
  /*Rprintf("%d %d\n",n_betas,betaincTau->size);*/
  double tau_rv,tau_resid, copyBeta=0.0;

  /** need to replace variable fixed_beta with x **/
  copyBeta=gsl_vector_get(betaincTau,fixed_beta);/** store value so can reset later */
  gsl_vector_set(betaincTau,fixed_beta,x);

  tau_rv=gsl_vector_get(betaincTau,betaincTau->size-2);/** extract the tau-precision from *beta - last entry */
  /*Rprintf("g_outer_rv tau=%f\n",tau_rv);*/
  tau_resid=gsl_vector_get(betaincTau,betaincTau->size-1);/** extract the tau-precision from *beta - last entry */
  /*Rprintf("g_outer_resid tau=%f\n",tau_resid);*/

  if(tau_rv<=0.0){/*Rprintf("tau_rv negative=%e in g_outer_gaus_single!\n",tau_rv);*/
                 /** aborting so re-copy value of beta changed back to what it was since passed by memory **/
		 /** this might legitimately occur when trying to use a central difference - which is caught by testing for isnan */
                 gsl_vector_set(betaincTau,fixed_beta,copyBeta);
                 return(GSL_NAN);
		 /*error("");*/}

  if(tau_resid<=0.0){/*Rprintf("tau_resid negative=%e in g_outer_gaus_single!\n",tau_resid);*/
                    /** aborting so re-copy value of beta changed back to what it was since passed by memory **/
		    /** this might legitimately occur when trying to use a central difference - which is caught by testing for isnan */
                    gsl_vector_set(betaincTau,fixed_beta,copyBeta);
                    return(GSL_NAN);
                    /*error("");*/}

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

	  term1+= g_inner_gaus(betaincTau,designdata,j, epsabs_inner,maxiters_inner,verbose);
      }

/*Rprintf("term1 in g_outer=%f\n",term1);*/
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

  /** finally re-copy value of beta changed back to what it was since passed by memory **/
  gsl_vector_set(betaincTau,fixed_beta,copyBeta);

  if(gsl_isnan(gval)){error("g_outer_gaus_single\n");}

/*Rprintf("g_outer_final=%f term1=%f term2=%f term3=%f term4=%f term5=%f total=%f %d\n",gval,term1,term2,term3,term4,term5,term1+term2+term3+term4,n);*/
	return(gval);/** negative since its a minimiser */
}
