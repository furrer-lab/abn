#include <R.h>
#include <Rdefines.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "structs.h"
#include "utility.h"
#include "fit_single_node.h"
#include "node_binomial.h"
#include "node_gaussian.h"
#include "node_poisson.h"
#include "node_binomial_rv.h"
#include "node_gaussian_rv.h"
#include "node_poisson_rv.h"
#include <time.h>

#define DEBUG_12

/** ******************************************************************************************************************************************************/
/** compute all the node scores to create a cache ********************************************************************************************************/
/** ******************************************************************************************************************************************************/
SEXP fit_single_node(SEXP R_obsdata, SEXP R_child, SEXP R_parents, SEXP R_numVars, SEXP R_vartype, SEXP R_maxparents,
		     SEXP R_priors_mean, SEXP R_priors_sd,SEXP R_priors_gamshape,SEXP R_priors_gamscale,
		     SEXP R_maxiters, SEXP R_epsabs, SEXP R_verbose, SEXP R_errorverbose, SEXP R_trace,
		     SEXP R_groupedvars, SEXP R_groupids, SEXP R_epsabs_inner,SEXP R_maxiters_inner,
		     SEXP R_finitestepsize, SEXP R_hparams, SEXP R_maxiters_hessian, SEXP R_ModesONLY,
		     SEXP R_max_hessian_error,SEXP R_myfactor_brent, SEXP R_maxiters_hessian_brent, SEXP R_num_intervals_brent)
{

/** ****************/
/** declarations **/
unsigned int i,k;
 int errverbose,verbose,trace;
datamatrix data,designmatrix;
const double priormean=asReal(R_priors_mean);/*Rprintf("priormean=%f %f\n",priormean[0],priormean[5]);*/
const double priorsd=asReal(R_priors_sd);/*Rprintf("priorsd=%f %f\n",priorsd[0],priorsd[5]);*/
const double priorgamshape=asReal(R_priors_gamshape);  /*Rprintf("priorgamshape=%f %f\n",priorgamshape[0],priorgamshape[1]);*/
const double priorgamscale=asReal(R_priors_gamscale);  /*Rprintf("priorgamscale=%f %f\n",priorgamscale[0],priorgamscale[1]);*/
int numVars=asInteger(R_numVars);
int ModesONLY=asInteger(R_ModesONLY);
int storeModes=1;/** store the modes **/
double max_hessian_error=asReal(R_max_hessian_error);
double myfactor_brent=asReal(R_myfactor_brent);
int maxiters_hessian_brent=asInteger(R_maxiters_hessian_brent);
double num_intervals_brent=asReal(R_num_intervals_brent);
/*Rprintf("vartype: ");for(i=0;i<LENGTH(R_vartype);i++){Rprintf("%u ",vartype[i]);} Rprintf("\n");*/
network dag;
/*storage nodescore;*/
SEXP listresults;
SEXP tmplistentry;
int curnode=0;

/** end of declarations*/
/** *******************/
/** ********************************************/
/** parse function arguments - R data structs **/
/** get the data we require from the passed arguments e.g. how many variables/nodes to we have **/
const int maxiters=asInteger(R_maxiters);
const double epsabs=asReal(R_epsabs);
int maxiters_inner=asInteger(R_maxiters_inner);
int maxiters_hessian=asInteger(R_maxiters_hessian);
double epsabs_inner=asReal(R_epsabs_inner);
const int maxparents=asInteger(R_maxparents);
double finitestepsize=asReal(R_finitestepsize);
/*double h_lowerend=REAL(R_hparams)[0];
double h_upperend=REAL(R_hparams)[1];*/
double h_guess=REAL(R_hparams)[0];
double h_epsabs=REAL(R_hparams)[1];
verbose=asInteger(R_verbose);
errverbose=asInteger(R_errorverbose);
trace=asInteger(R_trace);

// If R_groupids is not a vector, then fail gracefully
if (!isVector(R_groupids)) {
  error("R_groupids is not a vector");
} else if (LENGTH(R_groupids) != LENGTH(VECTOR_ELT(R_obsdata,0))) {
  // If R_groupids is a vector, then check that it is of the correct length
  error("R_groupids is not the same length as R_obsdata");
} else {
  // all is well
  // Rprintf("R_groupids is a vector of the correct length\n");
}
/** end of argument parsing **/
/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                                */
/** generic code to create a list comprising of vectors of type double
   - currently overkill but useful template **/

PROTECT(listresults = allocVector(VECSXP, 1));
for(i=0;i<1;i++){
				PROTECT(tmplistentry=NEW_NUMERIC(numVars+2+4));
				SET_VECTOR_ELT(listresults, i, tmplistentry);
                                UNPROTECT(1);
				}
/** *******************************************************************************
***********************************************************************************
 STEP 1. read in combinations of parents into cache - scores all zero'd */

/** create the observed data */
make_data(R_obsdata,&data,R_groupids);
// // print out data. 1: simply the data; 2: data and groupIDs
// printDATA(&data,2);

make_dag(&dag, numVars,(SEXP)NULL,1,R_vartype,&maxparents,R_groupedvars);/** create an empty network but with max.parents set **/
// print out DAG
// printDAG(&dag,2);



   curnode=asInteger(R_child)-1;/** -1 since R indexes start at unity and not zero **/
                           /** important NOTE: the nodecache is indexed in terms of i and NOT curnode e.g. it runs from 0,,,numVarsinCache-1 but the actual nodeid are in whichnodes */

/*index=0;*/
         /** copy the parent combination in cache into dag  **/
	 for(k=0;k<numVars;k++){dag.defn[curnode][k]=INTEGER(R_parents)[k];}
/*Rprintf("child=%d\n",curnode+1);
for(k=0;k<numVars;k++){Rprintf("|%d|",dag.defn[curnode][k]);}Rprintf("\n");*/

                         switch(dag.varType[curnode])  /** choose which type of node we have */
                         {

			   case 1:{ /** binary/categorical node */
			           if(dag.groupedVars[curnode]){/** have grouped binary variable so node is a glmm */
				     calc_node_Score_binary_rv_R(&dag,&data,curnode,errverbose,trace, &designmatrix, priormean, priorsd,priorgamshape,priorgamscale,maxiters,epsabs,
								  storeModes,epsabs_inner,maxiters_inner,finitestepsize, verbose,
								  h_guess,h_epsabs,maxiters_hessian,ModesONLY,
								  max_hessian_error,myfactor_brent, maxiters_hessian_brent, num_intervals_brent);
				    } else {/** not grouped so node is a glm **/
                                    calc_node_Score_binary(&dag,&data,curnode,errverbose, &designmatrix, priormean, priorsd,maxiters,epsabs,storeModes); }
                                    /** results are in dag->nodeScores and dag->modes (if storeModes=TRUE) */
				    storeNodeResults(listresults,&dag,1,curnode,dag.varType[curnode]);
                                    break;
                                   }

                           case 2:{ /** gaussian node */
			           if(dag.groupedVars[curnode]){/** have grouped binary variable so node is a glmm */
				     calc_node_Score_gaus_rv_R(&dag,&data,curnode,errverbose,trace, &designmatrix, priormean, priorsd,priorgamshape,priorgamscale,maxiters,epsabs,
								  storeModes,epsabs_inner,maxiters_inner,finitestepsize, verbose,
								  h_guess,h_epsabs,maxiters_hessian,ModesONLY,
								  max_hessian_error,myfactor_brent, maxiters_hessian_brent, num_intervals_brent);
				    } else {/** not grouped so node is a glm **/
                                    calc_node_Score_gaus(&dag,&data,curnode,errverbose, &designmatrix, priormean, priorsd,priorgamshape,priorgamscale,maxiters,epsabs, storeModes);}
				    storeNodeResults(listresults,&dag,1,curnode,dag.varType[curnode]);
                                    break;
                                   }

			    case 3:{ /** poisson node */
			            if(dag.groupedVars[curnode]){/** have grouped poisson variable so node is a glmm */
				      calc_node_Score_pois_rv_R(&dag,&data,curnode,errverbose,trace, &designmatrix, priormean, priorsd,priorgamshape,priorgamscale,maxiters,epsabs,
								  storeModes,epsabs_inner,maxiters_inner,finitestepsize, verbose,
								  h_guess,h_epsabs,maxiters_hessian,ModesONLY,
								  max_hessian_error,myfactor_brent, maxiters_hessian_brent, num_intervals_brent);
				    } else {/** not grouped so node is a glm **/
				      calc_node_Score_pois(&dag,&data,curnode,errverbose, &designmatrix, priormean, priorsd,maxiters,epsabs, storeModes);}
				    /** results are in dag->nodeScores and dag->modes (if storeModes=TRUE) */
				    storeNodeResults(listresults,&dag,1,curnode,dag.varType[curnode]);
				    break;
                                   }

                           default: {Rprintf("dag.varType[i]=%d\n",dag.varType[curnode]);error("in default switch - should never get here!");}
                         }


                         R_CheckUserInterrupt();/** allow an interupt from R console */


/*gsl_set_error_handler (NULL);*//** restore the error handler*/
free_dag(&dag);/** free any memory not allocated using Ralloc() */
UNPROTECT(1);

return(listresults);

}

/** *************************************************************************************************/
/** *************************************************************************************************/
/** copy results from C into R SEXP *****************************************************************/
void storeNodeResults(SEXP results,network *dag,int storeModes, int nodeid, int vartype)
{
  int i;
/** format is results has one vector for each node of length numNodes+3, the format is
    vec= node score, then each of the parameter modes - including entries DBL_MAX if node not in model**/
  REAL(VECTOR_ELT(results,0))[0]=dag->nodeScores[nodeid]; /** first entry is the mlik */
  REAL(VECTOR_ELT(results,0))[1]=dag->nodeScoresErrCode[nodeid];
  REAL(VECTOR_ELT(results,0))[2]=dag->hessianError[nodeid];
  /** next get all the parameters including intercept and extra term for gaussian precision parameter etc**/
  /** +4 = mlik+intercept+gauss precision + group precision**/
  for(i=1+2;i<dag->numNodes+2+4;i++){REAL(VECTOR_ELT(results,0))[i]=gsl_matrix_get(dag->modes,nodeid,i-1-2);}

}
