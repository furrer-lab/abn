#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>
#include "structs.h"
#include "utility.h" 
#include "cycles.h"
#include "buildcachematrix.h"


/** MAIN FUNCTION **/
SEXP buildcachematrix(SEXP R_dagdim, SEXP R_dagbanned, SEXP R_dagretained, SEXP R_maxparents, SEXP R_whichnodes)
{

int i,j,k,total,drop;
network dag,dagban,dagretain;
const int numNodes=asInteger(R_dagdim);
const int *maxparents=INTEGER(R_maxparents);
gsl_combination *c;
int curnode=0;
int curparentindex=0;
SEXP listresults,nodeidvector,parentsvector;
int *ptr_parents,*ptr_nodeid;
int localindex;/** this indexes every individual entry (not just row) for all parent combinations */
/** setup some workspaces **/
make_dag(&dag,numNodes,R_dagbanned,0, (SEXP)NULL,(int*)NULL,(SEXP)NULL);/** just copy banned into dag - this is blanked out to zeros anyway later**/
make_dag(&dagban,numNodes,R_dagbanned,0,(SEXP)NULL,(int*)NULL,(SEXP)NULL);
make_dag(&dagretain,numNodes,R_dagretained,0,(SEXP)NULL,(int*)NULL,(SEXP)NULL);

/** ********************************************************************************/
/** STEP 0. - create R storage for sending results back => need total number of entries. To do this we iterate over **/
/** all possible combinations and check for banned arcs and retained arcs, dropping those combinations which do not **/
/**meet this criteria **/
/** ***************************************************************************************/

total=0;
drop=0;
for(i=0;i<LENGTH(R_whichnodes);i++){/** for each node passed get all network scores for all possible subsets subject to cardinality maxparents**/
  curnode=INTEGER(R_whichnodes)[i]-1;/** -1 since R indexes start at unity and not zero **/ 
  /** get all possible subsets **/
   for(j=0; j<=maxparents[curnode]; j++){/** allow each node to have different numbers of max. parents */  
    c = gsl_combination_calloc (numNodes, j);/** setup for all combinations of (maxparents choose N) 
                                                 NOTE: use N and not N-1 since it means not having to remap the indices but we then need to discard
                                                 any combination where the parent=1 and the parent=curnode **/
    do
       {drop=0;/** a flag */
       for(k=0;k<numNodes;k++){dag.defn[curnode][k]=0;}/** reset the node parents to independence **/
       for(k=0;k<gsl_combination_k(c);k++){curparentindex=gsl_combination_get (c, k);
	                                   if( curnode==curparentindex ||  /** invalid combination as parent is also the child */
                                               dagban.defn[curnode][curparentindex]==1 ){ /**this arc is banned so invalid combination */ 
                                               /*Rprintf("banned\n");*/
					       drop=1;break;
                                            }/** so skip to next combination **/
                                           dag.defn[curnode][curparentindex]=1;/** store combination */
                                        }
        /** now check for retained links **/
       if(!drop){/** do if not aleady set drop to 1 TRUE */
	 for(k=0;k<numNodes;k++){if(isretained(curnode,dag.defn,dagretain.defn,numNodes)){ } else {/*Rprintf("not retained\n");*/
	                                                                                       drop=1;break;}}
	 /** if not retained then skip to next iteration */
       }
       /** if get to here must have a valid combination */
       if(!drop){/** have a valid combination to increment counter */
	 total++;}
      
      }
       while (gsl_combination_next (c) == GSL_SUCCESS);
       gsl_combination_free (c);
    
    }
}
/*Rprintf("have %d combinations\n",total);*/
/** ********************************************************************************/
/** We now simply repeat above but this time store the results                     */
/** ********************************************************************************/
/** create storage for sending back to R **/
PROTECT(listresults = allocVector(VECSXP, 2));/** number of elements in the outer most list**/
         PROTECT(nodeidvector=NEW_INTEGER(total));/** a single vector containing the node id from 1-n */
	 SET_VECTOR_ELT(listresults, 0, nodeidvector);
	 PROTECT(parentsvector = NEW_INTEGER(total*numNodes));/** one long vector since the matrix is filled by cols not rows anyway so have to transform (back in R) */
	 SET_VECTOR_ELT(listresults, 1, parentsvector);
	 UNPROTECT(2);
	 
ptr_parents = INTEGER(parentsvector);/** point to the location of this matrix - actually one long vector **/
ptr_nodeid = INTEGER(nodeidvector);
/** ********************************************************************************/
/** ********************************************************************************/
localindex=0;	 
total=0;
drop=0;
for(i=0;i<LENGTH(R_whichnodes);i++){/** for each node passed get all network scores for all possible subsets subject to cardinality maxparents**/
  curnode=INTEGER(R_whichnodes)[i]-1;/** -1 since R indexes start at unity and not zero **/ 
  /** get all possible subsets **/
   for(j=0; j<=maxparents[curnode]; j++){/** allow each node to have different numbers of max. parents */  
    c = gsl_combination_calloc (numNodes, j);/** setup for all combinations of (maxparents choose N) 
                                                 NOTE: use N and not N-1 since it means not having to remap the indices but we then need to discard
                                                 any combination where the parent=1 and the parent=curnode **/
    do
       {drop=0;/** a flag */
       for(k=0;k<numNodes;k++){dag.defn[curnode][k]=0;}/** reset the node parents to independence **/
       for(k=0;k<gsl_combination_k(c);k++){curparentindex=gsl_combination_get (c, k);
	                                   if( curnode==curparentindex ||  /** invalid combination as parent is also the child */
                                               dagban.defn[curnode][curparentindex]==1 ){ /**this arc is banned so invalid combination */ 
                                               /*Rprintf("banned\n");*/
					       drop=1;break;
                                            }/** so skip to next combination **/
                                           dag.defn[curnode][curparentindex]=1;/** store combination */
                                        }
        /** now check for retained links **/
       if(!drop){/** do if not aleady set drop to 1 TRUE */
	 for(k=0;k<numNodes;k++){if(isretained(curnode,dag.defn,dagretain.defn,numNodes)){ } else {/*Rprintf("not retained\n");*/
	                                                                                       drop=1;break;}}
	 /** if not retained then skip to next iteration */
       }
       /** if get to here must have a valid combination */
       if(!drop){/** have a valid combination to increment counter */
	         /** now store the model defin into ***************************/
		 ptr_nodeid[total]=curnode+1;/** back to indexing from 1 **/
		 for(k=0;k<numNodes;k++){ptr_parents[localindex++]=dag.defn[curnode][k];}/** copy to store **/
	 total++;}
      
      }
       while (gsl_combination_next (c) == GSL_SUCCESS);
       gsl_combination_free (c);
    
    }
}	 
	 
/** test for free any memory not used */
     free_dag(&dag);
     free_dag(&dagban);
     free_dag(&dagretain);
     
	 
UNPROTECT(1);

return(listresults);	 
}



/** END OF MAIN **/

/** **************************************************************************************************/   
/** FUNCTIONS DOWN HERE ******************************************************************************/
/** **************************************************************************************************/                                                                                               
int isretained(int curnode, int **dagdefn,int **retained_m, int numNodes){

  /** check if dagdefn[curnode][..] contains all arcs - as a subset perhaps - which should be retained */
 int j; 
   for(j=0;j<numNodes;j++){/** for each col in retained */
           if(retained_m[curnode][j]==1){/** want to retain this arc **/
	      if(dagdefn[curnode][j]==1){/** the jth arc is retained */
	      } else {return(0);} /** we are finished as found an arc which was not retained */
	   }
  }
/** to get to here must have retained all the arcs */  
return(1);
}  
  

