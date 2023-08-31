#include <stdio.h>
#include <stdlib.h>
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
#include "cycles.h"
#include "searchhill.h"

#define DEBUG_swap1
#define DEBUGimpscore1

/** pass a DAG definition from R and this returns an error if it contains a cycle ***/
/** MAIN FUNCTION **/
SEXP searchhill(SEXP R_children, SEXP R_nodecache, SEXP R_nodescores, SEXP R_numVars, SEXP R_numRows, SEXP R_numparspernode, SEXP R_seed, 
		SEXP R_verbose, SEXP R_timing, SEXP R_usestartdag, SEXP R_startdag, SEXP R_numsearches,SEXP R_dagretained)
{

cache nodecache;
network dag_scratch,dag_cur,dag_best,dagretain;
cycle cyclestore;
const gsl_rng_type *T;
gsl_rng *r;
clock_t start=0; clock_t end=0;
double elapsed;
unsigned int first=1;
int networkindex=0;
unsigned long int seed=asInteger(R_seed);
int numVars=asInteger(R_numVars);
int numRows=asInteger(R_numRows);
int verbose=asInteger(R_verbose);
int timingon=asInteger(R_timing);
int usestartdag=asInteger(R_usestartdag);
int nosearches=asInteger(R_numsearches);
int *store=(int *)R_alloc(numRows,sizeof(int));
/*int chosenindex=0;*/
int betterscore;/** if true then hill iter found a better DAG */
int searchNum;
SEXP listresults, scorevector,ans;
/*int storecounter=0;*/
/*make_dag(&dag_cur, numVars,(SEXP)NULL,1,(SEXP)NULL,(int*)NULL);*//** create an empty network **/
make_dag(&dag_scratch, numVars,(SEXP)NULL,1,(SEXP)NULL,(int*)NULL,(SEXP)NULL);/** create an empty network **/
make_dag(&dag_best, numVars,(SEXP)NULL,1,(SEXP)NULL,(int*)NULL,(SEXP)NULL);/** create an empty network **/
make_dag(&dagretain,numVars,R_dagretained,0,(SEXP)NULL,(int*)NULL,(SEXP)NULL);
init_hascycle(&cyclestore,&dag_scratch); /** initialise storage only need numnodes from dag object */
make_nodecache(&nodecache, numVars,numVars, numRows,R_numparspernode, R_children, R_nodecache, R_nodescores);/** this rolls data from R into local C arrays **/

/** ********************************************************************************/
/** SEXP Storage                                                                   */
/** ********************************************************************************/
/** create storage for sending back to R **/
 PROTECT(listresults = allocVector(VECSXP, nosearches*2+1));/** number of elements in the outer most list 
                                                                        two matrices for each search plus one score vector**/
         PROTECT(scorevector=NEW_NUMERIC(nosearches*2));/** a single vector containing the network score for each step of the search */
         SET_VECTOR_ELT(listresults, 0, scorevector);/** assign the above vector to the first entry in the R list */
         UNPROTECT(1);
/** ********************************************************************************/
/** ********************************************************************************/
/** create a generator - use default which is likely mersenne-twister - any gsl generator is sensible */
gsl_rng_env_setup();
T = gsl_rng_default;
r = gsl_rng_alloc (T);
gsl_rng_set (r, seed);

if(usestartdag){make_dag(&dag_cur, numVars,R_startdag,0,(SEXP)NULL,(int*)NULL,(SEXP)NULL);/** user supplied DAG **/
 } else {/** use a randomly chosen DAG so create empty dag_cur **/
         make_dag(&dag_cur, numVars,(SEXP)NULL,1,(SEXP)NULL,(int*)NULL,(SEXP)NULL);/** create an empty network **/
          }
          
searchNum=0;

while(searchNum<nosearches){

if(timingon){start = clock();}
 
 /** ********* GENERATE INITIAL NETWORK ***********************/
 /** Step 1. generate a new DAG which is a random choice from the node cache - each parent combination have same prob of being
    chosen and check that it is acyclic, and also get its network score **/
 /** use a randomly chosen DAG **/
  if(!usestartdag){
         generateRandomDAG(r,&dag_cur,&nodecache,store,&cyclestore,&dagretain,verbose);
	 if(hascycle(&cyclestore,&dag_cur)){error("initial DAG definition is not acyclic!");}}
 
 lookupscores(&dag_cur,&nodecache);/** lookup DAG in cache to find network score and also locations in cache which comprise the DAG **/
  
 if(verbose){Rprintf("init DAG\n");printDAG(&dag_cur,2);}
 
  copynetwork(&dag_cur,&dag_scratch);/** make a copy of the initial network in case new network is worse: dag into dag_cur*/
  copynetwork(&dag_cur,&dag_best);/** set all networks to the same */ 
  
  PROTECT(ans = allocMatrix(INTSXP, dag_best.numNodes, dag_best.numNodes));
         store_results(listresults,&dag_best,networkindex++, ans,0);
	 UNPROTECT(1);
   first=1;  
   

   while(first || dag_best.networkScore>dag_cur.networkScore){
   if(!first){/** skip this on first run **/
              if(verbose){Rprintf("DAG\n");printDAG(&dag_best,2);}
              /*Rprintf("DAG\n");Rprintf("best\n");printDAG(&dag_best,2);Rprintf("cur\n");printDAG(&dag_cur,2);Rprintf("scratch\n");printDAG(&dag_scratch,2);*/
              copynetwork(&dag_best,&dag_scratch);/** make a copy of the initial network in case new network is worse: dag into dag_cur*/
              copynetwork(&dag_best,&dag_cur);/** set all networks to the same */
   }
   /** Step 3. now iterate over the entire parent cache - ignoring those "rows" in the current DAG and in each case build a new DAG
   check it is acyclic and then compute its score **/
  /** dag is the network which will be tinkered with, dag_cur is the "baseline" at this step, and dag_best is to hold any new better network */
  hillSingleIteration(&dag_scratch,&nodecache,&betterscore,&dag_cur,&cyclestore,&dag_best, verbose);
  first=0;
  /*tmp++;*/
 /* if(tmp>4){break;}*/
 }
  
  PROTECT(ans = allocMatrix(INTSXP, dag_best.numNodes, dag_best.numNodes));
         store_results(listresults,&dag_best,networkindex++, ans,0);/** first 0 is the index of networks to save ignoring first scorevector */
	 UNPROTECT(1);
	 
  if(timingon){end = clock();
 elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
 if(verbose){Rprintf("CPU time:%10.6f\n",elapsed);};} 

 searchNum++;
 if(searchNum%100==0){if(verbose){Rprintf("%d out of %d\n",searchNum,nosearches);};};
}

gsl_rng_free (r);

/** test for free any memory not used */
     free_dag(&dag_scratch);
     free_dag(&dag_best);
     free_dag(&dagretain);
    
    /** TEST GK 20.03.2021*/
    free_dag(&dag_cur);
    
UNPROTECT(1);

return(listresults);

}

/** ********************************************************************************************************************/
/** ********************************************************************************************************************/
int generateRandomDAG(gsl_rng *r,network *dag,cache *nodecache,int *store,cycle *cyclestore, network *dagretain, int verbose)
{
int i,j;
int chosenindex=0;
/*do{*/
 /** choose a random set of parents from each node **/
 /*if(verbose){Rprintf("attempt %d\n",curtry+1);}*/ 
 
 for(i=0;i<dag->numNodes;i++){/** for each node find a parent combination in the cache of allowed combinations **/
   /** need to create fill temp array first **/
  for(j=0;j<nodecache->numparcombs[i];j++){store[j]=j;} /** arrayx with indices from 0...to M-1 with M total number of parent combs for node i **/
  gsl_ran_choose (r, &chosenindex, 1, store,nodecache->numparcombs[i], sizeof (int));
  /** now for node i copy its parent combination with index "chosenindex" into DAG dag **/
  /*Rprintf("node=%d chosen index=%d\n",i+1,chosenindex+1);*/
  /*dag->defn[i]=nodecache->defn[i][chosenindex];*//** just copy address - DO NOT DO THIS - nonsense as it can change the cache!**/
  for(j=0;j<dag->numNodes;j++){dag->defn[i][j]=nodecache->defn[i][chosenindex][j];}
 }

 /*for(i=0;i<numVars;i++){for(j=0;j<numVars;j++){Rprintf("%d ",dag->defn[i][j]);}Rprintf("\n");} */

 /** we now have the new randomly generated DAG so now check that it is acyclic **/
 checkandfixcycle(cyclestore,dag,r,dagretain,verbose); /** checks for a cycle and fixes it if it finds one **/
 /*yescycle=hascycle(cyclestore,dag);*//** check supposedly fixed model */  
 
 /*curtry++;
 R_CheckUserInterrupt();*//** allow an interupt from R console */ 
/*} while( yescycle && curtry<maxtries); */ /** continue searching is have cycle or exceeded try limit or DAG has too many arcs **/

return(1); /** return true if either of these are positive */

}

/** ********************************************************************************************************************/
/** ********************************************************************************************************************/
/** *************************************************************************************************/
/** try and make this inline since its a v.small function but helps clarity COPY DAG to DAG_SCRATCH */
/** *************************************************************************************************/
void copynetwork(network *dagsrc,network *dagdest)
{
 int i,j; 
 dagdest->networkScore=0.0;
 
 dagdest->numNodes=dagsrc->numNodes;/** copy number of nodes **/
 /*dagdest->networkScore=dagsrc->networkScore;*//** copy network score **/
 
 for(i=0;i<dagsrc->numNodes;i++){ 
   dagdest->nodeScores[i]=dagsrc->nodeScores[i];/** copy node scores **/
   dagdest->networkScore+=dagdest->nodeScores[i];
   dagdest->locationInCache[i]=dagsrc->locationInCache[i]; /** copy position of parents in cache **/
   for(j=0;j<dagsrc->numNodes;j++){dagdest->defn[i][j]=dagsrc->defn[i][j];}} /** copy DAG defn **/
 
}
/** ********************************************************************************************************************/
/** ********************************************************************************************************************/
/** perform a single step of a hill climber search *********************************************************************/
void hillSingleIteration(network *dag_scratch,cache *nodecache,int *betterscore,network *dag_cur, cycle *cyclestore, network *dag_best, int verbose){
  
  /**note dag_scratch, dag_cur and dag_best all start out identical **/
  int i,j;
  /** for each member of the cache of possible parents swap-in one at a time to see if there is a better DAG **/
  for(i=0;i<nodecache->numVars;i++){ /** for each variable in the current DAG */
    /*Rprintf("node=%d ",i);*/
      for(j=0;j<nodecache->numparcombs[i];j++){ /** for all possible allowed parent combinations for that node */
	/*Rprintf("parentcomb=%d \n",j);*/
             /** has the current DAG - "dag" got this parent combination, if not **/
	     if(dag_scratch->locationInCache[i]!=j){ /** parent combination for node i is not held at position j in nodecache so is NOT in the current DAG **/
	         /** swap the current parent combination with new one **/
	         swapnode(i,j,dag_scratch,nodecache);/** change parent combination for node i for cache[i][combination j], also update its score */
	         if(!hascycle(cyclestore,dag_scratch) && improvedscore(dag_scratch,dag_best)){/** check for cycle in new DAG **/
	             /** if does not exceed arc limit, and has no cycle, then compute network score and compare with current DAG score **/       
	             /*Rprintf("found better\n");*/ /** then our new dag can overwrite the existing one **/
		     /** 1. found a better model to copy dag to new best model **/
		     copynetwork(dag_scratch,dag_best);
		     if(verbose){printDAG(dag_best,2);}
	        } /*else {Rprintf("have cycle\n");*/
		        
           /** have changed the original dag for better or worse so reset for another try */
		     copynetwork(dag_cur,dag_scratch);
              } /** end of same combination check */
       
       } /** end of all combinations for node i */
   } /** end of all nodes */

 /** dag_best is the finally model */ 
 
}
/** ********************************************************************************************************************/
/** ********************************************************************************************************************/
/** change a parent combination for a single node *********************************************************************/
void swapnode(int nodeindex,int combinationindex,network *dag,cache *nodecache)
{

int i;

#ifdef DEBUG_swap
for(i=0;i<dag->numNodes;i++){
   for(j=0;j<dag->numNodes;j++){
     Rprintf("%d ",dag->defn[i][j]);
   }Rprintf("\n");
}
Rprintf("swap out node %d\n",nodeindex+1);
for(i=0;i<dag->numNodes;i++){Rprintf("%d ",nodecache->defn[nodeindex][combinationindex][i]);}Rprintf("\n");

#endif

/** swap out a single node combination */
for(i=0;i<dag->numNodes;i++){ /** overwrite the parent combination in dag with that in nodecache **/
             dag->defn[nodeindex][i]=nodecache->defn[nodeindex][combinationindex][i];
}

/** also update the location in cache - **needed** as if this is a good swap then will copy entire model**/
dag->locationInCache[nodeindex]=combinationindex;

/** also update the score for the node being changed **/
dag->nodeScores[nodeindex]=nodecache->nodeScores[nodeindex][combinationindex];

/*Rprintf("curNS=%f\n",dag->networkScore);
dag->networkScore=0.0;
for(i=0;i<dag->numNodes;i++){dag->networkScore+=dag->nodeScores[i];}
Rprintf("newNS=%f\n",dag->networkScore);
*/

#ifdef DEBUG_swap
Rprintf("new DAG\n");
for(i=0;i<dag->numNodes;i++){
   for(j=0;j<dag->numNodes;j++){
     Rprintf("%d ",dag->defn[i][j]);
   }Rprintf("\n");
}Rprintf("\n");


#endif



}

/** ********************************************************************************************************************/
/** ********************************************************************************************************************/
/** compare two DAGs and compare network scores    *********************************************************************/
int improvedscore(network *dag, network *dag_orig)
{

  int i;
  double netscore_dag=0.0;
  double netscore_dagorig=0.0;
  
  for(i=0;i<dag->numNodes;i++){
                               netscore_dag+=dag->nodeScores[i];
                               netscore_dagorig+=dag_orig->nodeScores[i];
  } 
  
  /*dag->networkScore=netscore_dag;*//** save new better score */
  
  #ifdef DEBUGimpscore
  Rprintf("old score=%f new score=%f\n",netscore_dagorig,netscore_dag);
  #endif
  if(netscore_dag>netscore_dagorig){/** new DAG has better score than previous one **/
     return(1);
  } else {return(0);}   

}
/** ********************************************************************************************************************/
/** ********************************************************************************************************************/
