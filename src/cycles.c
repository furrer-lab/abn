#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include "structs.h"
#include "utility.h" 
#include "cycles.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h> 
/** pass a DAG definition from R and this returns the value one if it contains a cycle ***/
/** MAIN FUNCTION **/
SEXP checkforcycles(SEXP R_dag, SEXP R_dagdim)
{

  network dag;
  cycle cyclestore;
  SEXP err = PROTECT(allocVector(INTSXP, 1));
  
  INTEGER(err)[0] = 0;
  make_dag(&dag, asInteger(R_dagdim), R_dag, 0, (SEXP)NULL, (int*)NULL,(SEXP)NULL);
  init_hascycle(&cyclestore, &dag); /** initialise storage but needs to be passed down through generate_random_dag etc */
  if (hascycle(&cyclestore, &dag)) {
    /*       Rprintf("well then\n");  */
    INTEGER(err)[0] = 1;
    /* using `error(msg)` does not work as variables will not be properly freed/cleaned. */
  }
  free_dag(&dag);
  UNPROTECT(1);
  return(err);
}
/** END OF MAIN **/

/** **************************************************************************************************/
/** create storage for doing cycle check - to replace previous static allocation calls */
/** **************************************************************************************************/
void init_hascycle(cycle *cyclestore,network *dag){

unsigned int i;
unsigned int numnodes=dag->numNodes;
unsigned int *isactive, *incomingedges;
unsigned int **graph,*graphtmp;
          
   isactive=(unsigned int *)R_alloc(numnodes,sizeof(unsigned  int));
   /*isactive_scratch=(unsigned int *)R_alloc(numnodes,sizeof(unsigned  int));*/
   incomingedges=(unsigned  int *)R_alloc(numnodes,sizeof(unsigned  int));
   graph=(unsigned  int **)R_alloc(numnodes,sizeof(unsigned  int*));/** create storage for a copy of the dag->defn[][] */
   for(i=0;i<numnodes;i++){graphtmp=(unsigned  int *)R_alloc( numnodes,sizeof(unsigned  int)); graph[i]=graphtmp;}

cyclestore->isactive=isactive;
cyclestore->incomingedges=incomingedges;
cyclestore->graph=graph;


}

/** *************************************************************************************************/
/** check for cycle in graph  - do this by checking for a topolgical ordering */
/** *************************************************************************************************/
unsigned int hascycle(cycle *cyclestore,network *dag){

unsigned int numnodes=dag->numNodes;
unsigned int i,j, nodesexamined,success;
unsigned int *isactive=cyclestore->isactive;
unsigned int *incomingedges=cyclestore->incomingedges;
unsigned int **graph=cyclestore->graph;

for(i=0;i<numnodes;i++){isactive[i]=1;} /** all nodes initially active */

/** copy the current network definition into graph[][]*/
for(i=0;i<numnodes;i++){for(j=0;j<numnodes;j++){graph[i][j]=dag->defn[i][j];}}

/** calc number of incoming edges for each child  and put into incomingedges**/
 get_numincomingedges(incomingedges,graph, numnodes);


/** find a node with no incoming edges - see lecture11.pdf in ../articles folder **/
nodesexamined=0;
success=1;
while(success){
        success=0;
	for(i=0;i<numnodes;i++){
        	  if(isactive[i] && !incomingedges[i]){/** if not yet examined and no incoming edges */
                	isactive[i]=0; /** set this child to inactive*/
                	
                        /** remove all OUTGOING links from node i, e.g. wherever i is a parent */
                	for(j=0;j<numnodes;j++){graph[j][i]=0;}
                	get_numincomingedges(incomingedges,graph,numnodes);
           	        success=1; nodesexamined++;}
			}
           }
         

      if(nodesexamined==numnodes){return(0);/** no cycle */                               
      } else {/*Rprintf("=>%d %d\n",nodesexamined, numnodes);*/
              return(1);} /** found a cycle */
 

}

/** *************************************************************************************************/
/** v.small function but helps clarity in hascycle() */
/** *************************************************************************************************/
void get_numincomingedges(unsigned int *incomingedges,unsigned int **graph, unsigned int numnodes)
{ /** count up how many parents each child has **/
unsigned int i,j;
unsigned int numincomedge;
for(i=0;i<numnodes;i++){/** for each child */
        numincomedge=0;
 	for(j=0;j<numnodes;j++){numincomedge+=graph[i][j];
                                }
        incomingedges[i]=numincomedge;
        } 	

}
/** ****************************************************************************************************/
/** check for cycle in graph  - do this by checking for a topolgical ordering, and if a cycle is found */
/** then remove it by dropping all outcoming arcs from a node and repeating                            */
/** ****************************************************************************************************/
void checkandfixcycle(cycle *cyclestore,network *dag, gsl_rng *r, network *dagretain, int verbose){

unsigned int numnodes=dag->numNodes;
unsigned int i,j, nodesexamined,success;
unsigned int *isactive=cyclestore->isactive;
unsigned int *incomingedges=cyclestore->incomingedges;
unsigned int **graph=cyclestore->graph;
gsl_permutation *p = gsl_permutation_alloc (numnodes);
unsigned int curnode=0;
int editNode=-1;
for(i=0;i<numnodes;i++){isactive[i]=1;} /** all nodes initially active */
/** copy the current network definition into graph[][]*/
for(i=0;i<numnodes;i++){for(j=0;j<numnodes;j++){graph[i][j]=dag->defn[i][j];}}

/** calc number of incoming edges for each child  and put into incomingedges**/
 get_numincomingedges(incomingedges,graph, numnodes);
if(verbose){Rprintf("start DAG\n");printDAG(dag,1);}
/** need to randomise which node we check first i.e. which cycle we find first */
gsl_permutation_init (p);
gsl_ran_shuffle (r, p->data, numnodes, sizeof(size_t));

/** find a node with no incoming edges - see lecture11.pdf in ../articles folder **/
nodesexamined=0;
success=1;
while(success){
        success=0;
	for(i=0;i<numnodes;i++){/** iterate over each node and find a node with no incoming edges */
	  curnode=gsl_permutation_get(p,i);/** choose nodes to check in a RANDOM order */
        	  if(isactive[curnode] && !incomingedges[curnode]){/** if not yet examined and no incoming edges */
                	isactive[curnode]=0; if(verbose){Rprintf("node %d is now inactive\n",curnode);}/** set this child to inactive*/
                        /** remove all OUTGOING links from node i, e.g. wherever i is a parent */
                	for(j=0;j<numnodes;j++){graph[j][curnode]=0;}
                	get_numincomingedges(incomingedges,graph,numnodes);
           	        success=1; nodesexamined++;}
			}
			
       if(success==0 && nodesexamined!=numnodes){
         /**to get to here means that all of the active nodes have at least one child and so we grab the first of these and remove the children and then repeat above */
         /** now remove all OUTGOING links from node = editNode, e.g. wherever i is a parent: NOTE - make changes to original dag here */
             for(i=0;i<numnodes;i++){if(isactive[gsl_permutation_get(p,i)]){editNode=gsl_permutation_get(p,i);break;}} /** get first active node */
	     droplinks(dag,dagretain->defn,editNode);/*for(j=0;j<numnodes;j++){dag->defn[j][editNode]=0;}*/ if(verbose){Rprintf("dropping links from node %d\n",editNode);} /** remove links **/
             for(i=0;i<numnodes;i++){for(j=0;j<numnodes;j++){graph[i][j]=dag->defn[i][j];}} /** copy into graph **/
	     if(verbose){printDAG(dag,1);}
	     isactive[editNode]=0;   
             get_numincomingedges(incomingedges,graph,numnodes); /** get new incoming links in graph */
	     nodesexamined++;
         success=1;/** repeat cycle check within while loop **/}
         
 }
 gsl_permutation_free(p);
         
}

/** *****************************************************************************************************/
/** to fix a cycle when creating a random network we remove arcs - but we should NOT remove arcs which  */
/** are present in a retain matrix. This function removes arcs while checking for this                  */
/** *****************************************************************************************************/
void droplinks(network *dag, int **retaingraph, unsigned int editnode)
{
 unsigned int numnodes=dag->numNodes; 
 unsigned int i;
 for(i=0;i<numnodes;i++){
    if(!retaingraph[i][editnode]){/** if outgoing arc does not need to be retained then drop **/
                                  dag->defn[i][editnode]=0;} /*else {Rprintf("found a retained arc child=%d parent=%d\n",i,editnode);}*/
 }

}
