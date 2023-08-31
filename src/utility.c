#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include "structs.h" 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#define DEBUG_groupdefn1
#define DEBUG_cachedefn1
/** **************************************************************************************************/   
/** ONLY FUNCTIONS HERE ******************************************************************************/
/** **************************************************************************************************/
/** copy a dag from R into C struct network **********************************************************/

void make_dag(network *dag, int dim, SEXP R_dag, int empty, SEXP R_vartype, const int *maxparents, SEXP R_groupedvars){

int **model_defn,*model_defnA;
int i,j;
dag->numNodes=dim;/** store number of nodes - e.g. observed variables in the DAG **/

if(maxparents){dag->maxparents=*maxparents;}

/** create storage **/
model_defn=(int **)R_alloc(dim,sizeof(int*));/** each row is a variable **/ 
	for(i=0;i<dim;i++){model_defnA=(int *)R_alloc(dim,sizeof(int)); model_defn[i]=model_defnA;}
dag->defn=model_defn;

dag->modes = gsl_matrix_alloc (dim,dim+3);/** this will hold a matrix of modes for each parameter 
                                                only used in posterior density estimates, +2 as +1 for intercept term
                                                and then some nodes might be gaussian and need another +1 and then +1 for group precision */

if(!empty){/** fill up new dag with contents of R_dag **/
/** copy contents of R_dag into array - NOTE that as R_dag is an R MATRIX it is single dimension and just needs to be unrolled */
for(j=0;j<dag->numNodes;j++){for(i=0;i<dag->numNodes;i++){dag->defn[i][j]=INTEGER(R_dag)[i+j*dag->numNodes];}} 
} else { /** fill with zeros **/
         for(j=0;j<dag->numNodes;j++){for(i=0;i<dag->numNodes;i++){dag->defn[i][j]=0;}}
}

#ifdef DEBUG_dagdefn   
for(i=0;i<dag->numNodes;i++){for(j=0;j<dag->numNodes;j++){Rprintf("%u ",dag->defn[i][j]);} Rprintf("\n");}    
#endif


/** also create storage for 1-d array which will hold which "row" in the cache this DAG is built from **/
dag->locationInCache=(int *)R_alloc(dim,sizeof(int));/** e.g. locationInCache[2]=30 means 
                                                               third node has parentcombination 30 in cache[2][...] **/
for(i=0;i<dag->numNodes;i++){dag->locationInCache[i]=0;} /** set to empty */                                                               

dag->nodeScores=(double *)R_alloc(dim,sizeof(double));/** will hold individual node scores  **/

for(i=0;i<dag->numNodes;i++){dag->nodeScores[i]=0.0;} /** set to */
 
dag->nodeScoresErrCode=(int *)R_alloc(dim,sizeof(int));/** will hold individual node scores  **/

for(i=0;i<dag->numNodes;i++){dag->nodeScoresErrCode[i]=0;} /** set to 0 - no error */ 

dag->hessianError=(double *)R_alloc(dim,sizeof(double));/** will hold individual node scores  **/

for(i=0;i<dag->numNodes;i++){dag->hessianError[i]=0.0;} /** set to 0 - no error */ 

dag->networkScore=0.0;/** dummy **/  

/** store vartype - if actually passed **/
if(R_vartype){
dag->varType=(int *)R_alloc(dim,sizeof(int));
for(i=0;i<dag->numNodes;i++){dag->varType[i]=INTEGER(R_vartype)[i];}
}

/** get the indexes of the nodes who need adjustment for correlation e.g. rv's */
if(R_groupedvars){
  dag->groupedVars=(int *)R_alloc(dim,sizeof(int));
  for(i=0;i<dag->numNodes;i++){dag->groupedVars[i]=0;} /** set all to zero */
  for(i=0;i<dag->numNodes;i++){/** for each node */
    
    for(j=0;j<LENGTH(R_groupedvars);j++){/** for each grouped variable */
                      if(i==INTEGER(R_groupedvars)[j]){/** this node need grouping adjustment **/
			                              dag->groupedVars[i]=1;break;}
    }
  }
  #ifdef DEBUG_groupdefn   
Rprintf("groupvars-\n");
for(i=0;i<dag->numNodes;i++){Rprintf("%d ",dag->groupedVars[i]);} Rprintf("\n");   
#endif
}  
  
#ifdef DEBUG_dagdefn   
Rprintf("vartypes\n");
for(i=0;i<dag->numNodes;i++){Rprintf("%d ",dag->varType[i]);} Rprintf("\n");   
#endif


}
/** *************************************************************************************************/
/** free up GSL memory ******************************************************************************/
void free_dag(network *dag)
{
    
gsl_matrix_free(dag->modes);  
  
}
/** **************************************************************************************************/
/** **************************************************************************************************/
/** copy data.frame into a C struct datamatrix *******************************************************/
/** **************************************************************************************************/
void make_data(SEXP R_obsdata,datamatrix *obsdata, SEXP R_groupids)
{

int numDataPts,numVars,i,j;
double **data, *tmpdata;
numVars=LENGTH(R_obsdata);/** number of columns in data.frame */
numDataPts=LENGTH(VECTOR_ELT(R_obsdata,0));


/** create a copy of R_data.frame into 2-d C array of ints - note: each CASE is an array NOT each variable */ 
data=(double **)R_alloc( (numDataPts),sizeof(double*));/** number of ROWS*/
	for(i=0;i<numDataPts;i++){tmpdata=(double *)R_alloc( numVars,sizeof(double)); data[i]=tmpdata;} 

  for(i=0;i<numDataPts;i++){/** for each CASE/observation **/
     for(j=0;j<numVars;j++){/** for each variable **/
                             data[i][j]=REAL(VECTOR_ELT(R_obsdata,j))[i];/** copy data.frame cell entry into C 2-d array entry **/
                             
                              }
     }
   
/** store vartype - if actually passed **/
if(R_groupids){
obsdata->groupIDs=(int *)R_alloc(numDataPts,sizeof(int));
for(i=0;i<numDataPts;i++){obsdata->groupIDs[i]=INTEGER(R_groupids)[i];}
#ifdef DEBUG_groupdefn   
Rprintf("groupids\n");
for(i=0;i<numDataPts;i++){Rprintf("%d ",obsdata->groupIDs[i]);} Rprintf("\n");   
#endif
}

/*for(j=0;j<numVars;j++){Rprintf("%d ",myvartype[j]);}Rprintf("\n"); */  
obsdata->defn=data;/** original observed data */
obsdata->numVars=numVars;/** total number of variables */
obsdata->numDataPts=numDataPts;/** total number of case/observations */

} 

/** **************************************************************************************************/
/** copy a dag from R into C struct network **********************************************************/
void make_nodecache(cache *nodecache, int vars, int cols, int rows, SEXP R_numparcombs, SEXP R_children, SEXP R_cachedefn, SEXP R_nodescores)
{
/** vars is the number of variables which the cache is wanted, cols is the number of variables in the DAG */
/** whichnodes is the node indexes from 1 **/

int ***defn,**defnA,*defnB,*numParentcombs,**tmpbig,*tmpbigA;
double **scores,*scoresA;
int i,j,k,index;
nodecache->numVars=vars;/** store number of nodes - e.g. observed variables in the DAG **/
nodecache->numRows=rows;/** store number of nodes - e.g. observed variables in the DAG **/

/** create storage **/
numParentcombs=(int *)R_alloc(vars,sizeof(int)); /** will hold number of parent combinations per node **/
for(i=0;i<vars;i++){numParentcombs[i]=INTEGER(R_numparcombs)[i];} 

nodecache->numparcombs=numParentcombs;

/** create storage 3-d array, first dim is nodeid, second is the matrix of parent combs for that node,
 e.g. defn[i][j] is the jth parent combination for node i => a vector **/
defn=(int ***)R_alloc(vars,sizeof(int**));/** for each node **/ 
	for(i=0;i<vars;i++){defnA=(int **)R_alloc(nodecache->numparcombs[i],sizeof(int*)); defn[i]=defnA;
	                    for(j=0;j<nodecache->numparcombs[i];j++){defnB=(int *)R_alloc(cols,sizeof(int)); defnA[j]=defnB; /** note cols not vars **/
	                    }
	}

nodecache->defn=defn;/** point storage **/

/** now create a TEMPORARY 2-d MATRIX - which the R_cachedefn is copied into and then this is rolled into the 3-d array*/
/** tmpbig is only used to make rolling the data into the 3-d array easy, e.h. go from a 2-d to 3-d not 1-d to 3-d **/
tmpbig=(int **)R_alloc(rows,sizeof(int*));
	for(i=0;i<rows;i++){tmpbigA=(int *)R_alloc(cols,sizeof(int)); tmpbig[i]=tmpbigA;}
for(j=0;j<cols;j++){for(i=0;i<rows;i++){tmpbig[i][j]=INTEGER(R_cachedefn)[i+j*rows];}} 

scores=(double **)R_alloc(vars,sizeof(double*));
	for(i=0;i<vars;i++){scoresA=(double *)R_alloc(nodecache->numparcombs[i],sizeof(double)); scores[i]=scoresA;}

nodecache->nodeScores=scores;

/** now roll this into the 3-d array **/
index=0;/** this is the "row index" in the 2-d  array */
for(i=0;i<vars;i++){/** for each node in DAG **/
   for(j=0;j<nodecache->numparcombs[i];j++){/** for each of its parent combinations **/
      if(R_nodescores){nodecache->nodeScores[i][j]=REAL(R_nodescores)[index];/** copy scores from 1-d to 2-d **/
	               if(ISNAN(nodecache->nodeScores[i][j])){nodecache->nodeScores[i][j]= -DBL_MAX;} /** if missing then assign to worst possible value */
      } else {/*error("null R_nodescores");*/nodecache->nodeScores[i][j]=0.0;} /** is R_nodescores is null - which will be the case when making the buildscorecache() then fill with zeros **/	
      for(k=0;k<cols;k++){nodecache->defn[i][j][k]=tmpbig[index][k];}
      index++;}}

#ifdef DEBUG_cachedefn 
for(i=0;i<rows;i++){for(j=0;j<cols;j++){Rprintf("%u ",tmpbig[i][j]);} Rprintf("\n");}    
Rprintf("---\n\n");
for(i=0;i<cols;i++){Rprintf("%d |\n",i+1);
   for(j=0;j<nodecache->numparcombs[i];j++){
      for(k=0;k<cols;k++){Rprintf("%u ",nodecache->defn[i][j][k]);}Rprintf("score=%f\n",nodecache->nodeScores[i][j]);}}
      
#endif

}
/** **************************************************************************************************/
/** **************************************************************************************************/
/** given a DAG then lookup in cache to find out rwos and get node scores    *************************/
/** **************************************************************************************************/
int lookupscores(network *dag,cache *nodecache)
{
  
int i,j,k;
int ismatch=0;
int foundmatch=0;
/*double myscore=0.0;*/
dag->networkScore=0.0;

for(i=0;i<dag->numNodes;i++){/** iterate through each node **/
  foundmatch=0;
  for(j=0;j<nodecache->numparcombs[i];j++){/** iterate through each parent combination **/
    ismatch=1;/** set to true for a parent match **/
    for(k=0;k<dag->numNodes;k++){/** iterate through members of each parent combination **/
                                 if(dag->defn[i][k]!=nodecache->defn[i][j][k]){ismatch=0;break;}
    }
    if(ismatch){/** is this is still true then must have an exact match so get node score and then skip to next node **/
                dag->nodeScores[i]=nodecache->nodeScores[i][j];
		dag->networkScore+=dag->nodeScores[i];/** save total score */
                foundmatch=1;
		dag->locationInCache[i]=j;/** save position in cache */
		break;/** skip to next node */
               }
    }
 if(!foundmatch){/** if get to here and foundmatch is still negative then that says the DAG is invalid e.g. requested parent comb not exist in cache**/
   error("DAG not found in cache!");
   /*dag->networkScore= -DBL_MAX;*//** e.g. just a dummy value saying model is invalid **/
   /*return(foundmatch);*/
 }
 /** if foundmatch is true then continue to next node **/

}

/*dag->networkScore=myscore;*/

return(1); 
  
}

/** **************************************************************************************************/
/** convenience function to print a dag                                      *************************/
/** **************************************************************************************************/
void printDAG(network *dag,int what)
{
int i,j;
double score=0.0;

  switch(what){
    case 1: { /** just print DAG definition */
                  for(i=0;i<dag->numNodes;i++){Rprintf("--");}Rprintf("\n");
		  for(i=0;i<dag->numNodes;i++){
                       for(j=0;j<dag->numNodes;j++){
                              Rprintf("%d ",dag->defn[i][j]);}Rprintf("\n");
                                               }
             break;} 
    
    case 2: { /** DAG definition and total network score **/
                for(i=0;i<dag->numNodes;i++){Rprintf("--");}Rprintf("\n");
		  for(i=0;i<dag->numNodes;i++){
                       for(j=0;j<dag->numNodes;j++){
                              Rprintf("%d ",dag->defn[i][j]);}Rprintf("\n");
                                               }
                  for(j=0;j<dag->numNodes;j++){Rprintf("nodescore=%f\n",dag->nodeScores[j]);score+=dag->nodeScores[j];}                            
                  Rprintf("-- log mlik for DAG: %f --\n",score);
		  for(j=0;j<dag->numNodes;j++){Rprintf("--");}Rprintf("\n");
		  
		  
             break;}
      
    default: error("printDAG - should never get here!");
    
    
    
    
  }


}
/** **************************************************************************************************/
/** convenience function to print a data set                                 *************************/
/** **************************************************************************************************/
void printDATA(datamatrix *data,int what)
{
int i,j;

  switch(what){
    case 1: { /** just print DATA definition */
                  for(i=0;i<data->numVars;i++){Rprintf("--");}Rprintf("\n");
		  for(i=0;i<data->numDataPts;i++){
                       for(j=0;j<data->numVars;j++){
                              Rprintf("%f ",data->defn[i][j]);}Rprintf("\n");
                                               }
             break;} 
      
    default: error("printDATA - should never get here!");
    
  }


}
/** **************************************************************************************************/
/** convenience function to print a node cache set                                 *************************/
/** **************************************************************************************************/
void printCACHE(cache *nodecache,int what)
{
int i,j,k;
int **nodecombs;
  
  switch(what){
    case 1: { /** just print DATA definition */
              Rprintf("-----Parent Definitions-----\n");
		  for(i=0;i<nodecache->numVars;i++){
		    nodecombs=nodecache->defn[i];
                       for(j=0;j<nodecache->numparcombs[i];j++){
			 Rprintf("Node %d:\t",i+1);
			 for(k=0;k<nodecache->numVars;k++){
                              Rprintf("%d ",nodecombs[j][k]);}Rprintf("\tscore=%f\n",nodecache->nodeScores[i][j]);
                                               }
		  }
             break;} 
      
    default: error("printCACHE - should never get here!");
    
  }


}
/** **************************************************************************************************/
/** **************************************************************************************************/
void store_results(SEXP R_listresults,network *dag, int iter, SEXP ans, int verbose)
{
  
int *rans;
int i,j;
double score=0;
/** store the network score **/
for(i=0;i<dag->numNodes;i++){score+=dag->nodeScores[i];}

REAL(VECTOR_ELT(R_listresults,0))[iter]=score;
/** store the network as a matrix **/
       rans = INTEGER(ans);
       for(i = 0; i < dag->numNodes; i++) {/** fill by row */
         for(j = 0; j < dag->numNodes; j++){
           rans[i + (dag->numNodes)*j] = dag->defn[i][j];}
       }
SET_VECTOR_ELT(R_listresults, iter+1, ans);

if(verbose){
for(i = 0; i < dag->numNodes; i++) {
         for(j = 0; j < dag->numNodes; j++){Rprintf("%d|",dag->defn[i][j]);}Rprintf("\n");}Rprintf("\n");
}

}

/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/       
void print_state (int iter, gsl_multiroot_fdfsolver * s)
{
 unsigned int i=0;
    Rprintf ("iter = %3u\n",iter);
    
    for(i=0;i< (s->x)->size-1; i++){
          Rprintf ("x=%5.10f ",gsl_vector_get (s->x, i));}
          Rprintf ("x=%5.10f\n",gsl_vector_get (s->x, (s->x)->size-1));
	  
    for(i=0;i< (s->x)->size-1; i++){
          Rprintf ("f(x)=%5.10f ",gsl_vector_get (s->f, i));}
          Rprintf ("f(x)=%5.10f\n",gsl_vector_get (s->f, (s->x)->size-1));   
    
  
}   
   
/** *************************************************************************************
*****************************************************************************************
*****************************************************************************************/
void print_state_min (int iter, gsl_multimin_fdfminimizer  * s)
{
 unsigned int i=0;
    Rprintf ("iter = %3u\n",iter);
    
    for(i=0;i< (s->x)->size-1; i++){
          Rprintf ("x=%5.10f ",gsl_vector_get (s->x, i));}
          Rprintf ("x=%5.10f",gsl_vector_get (s->x, (s->x)->size-1));
	  
    
          Rprintf ("f(x)=%5.10f\n",s->f);
          
    
  
}   
