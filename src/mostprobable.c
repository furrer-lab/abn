#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>
#include "mobius.h"
#include <time.h>

#define PRINTMEee
#define NotLOGa
#define UseLOG
#define JUNKa
#define PRINT_RES

SEXP mostprobable_C(SEXP R_localscoreslist, SEXP R_numnodes, SEXP R_indexstart, SEXP R_indexend, SEXP R_priorchoice, SEXP R_verbose)
{
/** ****************/
/** declarations **/
unsigned int numRows,i,j,k,l;
unsigned int numNodes;
/*double *f_hat_0,*f_hat_1;*/
int priorchoice;
SEXP listresults,tmplistentry;
int *ptr_nodeid=INTEGER(VECTOR_ELT(R_localscoreslist,0));/** vector of node ids from zero **/
double *ptr_score=REAL(VECTOR_ELT(R_localscoreslist,2)); /** vector of node scores */
clock_t startclock=0; clock_t endclock=0;double elapsed;int timingon=1;
int **parents,*parentstmp;
double tmp=0.0;
double *alphalong;
int *alphalonglookup,**final_results,*final_resultstmp;
int *parents_numparents,*start,*end;
int tot=0;
/*int myindex=0;*/
/* int **banned,*bannedtmp,**retained,*retainedtmp;*/
int *S,*Smax;
int cardinality;
double totaldbl;int total;double *g;
int **g_lookup,*g_tmp,*g_cardinality,*g_maxlookup,*g_mynodemax;
gsl_combination *c;
int g_index;
int mynode=0;
int mynodemax=0;
double myalphamax,myalpha;
int myalphaparentindex;
int myalphaparentindexmax=0;
int minindex;
int curindex=0;
int Match=0;
int *Sprime;
/*int Sprimecard=0;*/
int update=-1;
int lastindex;
int startindex=0;
int *curS;
int got=0;
priorchoice=asInteger(R_priorchoice);
numRows=LENGTH(VECTOR_ELT(R_localscoreslist,0));/** number of local node scores */ /*Rprintf("no. of score=%d\n",numRows);*/
numNodes=asInteger(R_numnodes);/** number of nodes in the network */

int verbose=asInteger(R_verbose);

/*gsl_permutation *p;*/
/*Rprintf("note: using logs - experimental...\n");*/
/** *******************************************************************************
***********************************************************************************
STEP 0. - create R storage for sending results back                              */
/** ******************************************************************************/
/** generic code to create a list comprising of vectors
   - currently overkill but useful template **/
PROTECT(listresults = allocVector(VECSXP, 1));
for(i=0;i<1;i++){
				PROTECT(tmplistentry=NEW_INTEGER(numNodes*numNodes));
				SET_VECTOR_ELT(listresults, i, tmplistentry);
                                UNPROTECT(1);
				}

/** *******************************************************************************/
/** *******************************************************************************
***********************************************************************************
STEP 1. memory allocations for housekeeping stuff                                */
/** *******************************************************************************/
/**  R matrix is really a long vector - going to create a 2-d array of pointers to this as it will produce easier to read code */
parents=(int **)R_alloc( (numRows),sizeof(int*));/** number of ROWS*/
	for(i=0;i<numRows;i++){parentstmp=(int *)R_alloc( numNodes,sizeof(int));parents[i]=parentstmp;}
/** also need later the number of parents in each parent combination  **/
parents_numparents=(int *)R_alloc( (numRows),sizeof(int));/** number of ROWS*/

/** create a matrix which will contain the parent combinations for a single node - currently
 *  these are the same pattern for all nodes e.g. same max.parent limit and no ban or retain lists */
/*numsubsets=(int *)R_alloc(numNodes,sizeof(int));*//** number of subsets is the number of parent combinations per node */
/*for(i=0;i<numNodes;i++){numsubsets[i]=INTEGER(R_numsubsets)[i];}*/

start=(int *)R_alloc(numNodes,sizeof(int));;/** number of subsets is the number of parent combinations per node */
for(i=0;i<numNodes;i++){start[i]=INTEGER(R_indexstart)[i];}

end=(int *)R_alloc(numNodes,sizeof(int));;/** number of subsets is the number of parent combinations per node */
for(i=0;i<numNodes;i++){end[i]=INTEGER(R_indexend)[i];}

/** need an alpha store e.g. alpha[i][] is alpha_i and alpha[i][j] is the jth member of alpha_i **/
alphalong=(double *)R_alloc( numRows,sizeof(double));/** one for each node*/

/*alpha=(double **)R_alloc( numNodes,sizeof(double*));*//** one for each node*/
	/*for(i=0;i<numNodes;i++){alphatmp=(double *)R_alloc( numsubsets[i],sizeof(double)); alpha[i]=alphatmp;}*/ /** one for each subset **/

alphalonglookup=(int *)R_alloc( (numRows),sizeof(int));/** one for each node*/

/*alphalookup=(int **)R_alloc( (numNodes),sizeof(int*));*//** one for each node*/
	/*for(i=0;i<numNodes;i++){alphalookuptmp=(int *)R_alloc( numsubsets[i],sizeof(int)); alphalookup[i]=alphalookuptmp;}*/ /** one for each subset **/

final_results=(int **)R_alloc( (numNodes),sizeof(int*));/** one for each node*/
	for(i=0;i<numNodes;i++){final_resultstmp=(int *)R_alloc( numNodes,sizeof(int)); final_results[i]=final_resultstmp;} /** one for each subset **/
#ifdef JUNK
banned=(int **)R_alloc( (numNodes),sizeof(int*));/** one for each node*/
	for(i=0;i<numNodes;i++){bannedtmp=(int *)R_alloc( numNodes,sizeof(int)); banned[i]=bannedtmp;} /** banned matrix**/
retained=(int **)R_alloc( (numNodes),sizeof(int*));/** one for each node*/
	for(i=0;i<numNodes;i++){retainedtmp=(int *)R_alloc( numNodes,sizeof(int)); retained[i]=retainedtmp;} /** banned matrix**/
#endif
/** *******************************************************************************/
/** *******************************************************************************
***********************************************************************************
STEP 2. housekeeping mostly for housekeeping stuff                                */
/** *******************************************************************************/
/** "copy" R "matrix" containing the parent combinations into the 2d C array n.b. NOT COPYING BUT USING ADDRESS LOOKUPS*/
/** this takes up just as much space.....! */

for(i=0;i<numRows;i++){for(j=0;j<numNodes;j++){parents[i][j]=INTEGER(VECTOR_ELT(R_localscoreslist,1))[i+numRows*j];}}
#ifdef JUNK
/** banned and retain matrices - unroll from R into 2d C array*/
for(i=0;i<numNodes;i++){
   for(j=0;j<numNodes;j++){banned[i][j]=INTEGER(VECTOR_ELT(R_localscoreslist,3))[i+numNodes*j];
                           retained[i][j]=INTEGER(VECTOR_ELT(R_localscoreslist,4))[i+numNodes*j];
  }
}
#endif
/** now iterate over all parent combinations for all nodes and calculate the prior value for each qprime() and
    also calculate the number of parents in each parent combination. The log score for each node combination is
    overwritten by adding in the prior value - save allocating more memory */
for(i=0;i<numRows;i++){/** for every combination of parents at each node passed */

   /** calculate the beta_i(G_i) values p557 in Koivisto - needed for each part of the later alpha calc**/
   tot=0;for(j=0;j<numNodes;j++){if(parents[i][j]==1){tot++;}} parents_numparents[i]=tot;/** calculate no.parents per node which we are looping */
   tmp= qprime_prior(numNodes,parents_numparents[i],priorchoice);/** unlike in getposterior_features() we now have feature=1 since consider all features  */
   /** note - using logs for prior*/
   ptr_score[i]+= gsl_sf_log(tmp);/** overwrite existing log network score to logscore+prior */
   }

for(i=0;i<numNodes;i++){for(j=0;j<numNodes;j++){final_results[i][j]= -1;}} /** initialise */
/** *******************************************************************************/
/** *******************************************************************************
***********************************************************************************
STEP 3. we now have a vector ptr_score[] with all the values needed for the alpha_i(S) calc **/
/** for each parent combination we need to find the subset of that combination with the max ptr_score **/
/** each subset must meet the ban and retain constrainsts                  */

 getAlphaMax(ptr_score,ptr_nodeid,parents,alphalong,alphalonglookup,numNodes,numRows);

  R_CheckUserInterrupt();/** allow an interupt from R console */

if(verbose){Rprintf("Step1. completed max alpha_i(S) for all i and S\n");};

/** for this to work we need alphalong and alphalonglookup */
/** *************************************************************************************************************************************/
/** **** We now have alpha_i(S) for every i and S where this contains the maximum value and a lookup of the relevant parent combination */
/** **** in alpha[i][s] and alphalookup[i][s] respectively ******************************************************************************/
/** alpha is the form alpha[i][j] is the ith node (indexing from 0) and then the jth parent combination (indexing from 0) for that node */
/** We now find g(V) but where we want to find the MAX of the sum at each stage in the recursive g() calc. this is a pain but seems to **/
/** work, however it is probably inefficient. The hardest part is finding out which parents are those for each node and this is current */
/** implemented by re-running the recursion for successively smaller sets and each times gives us the optimal combination for one node  */
/** *************************************************************************************************************************************/

/** *************************************************************************************************************************************/
/** *************************************************************************************************************************************/
/** The dynamic programming part of Koivisto.
 *  We want evaluate g(S) for S={0,...n-1} which involves computing g(S) for all subsets of S and this is done by breaking this into    */
/** nested blocks in the way of dynamic programming - this is key as the simpler way is to use the sum over all possible permutations   */
/** see Koivisto top of p.557 but this involves much more computation as there are vastly more permutations than subsets it iterate     */
/** over. We first compute g(empty set) then g(sets of cardinality=1 and then iterate upwards from their                                */
/** *************************************************************************************************************************************/

/*if(maxparents< numNodes-1){error("must use max. possible number of parents\n");}*//** fix this later */

/** we are going to use a FIXED size array to represent each subset, where a 1 gives whether the current col/parent is a member or not  */
S=(int *)R_alloc(numNodes,sizeof(int));
/*int *Sprime_copy=(int *)R_alloc(numNodes,sizeof(int));
int *Sprime_indexes=(int *)R_alloc(numNodes,sizeof(int));*/
Smax=(int *)R_alloc(numNodes,sizeof(int));
/*int *Sprime_short=(int *)R_alloc(numNodes-1,sizeof(int));*/
for(i=0;i<numNodes;i++){S[i]=0;Smax[i]=0;} /** set up an array to hold parent combinations **/
cardinality=0;/** cardinality */
totaldbl=0.0;
/** ******************************/
/** now separate case where there is no parent limit e.g. maxparents=numNodes-1 from where there is a parent limit
 *  the big difference here is that in the former we can calculate g(S) for S=(0,1,...,numNodes-1) direct from the all the g(S) previously
 *  computed where |S|<=numNodes-1, however in the latter more steps are subsequentl required *******************************************/
/*if(maxparents==numNodes-1){max_S_cardinality=numNodes;
} else {max_S_cardinality= maxparents+1;}
*/
for(i=0;i<=numNodes;i++){totaldbl+=gsl_sf_choose(numNodes,i);} /** sum of numNodes choose 0,...,numNodes choose maxparents **/
total=rint(totaldbl);/** this is NEEDED - as gsl_sf_choose() returns double and just casting to int does not work - e.g. may be out by 1 rint() seems to work fine **/
if(verbose){Rprintf("Total sets g(S) to be evaluated over: %d\n",total);};
/** total is the overall number of different subsets, including empty and full sets */
g=(double *)R_alloc(total,sizeof(double));/** will hold all values of g() */
for(i=0;i<total;i++){g[i]=0.0;} /** initalise */

g_cardinality=(int *)R_alloc(total,sizeof(int));
g_mynodemax=(int *)R_alloc(total,sizeof(int));
g_lookup=(int **)R_alloc(total,sizeof(int*));/** a row for each subset, and the cols define that subset **/
	for(i=0;i<total;i++){g_tmp=(int *)R_alloc(numNodes,sizeof(int)); g_lookup[i]=g_tmp;}
g_maxlookup=(int *)R_alloc(total,sizeof(int));/** a row for each subset, and the cols define that subset **/
	/*for(i=0;i<total;i++){g_tmp=(int *)R_alloc(numNodes,sizeof(int)); g_maxlookup[i]=g_tmp;}*/

/** go through and allocate all possible parent sets needed from N choose 0, through to N choose N **/
/** g_lookup will contain every possible subset expressed as 0001 0101 etc where the parents are 1's - a fixed width representation **/
/** we do NOT restrict these due to banned or retain lists as g() does not refer to a single node but multiple nodes */
g_index=0;
for(cardinality=0;cardinality<=numNodes;cardinality++){/** for each size of subset including empty and up to full set restricted by parent number*/
c = gsl_combination_calloc (numNodes, cardinality);/** setup for all combinations of (lenS choose card) **/
    do
       {
       for(k=0;k<numNodes;k++){g_lookup[g_index][k]=0;} /** initialise to empty set - all zeros**/
       for(k=0;k<gsl_combination_k(c);k++){g_lookup[g_index][gsl_combination_get (c, k)]=1;} /** set the parents for this combination */
       g_cardinality[g_index]=cardinality;/** store the size of this current set */
       /*for(k=0;k<numNodes;k++){Rprintf("%d",g_lookup[g_index][k]);}Rprintf("\tcardinality=%d\n",g_cardinality[g_index]);*/
       g_index++; /** increment to next entry in g_lookup/g_cardinality */
      }
       while (gsl_combination_next (c) == GSL_SUCCESS);
       gsl_combination_free (c);

} /** end of cardinality loop */

/** *************************************************************************************************************************************/
/** so now g_lookup contains all the possible subsets S of g(S) needed - we now need to calc the actual values for g(s) **/
/** *************************************************************************************************************************************/
mynode=0;
mynodemax=0;
myalphaparentindexmax=0;
curindex=0;
Match=0;
/*int Sprimecard=0;*/
update=-1;
g_index=0;/** will index each new entry into g() */
/** **************** CARDINALITY ZERO ***************************************************************************************************/
g[g_index]=1.0;/** g(empty set) **/
/** no corresponding g_maxlookup as there is no node here **/
g_maxlookup[g_index]= -1;
g_mynodemax[g_index]=-1;/** impossible since boundary */
 /*  Rprintf("-----------------------------------------------------\n");
   Rprintf("g(S); [fixed width] S=| ");
   for(j=0;j<numNodes;j++){Rprintf("%d ",g_lookup[g_index][j]);}Rprintf("|\n");
   Rprintf("MAX: node=EMPTY");Rprintf("\n");
   Rprintf("-----------------------------------------------------\n");*/
   g_index++;

/** **************** All S with CARDINALITY=1 ********************************************************************************************/
/** need to adjust for the retained arcs here since that means that e.g. g(1)=alpha_1(empty) is not allowed when node one has a retained arc */
  for(j=0;j<numNodes;j++){
  update=0;
  for(i=0;i<numRows;i++){
			   if(ptr_nodeid[i]==j){/**have node j **/
			     Match=1;
			     for(k=0;k<numNodes;k++){/** for each parent **/
			       if(parents[i][k]!=0){Match=0;break;}
			     }
                             if(Match==1){ /** have the empty set for that node - this is present in the node sets **/
			       g[g_index]=alphalong[i];
                               g_mynodemax[g_index]=j;
			       g_maxlookup[g_index]=alphalonglookup[i];/** this holds the index in g_lookup for the parents of Sprime - see later */
			       update=1;
			       break;} /** break since we are finished with this node */
			    /** set flag to say this node has been processed */
			   }
  }
 if(!update){/** did not find the empty set for this node so it must have at least one arc to be retained **/
                               g[g_index]= -DBL_MAX;
                               g_mynodemax[g_index]= j;
			       g_maxlookup[g_index]= -1; /** since their is NO parent combination subset of this which is valid */
 }
 g_index++;
}
   /*Rprintf("-----------------------------------------------------\n");
   Rprintf("g(S); [fixed width] S=| ");
   for(j=0;j<numNodes;j++){Rprintf("%d ",g_lookup[g_index][j]);}Rprintf("|\n");
   Rprintf("MAX: node=%d parents=",g_mynodemax[g_index]);for(j=0;j<numNodes;j++){Rprintf("%d ",parents[ g_maxlookup[g_index] ][j]);}Rprintf("\n");
   Rprintf("-----------------------------------------------------\n");
   */
 /** since know that first parent combination in each node is empty set */
/** above is g(0)=alpha_0(empty), g(1)=alpha_1(empty).....g(n-1)=alpha_{n-1}(empty), also this matches the first combination in g_lookup */
/** **************************************************************************************************************************************/

/** **************** NOW FOR EACH S at a time and all S have cardinality>1*****************************************************************/
/** **************************************************************************************************************************************/
/** calculate a single g(S) where S = the set in g_lookup[g_index] and has cardinality g_cardinality[g_index] **/
/** NOTE: S is in form fixed 00110 representation e.g. {2,3} */
/** ***************************************************************************************************************************************/
while(g_index<total){
  if(timingon && g_index%10000==0){startclock = clock();}
  if(g_index%10000==0){if(verbose){Rprintf("processed set %d of %d\n",g_index,total);R_CheckUserInterrupt();};/** allow an interupt from R console */}
myalphamax= -DBL_MAX;/** will hold the value of g(S) **/
minindex=0;/** to enable alpha_index - to get correct index must be reset for each new S */
for(i=0;i<g_cardinality[g_index];i++){/** for each member of current S defined in g_lookup[g_index] */

        /** 1. Find NODE: need alpha_mynode (S^prime) where mynode is the ith member of set S, and set S is stored in g_lookup[g_index] **/
        for(k=minindex;k<numNodes;k++){/** starts at minindex since we want the first "1" and then the second "1" in S*/
	                               if(g_lookup[g_index][k]==1){/** if Kth member of S is a parent then */
					                           mynode=k;
	                                                           minindex=k+1;/** update so on next iteration look for ith node */
								   break;}} /** find the index for alpha e.g. alpha[mynode] **/
	/** we now want to find the set S^prime which is S\{mynode} **/
	for(k=0;k<numNodes;k++){S[k]=g_lookup[g_index][k];} /** copy full set to S **/
        S[mynode]=0;/** remove the alpha index - mynode - from the set S, S is now S^prime **/
        Sprime=S;/** point Sprime to this to avoid confusion */
       /* Sprimecard=g_cardinality[g_index]-1; */
/*	j=0;for(k=0;k<numNodes;k++){if(Sprime[k]==1){Sprime_indexes[j++]=k;}}*/ /** Sprime_indexes contains the col indexes of the parents */
        /** now also create a version of Sprime with child removed not just set to zero - Sprime_short - so length is one less */
        /*j=0;for(k=0;k<numNodes;k++){if(k!=mynode){Sprime_short[j++]=Sprime[k];}}*/

      /*Rprintf("Sprime: node=%d |",mynode);for(j=0;j<numNodes;j++){Rprintf("%d ",Sprime[j]);}Rprintf("\n");*/
       /*for(k=0;k<numRows;k++){for(j=0;j<numNodes;j++){Rprintf("%d ",parents[k][j]);}Rprintf("\n");}*/

	/** 2. Find ALPHA(Sprime): ****************************************************/
	/** need to search over all entries in parents[mynode] to find the best subset of Sprime. Note by construction all combinations in parents
	 * have cardinality <=|maxparents| and satify banned and retain - hence use parents as the source of all combinations and not create
	   combinations to check against parents */

	getAlphaMaxSingle(Sprime,mynode,ptr_score,ptr_nodeid,parents,alphalong,alphalonglookup,numNodes,numRows,&myalpha,&myalphaparentindex,start,end);

	/** 3. FIND g(Sprime) *******************************/
	/** we now want to find the index in g_lookup which matches Sprime and then get the corresponding value for g() - we should have already calculated these **/
	for(k=0;k<g_index;k++){/** for all the g()'s of lower cardinality */
	    curindex=k; /** store for later in case we have a hit */
	    if(g_cardinality[k]== (g_cardinality[g_index]-1)){/** only look in g_lookup[k] where no. parents must be one less than current S */
	    Match=1;
	    for(j=0;j<numNodes;j++){if(g_lookup[k][j]!=Sprime[j]){Match=0;break;}}/** a mismatch to set to false*/
	    if(Match==1){break;}/** found a match so can finish searching and curindex is what we want**/
	    /*lookupindex=k;*//** row in g_lookup which has our parents in it */
	}
	}
	if(Match==0){error("no match for Sprime\n");}

	if(myalpha> -DBL_MAX){
	   myalpha+=g[curindex];}/** calculate alpha_mynode(Sprime)*g(Sprime): using + since log scale **/

	/*Rprintf("alpha[%d]*g=%f\n",mynode,myalpha);*/

	/** 4. get MAX **************************************/
	/*if(myalpha==myalphamax && myalpha!= -DBL_MAX){Rprintf("-high mem--- Note: have multiple possible max-----\n");}*/
	if(myalpha>myalphamax){/** got a new better "parent set" so update score and also save nodeid and Sprime **/
	                       myalphamax=myalpha;/** best value of alpha_mynode(Sprime)*g(Sprime) for the set S - so far**/
			       mynodemax=mynode; /** the node mynode in the current best value **/
	                       myalphaparentindexmax=myalphaparentindex;/** the best SUBSET of Sprime so far **/
	                      /* Rprintf("mymaxbestalpha %d\n",myalphaparentindexmax);*/
	}

}/** for loop over all members of current set S *************/


   /** now store the best term found in g(S) i.e. which of the various alpha_i(Sprime)g(Sprime) terms was MAX */
   g[g_index]=myalphamax;/** store best value - only really useful for debugging */
   g_mynodemax[g_index]=mynodemax;/** get the node */
   g_maxlookup[g_index]=myalphaparentindexmax;/** get the parents of the node - not Sprime but the best subset of Sprime */

  /* Rprintf("-----------------------------------------------------\n");
   Rprintf("g(S); [fixed width] S=| ");
   for(j=0;j<numNodes;j++){Rprintf("%d ",g_lookup[g_index][j]);}Rprintf("|\n");
   if(g_maxlookup[g_index]!= -1 ){
   Rprintf("MAX: node=%d parents=",g_mynodemax[g_index]);for(j=0;j<numNodes;j++){Rprintf("%d ",parents[ g_maxlookup[g_index] ][j]);}Rprintf("\n");
   } else {Rprintf("MAX: node=%d parents=: no valid subsets\n",g_mynodemax[g_index]);}
   Rprintf("-----------------------------------------------------\n");
  */
 g_index++; /** increment to next set S **/

 if(timingon && g_index%10000==0 ){endclock = clock();
 elapsed = ((double) (endclock - startclock)) / CLOCKS_PER_SEC;
 if(verbose){Rprintf("CPU time: %10.6f secs\n",elapsed);};};


} /** end of loop over all subsets of S=0,....max_parents */

/** **********************************************************************************************************/
/** NOW TO ACTUALLY INTERPRET THE RESULTS - we have all of g(S) for all subsets S of 0,...,n-1 but need to   */
/** unpick the results to actually find out what is the most probable DAG ***********************************/

   lastindex=total-1;
   startindex=0;
   curS=S;
   got=0;
   for(i=0;i<numNodes;i++){curS[i]=0;} /** reset */
   /** the first node found is an easy special case **/
   /*for(i=0;i<total;i++){for(j=0;j<numNodes;j++){Rprintf("%d ",g_lookup[i][j]);}Rprintf("| %d--\n",g_mynodemax[i]);}*/
  /* Rprintf("got node=%d [index=%d] score=%10.10f\t|parents=",g_mynodemax[lastindex],lastindex,g[lastindex]);*/
   for(k=0;k<numNodes;k++){final_results[ g_mynodemax[lastindex] ][k]=parents[ g_maxlookup[lastindex] ][k];}
   /** We start with the g(S); S={0,...,n-1} which is easy since this is the LAST term in g[] **/
   mynode=g_mynodemax[lastindex];/** we now have a node, next need its parents */
   /*for(k=0;k<numNodes;k++){Rprintf("%d ",parents[ g_maxlookup[lastindex] ][k]);}Rprintf("\n");*/
   curS[got++]=mynode;
   startindex=1;
   for(l=1;l<numNodes;l++){/** for each node - first node is done */
   /** now loop through g_lookup from the end backwards and find first subset of g() which does NOT contain g_mynodemax[total-0] **/
   for(i=startindex;i<total;i++){/** find node **/
                       /*Rprintf("=%d\n",lastindex-i);*/
		       Match=0;
		       for(j=0;j<got;j++){
                       if(g_lookup[lastindex-i][curS[j]]!=0){Match=1;}}/** have found the next value needed for g(S) */
			                                if(Match==0){
							  /*Rprintf("got node=%d [index=%d] score=%10.10f\t|parents=",g_mynodemax[lastindex-i],lastindex-i,g[lastindex-i]);*/
		                                          /*for(k=0;k<numNodes;k++){Rprintf("%d ",parents[ g_maxlookup[lastindex-i] ][k]);}Rprintf("\n");*/
							   for(k=0;k<numNodes;k++){final_results[ g_mynodemax[lastindex-i] ][k]=parents[ g_maxlookup[lastindex-i] ][k];} /** store final results */
							     mynode=g_mynodemax[lastindex -i];/** now get the node for this new value of g **/
		                                            curS[got++]=mynode;
		                                            startindex=i+1;/** now repeat but start search from next slot up **/
							    break;}
                         }
   }



/** roll final results: i.e. matrix defining the most probable DAG into a single long vector row by row **/
k=0;
for(i=0;i<numNodes;i++){
   for(j=0;j<numNodes;j++){/*Rprintf("%d ",final_results[i][j]);}Rprintf("\n");}*/
                  INTEGER(VECTOR_ELT(listresults,0))[k++]=final_results[i][j];
   }}

UNPROTECT(1);

return(listresults);
/*return(R_NilValue);*/

}

/** ********************************************************************************************************************/




