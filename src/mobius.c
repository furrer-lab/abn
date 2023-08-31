#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf.h>
#include "mobius.h"

/** *********************************************************************************************/
/** *********************************************************************************************/
double qprime_prior(int numNodes, int numparents, int priorchoice)
{
 /** two possible priors - the first priorchoice=1 is the simple uniform prior i.e. all sets of parents are equally likely **/
 /**                       the second prior priorchoice=2 is the one in Koivisto which say that all parent combinations of the same size are equally likely **/
 /**                       the FIRST prior seems better as it can be directly compared with the heuristic hill climber and also, an issue with the Koivisto prior **/
 /**                       is that is gives the same prior weight to no parents as having all possible parents and that can give daft results                    **/
 
  if(priorchoice==1){/** this is a uniform prior which is the same as no structural prior and comparable with random re-start heuristics **/
                    return(1.0);
                    }
 if(priorchoice==2){ /** this is the prior mentioned in Koivisto - all subsets of same size equally likely **/
                    return(1.0/gsl_sf_choose(numNodes-1,numparents));
                   }


 error("should never be here!\n");
 return(-1);
}

/** *********************************************************************************************/
/** *********************************************************************************************/
/** search all subsets of combinations for each node and then find the best subject to ban and retain constraints   */
void getAlphaMax(double *ptr_score,int *ptr_nodeid,int **parents,double *alpha,int *alphalookup, int numNodes, int numRows)
{
 int i,j,mynode;
 int mybestindex=-1;
 int curcomb;
 double myscore,mybestscore;
 
/* for(i=0;i<numRows;i++){
   for(j=0;j<numNodes;j++){
     Rprintf("%d ",(parents[i][j]));
   }Rprintf("\n");}
  */
 /** go through each row - a single parent combination for each node and search for all subsets of that combination 
  *  within the same node, and if any has a better score then select it */
 for(i=0;i<numRows;i++){ /** for the current parent combination */
    /** the current combination */
    curcomb=i;
    /** get its node **/
    mynode=ptr_nodeid[i];
    /** get its score */
    myscore=ptr_score[i];
    /** set the best score over all subsets - initially this is for the set S, where T are all the subsets of S*/
    mybestscore= myscore;
    mybestindex= i;
    for(j=0;j<numRows;j++){/** for each parent combination - all nodes**/
                         if(i!=j && ptr_nodeid[j]==mynode){ /** only interested in those for node mynode and don't check against itself*/
			   /** check if j combination is a subset */
			   if(issubset(parents,curcomb,j,numNodes)){
			      /** if combination j is a subset, and comb j meets ban and retain criteria then is the score better than the best score found so far? */
                                  if(ptr_score[j]>mybestscore){mybestscore=ptr_score[j];/** have a better score **/
				                               mybestindex=j;/** save the location of the best parent combination */
				  }
                           } else {continue;} /** combination j is not a subset of curcomb so skip to next choice */ 
			 }
    }
  /** by here we should have the best score and location of best subset for combination curcomb **/
  alpha[i]=mybestscore;
  alphalookup[i]=mybestindex;    
  }
  
}
/** *********************************************************************************************/
/** *********************************************************************************************/
/** search all subsets of combinations for each node and then find the best subject to ban and retain constraints   */
/** note we do NOT need to check for ban or retain constraints since very combination in parents must satisfy this */
void getAlphaMaxSingle(int *Sprime, int mynode, double *ptr_score,int *ptr_nodeid,int **parents,double *alpha,int *alphalookup,
		       int numNodes, int numRows, double *myalpha, int *myalphaparentindex, int *start, int *end)
{
 int j;
 int mybestindex=-1;
 double mybestscore;
  
 /** for the a single parent combination Sprime search for all subsets of that combination 
  *  within the same node, and if any has a better score then select it */
    /** set the best score over all subsets - initially this is for the set S, where T are all the subsets of S*/
    mybestscore= -DBL_MAX;
    mybestindex= -1; 
    /** we now that for node mynode its combinations start at start[mynode] and end at end[mynode] */
         for(j=start[mynode];j<=end[mynode];j++){/** for each parent combination for node mynode **/
                         if(ptr_nodeid[j]==mynode){ /** only interested in those for node mynode and don't check against itself*/
			    /** we know how many combinations each node has */
			   /** check if j combination is a subset */
			    if(issubset1D(Sprime,parents[j],numNodes)){
			      /** if combination j is a subset, and comb j meets ban and retain criteria then is the score better than the best score found so far? */   
			      if(ptr_score[j]>mybestscore){mybestscore=ptr_score[j];/** have a better score **/
				                           mybestindex=j;/** save the location of the best parent combination */
			/*Rprintf("parents=");for(i=0;i<numNodes;i++){Rprintf("%d ",parents[j][i]);}Rprintf("\n");*/                               
                               }  /** combination j is not a subset of curcomb so skip to next choice */ 
			    }
			}
    }
  /*if(mybestscore== -DBL_MAX){error("no valid subset found in getAlphaMaxSingle!\n");}*/
  /** if we don't find a combination which is subset of the Sprime needed - which is possible e.g. alpha_0(1) if 2->0 is in the retain list for node 0 */
  /** is no match found then score is -DBL_MAX and parentindex = -1  **/
  *myalpha=mybestscore;
  *myalphaparentindex=mybestindex;    
  
  
}
/** *********************************************************************************************/
/** *********************************************************************************************/  
int issubset(int **parents,int curcomb,int newcomb,int numNodes)
{
 /** compare two vectors vec1=parents[curcomb][j] for j=0,..numNodes-1 and vec2=parents[newcomb][j] 
  and find out whether newcomb is a subset of curcomb.*/ 
 int i;
 for(i=0;i<numNodes;i++){/** for each parent */
                 if( parents[newcomb][i] != parents[curcomb][i] ){ /** if any col is different then
		                                                        it must be curcomb=1 and newcomb=0 */
		    if( parents[newcomb][i]==1 && parents[curcomb][i]==0 ){
		    return(0);/** do not have a subset so are done*/ 
		 }  
		}
 }
  /** if get to here then have a subset */ 
     return(1);
   
}  
 /** *********************************************************************************************/
/** *********************************************************************************************/  
int issubset1D(int *Sprime,int *parents,int numNodes)
{
 /** compare two vectors: is parents a subset of Sprime?*/ 
 int i;
 for(i=0;i<numNodes;i++){/** for each parent */
                 if( parents[i] != Sprime[i] ){ /** if any col is different then
		                                                        it must be curcomb=1 and newcomb=0 */
		    if( parents[i]==1 && Sprime[i]==0 ){
		    return(0);/** do not have a subset so are done*/ 
		 }  
		}
 }
  /** if get to here then have a subset */ 
     return(1);
   
}

