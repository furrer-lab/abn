/***************************************************************************
 *   Copyright (C) 2006 by F. I. Lewis   
 *   fraser.lewis@ed.ac.uk   
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/** *********************************************************************** 
 * definitions of structures
 *
 ***************************************************************************/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
/** designed to hold a network definition **/
struct network_struct {
      int **defn;/** each row a variable and each col the parents of the variable, indexes from 0*/
      int *locationInCache;/** gives the "rows" in the cache which make up this DAG **/
      unsigned int numNodes;/** total number of variables/nodes **/
      double *nodeScores;/** array of individual scores **/
      int *nodeScoresErrCode;
      double *hessianError;
      double networkScore;/** log marginal likelihood */
      int *varType;/** holds the type of distribution **/
      int maxparents;/** **/
      gsl_matrix *modes;/** will hold the parameter points estimates from the laplace approx **/
      int *groupedVars;/** indexes of nodes with random effects */
};

typedef struct network_struct network;

/** designed to hold a network definition **/
struct cycle_struct {
     unsigned int *isactive;
     unsigned int *incomingedges;
     unsigned int **graph;
     
};

typedef struct cycle_struct cycle;

/** designed to hold a network definition **/
struct cache_struct {
      int ***defn;/** each row a variable and each col the parents of the variable, indexes from 0*/
      unsigned int numVars;/** total number of variables in a complete DAG**/
      unsigned int numRows;/** total number of different parent combinations **/
      int *numparcombs;/** number of parent combinations per node - number of rows e.g. j in  defn[i][j][k] */
      double **nodeScores;/** hold the score for each parent combination for each node */
};

typedef struct cache_struct cache;




typedef struct diskdatamatrix_struct diskdatamatrix;

/** designed to hold integer data and column names **/
struct datamatrix_struct {
      double **defn;
      gsl_matrix *datamatrix;
      int numDataPts;/** total number of datapoints**/
      int numVars;/** total number of variables/nodes **/
     /* char **namesVars;*//** array of strings denoting node/variable names, in order of data columns */
      int *numVarlevels;/** number of unique categories per variable */
     /*double *weights;*//** hold a double value which is the weight of each observed data point */
     /* double *xmin;
      double *xmax;*/
      unsigned int numparams;
      gsl_vector *priormean;
      gsl_vector *priorsd;
      gsl_vector *priorgamshape;
	    gsl_vector *priorgamscale;
      gsl_vector *gslvec1;
      gsl_vector *gslvec2;
      /*double relerr;*/
      gsl_vector *Y;
      int *groupIDs;
      int numUnqGrps;
      gsl_matrix **array_of_designs;
      gsl_vector **array_of_Y;
      gsl_matrix *datamatrix_noRV;
  
};

typedef struct datamatrix_struct datamatrix;

struct fnparams
       {
         gsl_vector *Y;
	 gsl_vector *vectmp1;
	 gsl_vector *vectmp2;
	 gsl_vector *vectmp1long;
	 gsl_vector *vectmp2long;
	  gsl_vector *vectmp3long;
	 gsl_vector *term1;
	 gsl_vector *term2;
	 gsl_vector *term3;
       gsl_matrix *X;
  gsl_matrix *mattmp1;
  gsl_matrix *mattmp2;
  gsl_matrix *mattmp3;
  gsl_matrix *mattmp4;
	 gsl_vector *priormean;
	 gsl_vector *priorsd;
	 gsl_vector *priorgamshape;
	 gsl_vector *priorgamscale;
	 gsl_vector *betafull;
	 gsl_vector *dgvaluesfull;
	 double betafixed;
	 int betaindex;
	 gsl_vector *dgvalues;
	 gsl_matrix *hessgvalues;
	 gsl_matrix *hessgvalues3pt;
	 gsl_vector *beta;
	 gsl_vector *betastatic;
	 gsl_permutation *perm;
	 double inits_adj;
	 double gvalue;
	 datamatrix *designdata;
	 gsl_vector *betaincTau;
	 gsl_vector *betaincTau_bu;
	 int fixed_index;
	 int fixed_index2;
	 double epsabs_inner;
	 int maxiters_inner;
	 int verbose;
	 double finitestepsize;
	 int nDim;
	 int mDim;
	 double logscore;
	 double logscore3pt;
	
	};
