#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

void make_dag(network *dag, int dim, SEXP R_dag, int empty, SEXP R_vartype, const int *maxparents, SEXP R_groupedvars);
void make_nodecache(cache *cache, int vars,  int cols, int rows, SEXP R_numparcombs, SEXP R_children, SEXP R_cachedefn, SEXP R_nodescores);
int lookupscores(network *dag,cache *nodecache);
void printDAG(network *dag,int what);
void store_results(SEXP R_listresults,network *dag, int iter, SEXP ans, int verbose);
void printDATA(datamatrix *data,int what);
void make_data(SEXP R_obsdata,datamatrix *obsdata, SEXP R_groupids);
void printCACHE(cache *nodecache,int what);
void free_dag(network *dag);
void print_state (int iter, gsl_multiroot_fdfsolver * s);
